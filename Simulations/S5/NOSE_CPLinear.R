# library(doParallel)
# library(foreach)
# library(nimble)
# library(dplyr)
# library(purrr)
# library(MCMCglmm)
# library(coda)


NOSE_LinearCP <- function(y, ### response
                          X, ### covariate
                          t, ### states
                          minD = 2, ### minimum distance
                          L = 25, ### truncation number
                        nchains = 4, ### number independent MCMC chains
                        niter = 12000, ### length of each MCMC chain
                        nburnin = 2000, ### burn in 
                        thin = 10 ### MCMC thining
                        ){
  N <- length(y)
  
  #### nimble code
  NimbleCode <- nimbleCode({
    for (i in 1:N) {
      y[i] ~ dnorm(mu_y[i], psi)
      mu_y[i] <- beta1[t[i]]*X[i] + beta0
    }  
    
    for (l in 1:L) {
      theta[l] ~ dunif(0, TN) 
      p[l] ~ dbeta(alpha, 1)
      z[l] ~ dbern(prod(p[1:l]))
      v0[l] ~ dt(mu = 0, tau = 1, df = 1)
      h0[l] <-  v0[l] * z[l] 

    }
    
    for(q in 1:TN){
      for(l in 1:L){
        res0[l, q] <-  (theta[l] <= ut[q]) * h0[l]
      }
      beta1[q] <- sum(res0[1:L, q])
    }
    psi ~ dgamma(1, 1) 
    alpha ~ dgamma(.0001, .001)
    beta0 ~ dnorm(0, sd = 100)
    
  })
  
  
  para <- list(L = L, niter = niter, nburnin = nburnin, thin = thin)
  
  
  ##### discrimination step
  CP_detect <- function(mcmc_res, minD){
    
    mymodel.samples <- as.data.frame(mcmc_res)
    beta1.samples <- mymodel.samples[, grep("^beta1",paste0("^",colnames(mymodel.samples), sep=""), fixed=T)]
    
    ess <- mean(apply(beta1.samples, 2, effectiveSize)) ### extrcat ess
  
    res.sim <- apply(beta1.samples, 2, posterior.mode) ### pointwise MAP
    
    step.0 <- diff(res.sim)
    
   #sd_diff <- sd(step.0)
    sd_diff <- sd(step.0[which(step.0 > quantile(step.0, .0005) & step.0 < quantile(step.0, .9995))])
    step.1 <- ifelse(abs(step.0)  > (3*(sd_diff)), step.0, 0)  
    loc.1 <- which(step.1 != 0)
    
    res.1 <- res.sim[loc.1+1]
    
    count <- length(res.1)
    
    cp_loc <- loc.1
    
    knot <- c(1, cp_loc, N)
    knot_temp <- knot
    threshold <- minD
    for(j in 2:length(knot)){
      if((knot[j] - knot[(j-1)]) <= threshold){
        knot_temp[j-1] <- NA
      }
    }
    knot <- c(na.omit(knot_temp))
    if(max(knot) == N){
      cp_loc <- knot[2:(length(knot)-1)]
    }
    else{
      cp_loc <- knot[-1]
    }
    cp_sim <- c()
    for(j in 2:length(knot)){
      cp_sim  <- c(cp_sim, mean(res.sim[knot[j-1] : knot[j]]))
    }
    
    count.1 <- length(cp_loc) ### number of change-points
    
    
    #### nuisance parameter estiamtion
    psi.samples <- mymodel.samples[, grep("^psi",paste0("^",colnames(mymodel.samples), sep=""), fixed=T)]
    psi.est <- mean(psi.samples)
    
    beta0.samples <- mymodel.samples[, grep("^beta0",paste0("^",colnames(mymodel.samples), sep=""), fixed=T)]
    beta0.est <- mean(beta0.samples)
    
    
    result.cp <- list("Number of CPs" = count.1, 
                      "location" = cp_loc, "beta1" = cp_sim,  
                      'beta0' = beta0.est, 
                      "variance estimate" = 1/psi.est)
    z <- list(result = result.cp, 
              samples = mymodel.samples, 
              effectiveSize = ess)
    return(z)
  }
  
  NOSErun <- function(seed, data, t, code, para){
    y <- data
    N <- length(y)    #### data size
    TN <- length(unique(t))     ### length of t
    L <- para$L
    niter <- para$niter
    nburnin <- para$nburnin
    thin <- para$thin
    NimbleCode <- NimbleCode
    NimbleData <- list(y = y, X=X, ut = c(1:TN), t = t)
    Nimbleconsts <- list(L = L, N = N, TN = TN)
    NimbleInits <- list(psi = 1, alpha = 1, 
                        p = rbeta(L, 1, 1), 
                        v0 = rcauchy(L, 0, 1), beta0 = rnorm(1))
    NOSEmean <- nimbleModel(code = NimbleCode, name = "NOSE", constants = Nimbleconsts,
                            data = NimbleData, inits = NimbleInits)
    compNOSE <- compileNimble(NOSEmean)
    NOSEconf <- configureMCMC(NOSEmean)
    
    ##### RJMCMC
    NOSEconf$addMonitors(c("psi", "beta1", 'beta0', 'z', 'v0'))
    configureRJ(NOSEconf, targetNodes = 'v0', indicatorNodes = 'z', control = list(mean = 0, scale = .2))
    
    NOSEmcmc <- buildMCMC(NOSEconf)
    compNOSEmcmc <- compileNimble(NOSEmcmc, project = NOSEmean)
    compNOSE$setInits(NimbleInits)
    start <- Sys.time()
    mcmc.out <- runMCMC(compNOSEmcmc, niter = niter, setSeed = seed, thin = thin, nburnin = nburnin, nchains = 1)  
    end <- Sys.time()
    time <- difftime(end, start, units = 'secs')[[1]]
    time <- rep(time, (niter - nburnin)/thin)
    mcmc.out <- cbind(mcmc.out, time)
    return(mcmc.out)
  }
  
  
  
  cl <- makeCluster(4)
  registerDoParallel(cl)
  
  #### paralell sampling
  mcmc_res <- foreach(i = 1:nchains,
                      .combine = 'rbind',
                      .packages = c('nimble')) %dopar% {NOSErun(seed = i, data = y, t=t, 
                                                                code = NimbleCode, para = para)}
  stopCluster(cl)
  
  rm(cl)
  
  time <- mean(mcmc_res[, which(colnames(mcmc_res) == 'time')])
  
  res_NOSE <- CP_detect(mcmc_res, minD = minD)
  
  PACE <- time/res_NOSE$effectiveSize
  
  return(list(CP_detect = res_NOSE$result, time = time, PACE = PACE))
  
}
