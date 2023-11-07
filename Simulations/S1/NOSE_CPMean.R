

NOSE_meanCP <- function(y, ### data
                        minD = 2, ### minimum distance
                        L = 25, ### turncation number
                        nchains = 4, ### num of MCMC chains
                        niter = 12000, ### MCMC runs
                        nburnin = 2000, ### burnin
                        thin = 10 ### MCMC thining
                        ){
  
  N <- length(y)
  ### nimble code
  NimbleCode <- nimbleCode({
    for (i in 1:N) {
      y[i] ~ dnorm(mu_y[i], psi)
      mu_y[i] <- res[i]
      res[i] <- sum(res0[1:L, i])
    }  
    
    ### IBP weighted spike-and-slab
    for (l in 1:L) {
      theta[l] ~ dunif(0, N) 
      p[l] ~ dbeta(alpha, 1)
      z[l] ~ dbern(prod(p[1:l]))
      v0[l] ~ dt(mu = 0, tau = 1, df = 1) ### cauchy slab
      h0[l] <-  v0[l] * z[l] 
    }
    
    for(i in 1:N){
      for(l in 1:L){
        res0[l, i] <-  (theta[l] <= i) * h0[l]
      }
    }
    psi ~ dgamma(1, 1) 
    alpha ~ dgamma(.0001, .001) ### gamma hyperprior
  })
  
  
  para <- list(L = L, niter = niter, nburnin = nburnin, thin = thin)
  
  CP_detect <- function(mcmc_res, minD){
    
    mymodel.samples <- as.data.frame(mcmc_res)
    res.samples <- mymodel.samples[, grep("^res",paste0("^",colnames(mymodel.samples), sep=""), fixed=T)]
    ess <- mean(apply(res.samples, 2, effectiveSize)) ### extract ess
    res.sim <- apply(res.samples, 2, posterior.mode) ### pointwise MAPs
    
    
    step.0 <- diff(res.sim)
    
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
    
    count.1 <- length(cp_loc) ### number of CPs
 
    #### nuisance parameters estimation 
    psi.samples <- mymodel.samples[, grep("^psi",paste0("^",colnames(mymodel.samples), sep=""), fixed=T)]
    psi.est <- mean(psi.samples)
    
    result.cp <- list("Number of CPs" = count.1, 
                      "location" = cp_loc, "piecewise mean" = cp_sim,  
                      "variance estimate" = 1/psi.est)
    z <- list(result = result.cp, 
              samples = mymodel.samples, 
              effectiveSize = ess)
    return(z)
  }
  
  NOSErun <- function(seed, data, code, para){
    y <- data
    N <- length(y)
    L <- para$L
    niter <- para$niter
    nburnin <- para$nburnin
    thin <- para$thin
    NimbleCode <- NimbleCode
    NimbleData <- list(y = y)
    Nimbleconsts <- list(L = L, N = N)
    NimbleInits <- list(psi = 1, alpha = 1, p = rbeta(L, 1, 1), 
                        v0 = rcauchy(L, 0, 1))
    NOSEmean <- nimbleModel(code = NimbleCode, name = "NOSE", constants = Nimbleconsts,
                      data = NimbleData, inits = NimbleInits)
    compNOSE <- compileNimble(NOSEmean)
    NOSEconf <- configureMCMC(NOSEmean)
    
    #### RJMCMC
    NOSEconf$addMonitors(c("psi", "res", 'z', 'v0'))
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
  
  #### parallel MCMC
  mcmc_res <- foreach(i = 1:nchains,
                .combine = 'rbind',
                .packages = c('nimble')) %dopar% {NOSErun(seed = i, data = y, 
                                                          code = NimbleCode, para = para)}
  stopCluster(cl)
  
  rm(cl)
  gc()
  
  time <- mean(mcmc_res[, which(colnames(mcmc_res) == 'time')])
  
  res_NOSE <- CP_detect(mcmc_res, minD = minD)
  
  PACE <- time/res_NOSE$effectiveSize
  
  return(list(CP_detect = res_NOSE$result, time = time, PACE = PACE))

}
