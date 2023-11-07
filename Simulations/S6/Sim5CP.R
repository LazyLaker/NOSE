
rm(list = ls())
library(parallel)
library(doParallel)
library(foreach)
library(nimble)
library(coda)
library(dplyr)
library(purrr)
library(MCMCglmm)
library(flexclust)
load('Reg5CPy.RData')
load('Reg5CPx.RData')
source('NOSE_CPlinear.R')


tau <- c(40, 80, 120, 160, 200)
thresh <- 15

Hausdorff <- function(a, b, n){
  
  s1 <- sapply(a, FUN = function(x)min(abs(x-b)))
  s2 <- sapply(b, FUN = function(x)min(abs(b-x)))
  return(max(as.numeric(s1), as.numeric(s2))/n)
}

res <- c()

for(m in 1:300){
  
  y <- as.numeric(yDat[m, ])
  X <- as.numeric(xDat[m, ])
  n <- 480
  t <- c()
  for(i in 1:n){
    t <- c(t, ceiling(i/2))
  }
  
  res_NOSE <- NOSE_LinearCP(y, X, t, minD = 10, niter = 14000)
  
  Num_NOSE <- res_NOSE$CP_detect$`Number of CPs`
  
  etau_NOSE <- res_NOSE$CP_detect$location
  
  
  resNOSE_TP <- sapply(etau_NOSE , FUN = function(x) ifelse(abs(x - tau)<=thresh, 1, 0))
  
  NOSE_TP <- sum(apply(resNOSE_TP, 1, FUN = function(x) ifelse(sum(x) > 0, 1, 0)))
  
  NOSE_FP <- max(0, length(etau_NOSE)- NOSE_TP)
  
  Hausdorff_NOSE <- Hausdorff(etau_NOSE, tau, n)
  
  time_NOSE <- res_NOSE$time
  PACE_NOSE <- res_NOSE$PACE
  
  
  end <- Sys.time()
  cat(m, '\n')
  
  Num <- c(Num_NOSE)
  Prec <- c(NOSE_TP/Num_NOSE)
  recall <- c(NOSE_TP/5)
  Haus <- c(Hausdorff_NOSE)
  
  res <- rbind(res, c(Num, Prec, recall, Haus, time_NOSE, PACE_NOSE))
  res <- data.frame(res)  
  colnames(res) <- c('Num_NOSE', 
                     'Prec_NOSE',
                     'Rec_NOSE', 
                     'Haus_NOSE')

}

