rm(list = ls())

library(parallel)
library(doParallel)
library(nimble)
library(not)
library(dplyr)
library(purrr)
library(MCMCglmm)
library(stepR)
library(flexclust)
library(changepoint)
library(gfpop)
load('PoisTeeth7CP.RData')

source('NOSE_CPPois.R')

tau <- c(64, 128, 192, 256, 320, 384, 448)
thresh <- 15
n <- 512
res <- c()

Hausdorff <- function(a, b, n){
  
  s1 <- sapply(a, FUN = function(x)min(abs(x-b)))
  s2 <- sapply(b, FUN = function(x)min(abs(b-x)))
  return(max(s1, s2)/n)
}



for(m in 1:nrow(Dat)){
  start <- Sys.time()
  y <- as.numeric(Dat[m, ])
  res_NOSE <- NOSE_PoisCP(y, minD = 10)
  res_not <- features(not(y, contrast = "pcwsConstMeanVar"))
  res_PELT <- cpt.meanvar(y, penalty = 'MBIC', method = 'PELT')
  res_SMUCE <- smuceR(y, family = 'poisson')
  myGraph <- graph(penalty = 2 * sdDiff(y)^2 * log(n), type = "updown")
  res_fpop <- gfpop(data = y, mygraph = myGraph, type = "poisson")
  
  
  
  Num_NOSE <- res_NOSE$CP_detect$`Number of CPs`
  Num_not <- length(res_not$cpt)
  Num_PELT <- length(res_PELT@cpts) - 1
  Num_SMUCE <- length(res_SMUCE$leftIndex[-1])
  Num_fpop <- length(res_fpop$changepoints[-1])
  

  etau_NOSE <- res_NOSE$CP_detect$location
  etau_not <- res_not$cpt
  etau_PELT <- res_PELT@cpts[-length(res_PELT@cpts)]
  etau_SMUCE <- res_SMUCE$leftEnd[-1]
  etau_fpop <- sort(res_fpop$changepoints, decreasing = T)[-1]
  
  
  resNOSE_TP <- sapply(etau_NOSE , FUN = function(x) ifelse(abs(x - tau)<=thresh, 1, 0))
  
  NOSE_TP <- sum(apply(resNOSE_TP, 1, FUN = function(x) ifelse(sum(x) > 0, 1, 0)))
  
  NOSE_FP <- max(0, length(etau_NOSE)- NOSE_TP)
  
  Hausdorff_NOSE <- Hausdorff(etau_NOSE, tau, n)
  
  time_NOSE <- res_NOSE$time
  PACE_NOSE <- res_NOSE$PACE
  
  
  resnot_TP <- sapply(etau_not , FUN = function(x) ifelse(abs(x - tau)<=thresh, 1, 0))
  
  not_TP <- sum(apply(resnot_TP, 1, FUN = function(x) ifelse(sum(x) > 0, 1, 0)))
  
  not_FP <- max(0, length(etau_not)- not_TP)
  
  Hausdorff_not <- Hausdorff(etau_not, tau, n)
  
  if(Num_PELT == 0){
    PELT_TP <- 0
    PELT_FP <- 0
    Hausdorff_PELT <- 0
  }
  
  else{
    resPELT_TP <- sapply(etau_PELT , FUN = function(x) ifelse(abs(x - tau)<=thresh, 1, 0))
    
    PELT_TP <- sum(apply(resPELT_TP, 1, FUN = function(x) ifelse(sum(x) > 0, 1, 0)))
    
    PELT_FP <- max(0, length(etau_PELT)- PELT_TP)
    
    Hausdorff_PELT <- Hausdorff(etau_PELT, tau, n)
  }
  
  
  resSMUCE_TP <- sapply(etau_SMUCE , FUN = function(x) ifelse(abs(x - tau)<=thresh, 1, 0))
  
  SMUCE_TP <- sum(apply(resSMUCE_TP, 1, FUN = function(x) ifelse(sum(x) > 0, 1, 0)))
  
  SMUCE_FP <- max(0, length(etau_SMUCE)- SMUCE_TP)
  
  Hausdorff_SMUCE <- Hausdorff(etau_SMUCE, tau, n)
  
  
  resfpop_TP <- sapply(etau_fpop , FUN = function(x) ifelse(abs(x - tau)<=thresh, 1, 0))
  
  fpop_TP <- sum(apply(resfpop_TP, 1, FUN = function(x) ifelse(sum(x) > 0, 1, 0)))
  
  fpop_FP <- max(0, length(etau_fpop)- fpop_TP)
  
  Hausdorff_fpop <- Hausdorff(etau_fpop, tau, n)
  
  
  
  
  end <- Sys.time()
  cat(m, '\n')
  cat(end-start, '\n')
  
  
  Num <- c(Num_NOSE, Num_not, Num_PELT, Num_SMUCE, Num_fpop)
  Prec <- c(NOSE_TP/(Num_NOSE), 
            not_TP/(Num_not),
            PELT_TP/(Num_PELT), SMUCE_TP/Num_SMUCE, fpop_TP/Num_SMUCE)
  recall <- c(NOSE_TP/7, not_TP/7, PELT_TP/7, SMUCE_TP/7, fpop_TP/7)
  
  
  Haus <- c(Hausdorff_NOSE, Hausdorff_not, 
            Hausdorff_PELT, Hausdorff_SMUCE, Hausdorff_fpop)
  
  res <- rbind(res, c(Num, Prec, recall, Haus, time_NOSE, PACE_NOSE))
  res <- data.frame(res)
  
  names(res) <- c('Num_NOSE', 'Numnot', 
                  'Num_PELT', 'Num_SMUCE', 'Num_fpop',
                  'Prec_NOSE', 'Prec_not', 
                  'Prec_PELT', 'Prec_SMUCE', 'Prec_fpop', 
                  'Recall_NOSE', 'Recall_not', 
                  'Recall_PELT', 'Recall_SMUCE', 'Recall_fpop',
                  'Haus_NOSE', 'Haus_not', 
                  'Haus_PELT', 'Haus_SMUCE', 'Haus_fpop', 
                  'time_NOSE', 'PACE_NOSE'
  )
  

  
}




