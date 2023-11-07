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
library(wbs)
library(mosum)
load('Teeth7CP.RData')

source('NOSE_CPMean.R')

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
  res_NOSE <- NOSE_meanCP(y, minD = 10)
  res_not <- features(not(y, contrast = "pcwsConstMean"))
  res_PELT <- cpt.mean(y, penalty = 'MBIC', method = 'PELT')
  res_mosum <- mosum(y, G = 40)
  res_wbs <- wbs(y)
  res_SMUCE <- smuceR(y)
  
  
  Num_NOSE <- res_NOSE$CP_detect$`Number of CPs`
  Num_not <- length(res_not$cpt)
  Num_SMUCE <- length(res_SMUCE$leftIndex[-1])
  Num_PELT <- length(res_PELT@cpts) - 1
  Num_wbs <- length(res_wbs$cpt$cpt.th[[1]])
  Num_mosum <- length(res_mosum$cpts.info$cpts)
  
  etau_NOSE <- res_NOSE$CP_detect$location
  etau_not <- res_not$cpt
  etau_SMUCE <- res_SMUCE$leftEnd[-1]
  etau_PELT <- res_PELT@cpts[-length(res_PELT@cpts)]
  etau_wbs <- res_wbs$cpt$cpt.th[[1]]
  etau_mosum <- res_mosum$cpts.info$cpts
  
  
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
  
  
  
  resPELT_TP <- sapply(etau_PELT , FUN = function(x) ifelse(abs(x - tau)<=thresh, 1, 0))
  
  PELT_TP <- sum(apply(resPELT_TP, 1, FUN = function(x) ifelse(sum(x) > 0, 1, 0)))
  
  PELT_FP <- max(0, length(etau_PELT)- PELT_TP)
  
  Hausdorff_PELT <- Hausdorff(etau_PELT, tau, n)
  
  
  
  
  resSMUCE_TP <- sapply(etau_SMUCE , FUN = function(x) ifelse(abs(x - tau)<=thresh, 1, 0))
  
  SMUCE_TP <- sum(apply(resSMUCE_TP, 1, FUN = function(x) ifelse(sum(x) > 0, 1, 0)))
  
  SMUCE_FP <- max(0, length(etau_SMUCE)- SMUCE_TP)
  
  Hausdorff_SMUCE <- Hausdorff(etau_SMUCE, tau, n)
  
  
  
  
  reswbs_TP <- sapply(etau_wbs , FUN = function(x) ifelse(abs(x - tau)<=thresh, 1, 0))
  
  wbs_TP <- sum(apply(reswbs_TP, 1, FUN = function(x) ifelse(sum(x) > 0, 1, 0)))
  
  wbs_FP <- max(0, length(etau_wbs)- wbs_TP)
  
  Hausdorff_wbs <- Hausdorff(etau_wbs, tau, n)
  
  
  
  
  resmosum_TP <- sapply(etau_mosum , FUN = function(x) ifelse(abs(x - tau)<=thresh, 1, 0))
  
  mosum_TP <- sum(apply(resmosum_TP, 1, FUN = function(x) ifelse(sum(x) > 0, 1, 0)))
  
  mosum_FP <- max(0, length(etau_mosum)- mosum_TP)
  
  Hausdorff_mosum <- Hausdorff(etau_mosum, tau, n)
  
  
  end <- Sys.time()
  cat(m, '\n')
  cat(end-start, '\n')
  
  
  Num <- c(Num_NOSE, Num_not, Num_SMUCE, Num_PELT, Num_wbs, Num_mosum)
  Prec <- c(NOSE_TP/(Num_NOSE), 
            not_TP/(Num_not),
            (SMUCE_TP)/(Num_SMUCE), PELT_TP/(Num_PELT), 
            wbs_TP/Num_wbs, mosum_TP/Num_mosum)
  recall <- c(NOSE_TP/7, not_TP/7, SMUCE_TP/7, PELT_TP/7, wbs_TP/7, mosum_TP/7)
  
  
  Haus <- c(Hausdorff_NOSE, Hausdorff_not, Hausdorff_SMUCE, 
            Hausdorff_PELT, Hausdorff_wbs, Hausdorff_mosum)
  
  res <- rbind(res, c(Num, Prec, recall, Haus, time_NOSE, PACE_NOSE))
  res <- data.frame(res)
  
  names(res) <- c('Num_NOSE', 'Numnot', 'Num_SMUCE', 
                  'Num_PELT', 'Num_wbs', 'Num_mosum',
                  'Prec_NOSE', 'Prec_not', 'Prec_SMUCE', 
                  'Prec_PELT', 'Prec_wbs', 'Prec_mosum',
                  'Recall_NOSE', 'Recall_not', 'Recall_SMUCE', 
                  'Recall_PELT', 'Recall_wbs', 'Recall_mosum',
                  'Haus_NOSE', 'Haus_not', 'Haus_SMUCE', 
                  'Haus_PELT', 'Haus_wbs', 'Haus_mosum', 
                  'time_NOSE', 'PACE_NOSE'
  )
  

  
}


  
  
