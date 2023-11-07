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
library(wbsts)
library(strucchange)
library(mosum)
library(SNSeg) #### kindly provided by the author jiangfy@fudan.edu.cn
load('ARData5CP.RData')

source('NOSE_CPLinear.R')

tau <- c(50, 100, 200, 300, 400)
thresh <- 10
n <- 450
res <- c()

Hausdorff <- function(a, b, n){
  
  s1 <- sapply(a, FUN = function(x)min(abs(x-b)))
  s2 <- sapply(b, FUN = function(x)min(abs(b-x)))
  return(max(s1, s2)/n)
}



for(m in 110:nrow(Dat)){
  start <- Sys.time()
  y <- as.numeric(Dat[m, ])
  X <- y[-n]
  y_res <- y[-1]
  N <- length(y_res)
  t <- 1:(n-1)
  res_NOSE <- NOSE_LinearCP(y_res, X, t, minD = 10)
  res_BP <-  breakpoints(y_res~X)
  #### SNseg result
  res_SNCP <- SNSeg_Uni(y, paras_to_test = 'acf', plot_SN = F)

  
  
  Num_NOSE <- res_NOSE$CP_detect$`Number of CPs`
  # Num_wbs <- length(res_wbs$cp.aft)
  Num_SNCP <- length(res_SNCP$est_cp)
  Num_BP <- length(res_BP$breakpoints)
  
  
  etau_NOSE <- res_NOSE$CP_detect$location
  # etau_wbs <- res_wbs$cp.aft
  etau_SNCP <- res_SNCP$est_cp
  etau_BP <- res_BP$breakpoints
  
  resNOSE_TP <- sapply(etau_NOSE , FUN = function(x) ifelse(abs(x - tau)<=thresh, 1, 0))
  
  NOSE_TP <- sum(apply(resNOSE_TP, 1, FUN = function(x) ifelse(sum(x) > 0, 1, 0)))
  
  NOSE_FP <- max(0, length(etau_NOSE)- NOSE_TP)
  
  Hausdorff_NOSE <- Hausdorff(etau_NOSE, tau, n)
  
  time_NOSE <- res_NOSE$time
  PACE_NOSE <- res_NOSE$PACE

  
  # if(is.null(etau_wbs)){
  #   wbs_TP <- 0 
  #   wbs_FP <- 0
  #   Hausdorff_wbs <- 0
  # }
  # else{
  #   reswbs_TP <- sapply(etau_wbs , FUN = function(x) ifelse(abs(x - tau)<=thresh, 1, 0))
  #   
  #   wbs_TP <- sum(apply(reswbs_TP, 1, FUN = function(x) ifelse(sum(x) > 0, 1, 0)))
  #   
  #   wbs_FP <- max(0, length(etau_wbs)- wbs_TP)
  #   
  #   Hausdorff_wbs <- Hausdorff(etau_wbs, tau, n)
  # }
  
  
  if(is.null(etau_SNCP)){
    SNCP_TP <- 0 
    SNCP_FP <- 0
    Hausdorff_SNCP <- 0
  }
  
  else{
  
  resSNCP_TP <- sapply(etau_SNCP , FUN = function(x) ifelse(abs(x - tau)<=thresh, 1, 0))
  
  SNCP_TP <- sum(apply(resSNCP_TP, 1, FUN = function(x) ifelse(sum(x) > 0, 1, 0)))
  
  SNCP_FP <- max(0, length(etau_SNCP)- SNCP_TP)
  
  Hausdorff_SNCP <- Hausdorff(etau_SNCP, tau, n)
  
  }
  
  if(is.null(etau_BP)){
    BP_TP <- 0 
    BP_FP <- 0
    Hausdorff_BP <- 0
  }
  else{
  
  resBP_TP <- sapply(etau_BP , FUN = function(x) ifelse(abs(x - tau)<=thresh, 1, 0))
  
  BP_TP <- sum(apply(resBP_TP, 1, FUN = function(x) ifelse(sum(x) > 0, 1, 0)))
  
  BP_FP <- max(0, length(etau_BP)- BP_TP)
  
  Hausdorff_BP <- Hausdorff(etau_BP, tau, n)
  }
  
  
  end <- Sys.time()
  cat(m, '\n')
  cat(end-start, '\n')
  
  
  Num <- c(Num_NOSE,  Num_SNCP, Num_BP)
  Prec <- c(NOSE_TP/(Num_NOSE), 
            SNCP_TP/Num_SNCP, BP_TP/Num_BP)
  recall <- c(NOSE_TP/5,  SNCP_TP/5, BP_TP/5)
  
  
  Haus <- c(Hausdorff_NOSE,  Hausdorff_SNCP, 
            Hausdorff_BP)
  
  res <- rbind(res, c(Num, Prec, recall, Haus, time_NOSE, PACE_NOSE))
  res <- data.frame(res)
  
  names(res) <- c('Num_NOSE',  'Num_SNCP', 
                  'Num_BP', 
                  'Prec_NOSE', 'Prec_SNCP', 
                  'Prec_BP', 
                  'Recall_NOSE', 'Recall_SNCP', 
                  'Recall_BP', 
                  'Haus_NOSE', 'Haus_SNCP', 
                  'Haus_BP', 
                  'time_NOSE', 'PACE_NOSE'
  )
  

  
}


  
  
