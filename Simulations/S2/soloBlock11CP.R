

library(solocp)
tau <- c(205, 267, 308, 472, 512, 820, 902, 1332, 1557, 1598, 1659)
thresh <- 15
n <- 2048
M <- 300
res <- c()

Hausdorff <- function(a, b, n){
  
  s1 <- sapply(a, FUN = function(x)min(abs(x-b)))
  s2 <- sapply(b, FUN = function(x)min(abs(b-x)))
  return(max(s1, s2)/n)
}

res <- c()
for(m in 1:M){
  y <- as.numeric(Dat[m, ])
  
  res_solo <- solocp_single(y, sigma = 10)
  
  etau_solo <- subset_changepoints(res_solo$ratio)
  
  num_solo <- length(etau_solo)
  
  res_soloTP <- sapply(etau_solo , FUN = function(x) ifelse(abs(x - tau)<=thresh, 1, 0))
  
  solo_TP <- sum(apply(res_soloTP, 1, FUN = function(x) ifelse(sum(x) > 0, 1, 0)))
  
  solo_FP <- max(0, length(etau_solo)- solo_TP)
  
  Hausdorff_solo <- Hausdorff(etau_solo, tau, n)
  
  res <- rbind(res, c(num_solo, solo_FP/num_solo, solo_TP/11, Hausdorff_solo))
  
}