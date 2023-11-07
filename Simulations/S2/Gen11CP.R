
rm(list = ls())
set.seed(123456)
M <- 300
n <- 2048
tau <- c(205, 267, 308, 472, 512, 820, 902, 1332, 1557, 1598, 1659)

seg_len <- diff(c(0, tau, n))

mu <- c(0, 14.64, -3.66, 7.32, -7.32, 10.98, -4.39, 3.29, 19.03, 7.68, 15.37, 0)

Dat <- matrix(0, M, n)
for(m in 1:M){
  y <- c()
  
  for(j in 1:12){
    y <- c(y, rnorm(seg_len[j], mu[j], 10))
  }
  
  Dat[m, ] <- y
}

Dat <- as.data.frame(Dat)

save(Dat, file = 'BlockMean11CP.RData')
