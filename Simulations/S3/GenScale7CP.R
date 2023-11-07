rm(list = ls())
set.seed(123456)
M <- 300


n <- 512 
tau <- c(64, 128, 192, 256, 320, 384, 448)
seg_len <- diff(c(0, tau, n))

Dat <- matrix(NA, M, n)

scale_t <- c(rep(1, 64), rep(.5, 64), rep(1, 64), rep(1.5, 64), rep(1, 64), 
             rep(.5, 64), rep(1, 64), rep(.5, 64))

for(m in 1:M){
  y <- rnorm(n, 1, scale_t)
  Dat[m, ] <- y
}



Dat <- as.data.frame(Dat)

save(Dat, file = 'ScaleTeeth7CP.RData')