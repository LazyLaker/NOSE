
rm(list = ls())
set.seed(123456)
M <- 300
n <- 512
Dat <- matrix(0, M, n)
for(m in 1:M){
  y <- c()
  y[1:64] <- rnorm(64, mean = 0, sd = sqrt(3))
  y[65:128] <- rnorm(64, mean = 1.5, sd = sqrt(3))
  y[129:192] <- rnorm(64, mean = 0, sd = sqrt(3))
  y[193:256] <- rnorm(64, mean = 1.5, sd = sqrt(3))
  y[257:320] <- rnorm(64, mean = 0, sd = sqrt(3))
  y[321:384] <- rnorm(64, mean = 1.5, sd = sqrt(3))
  y[385:448] <- rnorm(64, mean = 0, sd = sqrt(3))
  y[449:512] <- rnorm(64, mean = 1.5, sd = sqrt(3))
  Dat[m, ] <- y
}

Dat <- as.data.frame(Dat)

save(Dat, file = 'Teeth7CP.RData')
