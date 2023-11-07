rm(list = ls())
set.seed(1234567)
M <- 300
t <- c()
n <- 480
for(i in 1:n){
  t <- c(t, rep(i, 2))
}

beta1 <- c(rep(1, 80), rep(-1, 80), rep(.5, 80), rep(-.5, 80), rep(1, 80), rep(-1, 80))
beta0 <- .5
yDat <- matrix(0, M, n)
xDat <- matrix(0, M, n)

for(m in 1:M){
  x <- runif(n, -2, 2)
  epsilon <- rnorm(n, 0, 1)
  for(i in 1:n){
    yDat[m, i] <- beta0 + beta1[i] * x[i] + epsilon[i]
  }
  xDat[m, ] <- x
}

yDat <- data.frame(yDat)
xDat <- data.frame(xDat)

save(yDat, file = 'Reg5CPy.RData')
save(xDat, file = 'Reg5CPx.RData')





