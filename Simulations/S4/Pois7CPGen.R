rm(list = ls())
set.seed(123456)
M <- 300


n <- 512 
tau <- c(64, 128, 192, 256, 320, 384, 448)
seg_len <- diff(c(0, tau, n))

Dat <- matrix(NA, M, n)

pois_t <- c(rep(1, 64), rep(.25, 64), rep(2, 64), rep(1, 64), rep(3, 64), 
            rep(1.5, 64), rep(2.5, 64), rep(1, 64))

for(m in 1:M){
y <- rpois(n, pois_t)
Dat[m, ] <- y
}



Dat <- as.data.frame(Dat)

save(Dat, file = 'PoisTeeth7CP.RData')