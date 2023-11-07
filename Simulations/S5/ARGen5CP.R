
rm(list = ls())
set.seed(1234567)
M <- 300
n <- 450
Dat <- matrix(0, M, n)

for(m in 1:M){
  y <- numeric(n)
  y[1] <- rnorm(1)
  
  for(i in 2:n){
    if(i <= 50){
      y[i] <- 0.5*y[i-1] + rnorm(1)
    }
    else{
      if( i <= 100){
        y[i] <- -0.5*y[i-1] + rnorm(1)
      }
      else{
        if(i <= 200){
          y[i] <- 0.65*y[i-1] + rnorm(1)
        }
        else{
          if(i<= 300){
            y[i] <- -0.25*y[i-1] + rnorm(1)
          }
          else{
            if(i <= 400){
              y[i] <- -0.85*y[i-1] + rnorm(1)
            }
            else{
              y[i] <- 0.45*y[i-1] + rnorm(1)
            }
          }
        }
      }
    }
  }
  
  
  
  Dat[m, ] <- y
}

Dat <- as.data.frame(Dat)

save(Dat, file = 'ARData5CP.RData')