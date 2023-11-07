library(StepSignalMargiLike)

load('BlockMean11CP.RData')

n <- 2048
M <- 300
tau <- c(205, 267, 308, 472, 512, 820, 902, 1332, 1557, 1598, 1659)
thresh <- 15
minD <- 10

res <- c()


Hausdorff <- function(a, b, n){
  
  s1 <- sapply(a, FUN = function(x)min(abs(x-b)))
  s2 <- sapply(b, FUN = function(x)min(abs(b-x)))
  return(max(s1, s2)/n)
}

for(m in 1:M){

prior <- prior.norm.A(as.numeric(Dat[m, ]))
etau0 <- est.changepoints(as.numeric(Dat[m, ]), model = 'normal', prior = prior, max.segs = 25)


knot <- c(1, etau0, n)
knot_temp <- knot
threshold <- minD
for(j in 2:length(knot)){
  if((knot[j] - knot[(j-1)]) <= threshold){
    knot_temp[j-1] <- NA
  }
}
knot <- c(na.omit(knot_temp))

if(max(knot) == n){
  etau <- knot[2:(length(knot)-1)]
}
else{
  etau <- knot[-1]
}


num_step <- length(etau)

res_stepTP <- sapply(etau , FUN = function(x) ifelse(abs(x - tau)<=thresh, 1, 0))

step_TP <- sum(apply(res_stepTP, 1, FUN = function(x) ifelse(sum(x) > 0, 1, 0)))

step_FP <- max(0, length(etau)- step_TP)

Hausdorff_step <- Hausdorff(etau, tau, n)

res <- rbind(res, c(num_step, step_FP/num_step, step_TP/11, Hausdorff_step))

cat(m, '\n')

}

table(res[, 1])

apply(res, 2, mean)
