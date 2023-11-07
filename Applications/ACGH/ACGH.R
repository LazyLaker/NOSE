rm(list = ls())
set.seed(1234567)
library(nimble)
library(not)
library(foreach)
library(doParallel)
library(dplyr)
library(flexclust)
library(purrr)
library(MCMCglmm)
library(stepR)
library(wbs)
library(changepoint)
library(ecp)
library(DeCAFS)
library(ggplot2)
library(robseg) ### download from https://github.com/guillemr/robust-fpop
library(mosum)
source('NOSE_CPMean.R')

data("ACGH")


y <- ACGH$data[, 37]

n <- length(y)


res_NOSE <- NOSE_meanCP(y, L = 35, minD = 10,
                     nchains = 4, niter = 22000, nburnin = 2000, 
                     thin = 20)



res_not <- not(y, method = 'not', contrast = 'pcwsConstMean')

res_not <- features(res_not)

res_smuce <- smuceR(y)

res_PELT <- cpt.mean(y/mad(diff(y)/sqrt(2)), method = 'PELT', penalty = 'SIC', Q = 30)

res_wbs <- wbs(y)

res_mosum <- mosum(y, G = 40)

#res_ecp <- e.divisive(matrix(y, ncol = 1), min.size = 10)

res_fpop <- fpop_intern(y, pen.value = 1.345, test.stat = 'Huber')

plot(y, type = 'l')




dat <- data.frame(location = 1:length(y), number = y)

##### Original data

p <- ggplot(dat, aes(location, number)) + geom_point()
p <- p + geom_vline(xintercept  = c(1271, 1277), color = 'green')
p <- p + geom_vline(xintercept  = c(524, 583), color = 'blue', linetype = 'dashed')
p


#### NOSE
p <- ggplot(dat, aes(location, number)) + geom_point()

p <- p + geom_vline(xintercept  = res_NOSE$CP_detect$location, color = 'red')
p


#### SMUCE
p <- ggplot(dat, aes(location, number)) + geom_point()

p <- p + geom_vline(xintercept  = res_smuce$leftEnd[-1], color = 'red')
p


#### NOT
p <- ggplot(dat, aes(location, number)) + geom_point()
p <- p + geom_vline(xintercept  = res_not$cpt, color = 'red')
p <- p + geom_vline(xintercept  = c(1271, 1277), color = 'green')
p

#### FPOP
p <- ggplot(dat, aes(location, number)) + geom_point()
p <- p + geom_vline(xintercept  = sort(res_fpop$cpts, decreasing = T)[-1], color = 'red')
p <- p + geom_vline(xintercept = c(524, 582), color = 'blue', linetype = 'dashed')
p



####### remove outliers


win_width <- 10
resid <- numeric(n)

for(i in 1:n){

  resid[i] <- y[i] - median(y[which(abs((1:n)-i) <= win_width)])

}

outlier <- which(abs(resid) > 3*sd(resid))


y2 <- y
y2[outlier]  <- NA
y2 <- na.omit(y2)

res_not2<- features(not(y2))

res_fpop2 <- fpop_intern(y2, pen.value = 1.345, test.stat = 'Huber')



res_NOSE2 <- NOSE_meanCP(y2, L = 35, minD = 10,
                         nchains = 4, niter = 22000, nburnin = 2000,
                         thin = 20)


dat2 <- data.frame(number = y2, location = 1:length(y2))

#### Data_Outlier_removed
p <- ggplot(dat2, aes(location, number)) + geom_point()
p <- p + labs(title="Data without outliers")

p


##### NOSE Outlier removed
p <- ggplot(dat2, aes(location, number)) + geom_point()

p <- p + geom_vline(xintercept  = res_NOSE2$CP_detect$location, color = 'red')
p <- p + labs(title="NOSE")
p



##### NOT outlier removed
p <- ggplot(dat2, aes(location, number)) + geom_point()

p <- p + geom_vline(xintercept  = res_not2$cpt, color = 'red')
p <- p + labs(title="NOT")

p

#### R-FPOP outlier removed
p <- ggplot(dat2, aes(location, number)) + geom_point()

p <- p + geom_vline(xintercept = sort(res_fpop2$cpts, decreasing = T)[-1], color = 'red')
 p <- p + geom_vline(xintercept = c(518, 575), color = 'blue', linetype = 'dashed')
p <- p + labs(title="R-FPOP")

p

