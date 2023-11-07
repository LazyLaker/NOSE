rm(list = ls())
set.seed(1234567)
library(nimble)
library(not)
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
library(gfpop)
library(MCMCpack)
source('NOSE_CPScale.R')
load('agri.RData')

y <- as.numeric(agri)
n <- length(y)
res_NOSE <- NOSE_ScaleCP(y, minD = 10, niter = 22000, nburnin = 2000, thin = 20)


res_not <- res_not <- not(y, method = 'not', contrast = 'pcwsConstMeanVar')

res_not <- features(res_not)

res_pelt <- cpt.var(y, penalty = 'SIC', method = 'PELT')

 res_ecp <- e.divisive(matrix(y, ncol = 1), min.size = 10)


res_smuce <- smuceR(y-mean(y), family = 'gaussvar')


myGraph <- graph(type = "std", penalty = 2 * log(n), K = 10)
gfpop(data = y, mygraph = myGraph, type = "variance")



###### plot 5x2.5


dat <- data.frame(Days= 1:length(y), Returns = y)

#### NOSE result
p <- ggplot(dat, aes(Days, Returns)) + geom_line()

p <- p + geom_rect(aes(xmin = 0, xmax = 139, ymin = -5, ymax = 5), fill = NA, color = "blue")

p <- p + geom_rect(aes(xmin = 206, xmax = 425, ymin = -10, ymax = 10), 
                   fill = NA, color = "orange",lty = 2)

p <- p + geom_vline(xintercept  = res_NOSE$CP_detect$location, color = 'red')
p <- p + ggtitle('NOSE')

p


##### SMUCE Result

p <- ggplot(dat, aes(Days, Returns)) + geom_line() 

p <- p + geom_vline(xintercept  = res_smuce$leftEnd[-1], color = 'red') 

p <- p + ggtitle('SMUCE')
p


##### PELT Result
p <- ggplot(dat, aes(Days, Returns)) + geom_line()

p <- p + geom_vline(xintercept  = res_pelt@cpts[-which.max(res_pelt@cpts)], color = 'red')

p <- p + ggtitle('PELT')

p

##### NOT Result
p <- ggplot(dat, aes(Days, Returns)) + geom_line()

p <- p + geom_vline(xintercept  = res_not$cpt, color = 'red')

p <- p + ggtitle('NOT')

p


### original data
p <- ggplot(dat, aes(Days, Returns)) + geom_line()

p <- p + geom_rect(aes(xmin = 0, xmax = 137, ymin = -5, ymax = 5), fill = NA, color = "blue")

p <- p + geom_rect(aes(xmin = 200, xmax = 420, ymin = -12, ymax = 10), 
                   fill = NA, color = "orange",lty = 2)

p



#### centered data
dat2 <- data.frame(Days= 1:length(y), CentR = abs(y-mean(y)))

tau <- c(0, res_NOSE$CP_detect$location, length(y))

seg_len <- diff(tau)


sigma <- c()

for(j in 1:8){
  sigma <- c(sigma, rep(sd(y[(tau[j]+1) : tau[j+1]]), seg_len[j]))
}


#### centered DRAIP plot
p <- ggplot(dat2, aes(Days, CentR)) + geom_line()

p <- p + geom_rect(aes(xmin = 0, xmax = 139, ymin = 0, ymax = 5), fill = NA, color = "blue")

p <- p + geom_rect(aes(xmin = 206, xmax = 425, ymin = 0, ymax = 10), 
                   fill = NA, color = "orange")

p <- p + geom_line(aes(Days, sigma), color = 'red', size = 1)

p <- p + geom_vline(xintercept  = res_NOSE$CP_detect$location, color = 'red', linetype = 'dashed')

p <- p + labs(y = 'Absolute Centered Return')

p




####### plots of abnormal interval

df <- data.frame(y[205:344])
colnames(df) <- 'y'
ggplot(df, aes(sample=y)) +
  stat_qq(size = 2) + 
  stat_qq_line()


ggplot(df, aes(x=y)) + 
  geom_density(color="darkblue", fill="lightblue")











