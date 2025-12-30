set.seed(2213)
library(matrixStats)
library(invgamma)
library(LaplacesDemon)
library(momentfit)
library(gtools)
library(mvtnorm)
library(randomForest)
library(nnet)
library(neuralnet)
library(parallel)
library(glmnet)
options(mc.cores = 23)


### load functions
source("Neyman_EL_Ex2.R")  
source("Neyman_ETEL_Ex2.R") 
source("Neyman_HD_Ex2.R")  

### prior function
loglik.prior<-function(theta){
  sum(dnorm(theta,mean=mutheta,sd=sigtheta,log=TRUE))
}


#Set the prior parameters
mualp<-1;sigalp<-1
mutheta<-1;sigtheta<-1
### Random forest
### n = 40
### ETEL 
high_ETEL_Neyman_ex2_40<-mclapply(sample(c(1:50000000),1000),function(x) high_Neyman_ETEL_ex2(x,40))
para_high_ETEL_Neyman_ex2_40<-unlist(mclapply(high_ETEL_Neyman_ex2_40, '[[', "theta"))
coverge_high_ETEL_Neyman_ex2_40<-sum(unlist(mclapply(high_ETEL_Neyman_ex2_40, '[[', "ci")))/1000
coverge_high_ETEL_Neyman_ex2_40
mean(para_high_ETEL_Neyman_ex2_40)-1
sqrt((mean(para_high_ETEL_Neyman_ex2_40)-1)^2+var(para_high_ETEL_Neyman_ex2_40))

### EL
high_EL_Neyman_ex2_40<-mclapply(sample(c(1:50000000),1000),function(x) high_Neyman_EL_ex2(x,40))
para_high_EL_Neyman_ex2_40<-unlist(mclapply(high_EL_Neyman_ex2_40, '[[', "theta"))
coverge_high_EL_Neyman_ex2_40<-sum(unlist(mclapply(high_EL_Neyman_ex2_40, '[[', "ci")))/1000
coverge_high_EL_Neyman_ex2_40
mean(para_high_EL_Neyman_ex2_40)-1
sqrt((mean(para_high_EL_Neyman_ex2_40)-1)^2+var(para_high_EL_Neyman_ex2_40))

### HD
high_HD_Neyman_ex2_40<-mclapply(sample(c(1:50000000),1000),function(x) high_Neyman_HD_ex2(x,40))
para_high_HD_Neyman_ex2_40<-unlist(mclapply(high_HD_Neyman_ex2_40, '[[', "theta"))
coverge_high_HD_Neyman_ex2_40<-sum(unlist(mclapply(high_HD_Neyman_ex2_40, '[[', "ci")))/1000
coverge_high_HD_Neyman_ex2_40
mean(para_high_HD_Neyman_ex2_40)-1
sqrt((mean(para_high_HD_Neyman_ex2_40)-1)^2+var(para_high_HD_Neyman_ex2_40))


### LASSO
### n = 40
### ETEL 
lasso_ETEL_Neyman_ex2_40<-mclapply(sample(c(1:50000000),1000),function(x) lasso_Neyman_ETEL_ex2(x,40))
para_lasso_ETEL_Neyman_ex2_40<-unlist(mclapply(lasso_ETEL_Neyman_ex2_40, '[[', "theta"))
coverge_lasso_ETEL_Neyman_ex2_40<-sum(unlist(mclapply(lasso_ETEL_Neyman_ex2_40, '[[', "ci")))/1000
coverge_lasso_ETEL_Neyman_ex2_40
mean(para_lasso_ETEL_Neyman_ex2_40)-1
sqrt((mean(para_lasso_ETEL_Neyman_ex2_40)-1)^2+var(para_lasso_ETEL_Neyman_ex2_40))

### EL
lasso_EL_Neyman_ex2_40<-mclapply(sample(c(1:50000000),1000),function(x) lasso_Neyman_EL_ex2(x,40))
para_lasso_EL_Neyman_ex2_40<-unlist(mclapply(lasso_EL_Neyman_ex2_40, '[[', "theta"))
coverge_lasso_EL_Neyman_ex2_40<-sum(unlist(mclapply(lasso_EL_Neyman_ex2_40, '[[', "ci")))/1000
coverge_lasso_EL_Neyman_ex2_40
mean(para_lasso_EL_Neyman_ex2_40)-1
sqrt((mean(para_lasso_EL_Neyman_ex2_40)-1)^2+var(para_lasso_EL_Neyman_ex2_40))


### HD
lasso_HD_Neyman_ex2_40<-mclapply(sample(c(1:50000000),1000),function(x) lasso_Neyman_HD_ex2(x,40))
para_lasso_HD_Neyman_ex2_40<-unlist(mclapply(lasso_HD_Neyman_ex2_40, '[[', "theta"))
coverge_lasso_HD_Neyman_ex2_40<-sum(unlist(mclapply(lasso_HD_Neyman_ex2_40, '[[', "ci")))/1000
coverge_lasso_HD_Neyman_ex2_40
mean(para_lasso_HD_Neyman_ex2_40)-1
sqrt((mean(para_lasso_HD_Neyman_ex2_40)-1)^2+var(para_lasso_HD_Neyman_ex2_40))



### Neuro Network
### n = 40
### ETEL 
nnet_ETEL_Neyman_ex2_40<-mclapply(sample(c(1:50000000),1000),function(x) nnet_Neyman_ETEL_ex2(x,40))
para_nnet_ETEL_Neyman_ex2_40<-unlist(mclapply(nnet_ETEL_Neyman_ex2_40, '[[', "theta"))
coverge_nnet_ETEL_Neyman_ex2_40<-sum(unlist(mclapply(nnet_ETEL_Neyman_ex2_40, '[[', "ci")))/1000
coverge_nnet_ETEL_Neyman_ex2_40
mean(para_nnet_ETEL_Neyman_ex2_40)-1
sqrt((mean(para_nnet_ETEL_Neyman_ex2_40)-1)^2+var(para_nnet_ETEL_Neyman_ex2_40))

### EL
nnet_EL_Neyman_ex2_40<-mclapply(sample(c(1:50000000),1000),function(x) nnet_Neyman_EL_ex2(x,40))
para_nnet_EL_Neyman_ex2_40<-unlist(mclapply(nnet_EL_Neyman_ex2_40, '[[', "theta"))
coverge_nnet_EL_Neyman_ex2_40<-sum(unlist(mclapply(nnet_EL_Neyman_ex2_40, '[[', "ci")))/1000
coverge_nnet_EL_Neyman_ex2_40
mean(para_nnet_EL_Neyman_ex2_40)-1
sqrt((mean(para_nnet_EL_Neyman_ex2_40)-1)^2+var(para_nnet_EL_Neyman_ex2_40))

### HD
nnet_HD_Neyman_ex2_40<-mclapply(sample(c(1:50000000),1000),function(x) nnet_Neyman_HD_ex2(x,40))
para_nnet_HD_Neyman_ex2_40<-unlist(mclapply(nnet_HD_Neyman_ex2_40, '[[', "theta"))
coverge_nnet_HD_Neyman_ex2_40<-sum(unlist(mclapply(nnet_HD_Neyman_ex2_40, '[[', "ci")))/1000
coverge_nnet_HD_Neyman_ex2_40
mean(para_nnet_HD_Neyman_ex2_40)-1
sqrt((mean(para_nnet_HD_Neyman_ex2_40)-1)^2+var(para_nnet_HD_Neyman_ex2_40))


