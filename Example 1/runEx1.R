set.seed(213)
library(matrixStats)
library(invgamma)
library(LaplacesDemon)
library(momentfit)
library(gtools)
library(mvtnorm)
library(randomForest)
library(parallel)
options(mc.cores = 23)

### prior function
loglik.prior<-function(theta){
  sum(dnorm(theta,mean=mutheta,sd=sigtheta,log=TRUE))
}


#Set the prior parameters
mualp<-0;sigalp<-1000 
mutheta<-0;sigtheta<-1000
### Random forest
### n = 50
### ETEL 
high_ETEL_Neyman_50<-mclapply(sample(c(1:50000000),1000),function(x) high_Neyman_ETEL(x,50))
para_high_ETEL_Neyman_50<-unlist(mclapply(high_ETEL_Neyman_50, '[[', "theta"))
coverge_high_ETEL_Neyman_50<-sum(unlist(mclapply(high_ETEL_Neyman_50, '[[', "ci")))/1000
### EL
high_EL_Neyman_50<-mclapply(sample(c(1:50000000),1000),function(x) high_Neyman_EL(x,50))
para_high_EL_Neyman_50<-unlist(mclapply(high_EL_Neyman_50, '[[', "theta"))
coverge_high_EL_Neyman_50<-sum(unlist(mclapply(high_EL_Neyman_50, '[[', "ci")))/1000
### HD
high_HD_Neyman_50<-mclapply(sample(c(1:50000000),1000),function(x) high_Neyman_HD(x,50))
para_high_HD_Neyman_50<-unlist(mclapply(high_HD_Neyman_50, '[[', "theta"))
coverge_high_HD_Neyman_50<-sum(unlist(mclapply(high_HD_Neyman_50, '[[', "ci")))/1000


### n = 200
### ETEL 
high_ETEL_Neyman_200<-mclapply(sample(c(1:20000000),1000),function(x) high_Neyman_ETEL(x,200))
para_high_ETEL_Neyman_200<-unlist(mclapply(high_ETEL_Neyman_200, '[[', "theta"))
coverge_high_ETEL_Neyman_200<-sum(unlist(mclapply(high_ETEL_Neyman_200, '[[', "ci")))/1000
### EL
high_EL_Neyman_200<-mclapply(sample(c(1:20000000),1000),function(x) high_Neyman_EL(x,200))
para_high_EL_Neyman_200<-unlist(mclapply(high_EL_Neyman_200, '[[', "theta"))
coverge_high_EL_Neyman_200<-sum(unlist(mclapply(high_EL_Neyman_200, '[[', "ci")))/1000
### HD
high_HD_Neyman_200<-mclapply(sample(c(1:20000000),1000),function(x) high_Neyman_HD(x,200))
para_high_HD_Neyman_200<-unlist(mclapply(high_HD_Neyman_200, '[[', "theta"))
coverge_high_HD_Neyman_200<-sum(unlist(mclapply(high_HD_Neyman_200, '[[', "ci")))/1000



### LASSO
### n = 50
### ETEL 
lasso_ETEL_Neyman_50<-mclapply(sample(c(1:50000000),1000),function(x) high_Neyman_ETEL(x,50))
para_lasso_ETEL_Neyman_50<-unlist(mclapply(lasso_ETEL_Neyman_50, '[[', "theta"))
coverge_lasso_ETEL_Neyman_50<-sum(unlist(mclapply(lasso_ETEL_Neyman_50, '[[', "ci")))/1000
### EL
lasso_EL_Neyman_50<-mclapply(sample(c(1:50000000),1000),function(x) high_Neyman_EL(x,50))
para_lasso_EL_Neyman_50<-unlist(mclapply(lasso_EL_Neyman_50, '[[', "theta"))
coverge_lasso_EL_Neyman_50<-sum(unlist(mclapply(lasso_EL_Neyman_50, '[[', "ci")))/1000
### HD
lasso_HD_Neyman_50<-mclapply(sample(c(1:50000000),1000),function(x) lasso_Neyman_HD(x,50))
para_lasso_HD_Neyman_50<-unlist(mclapply(lasso_HD_Neyman_50, '[[', "theta"))
coverge_lasso_HD_Neyman_50<-sum(unlist(mclapply(lasso_HD_Neyman_50, '[[', "ci")))/1000


### n = 200
### ETEL 
lasso_ETEL_Neyman_200<-mclapply(sample(c(1:50000000),1000),function(x) high_Neyman_ETEL(x,200))
para_lasso_ETEL_Neyman_200<-unlist(mclapply(lasso_ETEL_Neyman_200, '[[', "theta"))
coverge_lasso_ETEL_Neyman_200<-sum(unlist(mclapply(lasso_ETEL_Neyman_200, '[[', "ci")))/1000
### EL
lasso_EL_Neyman_200<-mclapply(sample(c(1:20000000),1000),function(x) high_Neyman_EL(x,200))
para_lasso_EL_Neyman_200<-unlist(mclapply(lasso_EL_Neyman_200, '[[', "theta"))
coverge_lasso_EL_Neyman_200<-sum(unlist(mclapply(lasso_EL_Neyman_200, '[[', "ci")))/1000
### HD
lasso_HD_Neyman_200<-mclapply(sample(c(1:50000000),1000),function(x) lasso_Neyman_HD(x,200))
para_lasso_HD_Neyman_200<-unlist(mclapply(lasso_HD_Neyman_200, '[[', "theta"))
coverge_lasso_HD_Neyman_200<-sum(unlist(mclapply(lasso_HD_Neyman_200, '[[', "ci")))/1000


### Neuro Network
### n = 50
### ETEL 
nnet_ETEL_Neyman_50<-mclapply(sample(c(1:50000000),1000),function(x) nnet_Neyman_ETEL(x,50))
para_nnet_ETEL_Neyman_50<-unlist(mclapply(nnet_ETEL_Neyman_50, '[[', "theta"))
coverge_nnetETEL_Neyman_50<-sum(unlist(mclapply(nnet_ETEL_Neyman_50, '[[', "ci")))/1000
### EL
nnet_EL_Neyman_50<-mclapply(sample(c(1:50000000),1000),function(x) nnet_Neyman_EL(x,50))
para_nnet_EL_Neyman_50<-unlist(mclapply(nnet_EL_Neyman_50, '[[', "theta"))
coverge_nnet_EL_Neyman_50<-sum(unlist(mclapply(nnet_EL_Neyman_50, '[[', "ci")))/1000
### HD
nnet_HD_Neyman_50<-mclapply(sample(c(1:50000000),1000),function(x) nnet_Neyman_HD(x,50))
para_nnet_HD_Neyman_50<-unlist(mclapply(nnet_HD_Neyman_50, '[[', "theta"))
coverge_nnet_HD_Neyman_50<-sum(unlist(mclapply(nnet_HD_Neyman_50, '[[', "ci")))/1000


### n = 200
### ETEL 
nnet_ETEL_Neyman_200<-mclapply(sample(c(1:50000000),1000),function(x) nnet_Neyman_ETEL(x,200))
para_nnet_ETEL_Neyman_200<-unlist(mclapply(nnet_ETEL_Neyman_200, '[[', "theta"))
coverge_nnet_ETEL_Neyman_200<-sum(unlist(mclapply(nnet_ETEL_Neyman_200, '[[', "ci")))/1000
### EL
nnet_EL_Neyman_200<-mclapply(sample(c(1:50000000),1000),function(x) nnet_Neyman_EL(x,200))
para_nnet_EL_Neyman_200<-unlist(mclapply(nnet_EL_Neyman_200, '[[', "theta"))
coverge_nnet_EL_Neyman_200<-sum(unlist(mclapply(nnet_EL_Neyman_200, '[[', "ci")))/1000
### HD
nnet_HD_Neyman_200<-mclapply(sample(c(1:50000000),1000),function(x) nnet_Neyman_HD(x,200))
para_nnet_HD_Neyman_200<-unlist(mclapply(nnet_HD_Neyman_200, '[[', "theta"))
coverge_nnet_HD_Neyman_200<-sum(unlist(mclapply(nnet_HD_Neyman_200, '[[', "ci")))/1000


