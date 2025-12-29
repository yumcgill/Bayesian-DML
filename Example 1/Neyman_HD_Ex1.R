### Algorithm 2 with K=2, Bayesian HD
Neyman_HD<-function(Y_res,D_res){
  
  ## create the design matrix
  OR_pred <- model.matrix( ~ D_res - 1)
  ### calculate log likelihood for prior 
  
  ### estimating equations
  gf<-function(theta,deX){
    gf1<-apply(t(deX),1, function(t) t*((Y_res-deX%*%theta)))
    return(gf1)
  }
  
  ##initial theta
  old.theta<-lm.fit(OR_pred,Y_res)$coef
  emp_mod<-momentModel(gf, OR_pred, old.theta, grad=NULL, vcov="iid")
  
  ## solve the empirical likelihood to update theta
  old.theta<-solveGel(update(emp_mod, gelType="HD"), theta0=old.theta,tControl=list(method="Brent", lower = -10, upper = 10))$theta
  
  ### MCMC initialization 
  ###starting values
  
  nburn <- 1000
  nits<- 5000 + nburn
  theta.samp<-matrix(0,nrow=nits,ncol=ncol(OR_pred))
  theta.samp[1,]<-old.theta
  
  prop.sig<-vcov(lm(Y_res~OR_pred-1))
  
  fit <- gelFit(emp_mod,gelType="HD",theta0=old.theta,tControl=list(method="Brent", lower = -10, upper = 10))
  
  emp_log_lik_old<-sum(log(getImpProb(fit)$pt))
  
  log_prior_old<-loglik.prior(old.theta)
  
  for (i in 2:nits){
    
    ##symmetric proposal for theta
    para_p<-as.vector(rnorm(1,mean=old.theta,sd=sqrt(prop.sig)))
    ## log MH prob
    log_p<-dnorm((old.theta),mean=para_p,sd=sqrt(prop.sig),log=TRUE)-
      dnorm(para_p,mean=old.theta,sd=sqrt(prop.sig),log=TRUE)
    ## log prior density
    log_prior_new<-loglik.prior(para_p)  
    
    new.lam<-evalGel(emp_mod,para_p,gelType="HD",tControl=list(method="Brent", lower = -10, upper = 10))@lambda
    ## log emp likelihood value
    emp_log_lik_new<-sum(log(getImpProb(evalGel(emp_mod,para_p,new.lam,gelType="HD",tControl=list(method="Brent", lower = -10, upper = 10)))$pt))
    
    jmp<-min(log_p+emp_log_lik_new+log_prior_new-emp_log_lik_old-log_prior_old,0)
    
    
    if(log(runif(1)) <= jmp){ ### accept the move
      emp_log_lik_old<-emp_log_lik_new
      log_prior_old<-log_prior_new
      theta.samp[i,]<-old.theta<-para_p
      
    } else{## reject the move
      
      theta.samp[i,]<-old.theta
      
    }
    
  }
  
  post_ate<-theta.samp[(nburn+1):nits,]
  #post_ate<-theta.samp_fin
  return(list(post_bHD_pscov_ate = post_ate))
} 

### Simulation EXAMPLE 1 Random Forest 
high_Neyman_HD<-function(seed,n){
  p = 500
  sig = array(0.3,c(p,p))
  diag(sig) = 1
  expit = function(x) {1/(1+exp(-x))}
  x = matrix(rmvnorm(n,mean=rep(0, nrow(sig)), sigma=sig),n,p)
  d = rbinom(n, 1, p=expit(0.3*x[,1] + 0.2*x[,2] - 0.4*x[,5]))
  y = rnorm(n, mean=d + 0.5*x[,1] + x[,3] - 0.1*x[,4]-0.2*x[,7], sd=1)
  start_time<- Sys.time()
  samp <- sample(1:n, floor(n / 2))
  osamp <- setdiff(1:n, samp)
  rd<- rep(0,n)
  ry <- rep(0,n)
  fy1 <- randomForest(x[samp, ], y[samp], ntree = 500)
  yhat1 <- predict(fy1, x[osamp, ])
  ry[osamp] <- y[osamp] - yhat1
  
  fd1 <- randomForest(x[samp, ], d[samp], ntree = 500)
  dhat1 <- predict(fd1, x[osamp, ])
  rd[osamp] <- d[osamp] - dhat1
  
  fy2 <- randomForest(x[osamp, ], y[osamp], ntree = 500)
  yhat2 <- predict(fy2, x[samp, ])
  ry[samp] <- y[samp] - yhat2
  
  fd2 <- randomForest(x[osamp, ], d[osamp], ntree = 500)
  dhat2 <- predict(fd2, x[samp, ])
  rd[samp] <- d[samp] - dhat2
  
  res<-Neyman_HD(ry,rd)
  
  ena_time<- Sys.time()
  run_time <- ena_time - start_time
  
  ci<-0
  postmean_theta1 <- mean(res$post_bHD_pscov_ate)
  if(quantile(res$post_bHD_pscov_ate,prob=c(0.025))<1 && quantile(res$post_bHD_pscov_ate,prob=c(0.975))>1) {ci<-1}
  return(list(theta=postmean_theta1,ci=ci,run_time=run_time))
}

### Simulation EXAMPLE 1 LASSO

lasso_Neyman_HD<-function(seed,n){
  p = 500
  sig = array(0.3,c(p,p))
  diag(sig) = 1
  expit = function(x) {1/(1+exp(-x))}
  x = matrix(rmvnorm(n,mean=rep(0, nrow(sig)), sigma=sig),n,p)
  d = rbinom(n, 1, p=expit(0.3*x[,1] + 0.2*x[,2] - 0.4*x[,5]))
  y = rnorm(n, mean=d + 0.5*x[,1] + x[,3] - 0.1*x[,4]-0.2*x[,7], sd=1)
  start_time<- Sys.time()
  samp <- sample(1:n, floor(n / 2))
  osamp <- setdiff(1:n, samp)
  rd<- rep(0,n)
  ry <- rep(0,n)

  fy1 <- cv.glmnet(x[samp, ], y[samp], alpha = 1)
  yhat1 <- predict(fy1, x[osamp, ], s=fy1$lambda.min)
  ry[osamp] <- y[osamp] - yhat1
  
  fd1 <- cv.glmnet(x[samp, ], d[samp], family="binomial",alpha = 1)
  dhat1 <- predict(fd1, x[osamp, ],s=fy1$lambda.min,type="response")
  rd[osamp] <- d[osamp] - dhat1
  
  fy2 <- cv.glmnet(x[osamp, ], y[osamp], alpha = 1)
  yhat2 <- predict(fy2, x[samp, ], s=fy2$lambda.min)
  ry[samp] <- y[samp] - yhat2
  
  fd2 <- cv.glmnet(x[osamp, ], d[osamp],  family="binomial",alpha = 1)
  dhat2 <- predict(fd2, x[samp, ],s=fy1$lambda.min,type="response")
  rd[samp] <- d[samp] - dhat2
  
  res<-Neyman_HD(ry,rd)
  
  ena_time<- Sys.time()
  run_time <- ena_time - start_time
  
  ci<-0
  postmean_theta1 <- mean(res$post_bHD_pscov_ate)
  if(quantile(res$post_bHD_pscov_ate,prob=c(0.025))<1 && quantile(res$post_bHD_pscov_ate,prob=c(0.975))>1) {ci<-1}
  return(list(theta=postmean_theta1,ci=ci,run_time=run_time))
}


### Simulation EXAMPLE 1 neural network
nnet_Neyman_HD<-function(seed,n){
  p = 500
  sig = array(0.3,c(p,p))
  diag(sig) = 1
  expit = function(x) {1/(1+exp(-x))}
  x = matrix(rmvnorm(n,mean=rep(0, nrow(sig)), sigma=sig),n,p)
  d = rbinom(n, 1, p=expit(0.3*x[,1] + 0.2*x[,2] - 0.4*x[,5]))
  y = rnorm(n, mean=d + 0.5*x[,1] + x[,3] - 0.1*x[,4]-0.2*x[,7], sd=1)
  
  maxs <- apply(x, 2, max) 
  mins <- apply(x, 2, min)
  
  x.s <- as.data.frame(scale(x, center = mins, scale = maxs - mins))
  y.s <- scale(y, center = min(y), scale = max(y) - min(y))
 
  start_time<- Sys.time()
  samp <- sample(1:n, floor(n / 2))
  osamp <- setdiff(1:n, samp)
  rd<- rep(0,n)
  ry <- rep(0,n)
  fy1 <- nnet(x.s[samp,],y.s[samp], size=2,  maxit=1000, decay=0.01, MaxNWts=10000,  trace=FALSE,linout = TRUE)
  yhat1 <- predict(fy1, x.s[osamp, ])*(max(y) - min(y))+min(y)
  ry[osamp] <- y[osamp] - yhat1
  
  fd1 <- nnet(x.s[samp,],d[samp],size=2,  maxit=1000, decay=0.01, MaxNWts=10000,  trace=FALSE,linout = FALSE)
  dhat1 <- predict(fd1, x.s[osamp, ])
  rd[osamp] <- d[osamp] - dhat1
  
  fy2 <- nnet(x.s[osamp,],y.s[osamp], size=2,  maxit=1000, decay=0.01, MaxNWts=10000,  trace=FALSE,linout = TRUE)
  yhat2 <- predict(fy2, x.s[samp, ])*(max(y) - min(y))+min(y)
  ry[samp] <- y[samp] - yhat2
  
  fd2 <-  nnet(x.s[osamp,],d[osamp],size=2,  maxit=1000, decay=0.01, MaxNWts=10000,  trace=FALSE,linout = FALSE)
  dhat2 <- predict(fd2, x.s[samp, ])
  rd[samp] <- d[samp] - dhat2
  
  res<-Neyman_HD(ry,rd)
  
  ena_time<- Sys.time()
  run_time <- ena_time - start_time
  
  ci<-0
  postmean_theta1 <- mean(res$post_bHD_pscov_ate)
  if(quantile(res$post_bHD_pscov_ate,prob=c(0.025))<1 && quantile(res$post_bHD_pscov_ate,prob=c(0.975))>1) {ci<-1}
  return(list(theta=postmean_theta1,ci=ci,run_time=run_time))
}


