# ============================================================
# Simulation Study: Example 2
# ============================================================
# This script implements the Bayesian double machine-learning
# approach proposed in DiTraglia and Liu (2025).
#
# The Bayesian model is fitted using Stan, with posterior
# sampling performed using 6,000 iterations, including 1,000
# warm-up iterations.
# ============================================================


set.seed(2213)
library(rstan)
rstan_options(auto_write = TRUE)
library(MASS)
library(mvtnorm)
library(parallel)
options(mc.cores = 18)


BDML_ex2<-function(seed,n){
set.seed(seed)
p = 40
sig = array(0.05,c(p,p))
diag(sig) = 1
x = matrix(rmvnorm(n,mean=rep(0, nrow(sig)), sigma=sig),n,p)
d = rnorm(n, mean = 0.45*x[,1] + 0.9* x[, 2] - 0.4*x[,5]  ,sd = 1 )
y = rnorm(n, mean= d + 0.5*x[, 1] + x[,3] - 0.1*x[,4]-0.2*x[,7], sd=1)

stan_data <- list(N = n,P = p,X = x,D = d,Y = y)

start_time<- Sys.time()
fit <- stan(file = "BDML.stan",data = stan_data,iter = 6000,warmup = 1000,chains = 1,cores = 1)

ena_time<- Sys.time()
run_time <- ena_time - start_time

alpha_draws <- extract(fit, pars = "alpha")$alpha

# Posterior mean
alpha_mean <- mean(alpha_draws)

# 95% credible interval
alpha_ci <- quantile(alpha_draws, probs = c(0.025, 0.975))

ci<-0
if(alpha_ci[1]<1 && alpha_ci[2]>1) {ci<-1}
return(list(theta=alpha_mean,ci=ci,run_time=run_time))
}

BDML_ex2_40<-mclapply(sample(c(1:50000000),1000),function(x) BDML_ex2(x,40))
para_BDML_ex2_40<-unlist(mclapply(BDML_ex2_40, '[[', "theta"))
sqrt((mean(para_BDML_ex2_40)-1)^2 + var(para_BDML_ex2_40))
mean(para_BDML_ex2_100)-1
coverge_BDML_ex2_40<-sum(unlist(mclapply(BDML_ex2_40, '[[', "ci")))/1000


BDML_ex2_100<-mclapply(sample(c(1:50000000),1000),function(x) BDML_ex2(x,100))
mean(para_BDML_ex2_100)-1
sqrt((mean(para_BDML_ex2_100)-1)^2 + var(para_BDML_ex2_100))
coverge_BDML_ex2_100<-sum(unlist(mclapply(BDML_ex2_100, '[[', "ci")))/1000
