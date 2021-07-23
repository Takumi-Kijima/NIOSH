library(rstan)
library(bayesplot)
library(tcltk)
rstan_options(auto_write=TRUE)
options(mc.cores=parallel::detectCores())
max<-7
bayes_sim <- matrix(0,max-2,2)
N <- 100
beta1 <- -.9 
sigma <- 0.25
size<-3
set.seed(123)
pb <- txtProgressBar(min = 1, max = N, style = 3)
#setwd("")
for (i in  1:N) {  
  setTxtProgressBar(pb, i)
  Sys.sleep(0.01) 
  time=sample(c(10, 20, 40, 60, 240, 30, 120, 15, 480, 360, 5), size)  #runif(3,2,6) 
  logTime = log(time)
  logLC50 = 10+beta1*logTime+rnorm(size,0,sigma)
  data <- data.frame(logTime, logLC50)
  data_list<-list(N=nrow(data),LC=data$logLC50, time=data$logTime)
  mcmc_result<-stan(file="niosh.stan",data=data_list,seed=123)
  mcmc_sample<-rstan::extract(mcmc_result,permuted=FALSE)
  mcmc_combo(mcmc_sample,pars=c("Intercept","slope","sigma"))
  beta<-mcmc_sample[4001:8000]#extract betas
  n<- -1/beta  # ten berge
  bayes_sim[size-2,2]<- bayes_sim[size-2,2]+ifelse((as.numeric(quantile(n,.025))< -1/beta1 & as.numeric(quantile(n,.975))> -1/beta1),1,0)
}
