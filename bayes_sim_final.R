library(rstan)
library(bayesplot)
library(tcltk)
rstan_options(auto_write=TRUE)
options(mc.cores=1)
max<-10
bayes_sim <- matrix(0,max-2,5)
colnames(bayes_sim)<-c("Size","EDI cov","EDI width","HDI cov","HDI width")
N <- 1000
beta1 <- -.9 
sigma <- 0.25
k<-0

pb <- txtProgressBar(min = 1, max = N*(size-2), style = 3)
for(size in 3:max)
{

  set.seed(123)
  bayes_sim[size-2,1]<-size
  for (i in  1:N) {  
    k<-k+1
    setTxtProgressBar(pb, k)
    Sys.sleep(0.01) 
    time=sample(c(10, 20, 40, 60, 240, 30, 120, 15, 480, 360, 5), size)  #runif(3,2,6) 
    logTime = log(time)
    logLC50 = 10+beta1*logTime+rnorm(size,0,sigma)
    data <- data.frame(logTime, logLC50)
    data_list<-list(N=nrow(data),LC=data$logLC50, time=data$logTime)
    mcmc_result<-stan(file="niosh.stan",data=data_list,seed=123,refresh=0)
    mcmc_sample<-rstan::extract(mcmc_result,permuted=FALSE)
  
    beta<-mcmc_sample[4001:8000]#extract betas
    n<- -1/beta  # ten berge
    #ETI coverage
    bayes_sim[size-2,2]<- bayes_sim[size-2,2]+ifelse((as.numeric(quantile(n,.025))< -1/beta1 & as.numeric(quantile(n,.975))> -1/beta1),1,0)
    #sum of ETI width
    bayes_sim[size-2,3]<-bayes_sim[size-2,3]+as.numeric(quantile(n,.975)-quantile(n,.025))
    #HDI coverage
    bayes_sim[size-2,4]<- bayes_sim[size-2,4]+ifelse((as.numeric(hdi(n)[1])< -1/beta1 & as.numeric(hdi(n)[2])> -1/beta1),1,0)
    #sum of ETI width
    bayes_sim[size-2,5]<-bayes_sim[size-2,5]+as.numeric(hdi(n)[2]-hdi(n)[1])
    }
  
  bayes_sim[size-2,rep(2:5)]<- bayes_sim[size-2,rep(2:5)]/N #get average width and coverage prob
}
