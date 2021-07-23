library(rstan)
library(bayesplot)
library(tcltk)
#for faster compuation
rstan_options(auto_write=TRUE)
options(mc.cores=parallel::detectCores())


pb <- txtProgressBar(min = 1, max = length(models), style = 3) #show the progress
ci.bayes<-matrix(0,8,4)
colnames(ci.bayes)<-c("chemicals","n","Lower","Upper")
for (i in 1:length(models)) {  
  setTxtProgressBar(pb, i) 
  Sys.sleep(0.01)  
  data <- good_niosh %>% 
    filter(Chemical == (names(models)[i]))
  
  data_list<-list(N=nrow(data),LC=data$logLC50, time=data$logTime)
  mcmc_result<-stan(file="niosh.stan",data=data_list,seed=123)
  mcmc_result
  mcmc_sample<-rstan::extract(mcmc_result,permuted=FALSE)
  mcmc_combo(mcmc_sample,pars=c("Intercept","beta","sigma"))
  beta<-mcmc_sample[4001:8000]#extract betas
  n<--1/beta
  ci.bayes[i,3]<-as.numeric(round(quantile(n,.025),4))
  ci.bayes[i,4]<-as.numeric(round(quantile(n,.975),4))
  
}

for (i in 1:length(models)) {
  data <- good_niosh %>% 
    filter(Chemical == (names(models)[i]))
  mod<-summary(lm(data$logLC50~data$logTime))
  beta1<-as.numeric(mod$coefficients[2,1])
  ci.bayes[i,2] <- round(-1/beta1,4)
}

for (i in 1:length(models)) {
  ci.bayes[i,1] <- names(models)[i]
}

ci.bayes<-data.frame(ci.bayes)

ci.bayes
