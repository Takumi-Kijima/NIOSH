#Estimating ten Berge by jackknife method


jack_ci<-matrix(NA,8,3)#store results
jackknife <-function(data,sim_num)
{
  jack_result<-rep(NA,sim_num)
  for(i in 1:sim_num)
  {
    index<-sample(1:nrow(data),1) # pick one observation randomly
    x<-data$logTime[-index] #drop an observation
    y<-data$logLC50[-index]
    model<-lm(y~x) 
    jack_result[i]<--1/as.numeric(model$coefficients[2])
  }
  return(jack_result)
}

for (i in 1:length(models)) {  
  beta<-rep(NA,1000)
  data <- good_niosh_2 %>% 
    filter(Chemical == names(models)[i])
  
  results<-jackknife(data,1000)
  jack_ci[i,1]<-names(models)[i]
  jack_ci[i,2] <- quantile(results,.025) #lower bound of 95% CI. 
  jack_ci[i,3]  <- quantile(results,.975) #upper bound of 95% CI
}