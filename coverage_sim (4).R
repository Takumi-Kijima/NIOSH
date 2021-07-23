jackknife <-function(data,sim_num) {
  jack_result<-rep(NA,sim_num)
  for(i in 1:sim_num)
  {
    index<-sample(1:nrow(data),1)
    x<-data$logTime[-index]
    y<-data$logLC50[-index]
    model<-lm(y~x) 
    jack_result[i]<--1/as.numeric(model$coefficients[2])
  }
  return(jack_result)
}

tenBerge <- function(formula, data, indices) {
  d <- data[indices,] # allows boot to select sample
  mod <- lm(formula, data=d)
  
  beta1 <-mod$coefficients[2]
  
  return(-1/beta1)
}




coverage <- data.frame()

for (size in seq(3, 10, 1)) {
  
  boot_sim <- data.frame()
  set.seed(123)  #5
  N <- 1000
  beta1 <- -0.9 #-1.15
  sigma <- 0.25
  
  for (i in  1:N) {  #N
    
    time=sample(c(10, 20, 40, 60, 240, 30, 120, 15, 480, 360, 5), size)  #runif(3,2,6) 
    logTime = log(time)
    logLC50 = 10+beta1*logTime+rnorm(size,0,sigma)
    
    data <- data.frame(logTime, logLC50)

    results <- boot(data=data, statistic=tenBerge,
                    R=1000, formula=logLC50 ~ logTime)

    # get 95% confidence interval
    ci <- boot.ci(results, type="all")
    boot_sim[i,1] <- ci$bca[4]
    boot_sim[i,2] <- ci$bca[5]
    
    boot_sim[i,3] <- ci$perc[4]
    boot_sim[i,4] <- ci$perc[5]
    
    results_j<-jackknife(data,1000)
    boot_sim[i,5] <- quantile(results_j,c(.025))  #5
    boot_sim[i,6] <- quantile(results_j,c(.975))  #6
  }
  
  good_sim <- boot_sim %>% rename("bl"="V1", "bu"="V2", "pl"="V3", "pu"="V4", "jl"="V5", "ju"="V6")
  
  n <- -1/beta1
  good_sim2 <- good_sim %>% 
    mutate(isBCA = ifelse(bl<n&n<bu, 1, 0),
           isPerc = ifelse(pl<n&n<pu, 1, 0),
           isJack = ifelse(jl<n&n<ju, 1, 0))
  
  coverage[size, 1] = size
  coverage[size, 2] = mean(good_sim2$isBCA)
  coverage[size, 3] = mean(good_sim2$isPerc)
  coverage[size, 4] = mean(good_sim2$isJack)
}

