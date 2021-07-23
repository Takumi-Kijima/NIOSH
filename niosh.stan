data{
  int N; //sample size
  vector[N] time;//predictor
  vector[N] LC; //response
}

parameters{
  real Intercept;
  real slope;
  real<lower=0> sigma; // error sd
}



model{
    slope~normal(-0.7310924,0.2272379); //mu=(-1/.85-1/3.5)/2 #sd is determined so that P(-1/.85<slope<-1/3.5)=.95
    LC~normal( Intercept+slope*time,sigma);
}

