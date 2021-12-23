data {
  int<lower=1> N;
  int<lower=0,upper=1> low[N];
  vector[N] age;
  vector[N] lwt;
  real<lower=0,upper=1> z;
}

parameters {
  real b0;
  real b1;
  real b2;
}

model {
  b0 ~ normal(0,1000);
  b1 ~ normal(0,1000);
  b2 ~ normal(0,1000);
  low ~ bernoulli_logit(b0+b1*age+b2*lwt);
}


bugs_model <- function()
{
  z <- 0.2
  for(i in 1:N)
  {
    p[i] <-  1/ (1+exp(-(b0 + b1*age[i] + b2*lwt[i])))
    low[i]~dbern(p[i])

    nb_all[i] <- p[i]-(1-p[i])*z/(1-z)
    nb_model[i] <- step(pi[i]-z)*(p[i]-(1-p[i])*z/(1-z))
    nb_max[i] <- step(p[i]-z)*(p[i]-(1-p[i])*z/(1-z))
  }

  NB_all <- mean(nb_all[])
  NB_model <- mean(nb_model[])
  NB_max <- mean(nb_max[])

  b0~dnorm(0,0.00001)
  b1~dnorm(0,0.00001)
  b2~dnorm(0,0.00001)
}
