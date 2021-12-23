library(R2OpenBUGS)

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

model <- glm(low ~ age + lwt, family=binomial(link="logit"), data=data_set)
pi <- predict(model,type="response")
betahats <- unname(coefficients(model))

bugs_input<-list(N=dim(birthwt)[1],low=birthwt$low,age=birthwt$age,lwt=birthwt$lwt, pi=pi)

n_sim <- 20000 #More simulation runs because of autocorrelation in MCMC

out<-bugs(data=bugs_input, inits=list(list(b0=betahats[1], b1=betahats[2], b2=betahats[2])), model.file=bugs_model, parameters.to.save = c('NB_all','NB_model','NB_max'), n.chains = 1, n.burnin = 5000, n.iter = n_sim+5000, bugs.seed=1)

ENB_all <- mean(out$sims.list$NB_all); ENB_model <- mean(out$sims.list$NB_model); ENB_max <- mean(out$sims.list$NB_max)

EVPI<-ENB_max-max(0,ENB_model,ENB_all)
EVPI_r<-(ENB_max-max(0,ENB_all))/(ENB_model-max(0,ENB_all))
