library(MASS)
library(mvtnorm)
data_set <- birthwt

n <- dim(birthwt)[1]
n_sim <- 10000

z <- 0.2
set.seed(1234)
model <- glm(low ~ age + lwt, family=binomial(link="logit"), data=data_set)
pi <- predict(model, type="response") #Predicted risks

mu <- coefficients(model)
sigma <- vcov(model)

NB_all <- NB_model <- NB_max <- rep(0, n_sim)

betas <- mvtnorm::rmvnorm(n_sim,mu,sigma)

for(i in 1:n_sim)
{
  p <- 1 / (1+exp(-(betas[i,1] + betas[i,2]*data_set[,'age'] + betas[i,3]*data_set[,'lwt'])))

  NB_all[i] <- mean(p-(1-p)*z/(1-z)) #NB of treating all
  NB_model[i] <- mean((pi>z)*(p-(1-p)*z/(1-z))) #NB of using the model
  NB_max[i] <- mean((p>z)*(p-(1-p)*z/(1-z))) #NB of using the correct risks
}

ENB_all <- mean(NB_all); ENB_model <- mean(NB_model); ENB_max <- mean(NB_max)

EVPI<-ENB_max-max(0,ENB_model,ENB_all)
EVPI_r<-(ENB_max-max(0,ENB_all))/(ENB_model-max(0,ENB_all))

readr::write_rds(c(ENB_all,ENB_model,ENB_max,EVPI,EVPI_r),here("paperCode","supp_results","likelihood.rds"))

cat("Without any model, the expected NB of the best decisoin is ", ENB_base <- max(0,ENB_all))
cat("The expected incremental NB of the proposed model is ", INB_current <- max(0,ENB_all,ENB_model) - ENB_base)
cat("The expected incremental NB of the correct model is ", INB_perfect <- ENB_max - ENB_base)
