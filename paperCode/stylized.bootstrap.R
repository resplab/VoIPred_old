library(MASS)
data_set <- birthwt

n <- dim(birthwt)[1]
z <- 0.2 #This is the risk threshold
n_sim <- 1000

set.seed(1)

model <- glm(low ~ age + lwt, family=binomial(link="logit"), data=data_set) #Our risk
#1:
pi <- predict(model, type="response") #Predicted risks
#2:
NB_model <- NB_all <- NB_mx <- rep(0,n_sim)
for(i in 1:n_sim) #Looping over many simulations
{
  #2.1
  bs_data_set <- data_set[sample(1:n, n, replace = T),]
  bs_model <- glm(low ~ age + lwt, family=binomial(link="logit"), data=bs_data_set)
  #2.2
  p <- predict(bs_model, newdata = data_set, type="response")
  #2.3
  NB_all[i] <- mean(p-(1-p)*z/(1-z))
  NB_model[i] <- mean((pi>z)*(p-(1-p)*z/(1-z))) #NB of using the model
  NB_max[i] <- mean((p>z)*(p-(1-p)*z/(1-z))) #NB of using the correct risks
}
#3
ENB_all<-mean(NB_all);ENB_model<-mean(NB_model);ENB_max<-mean(NB_max)
#4
EVPI<-ENB_max-max(0,ENB_model,ENB_all)
EVPI_r<-(ENB_max-max(0,ENB_all))/(ENB_model-max(0,ENB_all))

