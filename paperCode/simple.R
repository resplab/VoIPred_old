library(MASS)

#birthw is an exemlary dataset containing the low birthweight status of newborns and some covairates
data_set <- birthwt
n <- dim(birthwt)[1]

z <- 0.2 #This is the risk threshold

model <- glm(low ~ age + lwt, family=binomial(link="logit"), data=data_set) #Our risk prediciton model

pi <- predict(model, type="response") #Predicted risks

#We hold simulation results in memory for exploration (not necessary)
out <- data.frame("NB_all"=double(),"NB_model"=double(),"NB_max"=double())

#Looping over 1000 simulations
for(i in 1:1000)
{
  bs_data_set <- data_set[sample(1:n, n, replace = T),] #Create a bootstrap sample and fit the model again
  bs_model <- glm(low ~ age + lwt, family=binomial(link="logit"), data=bs_data_set)
  #Predict risks from this model are applied to the original data
  p <- predict(bs_model, newdata = data_set, type="response")

  #Bayesiuan NB calculations. p are taken as random draws from the distribution of correct risks
  NB_all <- mean(p-(1-p)*z/(1-z)) #NB of treating all
  NB_model <- mean((pi>z)*(p-(1-p)*z/(1-z))) #NB of using the model
  NB_max <- mean((p>z)*(p-(1-p)*z/(1-z))) #NB of using the correct risks

  out[i,] <- c(NB_all,NB_model,NB_max)
}

EVPI <- mean(out$NB_max) - max(0,mean(out$NB_all),mean(out$NB_model))
cat("EVPI is ",EVPI)

#Some additional useful information
cat("Without any model, the expected NB of the best decisoin is ", NB_base <- max(0,mean(out$NB_all)))
cat("The expected incremental NB of the proposed model is ", INB_current <- max(0,mean(out$NB_all),mean(out$NB_model)) - NB_base)
cat("The expected incremental NB of the correct model is ", INB_perfect <- mean(out$NB_max) - NB_base)
cat("Relative EVPI is ", INB_perfect / INB_current)
