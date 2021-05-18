#' @import mvtnorm







# bootstrap: 0=no (parametric); 1=ordinary; 2=Bayesian
#' @export
evpp.glm<-function(reg_obj, n_sim=1000, bootstrap=0, lambdas=(1:99)/100)
{
  sample_size <- dim(reg_obj$data)[1]
  mu <- coefficients(reg_obj)
  sigma <- vcov(reg_obj)

  pi<-predict(reg_obj,type="response")

  NB_test <-  NB_all <-  NB_max  <- rep(0, length(lambdas))

  for(i in 1:n_sim)
  {
    if(bootstrap>0)
    {
      ws <- bootstrap(sample_size, bootstrap-1)
      bs_data <- cbind(reg_obj$data, ws=ws)
      bs_reg <- glm(data=bs_data, formula = reg_obj$formula, family=reg_obj$family, weights = ws)
      p <- as.vector(predict(bs_reg, newdata=reg_obj$data, type="response"))
    }
    else
    {
      new_betas <- rmvnorm(1, mean=mu, sigma=sigma)
      p <- as.vector(1/(1+exp(-(new_betas %*% t(model.matrix(reg_obj))))))
    }

    if(i==1) plot(pi,p)

    for(j in 1:length(lambdas))
    {
      NB_test[j] <- NB_test[j] + mean((p - (1 - p) * lambdas[j] / (1 - lambdas[j])) * (pi > lambdas[j]))

      NB_all[j] <- NB_all[j] + mean((p - (1 - p) * lambdas[j] / (1 - lambdas[j])) * 1)

      NB_max[j] <- NB_max[j] + mean((p - (1 - p) * lambdas[j] / (1 - lambdas[j])) * (p > lambdas[j]))
    }
    #cat('.')
  }

  EVPP <- (NB_max-pmax(0,NB_test,NB_all))/n_sim

  res <-cbind(lambda=lambdas, EVPP=EVPP, NB_all=NB_all/n_sim, NB_test=NB_test/n_sim, NB_max=NB_max/n_sim)

  return(res)
}





#' @export
#' @param reg_obj: any object that you can aplpy predict with new data to get predictions
#' @param x: The model matrix of predictors
#' @param y: The vector of responses
evpp.glmnet <- function(reg_obj, x, y, n_sim=1000, lambdas=(1:99)/100, Bayesian_bootstrap=F)
{
  sample_size <- dim(x)[1]

  pi<-predict(reg_obj, type="response", newx=x)

  NB_model <- rep(0, length(lambdas))
  NB_all <- NB_model
  NB_max <- NB_model

  for(i in 1:n_sim)
  {
    cat(".")
    weights <- bootstrap(sample_size, Bayesian = Bayesian_bootstrap)
    tmp <- cv.glmnet(x=x, y=y, family="binomial",weights = as.vector(weights))
    bs_reg <- glmnet(x=x, y=y, family="binomial", lambda=tmp$lambda.min, weights=as.vector(weights))
    p <- as.vector(predict(bs_reg, newx=x, type="response"))

    if(i==1) plot(pi,p)

    for(j in 1:length(lambdas))
    {
      NB_model[j] <- NB_model[j] + mean((p - (1 - p) * lambdas[j] / (1 - lambdas[j])) * (pi > lambdas[j]))
      NB_all[j] <- NB_all[j] + mean((p - (1 - p) * lambdas[j] / (1 - lambdas[j])) * 1)
      NB_max[j] <- NB_max[j] + mean((p - (1 - p) * lambdas[j] / (1 - lambdas[j])) * (p > lambdas[j]))
    }
    #cat('.')
  }

  EVPP <- (NB_max-pmax(0,NB_model,NB_all))/n_sim

  res <-cbind(lambda=lambdas, EVPP=EVPP, NB_all=NB_all/n_sim, NB_model=NB_model/n_sim, NB_max=NB_max/n_sim)

  return(res)
}












bootstrap <- function (n, Bayesian=F)
{
  if(Bayesian)
  {
    u <- c(0,sort(runif(n-1)),1)
    return((u[-1] - u[-length(u)])*n)
  }
  else
  {
    u <- rmultinom(1,n,rep(1/n,n))
    return(u)
  }
}








#' @export
process_results <- function(res, graphs=c("evpp","summit"), overlay=F, ...)
{
  out <- list()
  out$einb_current <- mean(pmax(0,res[,'NB_model'],res[,'NB_all'])-pmax(0,res[,'NB_all']))
  out$einb_perfect <- mean(res[,'NB_max']-pmax(0,res[,'NB_all']))
  out$uncertainty_index <- out$einb_perfect/out$einb_current
  out$eevpp <- out$einb_perfect-out$einb_current

  if(overlay)
  {
    fn <- lines
  }
  else
  {
    fn <- plot
  }

  if(!is.na(match("evpp",graphs)))
  {
    fn(res[,'lambda'], res[,'NB_max']-pmax(0,res[,'NB_model'],res[,'NB_all']),type='l', xlab="threshold", ylab="EVPP", ...)
  }

  if(!is.na(match("summit",graphs)))
  {
    fn(res[,'lambda'],pmax(0,res[,'NB_model'],res[,'NB_all'])-pmax(0,res[,'NB_all']),type='l', xlab='threshold', ylab='INB', ...)
    lines(res[,'lambda'],res[,'NB_max']-pmax(0,res[,'NB_all']),type='l',col="red")
  }

  return(out)
}





