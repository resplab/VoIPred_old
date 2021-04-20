#' @import mvtnorm


#' @export
generate_data <- function(n, betas)
{
  n_var <- length(betas)-1 #one is intercept
  X <- NULL

  for(i in 1:n_var)
  {
    X <- cbind(X, rnorm(n))
  }

  colnames(X)<-paste0("X",1:n_var)

  lin <- betas %*% t(cbind(1,X))

  Y <- rbinom(n, 1, 1/(1+exp(-lin)))

  message(paste0("P(Y=1)=",mean(Y)))

  return(as.data.frame(cbind(X, Y=Y)))
}






#' @export
stylized_sim <- function(betas=c(-2,0,1), sample_size=100, n_sim=10000, lambdas=rep(1:19)/20, bootstrap=F, reuse_data=F)
{
  formula=paste("Y~",paste0("X",1:(length(betas)-1),collapse = "+"))

  if(!reuse_data) sample_data <<-generate_data(sample_size, betas)

  reg <- glm(formula=formula, data=sample_data, family=binomial(link="logit"))

  mu <- coefficients(reg)
  sigma <- vcov(reg)

  pi <- predict(reg, type="response")
  sample_data[,'pi'] <- pi

  NB_test <- rep(0, length(lambdas))
  NB_all <- NB_test
  NB_max <- NB_test

  for(i in 1:n_sim)
  {
    if(bootstrap)
    {
      bs_data<-sample_data[ sample((1:sample_size),sample_size,replace = T),]
      bs_reg <- glm(data=bs_data, formula=Y~X1+X2, family = binomial(link="logit"))
      p <- as.vector(predict(bs_reg, newdata=sample_data, type="response"))
    }
    else
    {
      new_betas <- rmvnorm(1, mean=mu, sigma=sigma)
      p <- as.vector(1/(1+exp(-(new_betas %*% t(cbind(1,sample_data[,c('X1','X2')]))))))
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

  res <<- cbind(lambda=lambdas, NB_all=NB_all/n_sim, NB_test=NB_test/n_sim, NB_max=NB_max/n_sim)

  draw_graphs(res)

  return(res)
}



#' @export
draw_graphs <- function(res)
{
  plot(res[,'lambda'], res[,'NB_max']-pmax(0,res[,'NB_test'],res[,'NB_all']),type='l', xlab="z", ylab="EVPP", col='red')

  plot(res[,'lambda'],res[,'NB_max']-pmax(0,res[,'NB_all']),type='l', col='red', xlab='z', ylab='NB')
  lines(res[,'lambda'],res[,'NB_test']-pmax(0,res[,'NB_all']),type='l')

  plot(res[,'lambda'],res[,'NB_max'],type='l', col='red', xlab='z', ylab='NB')
  lines(res[,'lambda'],pmax(0,res[,'NB_test'],res[,'NB_all']),type='l')

  plot(res[,'lambda'], res[,'NB_max']/pmax(0,res[,'NB_test'],res[,'NB_all']) ,type='l',xlab='Threshold', ylab='Information Index', col='red')

  ii <- pmax(0,res[,'NB_test'],res[,'NB_all'])/res[,'NB_max']
 # ii[which(is.nan(ii))] <- 0
  plot(res[,'lambda'], ii,type='l',xlab='Threshold', ylab='Information Index 2', col='red')

}







#"gender","age10","oxygen","hosp1yr","sgrq10","fev1","nowsmk","LABA","LAMA"
#' @export
case_study<-function(covars=c("gender","age10","oxygen","sgrq10","fev1","nowsmk"), n_sim=1000, bootstrap=F, lambdas=(1:99)/100, only_severe=F)
{
  data(trials_data)
  trials_data[which(trials_data[,'trial']=="OPTIMAL"),'hosp1yr']<<-1

  #eclipse_data<-readRDS("validatedECLIPSE.RData")

  dev_data<-trials_data[which(trials_data[,'trial']=="OPTIMAL"),]
  dev_data<-as.data.frame(dev_data)
  dev_data<-dev_data[which(dev_data[,'period']==1),]
  missing_data<-which(is.na(rowSums(dev_data[,covars])))
  message("N removed due to missing data from development set:",length(missing_data))
  dev_data<-dev_data[-missing_data,]
  short_fu<-which(dev_data[,'event']==0 & dev_data[,'tte0']<0.5)
  message("N removed due to censoring from development set:",length(short_fu),"(",length(short_fu)/dim(dev_data)[1],")")
  dev_data<-dev_data[-short_fu,]
  dev_data['event_bin']<-(dev_data['event']>only_severe*1)*1


  formula<-as.formula(paste0("event_bin~",paste(covars,collapse="+"),""))
  reg<<-glm(data=dev_data,formula = formula, family=binomial(link="logit"))
  mu <- coefficients(reg)
  sigma <- vcov(reg)



  sample_data <<- dev_data #[1:500,]
  sample_size <- dim(sample_data)[1]
  pi<-predict(reg,type="response", newdata=sample_data)

  NB_test <- rep(0, length(lambdas))
  NB_all <- NB_test
  NB_max <- NB_test

  for(i in 1:n_sim)
  {
    if(bootstrap)
    {
      bs_data<-sample_data[sample((1:sample_size),sample_size,replace = T),]
      bs_reg <- glm(data=bs_data,formula = formula, family=binomial(link="logit"))
      p <- as.vector(predict(bs_reg, newdata=sample_data, type="response"))
    }
    else
    {
      new_betas <- rmvnorm(1, mean=mu, sigma=sigma)
      p <- as.vector(1/(1+exp(-(new_betas %*% t(cbind(1,sample_data[,c(covars)]))))))
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

  res <-cbind(lambda=lambdas, NB_all=NB_all/n_sim, NB_test=NB_test/n_sim, NB_max=NB_max/n_sim)

  draw_graphs(res)

  return(res)
}








#' @export
EVPP.glm<-function(reg_obj, n_sim=1000, bootstrap=F, lambdas=(1:99)/100)
{
  sample_size <- dim(reg$data)[1]
  mu <- coefficients(reg)
  sigma <- vcov(reg)

  pi<-predict(reg,type="response")

  NB_test <- rep(0, length(lambdas))
  NB_all <- NB_test
  NB_max <- NB_test

  for(i in 1:n_sim)
  {
    if(bootstrap)
    {
      bs_data<-reg$data[sample((1:sample_size),sample_size,replace = T),]
      bs_reg <- glm(data=bs_data,formula = reg$call$formula, family=reg$call$family)
      p <- as.vector(predict(bs_reg, newdata=reg$data, type="response"))
    }
    else
    {
      new_betas <- rmvnorm(1, mean=mu, sigma=sigma)
      p <- as.vector(1/(1+exp(-(new_betas %*% t(model.matrix(reg))))))
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
EVPP.glmnet <- function(reg_obj, x, y, n_sim=1000, lambdas=(1:99)/100)
{
  sample_size <- dim(x)[1]

  pi<-predict(reg_obj, type="response", newx=x)

  NB_test <- rep(0, length(lambdas))
  NB_all <- NB_test
  NB_max <- NB_test

  for(i in 1:n_sim)
  {
    cat(".")
    bs <- sample((1:sample_size),sample_size,replace = T)
    bs_x <- x[bs,]
    bs_y <- y[bs]
    tmp <- cv.glmnet(bs_x, bs_y, family="binomial")
    bs_reg <- glmnet(x=bs_x, y=bs_y, family="binomial", lambda=tmp$lambda.min)
    p <- as.vector(predict(bs_reg, newx=x, type="response"))

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
main_gusto <- function(n_sim=1000)
{
  data("gusto")
  master_formula <- day30 ~ sex + age + dia + miloc + pmi + htn +  smk + tx
  subsample <- sample(1:dim(gusto)[1],10000,T)
  x<-model.matrix(master_formula, data=gusto[subsample,])
  y<-gusto$day30[subsample]
  tmp <- cv.glmnet(x, y, family="binomial")
  regnet <<- glmnet(x,y,family="binomial",lambda=tmp$lambda.min)
  resnet <<- EVPP.glmnet(regnet,x,y,n_sim=n_sim)

  plot(resnet[,1:2],type='l')
}





#' @export
main_optimal <- function(n_sim=1000, only_severe=1)
{
  covars <- c("gender","age10","oxygen","hosp1yr","sgrq10","fev1","nowsmk","LABA","LAMA")

  data(trials_data)
  trials_data[which(trials_data[,'trial']=="OPTIMAL"),'hosp1yr']<<-1

  dev_data<-trials_data[which(trials_data[,'trial']=="OPTIMAL"),]
  dev_data<-as.data.frame(dev_data)
  dev_data<-dev_data[which(dev_data[,'period']==1),]
  missing_data<-which(is.na(rowSums(dev_data[,covars])))
  message("N removed due to missing data from development set:",length(missing_data))
  dev_data<-dev_data[-missing_data,]
  short_fu<-which(dev_data[,'event']==0 & dev_data[,'tte0']<0.5)
  message("N removed due to censoring from development set:",length(short_fu),"(",length(short_fu)/dim(dev_data)[1],")")
  dev_data<-dev_data[-short_fu,]
  dev_data['event_bin']<-(dev_data['event']>only_severe*1)*1

  master_formula<-as.formula(paste0("event_bin~",paste(covars,collapse="+"),""))

  subsample <- (1:dim(dev_data)[1]) #sample(1:dim(dev_data)[1],100,T)


  x<-model.matrix(master_formula, data=dev_data[subsample,])
  y<-dev_data[subsample,"event_bin"]

  tmp <- cv.glmnet(x, y, family="binomial")
  regnet <<- glmnet(x,y,family="binomial",lambda=tmp$lambda.min)
  resnet <<- EVPP.glmnet(regnet,x,y,n_sim=n_sim)

  plot(resnet[,1:2],type='l')
}
