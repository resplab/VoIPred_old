library(tidyverse)
library(glmnet)
library(pROC)


set.seed(3344)



case_study_OPTIMAL <-function(covars=c("age", "gender", "nowsmk", "oxygen", "indicated_statin", "fev1", "sgrq", "ster1yr", "BMI"),
                              n_sim=1000, bootstrap=T, lambdas=(1:99)/100, only_severe=F, alpha=1, Nfold=3, penalized=T,
                              min_FUT=0.1)
{
  data(trials_data)

  df_optimal <- trials_data %>%
    filter(trial=="OPTIMAL") %>%
    mutate(hosp1yr = 1) %>%
    group_by(ID) %>%
    mutate(num_events = sum(event > only_severe*1)) %>%
    mutate(FUT = max(tte0)) %>%
    ungroup() %>%
    filter(period==1) %>%
    select(-trial,-Days_In_Study,-hosp1yr,-period,-tte0,-event,-sgrq10,-age10,-MACRO,-OPTIMAL,-STATCOPE,
           -event_date,-event_month,-event_season,-event_seasonality,-event,
           -rand_month,-rand_season,-rand_seasonality,-bmi10,-randomized_azithromycin,
           -randomized_LAMA,-randomized_statin,-fev1pp100,-stgtotalc00)

  message("N of the development set:", nrow(df_optimal))

  missing_index <- (is.na(df_optimal) %>% rowSums() != 0)
  message("N removed due to missing data from development set: ", sum(missing_index))

  df_dev <- df_optimal[!missing_index,]

  message("N of the development set:", nrow(df_dev))

  short_fu_index <- (df_dev$FUT < min_FUT)

  message("N removed due to censoring from development set: ",sum(short_fu_index)," (",round(sum(short_fu_index)/nrow(df_dev)*100,1),"%)")

  df_dev <- df_dev %>%
    filter(FUT >= min_FUT)

  message("Final N of the development set:", nrow(df_dev))

  df_dev$trtgroup <- df_dev$trtgroup-1

  # sanity check: collinearity
  assertthat::are_equal(0,
                        sum((cor(df_dev[,-1]) * (abs(cor(df_dev[,-1]))>0.9) - diag(ncol(df_dev)-1))))

  # #checker: trtgroup spread over randomized_ICS and randomized_LABA
  # look <- df_dev %>%
  #   select(trtgroup,randomized_ICS,randomized_LABA) %>%
  #   mutate(check1 = as.numeric(trtgroup==1),
  #          check2 = as.numeric(trtgroup!=2),
  #          check3 =( randomized_ICS==check1) & (randomized_LABA == check2))
  # assertthat::are_equal(0,sum(look$check3==F))

  #df_dev <- df_dev %>%
  #  select(-trtgroup,-ID)

  # standardize X
  if(length(covars)==1){
    if(covars=="all"){
      covars <- c("age", "gender", "nowsmk",
                  "oxygen", "indicated_statin",
                  "fev1", "LAMA", "LABA", "ICS", "sgrq", "ster1yr",
                  "packyears", "fvc", "BMI")
    }
  }

  formula <- as.formula(paste0("num_events~",paste(c(covars,"randomized_ICS","randomized_LABA"),collapse="+"),""))

  y_dev <- df_dev$num_events
  x_dev <- model.matrix(formula, df_dev)
  offset_dev <- log(df_dev$FUT)
  x_dev_notx <- x_dev
  x_dev_notx[,'randomized_LABA'] <- 0
  x_dev_notx[,'randomized_ICS'] <- 0

  # check
  if(penalized){
    cv.reg <- cv.glmnet(x=x_dev, y=y_dev, family="poisson",offset=offset_dev,
                        alpha=alpha,nfolds = Nfold,type.measure = "mse")
    message(paste("lambda is"),cv.reg$lambda.min)
    reg <- glmnet(x=x_dev, y=y_dev, family="poisson",offset=offset_dev,
                        alpha=alpha,lambda =  cv.reg$lambda.min)
  } else{
    reg <- glmnet(x=x_dev, y=y_dev, family="poisson",offset=offset_dev,
                        alpha=alpha,lambda =  0)
  }

  #Record coefficients:
  print(coefficients(reg))
  write.table(as.matrix(coefficients(reg)),file="clipboard")
  #writeClipboard(coefficients(reg))

  pi <- 1 - exp(-as.vector(predict(reg,type="response", newx=x_dev_notx, newoffset = 0)))

  print(paste("AUC is ", roc(y_dev>0, pi)$auc))

  sample_data <- x_dev
  sample_size <- nrow(sample_data)

  NB_test <-  NB_all <-  NB_max  <- rep(0, length(lambdas))

  for(i in 1:n_sim)
  {
    if(i%%100==0) cat(i,",")
    if(bootstrap)
    {
      bs_index <- sample(1:length(y_dev),sample_size,replace = T)
      x_bs <- x_dev[bs_index,]
      y_bs <- y_dev[bs_index]
      offset_bs <- offset_dev[bs_index]
      if(penalized){
        cv_bs <- cv.glmnet(x=x_bs, y=y_bs, family="poisson",offset=offset_bs,
                         alpha=alpha,nfold=Nfold,type.measure = "mse")
        reg_bs <- glmnet(x=x_bs, y=y_bs, family="poisson",offset=offset_bs,
                         alpha=alpha,
                         lambda =  cv_bs$lambda.min)
      } else{
        reg_bs <- glmnet(x=x_bs, y=y_bs, family="poisson",offset=offset_bs,
                         alpha=alpha,
                         lambda =  0)
      }
      p <- 1 - exp(-(as.vector(predict(reg_bs,type="response", newx=x_dev_notx,newoffset = 0))))
    }
    else
    {
      stop("NOT IMPLEMENTED")
    }

    if(i==1) plot(pi,p)

    for(j in 1:length(lambdas))
    {
      NB_test[j] <- NB_test[j] + mean((p - (1 - p) * lambdas[j] / (1 - lambdas[j])) * (pi > lambdas[j]))

      NB_all[j] <- NB_all[j] + mean((p - (1 - p) * lambdas[j] / (1 - lambdas[j])) * 1)

      NB_max[j] <- NB_max[j] + mean((p - (1 - p) * lambdas[j] / (1 - lambdas[j])) * (p > lambdas[j]))
    }
  }

  res <-cbind(lambda=lambdas, NB_all=NB_all/n_sim, NB_test=NB_test/n_sim, NB_max=NB_max/n_sim)

  draw_graphs(res, col="red")

  return(res)
}

#look <- case_study_OPTIMAL(covars='all',n_sim = 50,alpha=1,Nfold=5)




#' @export
main_optimal <- function(n_sim=1000, only_severe=F, bootstrap=0)
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

  reg <- glm(master_formula, data=dev_data, family=binomial(link="log"))

  res <<- evpp.glm(reg,n_sim = n_sim, bootstrap = bootstrap)

  draw_graphs(res, col="red")
}







#"gender","age10","oxygen","hosp1yr","sgrq10","fev1","nowsmk","LABA","LAMA"
#' @export
case_study<-function(covars=c("gender","age10","oxygen","sgrq10","fev1","nowsmk"), n_sim=1000, bootstrap=0, lambdas=(1:99)/100, only_severe=F)
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
stylized_sim <- function(betas=c(-2,0,1), sample_sizes=c(100,250,1000), n_sim=10000, lambdas=rep(1:19)/20)
{
  formula=paste("Y~",paste0("X",1:(length(betas)-1),collapse = "+"))

  for(i in 1:sample_sizes)
  {
    for(j in 0:2)
    {
      set.seed(1234)
      res <- stylized_sim_internal(betas, sample_sizes[i], n_sim, lambdas, bootstrap=j)
      draw_graphs(res, col="red")
    }
  }
  return(res)
}



stylized_sim_internal <- function(betas, sample_size, n_sim, lambdas, bootstrap)
{
  formula=paste("Y~",paste0("X",1:(length(betas)-1),collapse = "+"))

  sample_data <<-generate_data(sample_size, betas)

  reg <- glm(formula=formula, data=sample_data, family=binomial(link="logit"))

  res <- evpp.glm(reg, n_sim, bootstrap = bootstrap)

  return(res)
}




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



