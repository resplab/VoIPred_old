library(tidyverse)
library(glmnet)

# Harry timesheet:
# 3 hours: April 11, 11am to 2pm
# 30 min : April 18, 2.30pm to 3pm
# : April 22, 9.30pm to 10.30pm

#' @export
draw_graphs <- function(res)
{
  plot(res[,'lambda'], res[,'NB_max']-pmax(0,res[,'NB_test'],res[,'NB_all']),type='l', xlab="z", ylab="EVPP", col='red')

  plot(res[,'lambda'],res[,'NB_max']-pmax(0,res[,'NB_all']),type='l', col='red', xlab='z', ylab='NB')
  lines(res[,'lambda'],(res[,'NB_test']-pmax(0,res[,'NB_all'])),type='l')

  plot(res[,'lambda'],res[,'NB_max'],type='l', col='red', xlab='z', ylab='NB')
  lines(res[,'lambda'],pmax(0,res[,'NB_test'],res[,'NB_all']),type='l')

  plot(res[,'lambda'], res[,'NB_max']/pmax(0,res[,'NB_test'],res[,'NB_all']) ,type='l',xlab='Threshold', ylab='Information Index', col='red')

  ii <- pmax(0,res[,'NB_test'],res[,'NB_all'])/res[,'NB_max']
  # ii[which(is.nan(ii))] <- 0
  plot(res[,'lambda'], ii,type='l',xlab='Threshold', ylab='Information Index 2', col='red')

}

case_study_OPTIMAL <-function(covars=c("gender","age","oxygen","sgrq","fev1","nowsmk","LAMA","LABA"),
                     n_sim=1000, bootstrap=T, lambdas=(1:99)/100, only_severe=F, alpha=1,Nfold=3,penalized=T,run_both=F,
                     min_FUT=0.05)
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

  df_dev <- df_dev %>%
    select(-trtgroup,-ID)

  # standardize X
  if(length(covars)==1){
    if(covars=="all"){
    covars <- c("age", "gender", "nowsmk",
                "oxygen", "fev1pp", "indicated_statin",
                "fev1", "LAMA", "LABA", "ICS", "randomized_LABA",
                "randomized_ICS", "sgrq", "Rand_Date", "ster1yr",
                "packyears", "fvc", "antib1yr", "BMI")
    }
  }

  formula <- as.formula(paste0("num_events~",paste(covars,collapse="+"),""))

  y_dev <- df_dev$num_events
  x_dev <- scale(df_dev %>%
    select(all_of(covars))) %>%
    as.matrix()
  offset_dev <- log(df_dev$FUT)


  # check
  if(penalized){
  cv.reg <- cv.glmnet(x=x_dev, y=y_dev, family="poisson",offset=offset_dev,
                      standardize = F,standardize.response = F,alpha=alpha,nfolds = Nfold)
  reg <- glmnet(x=x_dev, y=y_dev, family="poisson",offset=offset_dev,
                standardize = F,standardize.response = F,alpha=alpha,lambda =  cv.reg$lambda.min)
  } else{
    reg <- glmnet(x=x_dev, y=y_dev, family="poisson",offset=offset_dev,
                  standardize = F,standardize.response = F,alpha=alpha,lambda =  0)
  }
  pi <- 1 - exp(-as.vector(predict(reg,type="response", newx=x_dev,newoffset = 0)))


  sample_data <- x_dev
  sample_size <- nrow(sample_data)

  NB_test <-  NB_all <-  NB_max  <- rep(0, length(lambdas))

  for(i in 1:n_sim)
  {
    if(bootstrap)
    {
      bs_index <- sample(1:length(y_dev),sample_size,replace = T)
      x_bs <- x_dev[bs_index,]
      y_bs <- y_dev[bs_index]
      offset_bs <- offset_dev[bs_index]
      if(penalized){
        cv_bs <- cv.glmnet(x=x_bs, y=y_bs, family="poisson",offset=offset_bs,
                            standardize = F,standardize.response = F,alpha=alpha,nfold=Nfold)
        reg_bs <- glmnet(x=x_bs, y=y_bs, family="poisson",offset=offset_bs,
                      standardize = F,standardize.response = F,alpha=alpha,
                      lambda =  cv_bs$lambda.min)
      } else{
        reg_bs <- glmnet(x=x_bs, y=y_bs, family="poisson",offset=offset_bs,
                         standardize = F,standardize.response = F,alpha=alpha,
                         lambda =  0)
      }
      p <- 1 - exp(-(as.vector(predict(reg_bs,type="response", newx=x_dev,newoffset = 0))))
    }
    else
    {
      stop("NOT IMPLEMENTED")
      # new_betas <- rmvnorm(1, mean=mu, sigma=sigma)
      # p <- as.vector(1/(1+exp(-(new_betas %*% t(cbind(1,sample_data[,c(covars)]))))))
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

  draw_graphs(res)

  return(res)
}

N <- 10000
penalized <- case_study_OPTIMAL(covars='all',n_sim = N,alpha=1,Nfold=5,penalized = T)
not_penalized <- case_study_OPTIMAL(covars='all',n_sim = N,alpha=1,Nfold=5,penalized = F)

# save it
output_OPTIMAL <- rbind(penalized %>%
  as.data.frame() %>%
  mutate(type = "penalized"),
  (penalized %>%
     as.data.frame() %>%
     mutate(type = "unpenalized")))

library(here)

write_rds(output_OPTIMAL,here("results","output_OPTIMAL.rds"))

# coefs on the full data
return_full_model <-function(covars=c("gender","age","oxygen","sgrq","fev1","nowsmk","LAMA","LABA"),
                              n_sim=1000, bootstrap=T, lambdas=(1:99)/100, only_severe=F, alpha=1,Nfold=3,penalized=T,run_both=F,
                              min_FUT=0.05)
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

  # sanity check: collinearity
  assertthat::are_equal(0,
                        sum((cor(df_dev[,-1]) * (abs(cor(df_dev[,-1]))>0.9) - diag(ncol(df_dev)-1))))

  df_dev <- df_dev %>%
    select(-trtgroup,-ID)

  # standardize X
  if(length(covars)==1){
    if(covars=="all"){
      covars <- c("age", "gender", "nowsmk",
                  "oxygen", "fev1pp", "indicated_statin",
                  "fev1", "LAMA", "LABA", "ICS", "randomized_LABA",
                  "randomized_ICS", "sgrq", "Rand_Date", "ster1yr",
                  "packyears", "fvc", "antib1yr", "BMI")
    }
  }

  formula <- as.formula(paste0("num_events~",paste(covars,collapse="+"),""))

  y_dev <- df_dev$num_events
  x_dev <- scale(df_dev %>%
                   select(all_of(covars))) %>%
    as.matrix()
  offset_dev <- log(df_dev$FUT)


  # check
  if(penalized){
    cv.reg <- cv.glmnet(x=x_dev, y=y_dev, family="poisson",offset=offset_dev,
                        standardize = F,standardize.response = F,alpha=alpha,nfolds = Nfold)
    reg <- glmnet(x=x_dev, y=y_dev, family="poisson",offset=offset_dev,
                  standardize = F,standardize.response = F,alpha=alpha,lambda =  cv.reg$lambda.min)
  } else{
    reg <- glmnet(x=x_dev, y=y_dev, family="poisson",offset=offset_dev,
                  standardize = F,standardize.response = F,alpha=alpha,lambda =  0)
  }

  return(reg)
}

reg_p <- return_full_model(covars='all',n_sim = N,alpha=1,Nfold=5,penalized = T)
reg_up <- return_full_model(covars='all',n_sim = N,alpha=1,Nfold=5,penalized = F)
LASSO <- reg_p$beta[,1]
zero_index <- which(LASSO==0)
coef_table <- cbind(LASSO=format(round(LASSO,3),nsmall=3),
                    Unconstrained=format(round(reg_up$beta[,1],3),nsmall=3)) %>%
  data.frame()
coef_table$LASSO[zero_index] <- "."

write.csv(coef_table,here("results","coef_table_OPTIMAL.csv"),row.names=T)
