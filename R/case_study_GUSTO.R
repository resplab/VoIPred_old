library(glmnet)

set.seed(2222)

settings <- list()
settings$master_formula <- day30 ~ sex + age + dia + miloc + pmi + htn +  smk + tx

# bootstrap:0: parametric; 1:ordinary; 2:Bayesian
case_study_gusto <- function(n_sim=1000, subsample=1000)
{
  data("gusto")
  master_formula <- settings$master_formula
  if(!is.null(subsample))
  {
    sample <- gusto[sample(1:dim(gusto)[1],subsample,F),]
  }
  else
  {
    sample <- gusto
  }

  x <- model.matrix(master_formula,sample)
  y <- sample$day30

  cv_reg <- cv.glmnet(x, y, family="binomial")
  reg <- glmnet(x, y, family = "binomial", lambda = cv_reg$lambda.min)

  res0 <<- evpp.glmnet(reg, x, y, n_sim = n_sim, Bayesian_bootstrap = F)
  res1 <<- evpp.glmnet(reg, x, y, n_sim = n_sim, Bayesian_bootstrap = T)

  #tmp <- decision_curve(VoIPred:::aux$y, VoIPred:::aux$pi)
  plot(res0[,'lambda'],res0[,'dc_model']-res0[,'optimism'],type='l', xlab="Threshold", ylab="Net benefit")
  lines(res0[,'lambda'],res0[,'dc_model'],type='l',col="gray")
  lines(res0[,'lambda'],res0[,'NB_model'],type='l',col="green")
  lines(res1[,'lambda'],res1[,'NB_model'],type='l',col="blue")
  #lines(tmp[,'lambda'],tmp[,'NB_model'],col="red")

  #lines(res0[,'lambda'],res0[,'NB_max'],type='l',col="red")
  #lines(res1[,'lambda'],res1[,'NB_max'],type='l',col="orange")


  table_1 <- data.frame(
    cbind(
    covariate=colnames(VoIPred:::aux$coeffs),
    coeffs=format(t(unname(VoIPred:::aux$coeffs)),2,2),
    p_inc=format(unname(1-colMeans(VoIPred:::aux$bs_coeffs==0)),2,2),
    ci=paste(
      format(apply(X=VoIPred:::aux$bs_coeffs,2,FUN = quantile,0.025),2,2),
      format(apply(X=VoIPred:::aux$bs_coeffs,2,FUN = quantile,0.975),2,2),
      sep="–"
      )
    )
  )

  #rownames(table_1) <- colnames(VoIPred:::aux$coeffs)
  write.table(table_1,"clipboard",row.names = T)


  return(list(
    process_results(res0),
    process_results(res1),
    table_1=table_1
  ))
}




# Now change the sample size!
evpp_by_sample_size <- function(master_formula, n_sim=10, sample_sizes = c(250, 500, 1000, 2000, 8000, 32000, Inf))
{
  res <-NULL
  master_formula <- settings$master_formula

  for(i in 1:length(sample_sizes))
  {
    cat("sample size:",sample_sizes[i],"\n")
    if(is.infinite(sample_sizes[i]))
      sample <- gusto
    else
      sample <- gusto[sample(1:dim(gusto)[1],sample_sizes[i],F),]

    model_matrix <- model.matrix(master_formula,sample)
    cv_reg <- cv.glmnet(model_matrix, sample$day30, family="binomial")
    reg <- glmnet(model_matrix, sample$day30, family = "binomial", lambda = cv_reg$lambda.min)
    tmp <- process_results(evpp.glmnet(reg, model_matrix, sample$day30, n_sim = n_sim, Bayesian_bootstrap = F),graphs="")
    res <- rbind(res,c(sample_sizes[i],unlist(tmp)))
  }
  return(res)
}



calc_auc <- function(reg_obj, x, y, n_sim=1000)
{
  require(pROC)

  pi <- as.vector(predict(reg_obj,newx=x,type="response"))
  tmp <- pROC::roc(y,pi,quiet=T)
  plot(tmp)
  auc <- tmp$auc

  optimism <- 0

  n <-  length(y)
  for(i in 1:n_sim)
  {
    cat(".")

    bs <- sample(1:n,n,T)
    bs_x <- x[bs,]
    bs_y <- y[bs]

    tmp <- cv.glmnet(x=bs_x, y=bs_y, family="binomial")
    bs_reg <- glmnet(x=bs_x, y=bs_y, family="binomial", lambda=tmp$lambda.min)

    p_int <- as.vector(predict(bs_reg, newx=bs_x, type="response"))
    p_ext <- as.vector(predict(bs_reg, newx=x, type="response"))

    optimism <- optimism + pROC::roc(bs_y,p_int,quiet=T)$auc-pROC::roc(y,p_ext,quiet=T)$auc
  }

  return(c(auc=auc,optimism=optimism/n_sim))
}
