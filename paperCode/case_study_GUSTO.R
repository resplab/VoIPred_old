library(glmnet)
library(VoIPred)

set.seed(1234)

settings <- list()
settings$master_formula <- day30 ~ sex + age + dia + miloc + pmi + htn +smk + kill + tx
settings$default_th <- 0.1
settings$n_sim <- 100
settings$subsample <- 1000
settings$do_auc <- T   #If set to true, it calculates AUC with optimism correction with the same n_sim.
settings$do_sample_size <- T   #If set to true, it calculates EVCP as a function of sample size.


case_study_gusto <- function(load_file=NULL, save_file=NULL)
{
  results <<- list()

  data("gusto")
  gusto$kill <<- (as.numeric(gusto$Killip)>1)*1

  if(is.null(load_file))
  {
    master_formula <- settings$master_formula
    if(!is.null(settings$subsample))
    {
      sample <- gusto[sample(1:dim(gusto)[1],settings$subsample,F),]
    }
    else
    {
      sample <- gusto
    }

    x <- model.matrix(master_formula,sample)
    y <- sample$day30

    cv_reg <- cv.glmnet(x, y, family="binomial")
    reg <- glmnet(x, y, family = "binomial", lambda = cv_reg$lambda.min)

    results$reg_obj <<- reg

    if(settings$do_auc)
    {
      message("AUC calculations...")
      results$auc <<- calc_auc(reg, x, y, n_sim = settings$n_sim)
    }

    results$res0 <<- evcp.glmnet(reg, x, y, n_sim = settings$n_sim, Bayesian_bootstrap = F)

    results$coeffs <<- VoIPred:::aux$coeffs
    results$bs_coeffs <<- VoIPred:::aux$bs_coeffs

    results$res1 <<- evcp.glmnet(reg, x, y, n_sim = settings$n_sim, Bayesian_bootstrap = T)

    if(settings$do_sample_size)
    {
      message("\nEVCP by sample size calculations...\n")
      results$evcp_by_sample_size <<- evcp_by_sample_size(n_sim=1)
    }
  }
  else  #load_file is not null so load the results
  {
    message("Loading from ",load_file)
    tmp <- readRDS(file=load_file)
    .GlobalEnv$settings <- tmp$settings
    .GlobalEnv$results <- tmp$results
  }

  #Decision curve plot (Figure 1)
  plot(results$res0[,'lambda'],results$res0[,'dc_model']-results$res0[,'optimism'],type='l', xlab="Threshold", ylab="Net benefit", col="red", lwd=2)
  lines(results$res0[,'lambda'],results$res0[,'dc_all'],type='l',col="black")
  lines(results$res0[,'lambda'],results$res0[,'dc_all']*0,type='l',col="gray")
  lines(results$res0[,'lambda'],results$res0[,'NB_model'],type='l',col="orange" ,lwd=2)
  lines(results$res1[,'lambda'],results$res1[,'NB_model'],type='l',col="blue", lwd=2)
  #lines(tmp[,'lambda'],tmp[,'NB_model'],col="red")


  results$table_1 <<- data.frame(
    covariate=colnames(results$coeffs),
    point_estimate=format(t(unname(results$coeffs)),2,2),
    p_inc=format(unname(1-colMeans(results$bs_coeffs==0)),2,2),
    ci=paste(
      format(round(apply(X=results$bs_coeffs,2,FUN = quantile,0.025),3),nsmall=3),
      format(round(apply(X=results$bs_coeffs,2,FUN = quantile,0.975),3),nsmall=3),
      sep=","
      )
  )

  #rownames(table_1) <- colnames(VoIPred:::aux$coeffs)
  write.table(results$table_1,"clipboard",row.names = T)

  if(!is.null(save_file))
  {
    message("Saving to ",save_file)
    tmp <- list(settings=settings, results=results)
    saveRDS(tmp, file=save_file)
  }

  return(list(
    process_results(results$res0),
    process_results(results$res1),
    results$table1
  ))
}




# Now change the sample size!
evcp_by_sample_size <- function(n_sim, sample_sizes = c(250, 500, 1000, 2000, 4000, 8000, 16000, 32000, Inf))
{
  out <-NULL

  master_formula <- settings$master_formula

  for(i in 1:length(sample_sizes))
  {
    cat("\nsample size:",sample_sizes[i],"\n")
    evcp_th <- 0
    evcp_r <- 0

    for(j in 1:n_sim)
    {
      if(is.infinite(sample_sizes[i]))
      {
        sample <- gusto
        sample_sizes[i] <- dim(gusto)[1]
      }
      else
        sample <- gusto[sample(1:dim(gusto)[1],sample_sizes[i],replace=F),]

      model_matrix <- model.matrix(master_formula,sample)
      cv_reg <- cv.glmnet(model_matrix, sample$day30, family="binomial")
      reg <- glmnet(model_matrix, sample$day30, family = "binomial", lambda = cv_reg$lambda.min)

      res <- evcp.glmnet(reg, model_matrix, sample$day30, n_sim = settings$n_sim, Bayesian_bootstrap = F)

      index <- which(res[,'lambda']==settings$default_th)
      evcp_th <- evcp_th + res[index,'EVCP']/n_sim
      evcp_r <- evcp_r + process_results(res,graphs="")$evcp_r/n_sim
    }

    out <- rbind(out,c(sample_sizes[i], evcp_th=evcp_th, evcp_r=evcp_r))
  }
  return(out)
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
