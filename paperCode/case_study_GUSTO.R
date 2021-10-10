library(glmnet)
library(VoIPred)
library(rms)

settings <- list()
settings$master_formula <- day30 ~ age + miloc + pmi + kill + pmin(sysbp,100) + lsp(pulse,50) + htn + dia
settings$default_th <- 0.02
settings$n_sim <- 100 #if 0 wont do this part
settings$custom_th <- c(0.01,0.02,0.05,0.1)
settings$n_sim <- 100 #if 0 wont do this part
settings$subsample <- 1000
settings$auc_n_sim <- 0   #If set to 0, it will not calculate AUC with optimism correction with the same n_sim.
settings$sample_size_n_sim_outer <- 0 #if set to 0 will not do
settings$sample_size_n_sim_inner <- 100 #Voi calculations for each point within each iteration
settings$sample_sizes <- c(250, 500, 1000, 2000, 4000, 8000, 16000, 32000, Inf)

case_study_gusto <- function(load_file=NULL, save_file=NULL)
{
  #assign("last.warning", NULL, envir = baseenv())
  set.seed(1234)
  results <<- list()

  data("gusto")
  #gusto <<- gusto[which(gusto$tx=="SK"),]
  gusto$kill <<- (as.numeric(gusto$Killip)>1)*1

  if(is.null(load_file))
  {
    master_formula <- settings$master_formula
    if(!is.null(settings$subsample))
    {
      sample <- gusto[sample(1:dim(gusto)[1],settings$subsample,F),]
    }  else
    {
      sample <- gusto
    }

    x <- model.matrix(master_formula,sample)
    y <- sample$day30

    cv_reg <- cv.glmnet(x, y, family="binomial")
    reg <- glmnet(x, y, family = "binomial", lambda = cv_reg$lambda.min)

    results$reg_obj <<- reg

    if(settings$auc_n_sim)
    {
      message("AUC calculations...")
      results$auc <<- calc_auc(reg, x, y, n_sim = settings$auc_n_sim)
    }

    if(settings$n_sim>0)
    {
      message("VoI calculations...")
      results$res0 <<- voi.glmnet(reg, x, y, n_sim = settings$n_sim, Bayesian_bootstrap = F)

      results$coeffs <<- VoIPred:::aux$coeffs
      results$bs_coeffs <<- VoIPred:::aux$bs_coeffs

      results$res1 <<- voi.glmnet(reg, x, y, n_sim = settings$n_sim, Bayesian_bootstrap = T)
    }

    if(settings$sample_size_n_sim_outer>0)
    {
      message("\nvoi by sample size calculations...\n")
      results$voi_by_sample_size <<- voi_by_sample_size(n_sim=settings$sample_size_n_sim_outer, sample_sizes = settings$sample_sizes)
    }
  }
  else  #load_file is not null so load the results
  {
    message("Loading from ",load_file)
    tmp <- readRDS(file=load_file)
    .GlobalEnv$settings <- tmp$settings
    .GlobalEnv$results <- tmp$results
  }

  if(settings$n_sim>0)
  {
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
  }

  if(!is.null(save_file))
  {
    message("Saving to ",save_file)
    tmp <- list(settings=settings, results=results)
    saveRDS(tmp, file=save_file)
  }

  if(settings$n_sim>0)
  {
    return(list(
      process_results(results$res0),
      process_results(results$res1),
      results$table1
    ))
  }
  else
    return()
}




# Now change the sample size!
voi_by_sample_size <- function(n_sim, sample_sizes)
{
  out <-NULL

  master_formula <- settings$master_formula

  work_horse <- function(size)
  {
    if(is.infinite(size))
    {
      sample <- gusto
    }
    else
      sample <- gusto[sample(1:dim(gusto)[1],size,replace=F),]

    model_matrix <- model.matrix(master_formula,sample)
    res <- tryCatch(
    {
      cv_reg <- cv.glmnet(model_matrix, sample$day30, family="binomial")
      reg <- glmnet(model_matrix, sample$day30, family = "binomial", lambda = cv_reg$lambda.min)
      if(sum(as.numeric(reg$beta))<2)
      {
        warning("Degenerate model!")
      }
      else
      {
        voi.glmnet(reg, model_matrix, sample$day30, n_sim = settings$sample_size_n_sim_inner, Bayesian_bootstrap = F)
      }
    }, error=function(w)
    {
      message("ERROR:",w)
      return(work_horse(size))
    }
    , warning=function(w)
      {
        message("WARNING:",w)
        return(work_horse(size))
      }
    )

    return(res)
  }

  for(i in 1:length(sample_sizes))
  {
    cat("\nsample size:",sample_sizes[i],"\n")

    voi_th <- voi_r <- rep(0,length(settings$custom_th))

    for(j in 1:n_sim)
    {
      res <- work_horse(sample_sizes[i])

      if(is.infinite(sample_sizes[i])) sample_sizes[i] <- dim(gusto)[1]

      index <- which(res[,'lambda'] %in% settings$custom_th)
      voi_th <- voi_th + res[index,'voi']/n_sim
      voi_r <- voi_r + process_results(res,graphs="",th=settings$custom_th)$voi_r/n_sim
    }

    out <- rbind(out,c(sample_sizes[i], voi_th=voi_th, voi_r=voi_r))
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
