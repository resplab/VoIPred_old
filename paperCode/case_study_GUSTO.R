library(glmnet)
library(VoIPred)
library(rms)
library(pROC)

machine_id <- round(runif(1)*10^10)

settings <- list()
settings$master_formula <- day30 ~ age + miloc + pmi + kill + pmin(sysbp,100) + lsp(pulse,50) + htn + dia
settings$default_th <- 0.02
settings$custom_th <- c(0.01,0.02,0.05,0.1)
settings$n_sim <- 10 # if 0 wont do this part
settings$subsample <- 1000
settings$auc_n_sim <- 1   #If set to 0, it will not calculate AUC with optimism correction with the same n_sim.
settings$sample_size_n_sim_outer <- 0 #if set to 0 will not do
settings$sample_size_n_sim_inner <- 1000 #Voi calculations for each point within each iteration
settings$sample_sizes <- c(500, 1000, 2000, 4000, 8000, 16000, 32000, Inf)

case_study_gusto <- function(load_file=NULL, save_file=NULL, seed=1234)
{
  #assign("last.warning", NULL, envir = baseenv())
  set.seed(seed)
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

    pi <- predict(reg,newx=x, type="response")

    results$reg_obj <<- reg

    if(settings$auc_n_sim)
    {
      message("AUC calculations...")
      results$auc <<- calc_auc(reg, x, y, n_sim = settings$auc_n_sim)
    }

    if(settings$n_sim>0)
    {
      message("VoI calculations...")
      results$res0 <<- voi.glmnet2(master_formula, sample, pi, n_sim = settings$n_sim, Bayesian_bootstrap = F)

      results$coeffs <<- VoIPred:::aux$coeffs
      results$bs_coeffs <<- VoIPred:::aux$bs_coeffs

      results$res1 <<- voi.glmnet2(master_formula, sample, pi, n_sim = settings$n_sim, Bayesian_bootstrap = T)
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
      covariate=names(coefficients(reg)[-2,]),
      point_estimate=format(t(unname(coefficients(reg)[-2])),2,2),
      p_inc=format(unname(1-colMeans(results$bs_coeffs==0)),2,2),
      ci=paste(
        format(round(apply(X=results$bs_coeffs,2,FUN = quantile,0.025),3),nsmall=3),
        format(round(apply(X=results$bs_coeffs,2,FUN = quantile,0.975),3),nsmall=3),
        sep=","
        )
    )

    #rownames(table_1) <- colnames(VoIPred:::aux$coeffs)
    # write.table(results$table_1,"clipboard",row.names = T)
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
        require(glmnet)
        cv_reg <- cv.glmnet(model_matrix, sample$day30, family="binomial")
        reg <- glmnet(model_matrix, sample$day30, family = "binomial", lambda = cv_reg$lambda.min)
        if(sum(as.numeric(reg$beta))<2)
        {
          warning("Degenerate model!")
        }
        else
        {
          voi.glmnet2(master_formula, sample, n_sim = settings$sample_size_n_sim_inner, Bayesian_bootstrap = F)
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

    out <- rbind(out,c(machine_id=machine_id, sample_size=sample_sizes[i], voi_th=voi_th, voi_r=voi_r))
    # GRpush(out,T)
  }
  return(out)
}


# save each run
voi_by_sample_size_custom <- function(n_sim, sample_sizes,type)
{
  out <-NULL

  master_formula <- settings$master_formula

  work_horse <- function(size,Bayesian_bootstrap=F)
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
          # warning("Degenerate model!")
        }
        else
        {
          voi.glmnet(reg, model_matrix, sample$day30, n_sim = settings$sample_size_n_sim_inner, Bayesian_bootstrap = Bayesian_bootstrap)
        }
      }, error=function(w)
      {
        # message("ERROR:",w)
        return(work_horse(size))
      }
      , warning=function(w)
      {
        # message("WARNING:",w)
        return(work_horse(size))
      }
    )

    return(res)
  }

  for(i in 1:length(sample_sizes))
  {
    cat("\nsample size:",sample_sizes[i],"\n")

    # voi_th <- voi_r <- rep(0,length(settings$custom_th))


    if(type %in% c("both","regular")){
    valid_result <- TRUE

    while(valid_result)
    {
      if(is.infinite(sample_sizes[i])) sample_sizes[i] <- dim(gusto)[1]

      res <- work_horse(sample_sizes[i])

      index <- which(res[,'lambda'] %in% settings$custom_th)
      voi_th <- res[index,'voi']
      voi_r <- process_results(res,graphs="",th=settings$custom_th)$voi_r

      tmp_result <- c(sample_sizes[i], voi_th=voi_th,voi_r=voi_r,Bayesian=F)

      if(length(tmp_result) > 2*length(settings$custom_th)){
        print(tmp_result)
        write_rds(tmp_result,here("simulation_result","raw00",paste0("sim_",sample_sizes[i],"_",n_sim,".rds")))
        valid_result <- FALSE
      }
    }
    }

    if(type %in% c("both","Bayesian")){

    valid_result <- TRUE

    while(valid_result)
    {
      if(is.infinite(sample_sizes[i])) sample_sizes[i] <- dim(gusto)[1]

      res2 <- work_horse(sample_sizes[i],Bayesian_bootstrap=T)

      index <- which(res2[,'lambda'] %in% settings$custom_th)
      voi_th_Bayesian <- res2[index,'voi']
      voi_r_Bayesian <- process_results(res2,graphs="",th=settings$custom_th)$voi_r
      tmp_result <- c(sample_sizes[i], voi_th=voi_th_Bayesian,voi_r=voi_r_Bayesian,Bayesian=T)

      if(length(tmp_result) > 2*length(settings$custom_th)){
        print(tmp_result)
        write_rds(tmp_result,here("simulation_result","raw11",paste0("sim_",sample_sizes[i],"_",n_sim,".rds")))
        valid_result <- FALSE
      }
    }
    }
  }

  return("")
}

# for each sample, for each iteration, save the result and combine later
# 500

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



perturbed_scenario <- function(sample_size, event_p=NA, shrinkage_type="", shrinkage_factor=NA, bias_OR=NA, n_sim=200, seed=1234)
{
  set.seed(seed)

  sample <- gusto[sample((1:dim(gusto)[1]),sample_size),]

  if(!is.na(event_p))
  {
    cases <- which(sample$day30==1)
    controls <- which(sample$day30==0)
    p <- length(cases) / sample_size
    q <- event_p
    weights <- rep(NA, sample_size)
    weights[cases] <- q*(1-p)/p/(1-q)
    weights[controls] <- 1
    weights <- weights / sum(weights) * sample_size
  }
  else
  {
    weights <- rep(1,sample_size)
  }

  sample$kill <- (as.numeric(sample$Killip)>1)*1

  master_formula <- settings$master_formula

  x <- model.matrix(master_formula,sample)
  y <- sample$day30

  if(shrinkage_type=="auc")
  {
    cv_reg <- cv.glmnet(x, y, family="binomial", weights = weights, type.measure = "auc")
    best <- which.min(abs(cv_reg$cvm-shrinkage_factor))
    lambda <- cv_reg$lambda[best]
  }
  else if(shrinkage_type=="lambda")
  {
    lambda <- shrinkage_factor
  }
  else
  {
    cv_reg <- cv.glmnet(x, y, family="binomial", weights = weights)
    lambda <- cv_reg$lambda.min
  }

  reg <- glmnet(x, y, family = "binomial", lambda = lambda, weights = weights)

  pi <- as.vector(predict(reg, newx=x, type="response"))

#bias and intercept of the proposed model
  log_lins0 <- log(pi/(1-pi))
  tmp <- glm(y~log_lins0, family = binomial(link="logit"), weights = weights)
  b0 <- coefficients(tmp)[1]
  b1 <- coefficients(tmp)[2]

  auc <- roc(y,as.vector(predict(reg,newx=x)))$auc
  message("AUC is", auc)

  if(shrinkage_type=="noise")
  {
    log_lins <- log(pi/(1-pi)) + rnorm(length(pi),0,shrinkage_factor)
    #Make sure the model has the same intercept
    tmp <- glm(y~log_lins, family = binomial(link="logit"))
    B0 <- coefficients(tmp)[1]
    B1 <- coefficients(tmp)[2]
    log_lins1 <- -b0/b1 + (B0 + B1*log_lins)/b1
    pi <- as.vector(1/(1+exp(-log_lins1)))
  }

  if(!is.na(bias_OR))
  {
    reg$a0 <- reg$a0 + log(bias_OR)
    pi <- as.vector(predict(reg,newx=x,type="response"))
  }

  auc <- roc(y,pi)$auc
  message("AUC is", auc)

  res <- voi.glmnet2(master_formula, sample, pi, lambdas=settings$custom_th, n_sim = n_sim, weights = weights)

  c(res[,'voi'], auc=auc)
}


generte_perturbed_scenario <- function(seed=1234)
{
  out <- data.frame("sample_size"=integer(),"event_p"=double(), "noise_sd"=double(),"bias_OR"=double(),"voi_1"=double(), "voi_2"=double(),"voi_3"=double(),"voi_4"=double(),  "auc"=double())
  sample_sizes <- c(500,1000,2500,5000)
  event_ps <- c(NA,0.15,0.3,0.5)
  noise_sds <- c(1/3,2/3,1,3/2)
  bias_ORs <- c(1/2,3/4,4/3,2)
  for(sample_size in sample_sizes)
  {
    for(event_p in event_ps)
    {
      res <- perturbed_scenario(sample_size, event_p = event_p, seed=seed)
      out <- rbind(out,c(sample_size, event_p, NA, NA, res))
    }
    for(noise_sd in noise_sds)
    {
      res <- perturbed_scenario(sample_size, shrinkage_type="noise", shrinkage_factor = noise_sd, seed=seed)
      out <- rbind(out,c(sample_size, NA, noise_sd, NA, res))
    }
    for(bias_OR in bias_ORs)
    {
      res <- perturbed_scenario(sample_size, bias_OR = bias_OR, seed=seed)
      out <- rbind(out,c(sample_size, NA, NA, bias_OR, res))
    }
  }
  colnames(out) <- c("sample_size","event_p", "n_covar_remove","bias_OR","voi_1","voi_2","voi_3","voi_4", "auc")

  out
}

