# bootstrap:0: parametric; 1:ordinary; 2:Bayesian
#' @export
case_study_gusto <- function(n_sim=1000, subsample=NULL)
{
  data("gusto")
  master_formula <- day30 ~ sex + age + dia + miloc + pmi + htn +  smk + tx
  if(!is.null(subsample))
  {
    sample <- gusto[sample(1:dim(gusto)[1],subsample,F),]
  }
  else
  {
    sample <- gusto
  }

  model_matrix <- model.matrix(master_formula,sample)
  cv_reg <- cv.glmnet(model_matrix, sample$day30, family="binomial")
  reg <- glmnet(model_matrix, sample$day30, family = "binomial", lambda = cv_reg$lambda.min)

  res0 <<- evpp.glmnet(reg, model_matrix, sample$day30, n_sim = n_sim, Bayesian_bootstrap = F)
  res1 <<- evpp.glmnet(reg, model_matrix, sample$day30, n_sim = n_sim, Bayesian_bootstrap = T)


  #draw_graphs(res0)
  #draw_graphs(res1)
  #draw_graphs(res2)

  #plot(res0[,'lambda'], res0[,'NB_max']-pmax(0,res0[,'NB_model'],res0[,'NB_all']),type='l', xlab="z", ylab="EVPP", col='red')
  #lines(res1[,'lambda'], res1[,'NB_max']-pmax(0,res1[,'NB_model'],res1[,'NB_all']),type='l', xlab="z", ylab="EVPP", col='green')

  return(c(
    process_results(res0),
    process_results(res1)
  ))
}

