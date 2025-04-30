#' Calculates AIC and BIC from object fitted via LMMsolver
#'
#' @description `getAIC_BIC_LMM` is an utility function that takes an object of class
#'              `'LMMsolve'` fitted via `fit1ts()` or `fit2ts()` and calculates
#'              AIC, BIC and ED.
#'
#' @param fit An object of class `"LMMsolve"`
#' @param offset The vector of exposure times from dataLMM
#'
#' @return A list with:
#'          * `ED` effective dimension of the full model;
#'          * `EDbase` effective dimension of the baseline hazard only;
#'          * `Dev` deviance;
#'          * `AIC` the aic;
#'          * `BIC` the bic;
#'          * `n_beta` the number of estimated covariate parameters (if PH model).
#' @keywords internal
#'
getAIC_BIC_LMM <- function(fit, offset){
  n_beta <- length(fit$term.labels.f) - 2  # (intercept in baseline)

  ED <- sum(rev(fit$EDdf$Effective)[-1])
  ED_baseline <- ED - n_beta  #

  # deviance and AIC
  mu_hat <- offset * exp(fit$yhat)
  non_zero <- (fit$y >0)
  deviance <- 2 * sum(fit$y[non_zero] * log(fit$y[non_zero]/mu_hat[non_zero]))
  AIC <- deviance + 2*ED
  n_obs <- length(fit$y)
  BIC <- deviance + ED * log(n_obs)
  out <- data.frame(ED=ED, EDbase = ED_baseline, Dev=deviance, AIC=AIC, BIC=BIC,
                    n_beta=n_beta)
  return(out)
}
