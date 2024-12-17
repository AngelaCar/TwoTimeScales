#' Return the BIC of 1ts model
#'
#' @description `get_bic_fit_1d()` fits the 1ts model with or without individual
#' level covariates and it returns the BIC of the model.
#' See also `fit1tsmodel_ucminf()` and `fit1ts()`.
#'
#' @inheritParams get_aic_fit_1d
#'
#' @return the `bic` value of the fitted model.


get_bic_fit_1d <- function(lrho,
                           r, y,
                           Z = NULL,
                           Bs,
                           Ds,
                           Wprior = NULL,
                           control_algorithm = list(
                             maxiter = 20,
                             conv_crit = 1e-5,
                             verbose = FALSE,
                             monitor_ev = FALSE
                           )) {

  rho <- 10^lrho
  P <- rho * t(Ds) %*% Ds

  if (is.null(Z)) {
    fit <- iwls_1d(
      r = r, y = y,
      Bs = Bs,
      P = P,
      Wprior = Wprior,
      control_algorithm = control_algorithm
    )
  } else {
    fit <- GLAM_1d_covariates(
      R = r,
      Y = y,
      Z = Z,
      Bs = Bs,
      P = P,
      Wprior = Wprior,
      control_algorithm = control_algorithm
    )
  }
  bic <- fit$bic

  if (control_algorithm$monitor_ev) {
    print(paste0(
      "BIC for model with log_10(rho_s) = ",
      lrho,
      " : ", bic
    ))
  }

  return(bic)
}
