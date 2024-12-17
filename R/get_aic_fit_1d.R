#' Return the AIC of 1ts model
#'
#' @description `get_aic_fit_1d()` fits the 1ts model with or without individual
#' level covariates and it returns the AIC of the model.
#' See also `fit1tsmodel_ucminf()` and `fit1ts()`.
#'
#' @param r A vector of exposure times of length ns, or an array of dimension ns
#'   by n.
#' @param y A vector of event counts of length ns, or an array of dimension ns
#'   by n.
#' @param Z (optional) A regression matrix of covariates of dimension n by p.
#' @param lrho A starting value for \eqn{\log_{10}(\varrho_s)}. Default is 0.
#' @param Bs A matrix of B-splines for the time scale `s`.
#' @param Ds The difference matrix of the penalty.
#' @param Wprior An optional vector of a-priori weights.
#' @param control_algorithm A list with optional values for the parameters of
#'   the iterative processes:
#'   * `maxiter` The maximum number of iteration for the IWSL algorithm.
#'     Default is 20.
#'   * `conv_crit` The convergence criteria, expressed as difference between
#'     estimates at iteration i and i+1. Default is `1e-5`.
#'   * `verbose` A Boolean. Default is `FALSE`. If `TRUE` monitors the iteration
#'     process.
#'   * `monitor_ev` A Boolean. Default is `FALSE`. If `TRUE` monitors the
#'     evaluation of the model over the `log_10(rho_s)` values.
#'
#' @return The `aic` value of the fitted model.


get_aic_fit_1d <- function(lrho,
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

  aic <- fit$aic

  if (control_algorithm$monitor_ev) {
    print(paste0(
      "AIC for model with log_10(rho_s) = ",
      lrho,
      " : ", aic
    ))
  }

  return(aic)
}
