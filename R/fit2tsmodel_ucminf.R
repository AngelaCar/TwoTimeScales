#' Numerical optimization of the 2ts model
#'
#' @description `fit2tsmodel_ucminf()` performs a numerical optimization of the
#'   AIC or BIC of the two time scales model.
#'
#'   It finds the optimal values of `log_10(rho_u)` and `log_10(rho_s)`
#'   and returns the estimated optimal model.
#'   See also [ucminf::ucminf()].
#'
#' @inheritParams get_aic_fit_2d
#' @param optim_criterion The criterion to be used for optimization:
#' `"aic"` (default) or `"bic"`.
#'
#' @return An object of class `h2tsfit` with the following elements:
#'   * `optimal_model` A list containing the results of the optimal model.
#'   * `optimal_logrho` A vector with the optimal values of `log10(rho_u)` and
#'     `log10(rho_s)`.
#'   * `P_optimal` The optimal penalty matrix P.
#'
#' @import ucminf
#' @export

fit2tsmodel_ucminf <- function(Y, R,
                               Z = NULL,
                               optim_criterion = c("aic", "bic"),
                               lrho = c(0, 0),
                               Bu, Bs,
                               Iu, Is,
                               Du, Ds,
                               Wprior = NULL,
                               ridge = 0,
                               control_algorithm = list()) {

  # ---- Controls for iterative process ----
  con <- list(
    maxiter = 20,
    conv_crit = 1e-5,
    verbose = FALSE,
    monitor_ev = FALSE
  )
  Ncon <- names(con)
  namesCon <- names(control_algorithm)

  con[namesCon] <- control_algorithm
  if(length(namesCon[!namesCon %in% Ncon])>0) {
    warning("Undefined entries in control! Default settings are used.\n")
    warning("undefined keyword(s): ",
            paste(namesCon[! namesCon %in% Ncon], collapse = ", "))
  }

  # ---- Find optimal smoothing parameters ----
  if (optim_criterion == "aic"){
    op <- ucminf::ucminf(
      par = lrho,
      fn = get_aic_fit_2d,
      R = R, Y = Y, Z = Z,
      Bu = Bu, Bs = Bs,
      Iu = Iu, Is = Is,
      Du = Du, Ds = Ds,
      Wprior = Wprior,
      ridge = ridge,
      control_algorithm = con
    )
  }
  if(optim_criterion == "bic"){
    op <- ucminf::ucminf(
      par = lrho,
      fn = get_bic_fit_2d,
      R = R, Y = Y, Z = Z,
      Bu = Bu, Bs = Bs,
      Iu = Iu, Is = Is,
      Du = Du, Ds = Ds,
      Wprior = Wprior,
      ridge = ridge,
      control_algorithm = con
    )
  }

  # ---- With optimal smoothing parameters, calculate optimal model ----
  optim_lr <- op$par
  optim_r <- 10^optim_lr

  # ---- Construct penalty matrix P  ----
  P <- optim_r[1] * kronecker(Is, t(Du) %*% Du) +
    optim_r[2] * kronecker(t(Ds) %*% Ds, Iu)

  # model estimation
  if (is.null(Z)) { # no covariates
    mod <- GLAM_2d_no_covariates(
      R = R, Y = Y,
      Bu = Bu, Bs = Bs,
      Wprior = Wprior,
      P = P,
      ridge = ridge,
      control_algorithm = con
    )
  } else {
    mod <- GLAM_2d_covariates(
      R = R, Y = Y, # with covariates
      Bu = Bu, Bs = Bs,
      Z = Z,
      Wprior = Wprior,
      P = P,
      ridge = ridge,
      control_algorithm = con
    )
  }

  # ---- save results in a list ----
  results <- list(
    "optimal_model" = mod,
    "optimal_logrho" = optim_lr,
    "P_optim" = P)

  class(results) <- "haz2ts"

  return(results)
}
