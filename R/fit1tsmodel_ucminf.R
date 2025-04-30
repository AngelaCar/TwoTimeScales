#' Numerical optimization of the 1ts model
#'
#' @description `fit1tsmodel_ucminf()` performs a numerical optimization of the
#'   AIC or BIC of the one time scale model.
#'
#'   It finds the optimal values of \eqn{\log_{10}(\varrho_s)} and returns the estimated
#'   optimal model.
#'   See also [ucminf::ucminf()].
#'
#' @inheritParams get_aic_fit_1d
#' @param optim_criterion The criterion to be used for optimization:
#' `"aic"` (default) or `"bic"`.
#'
#' @return An object of class `haz1ts` with the following elements:
#'   * `optimal_model` A list containing the results of the optimal model.
#'   * `optimal_logrho` The optimal value of `log10(rho_s)`.
#'   * `P_optimal` The optimal penalty matrix P.
#'
#' @references Nielsen H, Mortensen S (2024).
#' _ucminf: General-Purpose Unconstrained Non-Linear Optimization_.
#' R package version 1.2.2, <https://CRAN.R-project.org/package=ucminf>
#'
#' @import ucminf
#' @keywords internal

fit1tsmodel_ucminf <- function(r, y,
                               Z = NULL,
                               lrho = 0,
                               Bs,
                               Ds,
                               Wprior = NULL,
                               optim_criterion = c("aic", "bic"),
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
      fn = get_aic_fit_1d,
      r = r, y = y,
      Z = Z,
      Bs = Bs,
      Ds = Ds,
      Wprior = Wprior,
      control_algorithm = con
    )
  }
  if(optim_criterion == "bic"){
    op <- ucminf::ucminf(
      par = lrho,
      fn = get_bic_fit_1d,
      r = r, y = y,
      Z = Z,
      Bs = Bs,
      Ds = Ds,
      Wprior = Wprior,
      control_algorithm = con
    )
  }

  # ---- With optimal smoothing parameters, calculate optimal model ----
  optim_lr <- op$par
  optim_r <- 10^optim_lr

  # ---- Construct penalty matrix P  ----
  P <- optim_r * t(Ds) %*% Ds

  if(is.null(Z)){
    opt_mod <- iwls_1d(
      r = r, y = y,
      Bs = Bs,
      P = P,
      control_algorithm = con
    )
  } else {
    opt_mod <- GLAM_1d_covariates(R = r,
                                  Y = y,
                                  Z = Z,
                                  Bs = Bs,
                                  P = P,
                                  Wprior = Wprior,
                                  control_algorithm = con)
  }

  # ---- save results in a list ----
  results <- list(
    "optimal_model" = opt_mod,
    "optimal_logrho" = optim_lr,
    "P_optim" = P)

  class(results) <- "haz1ts"
  return(results)
}
