#' Grid search for the optimal 1ts model
#'
#' @description `grid_search_1d()` performs a grid search for the minimum
#'   AIC or BIC of the one time scale model.
#'
#'   It finds the optimal values of `log_10(rho_s)` and returns the estimated
#'   optimal model.
#'
#' @param lrho A vector of `log_10(rho_s)` values.
#' @inheritParams get_aic_fit_1d
#' @param optim_criterion The criterion to be used for optimization:
#' `"aic"` (default) or `"bic"`.
#' @param par_gridsearch A list of parameters for the grid_search:
#'   * `plot_aic` A Boolean. Default is `FALSE`. If `TRUE`, plot the AIC values
#'     over the grid of `log_10(rhos)` values.
#'   * `plot_bic` A Boolean. Default is `FALSE`. If `TRUE`, plot the BIC values
#'     over the grid of `log_10(rhos)` values.
#'   * `return_aic` A Boolean. Default is `TRUE`. Return the AIC values.
#'   * `return_bic` A Boolean. Default is `TRUE`. Return the BIC values.
#'   * `mark_optimal` A Boolean. Default is `TRUE`. If the plot of the AIC or BIC
#'     values is returned, marks the optimal `log_10(rho_s)` in the plot.
#'   * `main_aic` The title of the AIC plot. Default is `"AIC grid"`.
#'   * `main_bic` The title of the BIC plot. Default is `"BIC grid"`.
#'
#' @return An object of class `h1tsfit` with the following elements:
#'   * `optimal_model` A list containing the results of the optimal model.
#'   * `optimal_logrho` The optimal value of `log10(rho_s)`.
#'   * `P_optimal` The optimal penalty matrix P.
#'   * `AIC` (if `par_gridsearch$return_aic == TRUE`) The vector of AIC values.
#'   * `BIC` (if `par_gridsearch$return_bic == TRUE`) The vector of BIC values.
#'
#' @importFrom graphics points
#' @export

grid_search_1d <- function(r, y,
                           Z = NULL,
                           lrho,
                           Bs,
                           Ds,
                           optim_criterion = c("aic", "bic"),
                           control_algorithm = list(
                             maxiter = 20,
                             conv_crit = 1e-5,
                             verbose = FALSE,
                             monitor_ev = FALSE
                           ),
                           par_gridsearch = list(
                             plot_aic = FALSE,
                             plot_bic = FALSE,
                             return_aic = TRUE,
                             return_bic = TRUE,
                             mark_optimal = TRUE,
                             main_aic = "AIC grid",
                             main_bic = "BIC grid"
                           )) {
  # ---- Controls for iterative process ----
  con <- control_algorithm
  gsp <- par_gridsearch

  # ---- Grid search ----
  aic <- lrho * 0
  bic <- lrho * 0

  for (i in 1:length(lrho)) {
    rho <- 10^lrho[i]
    P <- rho * t(Ds) %*% Ds

    if (is.null(Z)) {
      fit <- iwls_1d(
        r = r, y = y,
        Bs = Bs,
        P = P,
        control_algorithm = con
      )
    } else {
      fit <- GLAM_1d_covariates(
        R = r,
        Y = y,
        Z = Z,
        Bs = Bs,
        control_algorithm = con
      )
    }

    aic[i] <- fit$aic
    bic[i] <- fit$bic

    if (con$monitor_ev) {
      print(paste0(
        "AIC for model with log_10(rho_s) = ",
        lrho[i],
        " : ", aic[i]
      ))
      print(paste0(
        "BIC for model with log_10(rho_s) = ",
        lrho[i],
        " : ", bic[i]
      ))
    }
  }

  # ---- Optimal smoothing parameter ----
  which_opt_aic <- which(aic == min(aic))
  optim_lr_aic <- lrho[which_opt_aic]
  optim_r_aic <- 10^optim_lr_aic

  which_opt_bic <- which(bic == min(bic))
  optim_lr_bic <- lrho[which_opt_bic]
  optim_r_bic <- 10^optim_lr_bic

  if (gsp$plot_aic) {
    plot(lrho, aic,
      type = "l",
      main = gsp$main_aic
    )
    if (gsp$mark_optimal) {
      points(optim_lr_aic,
        aic[which_opt_aic],
        col = "red",
        pch = 15
      )
    }
  }
  if (gsp$plot_bic) {
    plot(lrho, bic,
      type = "l",
      main = gsp$main_bic
    )
    if (gsp$mark_optimal) {
      points(optim_lr_bic,
        bic[which_opt_bic],
        col = "red",
        pch = 15
      )
    }
  }

  # ---- Estimates optimal model ----
  if (optim_criterion == "aic") {
    optim_lr <- optim_lr_aic
    P <- optim_r_aic * t(Ds) %*% Ds
  } else {
    optim_lr <- optim_lr_bic
    P <- optim_r_bic * t(Ds) %*% Ds
  }

  # Model estimation
  if (is.null(Z)) {
    opt_mod <- iwls_1d(
      r = r, y = y,
      Bs = Bs,
      P = P,
      control_algorithm = con
    )
  } else {
    opt_mod <- GLAM_1d_covariates(
      R = r,
      Y = y,
      Z = Z,
      Bs = Bs,
      control_algorithm = con
    )
  }

  # ---- Save results in a list ----

  results <- list(
    "optimal_model" = opt_mod,
    "optimal_logrho" = optim_lr,
    "P_optimal" = P,
    "AIC" = aic,
    "BIC" = bic
  )
  if (gsp$return_aic == FALSE) {
    results$AIC <- NULL
  }
  if (gsp$return_bic == FALSE) {
    results$BIC <- NULL
  }

  class(results) <- "h1tsfit"
  return(results)
}
