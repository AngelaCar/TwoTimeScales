#' Grid search for the optimal 2ts model
#'
#' @description `grid_search_2d()` performs a grid search for the minimum
#'   AIC or BIC of the two time scales model.
#'
#'   It finds the optimal values of `log_10(rho_u)` and `log_10(rho_s)` and
#'   returns the estimated optimal model.
#'
#' @param lru A vector of `log_10(rho_u)` values.
#' @param lrs A vector of `log_10(rho_s)` values.
#' @inheritParams get_aic_fit_2d
#' @param optim_criterion The criterion to be used for optimization:
#' `"aic"` (default) or `"bic"`.
#' @param par_gridsearch A list of parameters for the grid_search:
#'   * `plot_aic` A Boolean. Default is `FALSE`. If `TRUE`, plot the AIC values
#'     over the grid of `log_10(rho_u)` and `log_10(rho_s)` values.
#'   * `plot_bic` A Boolean. Default is `FALSE`. If `TRUE`, plot the BIC values
#'     over the grid of `log_10(rho_u)` and `log_10(rho_s)` values.
#'   * `return_aic` A Boolean. Default is `TRUE`. Return the AIC values.
#'   * `return_bic` A Boolean. Default is `TRUE`. Return the BIC values.
#'   * `col` The color palette to be used for the AIC/BIC plot. Default is
#'     `grDevices::gray.colors(n=10)`.
#'   * `plot_contour` A Boolean. Default is `TRUE`. Adds white contour lines to
#'     the AIC/BIC plot.
#'   * `mark_optimal` A Boolean. Default is `TRUE`. If the plot of the AIC or BIC
#'     values is returned, marks the optimal combination of `log_10(rho_u)` and
#'     `log_10(rho_s)` in the plot.
#'   * `main_aic` The title of the AIC plot. Default is `"AIC grid"`.
#'   * `main_bic` The title of the BIC plot. Default is `"BIC grid"`.
#'
#' @return An object of class `h2tsfit` with the following elements:
#'   * `optimal_model` A list containing the results of the optimal model.
#'   * `optimal_logrho` The optimal couple of `log_10(rho_u)` and `log_10(rho_s)`
#'      values.
#'   * `P_optimal` The optimal penalty matrix P.
#'   * `AIC` (if `par_gridsearch$return_aic == TRUE`) The vector of AIC values.
#'   * `BIC` (if `par_gridsearch$return_bic == TRUE`) The vector of BIC values.
#'
#' @importFrom grDevices gray.colors
#' @importFrom graphics contour points
#' @export

grid_search_2d <- function(lru, lrs,
                           R, Y,
                           Bu, Bs,
                           Z = NULL,
                           Iu, Is,
                           Du, Ds,
                           Wprior = NULL,
                           ridge = 0,
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
                             col = grey.colors(n = 10),
                             plot_contour = FALSE,
                             mark_optimal = TRUE,
                             main_aic = "AIC grid",
                             main_bic = "BIC grid"
                           )) {
  # ---- Controls for iterative process and grid_search parameters ----
  con <- control_algorithm
  gsp <- par_gridsearch

  # ---- Grid search ----
  AIC <- matrix(
    0,
    length(lru),
    length(lrs)
  )
  BIC <- matrix(
    0,
    length(lru),
    length(lrs)
  )

  for (i in 1:length(lru)) {
    for (j in 1:length(lrs)) {
      ru <- 10^lru[i]
      rs <- 10^lrs[j]

      # ---- Construct the penalty matrix ----
      P <- ru * kronecker(Is, t(Du) %*% Du) + rs * kronecker(t(Ds) %*% Ds, Iu)

      # ---- Estimate the model with the correct function ----
      if (is.null(Z)) { # without covariates
        fit <- GLAM_2d_no_covariates(
          R = R, Y = Y,
          Bu = Bu, Bs = Bs,
          Wprior = Wprior,
          P = P,
          ridge = ridge,
          control_algorithm = con
        )
      } else { # with covariates
        fit <- GLAM_2d_covariates(
          R = R, Y = Y,
          Bu = Bu, Bs = Bs,
          Z = Z,
          Wprior = Wprior,
          P = P,
          ridge = ridge,
          control_algorithm = con
        )
      }

      if (con$monitor_ev) {
        print(paste0(
          "AIC for model with log_10(rho_u) =  ",
          lru[i],
          " and log_10(rho_s) = ",
          lrs[j],
          " : ", fit$aic
        ))
        print(paste0(
          "BIC for model with log_10(rho_u) =  ",
          lru[i],
          " and log_10(rho_s) = ",
          lrs[j],
          " : ", fit$bic
        ))
      }
      AIC[i, j] <- fit$aic
      BIC[i, j] <- fit$bic
    }
  }

  # ---- Optimal penalty ----
  which_opt_aic <- which(AIC == min(AIC), arr.ind = TRUE)
  optim_lr_aic <- c(lru[which_opt_aic[1]], lrs[which_opt_aic[2]])
  optim_r_aic <- 10^optim_lr_aic

  which_opt_bic <- which(BIC == min(BIC), arr.ind = TRUE)
  optim_lr_bic <- c(lru[which_opt_bic[1]], lrs[which_opt_bic[2]])
  optim_r_bic <- 10^optim_lr_bic

  if (gsp$plot_aic) {
    image.plot(lru, lrs, AIC,
      main = gsp$main_aic,
      col = gsp$col
    )
    if (gsp$plot_contour) {
      contour(lru, lrs, AIC,
        col = "white", add = TRUE
      )
    }
    if (gsp$mark_optimal) {
      points(optim_lr_aic[1], optim_lr_aic[2],
        col = "red",
        pch = 15
      )
    }
  }
  if (gsp$plot_bic) {
    image.plot(lru, lrs, BIC,
      main = gsp$main_bic,
      col = gsp$col
    )
    if (gsp$plot_contour) {
      contour(lru, lrs, BIC,
        col = "white", add = TRUE
      )
    }
    if (gsp$mark_optimal) {
      points(optim_lr_bic[1], optim_lr_bic[2],
        col = "red",
        pch = 15
      )
    }
  }

  # ---- Estimates optimal model ----
  if (optim_criterion == "aic") {
    optim_lr <- optim_lr_aic
    P <- optim_r_aic[1] * kronecker(Is, t(Du) %*% Du) +
      optim_r_aic[2] * kronecker(t(Ds) %*% Ds, Iu)
  } else {
    optim_lr <- optim_lr_bic
    P <- optim_r_bic[1] * kronecker(Is, t(Du) %*% Du) +
      optim_r_bic[2] * kronecker(t(Ds) %*% Ds, Iu)
  }

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
    mod <- GLAM_2d_covariates( # with covariates
      R = R, Y = Y,
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
    "P_optimal" = P,
    "AIC" = AIC,
    "BIC" = BIC
  )
  if (gsp$return_aic == FALSE) {
  results$AIC <- NULL
  }
  if (gsp$return_bic == FALSE) {
    results$BIC <- NULL
  }

  class(results) <- "haz2ts"
  return(results)

}
