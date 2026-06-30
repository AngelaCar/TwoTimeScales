#' Grid search for the optimal P-GAM model with two time scales
#'
#' @description `grid_search_pgam()` performs a grid search for the minimum
#'   AIC or BIC of the log-additive model with two time scales.
#'
#'   It finds the optimal values of `log_10(rho_u)` and `log_10(rho_s)` and
#'   returns the estimated optimal model.
#'
#' @param lru A vector of `log_10(rho_u)` values.
#' @param lrs A vector of `log_10(rho_s)` values.
#' @inheritParams get_aic_fit_pgam
#' @param optim_criterion The criterion to be used for optimization:
#' `"aic"` (default) or `"bic"`. BIC penalizes model complexity more strongly
#' than AIC, so that its usage is recommended when a smoother fit is preferable
#' (see also Camarda, 2012).
#'
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
#' @return An object of class `haz2tsPGAM` with the following elements:
#'   * `optimal_model` A list containing the results of the optimal model.
#'   * `optimal_logrho` The optimal couple of `log_10(rho_u)` and `log_10(rho_s)`
#'      values.
#'   * `P_optimal` The optimal penalty matrix P.
#'   * `AIC` (if `par_gridsearch$return_aic == TRUE`) The vector of AIC values.
#'   * `BIC` (if `par_gridsearch$return_bic == TRUE`) The vector of BIC values.
#'
#' @references Camarda, C. G. (2012). "MortalitySmooth: An R Package for
#'    Smoothing Poisson Counts with P-Splines."
#'    Journal of Statistical Software, 50(1), 1–24.
#'    https://doi.org/10.18637/jss.v050.i01
#'
#' @importFrom grDevices gray.colors
#' @importFrom graphics contour points
#' @keywords internal

grid_search_pgam <- function(lru, lrs,
                             r, y,
                             B, nb, nbu, nbs,
                             Z = NULL,
                             DutDu, DstDs,
                             ridge = 1e-6,
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
      Pu <- ru * DutDu
      Ps <- rs * DstDs
      P <- matrix(0, nb, nb)
      P[1:nbu, 1:nbu] <- Pu
      P[(nbu + 1):(nbs + nbu), (nbu + 1):(nbs + nbu)] <- Ps
      PpR <- P + ridge * diag(nb)

      if (is.null(Z)) {
        fit <- iwls_1d(r, y,
          Bs = B, P = PpR,
          control_algorithm = list(
            maxiter = con$maxiter,
            conv_crit = con$conv_crit,
            verbose = con$verbose
          )
        )
      } else {
        fit <- GLAM_1d_covariates(
          R = r,
          Y = y,
          Z = Z,
          ns = length(r),
          n = nrow(Z),
          Bs = B,
          P = PpR,
          Wprior = NULL,
          control_algorithm = list(
            maxiter = con$maxiter,
            conv_crit = con$conv_crit,
            verbose = con$verbose
          )
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
    Pu <- optim_r_aic[1] * DutDu
    Ps <- optim_r_aic[2] * DstDs
    P <- matrix(0, nb, nb)
    P[1:nbu, 1:nbu] <- Pu
    P[(nbu + 1):(nbs + nbu), (nbu + 1):(nbs + nbu)] <- Ps
    PpR <- P + ridge * diag(nb)

    if (is.null(Z)) {
      mod <- iwls_1d(r, y,
        Bs = B, P = PpR,
        control_algorithm = list(
          maxiter = con$maxiter,
          conv_crit = con$conv_crit,
          verbose = con$verbose
        )
      )
    } else {
      mod <- GLAM_1d_covariates(
        R = r,
        Y = y,
        Z = Z,
        ns = length(r),
        n = nrow(Z),
        Bs = B,
        P = PpR,
        Wprior = NULL,
        control_algorithm = list(
          maxiter = con$maxiter,
          conv_crit = con$conv_crit,
          verbose = con$verbose
        )
      )
    }
  } else {
    optim_lr <- optim_lr_bic
    Pu <- optim_r_bic[1] * DutDu
    Ps <- optim_r_bic[2] * DstDs
    P <- matrix(0, nb, nb)
    P[1:nbu, 1:nbu] <- Pu
    P[(nbu + 1):(nbs + nbu), (nbu + 1):(nbs + nbu)] <- Ps
    PpR <- P + ridge * diag(nb)

    if (is.null(Z)) {
      mod <- iwls_1d(r, y,
        Bs = B, P = PpR,
        control_algorithm = list(
          maxiter = con$maxiter,
          conv_crit = con$conv_crit,
          verbose = con$verbose
        )
      )
    } else {
      mod <- GLAM_1d_covariates(
        R = r,
        Y = y,
        Z = Z,
        ns = length(r),
        n = nrow(Z),
        Bs = B,
        P = PpR,
        Wprior = NULL,
        control_algorithm = list(
          maxiter = con$maxiter,
          conv_crit = con$conv_crit,
          verbose = con$verbose
        )
      )
    }

  }
  mu <- r * exp(B %*% mod$alpha)
  w <- c(mu)

  BWB <- t(B) %*% (w * B)
  BWBpP <- BWB + PpR
  K <- solve(BWBpP) %*% BWB
  K11 <- K[1:nbu, 1:nbu]
  K22 <- K[(nbu + 1):(nbu + nbs), (nbu + 1):(nbu + nbs)]
  EDu <- sum(diag(K11))
  EDs <- sum(diag(K22))

  # ---- save results in a list ----

  results <- list(
    "optimal_model" = mod,
    "optimal_logrho" = optim_lr,
    "P_optimal" = P,
    "EDu" = EDu,
    "EDs" = EDs,
    "AIC" = AIC,
    "BIC" = BIC
  )
  if (gsp$return_aic == FALSE) {
    results$AIC <- NULL
  }
  if (gsp$return_bic == FALSE) {
    results$BIC <- NULL
  }

  class(results) <- "haz2tsPGAM"
  return(results)
}

# ------------------------------------------------------------------------------
#' Numerical optimization of the P-GAM model with two time scales
#'
#' @description `fitpgammodel_ucminf()` performs a numerical optimization of the
#'   AIC or BIC of log-additive model with two time scales.
#'
#'   It finds the optimal values of `log_10(rho_u)` and `log_10(rho_s)`
#'   and returns the estimated optimal model.
#'   See also [ucminf::ucminf()].
#'
#' @inheritParams get_aic_fit_2d
#' @param optim_criterion The criterion to be used for optimization:
#' `"aic"` (default) or `"bic"`.
#'
#' @return An object of class `haz2tsPGAM` with the following elements:
#'   * `optimal_model` A list containing the results of the optimal model.
#'   * `optimal_logrho` A vector with the optimal values of `log10(rho_u)` and
#'     `log10(rho_s)`.
#'   * `P_optimal` The optimal penalty matrix P.
#' @references Nielsen H, Mortensen S (2024).
#' _ucminf: General-Purpose Unconstrained Non-Linear Optimization_.
#' R package version 1.2.2, <https://CRAN.R-project.org/package=ucminf>
#'
#' @import ucminf
#' @keywords internal

fitpgammodel_ucminf <- function(y, r,
                                Z = NULL,
                                optim_criterion = c("aic", "bic"),
                                lrho = c(0, 0),
                                B, nb, nbu, nbs,
                                DutDu, DstDs,
                                ridge = 1e-6,
                                control_algorithm = list()) {
  # ---- Controls for iterative process ----
  con <- list(
    maxiter = 20,
    conv_crit = 1e-5,
    verbose = FALSE,
    monitor_ev = FALSE,
    xtol = 1e-5
  )
  Ncon <- names(con)
  namesCon <- names(control_algorithm)

  con[namesCon] <- control_algorithm
  if (length(namesCon[!namesCon %in% Ncon]) > 0) {
    warning("Undefined entries in control! Default settings are used.\n")
    warning(
      "undefined keyword(s): ",
      paste(namesCon[!namesCon %in% Ncon], collapse = ", ")
    )
  }

  # ---- Find optimal smoothing parameters ----
  if (optim_criterion == "aic") {
    op <- ucminf::ucminf(
      par = lrho,
      fn = get_aic_fit_pgam,
      r = r, y = y, Z = Z,
      B = B, nb = nb, nbu = nbu, nbs = nbs,
      DutDu = DutDu, DstDs = DstDs,
      ridge = ridge,
      control_algorithm = con,
      control = list(xtol = con$xtol)
    )
  }
  if (optim_criterion == "bic") {
    op <- ucminf::ucminf(
      par = lrho,
      fn = get_bic_fit_pgam,
      r = r, y = y, Z = Z,
      B = B, nb = nb, nbu = nbu, nbs = nbs,
      DutDu = DutDu, DstDs = DstDs,
      ridge = ridge,
      control_algorithm = con,
      control = list(xtol = con$xtol)
    )
  }

  # ---- With optimal smoothing parameters, calculate optimal model ----
  optim_lr <- op$par
  optim_r <- 10^optim_lr

  # ---- Construct penalty matrix P  ----
  Pu <- optim_r[1] * DutDu
  Ps <- optim_r[2] * DstDs
  P <- matrix(0, nb, nb)
  P[1:nbu, 1:nbu] <- Pu
  P[(nbu + 1):(nbs + nbu), (nbu + 1):(nbs + nbu)] <- Ps
  PpR <- P + ridge * diag(nb)

  if (is.null(Z)) {
    mod <- iwls_1d(r, y,
      Bs = B, P = PpR,
      control_algorithm = list(
        maxiter = con$maxiter,
        conv_crit = con$conv_crit,
        verbose = con$verbose
      )
    )
  } else {
    mod <- GLAM_1d_covariates(
      R = r,
      Y = y,
      Z = Z,
      ns = length(r),
      n = nrow(Z),
      Bs = B,
      P = PpR,
      Wprior = NULL,
      control_algorithm = list(
        maxiter = con$maxiter,
        conv_crit = con$conv_crit,
        verbose = con$verbose
      )
    )
  }

  mu <- r * exp(B %*% mod$alpha)
  w <- c(mu)

  BWB <- t(B) %*% (w * B)
  BWBpP <- BWB + PpR
  K <- solve(BWBpP) %*% BWB
  K11 <- K[1:nbu, 1:nbu]
  K22 <- K[(nbu + 1):(nbu + nbs), (nbu + 1):(nbu + nbs)]
  EDu <- sum(diag(K11))
  EDs <- sum(diag(K22))

  # ---- save results in a list ----
  results <- list(
    "optimal_model" = mod,
    "optimal_logrho" = optim_lr,
    "P_optim" = P,
    "EDu" = EDu,
    "EDs" = EDs
  )

  class(results) <- "haz2tsPGAM"

  return(results)
}


#------------------------------------------------------------------------------
#' Return the AIC of P-GAM model
#'
#' @description `get_aic_fit_pgam()` fits the log-additive model with two time scales
#' with or without individual level covariates and it returns the AIC of the model.
#' See also `fitpgammodel_ucminf()` and `fitpgam()`.
#'
#' @param lrho A vector of two elements, the initial values for \eqn{\log_{10}(\varrho_u)}
#'   and \eqn{\log_{10}(\varrho_s)}.
#' @param r A vector of exposure times of length `nu` by `ns`.
#' @param y A vector of event counts of length `nu` by `ns`.
#' @param Z (optional) A regression matrix of covariates values of dimensions
#'   n by p.
#' @param B A matrix of B-splines.
#' @param nb The number of B-splines.
#' @param nbu The number of B-splines for the `u` dimension.
#' @param nbs The number of B-splines for the `s` dimension.
#' @param DutDu \eqn{D_u^TD_u} to multiply \eqn{\varrho_u}.
#' @param DstDs \eqn{D_s^TD_s} to multiply \eqn{\varrho_s}.
#' @param ridge A ridge penalty parameter: default is 0. This is useful when, in
#'    some cases the algorithm shows convergence problems. In this case, set to a small
#'    number, for example `1e-6`.
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
#' @keywords internal
#'
get_aic_fit_pgam <- function(lrho,
                             r, y,
                             Z = NULL,
                             B, nb, nbu, nbs,
                             DutDu, DstDs,
                             ridge = 1e-6,
                             control_algorithm = list(
                               maxiter = 20,
                               conv_crit = 1e-5,
                               verbose = FALSE,
                               monitor_ev = FALSE
                             )) {
  # ---- Construct the penalty matrix ----

  ru <- 10^lrho[1]
  rs <- 10^lrho[2]
  Pu <- ru * DutDu
  Ps <- rs * DstDs
  P <- matrix(0, nb, nb)
  P[1:nbu, 1:nbu] <- Pu
  P[(nbu + 1):(nbs + nbu), (nbu + 1):(nbs + nbu)] <- Ps
  PpR <- P + ridge * diag(nb)

  if (is.null(Z)) {
    fit <- iwls_1d(r, y,
      Bs = B, P = PpR,
      control_algorithm = list(
        maxiter = control_algorithm$maxiter,
        conv_crit = control_algorithm$conv_crit,
        verbose = control_algorithm$verbose
      )
    )
  } else {
    fit <- GLAM_1d_covariates(
      R = r,
      Y = y,
      Z = Z,
      ns = length(r),
      n = nrow(Z),
      Bs = B,
      P = PpR,
      Wprior = NULL,
      control_algorithm = list(
        maxiter = control_algorithm$maxiter,
        conv_crit = control_algorithm$conv_crit,
        verbose = control_algorithm$verbose
      )
    )
  }

  if (control_algorithm$monitor_ev) {
    print(paste0(
      "AIC for model with log_10(rho_u) =  ",
      lrho[1],
      " and log_10(rho_s) = ",
      lrho[2],
      " : ", fit$aic
    ))
  }

  # ---- Extract the aic and return it as result ----
  aic <- fit$aic
  return(aic)
}


#' Return the BIC of P-GAM model
#'
#' @description `get_bic_fit_pgam()` fits the log-additive model with two time scales
#' with or without individual level covariates and it returns the BIC of the model.
#' See also `fitpgammodel_ucminf()` and `fitpgam()`.
#'
#' @param lrho A vector of two elements, the initial values for \eqn{\log_{10}(\varrho_u)}
#'   and \eqn{\log_{10}(\varrho_s)}.
#' @param r A vector of exposure times of length `nu` by `ns`.
#' @param y A vector of event counts of length `nu` by `ns`.
#' @param Z (optional) A regression matrix of covariates values of dimensions
#'   n by p.
#' @param B A matrix of B-splines.
#' @param nb The number of B-splines.
#' @param nbu The number of B-splines for the `u` dimension.
#' @param nbs The number of B-splines for the `s` dimension.
#' @param DutDu \eqn{D_u^TD_u} to multiply \eqn{\varrho_u}.
#' @param DstDs \eqn{D_s^TD_s} to multiply \eqn{\varrho_s}.
#' @param ridge A ridge penalty parameter: default is 0. This is useful when, in
#'    some cases the algorithm shows convergence problems. In this case, set to a small
#'    number, for example `1e-6`.
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
#' @return The `bic` value of the fitted model.
#' @keywords internal
#'
get_bic_fit_pgam <- function(lrho,
                             r, y,
                             Z = NULL,
                             B, nb, nbu, nbs,
                             DutDu, DstDs,
                             ridge = 1e-6,
                             control_algorithm = list(
                               maxiter = 20,
                               conv_crit = 1e-5,
                               verbose = FALSE,
                               monitor_ev = FALSE
                             )) {
  # ---- Construct the penalty matrix ----
  ru <- 10^lrho[1]
  rs <- 10^lrho[2]
  Pu <- ru * DutDu
  Ps <- rs * DstDs
  P <- matrix(0, nb, nb)
  P[1:nbu, 1:nbu] <- Pu
  P[(nbu + 1):(nbs + nbu), (nbu + 1):(nbs + nbu)] <- Ps
  PpR <- P + ridge * diag(nb)

  if (is.null(Z)) {
    fit <- iwls_1d(r, y,
      Bs = B, P = PpR,
      control_algorithm = list(
        maxiter = control_algorithm$maxiter,
        conv_crit = control_algorithm$conv_crit,
        verbose = control_algorithm$verbose
      )
    )
  } else {
    fit <- GLAM_1d_covariates(
      R = r,
      Y = y,
      Z = Z,
      ns = length(r),
      n = nrow(Z),
      Bs = B,
      P = PpR,
      Wprior = NULL,
      control_algorithm = list(
        maxiter = control_algorithm$maxiter,
        conv_crit = control_algorithm$conv_crit,
        verbose = control_algorithm$verbose
      )
    )
  }

  if (control_algorithm$monitor_ev) {
    print(paste0(
      "BIC for model with log_10(rho_u) =  ",
      lrho[1],
      " and log_10(rho_s) = ",
      lrho[2],
      " : ", fit$bic
    ))
  }

  # ---- Extract the aic and return it as result ----
  bic <- fit$bic
  return(bic)
}
