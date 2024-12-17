#' Return the AIC of 2ts model
#'
#' @description `get_aic_fit_2d()` fits the 2ts model with or without individual
#' level covariates and it returns the AIC of the model.
#' See also `fit2tsmodel_ucminf()` and `fit2ts()`.
#'
#' @param lrho A vector of two elements, the initial values for \eqn{\log_{10}(\varrho_u)}
#'   and \eqn{\log_{10}(\varrho_s)}.
#' @param R A matrix (or 3d-array) of exposure times of dimension nu by ns
#'   (or nu by ns by n).
#' @param Y A matrix (or 3d-array) of event counts of dimension nu by ns
#'   (or nu by ns by n).
#' @param Z (optional) A regression matrix of covariates values of dimensions
#'   n by p.
#' @param Bu A matrix of B-splines for the `u` time scale of dimension nu by cu.
#' @param Bs A matrix of B-splines for the `s` time scale of dimension ns by cs.
#' @param Iu An identity matrix of dimension nbu by nbu.
#' @param Is An identity matrix of dimension nbs by nbs.
#' @param Du The difference matrix over `u`.
#' @param Ds The difference matrix over `s`.
#' @param Wprior An optional matrix of a-priori weights.
#' @param ridge A ridge penalty parameter: default is 0. This is useful when, in
#'    some cases the algorithm shows convergence problems. In this case, set to a small
#'    number, for example `1e-4`.
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


get_aic_fit_2d <- function(lrho,
                           R, Y,
                           Z = NULL,
                           Bu, Bs,
                           Iu, Is,
                           Du, Ds,
                           Wprior = NULL,
                           ridge = 0,
                           control_algorithm = list(
                             maxiter = 20,
                             conv_crit = 1e-5,
                             verbose = FALSE,
                             monitor_ev = FALSE
                           )
                           ) {
  # ---- Construct the penalty matrix ----
  ru <- 10^lrho[1]
  rs <- 10^lrho[2]

  P <- ru * kronecker(Is, t(Du) %*% Du) + rs * kronecker(t(Ds) %*% Ds, Iu)

  # ---- Estimate the model with the correct function ----
  if (is.null(Z)) { # without covariates
    fit <- GLAM_2d_no_covariates(
      R = R, Y = Y,
      Bu = Bu, Bs = Bs,
      Wprior = Wprior,
      P = P,
      ridge = ridge,
      control_algorithm = control_algorithm
    )
  } else { # with covariates
    fit <- GLAM_2d_covariates(
      R = R, Y = Y,
      Bu = Bu, Bs = Bs,
      Z = Z,
      Wprior = Wprior,
      P = P,
      ridge = ridge,
      control_algorithm = control_algorithm
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
