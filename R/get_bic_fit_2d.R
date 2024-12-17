#' Return the BIC of 2ts model
#'
#' @description `get_bic_fit_2d()` fits the 2ts model with or without individual
#' level covariates and it returns the BIC of the model.
#' See also `fit2tsmodel_ucminf()` and `fit2ts()`.
#'
#'
#' @inheritParams get_aic_fit_2d
#'
#' @return The `bic` value of the fitted model.


get_bic_fit_2d <- function(lrho,
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
      "BIC for model with log_10(rho_u) =  ",
      lrho[1],
      " and log_10(rho_s) = ",
      lrho[2],
      " : ", fit$aic
    ))
  }

  # ---- Extract the bic and return it as result ----
  bic <- fit$bic
  return(bic)
}
