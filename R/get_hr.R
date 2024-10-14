#' Get the Hazard Ratios with their Standard Errors
#'
#' @description `get_hr()` takes as input the results of a model with covariates
#' estimated by `fit2ts` or `fit1ts` and returns the estimated hazard ratios
#' together with their standard errors.
#'
#' @param fitted_model A list returned by the function `fit2ts` or `fit1ts`.
#'
#' @return A list with the following elements:
#'   * `HR` A vector of hazard ratios (calculated as \eqn{\exp(\hat\beta)}).
#'   * `SE_HR` A vector of Standard Errors for the hazard ratios calculated
#'     via the delta method.
#'   * `beta` A vector of the estimated \eqn{\hat\beta} coefficients.
#'   * `SE_beta` A vector of the Standard Errors for the beta coefficients.
#'
#' @export

get_hr <- function(fitted_model){
  if(is.null(fitted_model$optimal_model$beta)) stop("This model does not have covariates effects.")

  HR <- exp(fitted_model$optimal_model$beta)
  SE_HR <- HR * fitted_model$optimal_model$SE_beta

  results <- list(
    "beta" = fitted_model$optimal_model$beta,
    "SE_beta" = fitted_model$optimal_model$SE_beta,
    "HR" = HR,
    "SE_HR" = SE_HR
  )
}
