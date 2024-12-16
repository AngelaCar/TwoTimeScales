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
#' @examples
#' # Create some fake data - the bare minimum
#' id <- 1:20
#' u <- c(5.43, 3.25, 8.15, 5.53, 7.28, 6.61, 5.91, 4.94, 4.25, 3.86, 4.05, 6.86,
#'        4.94, 4.46, 2.14, 7.56, 5.55, 7.60, 6.46, 4.96)
#' s <- c(0.44, 4.89, 0.92, 1.81, 2.02, 1.55, 3.16, 6.36, 0.66, 2.02, 1.22, 3.96,
#'        7.07, 2.91, 3.38, 2.36, 1.74, 0.06, 5.76, 3.00)
#' ev <- c(1, 0, 0, 1, 0, 1, 0, 1, 1, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 1)
#' x1 <- c(0, 0, 0, 1, 0, 1, 1, 1, 1, 0, 0, 0, 1, 0, 1, 1, 0, 0, 1, 0)
#'
#' fakedata <- as.data.frame(cbind(id, u, s, ev, x1))
#' fakedata2ts <- prepare_data(data = fakedata,
#'                             u = "u",
#'                             s_out = "s",
#'                             ev = "ev",
#'                             ds = .5,
#'                             individual = TRUE,
#'                             covs = "x1")
#' # Fit a fake model - not optimal smoothing
#' fakemod <- fit2ts(fakedata2ts,
#'                   optim_method = "grid_search",
#'                   lrho = list(seq(1, 1.5, .5),
#'                               seq(1, 1.5, .5)))
#' get_hr(fakemod)


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
  return(results)
}
