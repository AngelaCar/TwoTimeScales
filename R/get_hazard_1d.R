#' Get estimated (log-)hazard curve with 1 time scale
#'
#' @description `get_hazard_1d()` takes as input the results of a model
#' estimated by `fit1ts` and it returns the estimated smooth log-hazard
#' and the smooth  hazard together with their standard errors.
#'
#' If the model includes covariates, then only the baseline (log-)hazard is returned.
#' It is possible to provide values that define a new grid for plotting.
#' If not specified, the plotting is done using the same B-splines basis as for
#' the estimation. The function will check if the parameters for the grid
#' provided by the user are compatible with those originally used to construct
#' the B-splines for estimating the model. If not, the grid will be adjusted
#' accordingly and a warning will be returned.
#'
#' @inheritParams plot.haz1ts
#' @param fitted_model is an object of class `"haz1ts"`, the output of the function `fit1ts()`.
#'
#' @return A list with the following elements:
#'   * `new_plot_grid` A list of new specifications of the grid for plotting.
#'   * `hazard` A vector containing the estimated hazard.
#'   * `loghazard` A vector containing the estimated log-hazard.
#'   * `log10hazard` A vector containing the estimated log10-hazard.
#'   * `SE_hazard` A vector containing the estimated SEs for the hazard.
#'   * `SE_loghazard` A vector containing the estimated SEs for the log-hazard.
#'   * `SE_log10hazard` A vector containing the estimated SEs for the log10-hazard
#' @export
#'
get_hazard_1d <- function(fitted_model, plot_grid = NULL) {
  Bbases <- fitted_model$optimal_model$Bbases

  if (!is.null(plot_grid) & length(plot_grid) != 3)
    stop ("Not enough elements provided in `plot_grid`.")
  if (is.null(plot_grid['smin'])){
    smin <- attributes(Bbases$Bs)$xl
    } else smin <- plot_grid['smin']

  if (is.null(plot_grid['smax'])) {
    smax <- attributes(Bbases$Bs)$xr
    } else smax <- plot_grid['smax']
  if (is.null(plot_grid['ds'])) {
    ds <- diff(attributes(Bbases$Bs)$x)[1]
    } else ds <- plot_grid['ds']


  if (smin < attributes(Bbases$Bs)$xl) {
    smin <- attributes(Bbases$Bs)$xl
    warning("`smin` is smaller than the lower limit of the domain of Bs. Left boundary adjusted to  =  ", attributes(Bbases$Bs)$xl)
  }
  if (smax > attributes(Bbases$Bs)$xr) {
    smax <- attributes(Bbases$Bs)$xr
    warning("`smax` is larger than the upper limit of the domain of Bs. Right boundary adjusted to  =  ", attributes(Bbases$Bs)$xr)
  }
  K <- ceiling((smax - smin) / ds)
  ints <- seq(smin, smin + K * ds, by = ds)

  # Evaluate old basis in new grid of points
  Bs <- JOPS::bbase(ints, nseg = attributes(Bbases$Bs)$nseg, bdeg = attributes(Bbases$Bs)$bdeg)
  new_grid <- list(
    "ints" = ints,
    "smin" = smin,
    "smax" = smax,
    "ds" = ds
  )

 # ---- Calculate (baseline) hazard ----
  eta <- Bs %*% fitted_model$optimal_model$alpha
  haz <- exp(eta)

  # ---- Calculate Standard Errors for the log-hazard ----
  se_eta <- Bs %*% fitted_model$optimal_model$SE_alpha

  # ---- Calculate the log10-hazard (baseline) ----
  log10haz <- log10(haz)

  # ---- Calculate Standard Errors for the hazard ----
  se_haz <- haz * se_eta

  # ---- Calculate Standard Errors for the log10-hazard ----
  const <- log(10)
  se_log10haz <- abs(1/(haz * const)) * se_haz

  # ---- Return results in a list ----
  results <- list(
    "new_plot_grid" = new_grid,
    "hazard" = haz,
    "loghazard" = eta,
    "log10hazard" = log10haz,
    "SE_hazard" = se_haz,
    "SE_loghazard" = se_eta,
    "SE_log10hazard" = se_log10haz
  )

  return(results)
}
