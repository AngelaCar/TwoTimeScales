#' Get estimated (log-)hazard values with 1 time scale
#'
#' @description `get_hazard_1d()` takes as input the results of a model
#' estimated by `fit1ts` and it returns the estimated values of the smooth log-hazard
#' and the smooth hazard together with their standard errors.
#'
#' If the model includes covariates, then only the baseline (log-)hazard is returned.
#' It is possible to provide values that define a new grid for evaluation of the
#' estimated hazard.
#' If not specified, the hazard is evaluated on the same grid used for the
#' binning of the data, and therefore the estimation of the model.
#' The function will check if the parameters for the new grid
#' provided by the user are compatible with those originally used to construct
#' the B-splines for estimating the model. If not, the grid will be adjusted
#' accordingly and a warning will be returned.
#'
#' @inheritParams plot.haz1ts
#' @param fitted_model is an object of class `"haz1ts"`, the output of the function `fit1ts()`.
#'
#' @return A list with the following elements:
#'   * `new_plot_grid` A list of parameters that specify the new grid, of the form
#'   list("ints", "smin", "smax", "ds")
#'   * `hazard` A vector containing the estimated hazard values.
#'   * `loghazard` A vector containing the estimated log-hazard values.
#'   * `log10hazard` A vector containing the estimated log10-hazard values.
#'   * `SE_hazard` A vector containing the estimated SEs for the hazard.
#'   * `SE_loghazard` A vector containing the estimated SEs for the log-hazard.
#'   * `SE_log10hazard` A vector containing the estimated SEs for the log10-hazard.
#' @export
#'
#' @examples
#' ## preparing data - no covariates
#' dt1ts <- prepare_data(
#'   data = reccolon2ts,
#'   s_in = "entrys",
#'   s_out = "timesr",
#'   events = "status",
#'   ds = 180
#' )
#'
#' ## fitting the model with fit1ts() - default options
#'
#' mod1 <- fit1ts(dt1ts)
#' # Obtain 1d hazard
#' get_hazard_1d(mod1)
#' # Change grid
#' get_hazard_1d(mod1,
#'   plot_grid = c(smin = 0, smax = 2730, ds = 30)
#' )
#'
get_hazard_1d <- function(fitted_model, plot_grid = NULL) {
  Bbases <- fitted_model$optimal_model$Bbases

  if (is.null(plot_grid)) {
    Bs <- Bbases$Bs
    new_grid <- list(
      "ints" = attributes(Bs)$x,
      "smin" = attributes(Bs)$xl,
      "smax" = attributes(Bs)$xr
    )
    new_grid$ds <- new_grid$ints[2] - new_grid$ints[1]
  } else {
    if (!is.null(plot_grid) & length(plot_grid) != 3) {
      stop("Not enough elements provided in `plot_grid`.")
    } else {
      smin <- plot_grid["smin"]
      smax <- plot_grid["smax"]
      ds <- plot_grid["ds"]
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
    }
  }


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
  se_log10haz <- abs(1 / (haz * const)) * se_haz

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
