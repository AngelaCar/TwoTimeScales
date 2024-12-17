#' Get estimated (log-)hazard values with 1 time scale
#'
#' @description `get_hazard_1d_LMM()` takes as input the results of a model
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
#' @param fitted_model is an object of class `"haz1tsLMM"`, the output of the function `fit1ts()`.
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
#' ## fitting the model with fit1ts()
#'
#' mod1 <- fit1ts(dt1ts,
#'   optim_method = "LMMsolver"
#' )
#' # Obtain 1d hazard
#' get_hazard_1d_LMM(mod1)
#' # Change grid
#' get_hazard_1d_LMM(mod1,
#'   plot_grid = c(smin = 0, smax = 2730, ds = 30)
#' )
#'
get_hazard_1d_LMM <- function(fitted_model, plot_grid = NULL) {
  if (is.null(plot_grid)) {
    ds <- attributes(fitted_model$optimal_model$splRes[[1]]$knots[[1]])$dx
    smin <- attributes(fitted_model$optimal_model$splRes[[1]]$knots[[1]])$xmin
    smax <- attributes(fitted_model$optimal_model$splRes[[1]]$knots[[1]])$xmax
    K <- ceiling((smax - smin) / ds)
    ints <- seq(smin, smin + K * ds, by = ds)
  } else {
    if (!is.null(plot_grid) & length(plot_grid) != 3) {
      stop("Not enough elements provided in `plot_grid`.")
    } else {
      if (!is.null(plot_grid) & length(plot_grid) != 3) {
        stop("Not enough elements provided in `plot_grid`.")
      } else {
        smin <- plot_grid["smin"]
        smax <- plot_grid["smax"]
        ds <- plot_grid["ds"]
        if (smin < attributes(fitted_model$optimal_model$splRes[[1]]$knots[[1]])$xmin) {
          smin <- attributes(fitted_model$optimal_model$splRes[[1]]$knots[[1]])$xmin
          warning(
            "`smin` is smaller than the lower limit of the domain of Bs. Left boundary adjusted to  =  ",
            attributes(fitted_model$optimal_model$splRes[[1]]$knots[[2]])$xmin
          )
        }
        if (smax > attributes(fitted_model$optimal_model$splRes[[1]]$knots[[1]])$xmax) {
          smax <- attributes(fitted_model$optimal_model$splRes[[1]]$knots[[1]])$xmax
          warning(
            "`smax` is larger than the upper limit of the domain of Bs. Right boundary adjusted to  =  ",
            attributes(fitted_model$optimal_model$splRes[[1]]$knots[[1]])$xmax
          )
        }
        K <- ceiling((smax - smin) / ds)
        ints <- seq(smin, smin + K * ds, by = ds)
      }
    }
  }


  new_grid <- data.frame("s" = ints)
  # ---- Get hazard ----
  trend1D <- obtainSmoothTrend(fitted_model$optimal_model,
    newdata = new_grid,
    includeIntercept = TRUE
  )
  haz <- trend1D$ypred
  se_haz <- trend1D$se
  eta <- log(haz)
  se_eta <- abs(1 / haz) * se_haz
  log10haz <- log10(haz)
  const <- log(10)
  se_log10haz <- abs(1 / (haz * const)) * se_haz

  new_grid <- list(
    "ints" = ints,
    "smin" = smin,
    "smax" = smax,
    "ds" = ds
  )

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
