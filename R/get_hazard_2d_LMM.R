#' Get estimated (log-)hazard surface with 2 time scales
#'
#' @description `get_hazard_2d_LMM()` takes as input an object of class `'haz2tsLMM'`
#'  and it returns the estimated smooth log-hazard, the log10-hazard and the
#'  hazard surface together with their standard errors.
#'
#' It is possible to provide values that define a new grid for evaluation of the
#' estimated hazard.
#' If not specified, the hazard is evaluated on the same grid used for the
#' binning of the data, and therefore the estimation of the model.
#' The function will check if the parameters for the new grid
#' provided by the user are compatible with those originally used to construct
#' the B-splines for estimating the model. If not, the grid will be adjusted
#' accordingly and a warning will be returned.
#'
#' @inheritParams plot.haz2tsLMM
#' @param fitted_model is an object of class `'haz2tsLMM'`
#' the output of the function `fit2ts()`.
#' @param tmax The maximum value of `t` that should be plotted.
#'
#' @return A list with the following elements:
#'   * `new_plot_grid` A list of parameters that specify the new grid, of the form
#'   list("intu", "umin", "umax", "du", "ints", "smin", "smax", "ds")
#'   * `hazard` A matrix containing the estimated hazard values.
#'   * `loghazard` A matrix containing the estimated log-hazard values.
#'   * `log10hazard`A matrix containing the estimated log10-hazard values.
#'   * `SE_hazard` A matrix containing the estimated SEs for the hazard.
#'   * `SE_loghazard` A matrix containing the estimated SEs for the log-hazard.
#'   * `SE_log10haz` A matrix containing the estimated SEs for the log10-hazard.
#' @export
#'
#' @examples
#' # Create some fake data - the bare minimum
#' id <- 1:20
#' u <- c(
#'   5.43, 3.25, 8.15, 5.53, 7.28, 6.61, 5.91, 4.94, 4.25, 3.86, 4.05, 6.86,
#'   4.94, 4.46, 2.14, 7.56, 5.55, 7.60, 6.46, 4.96
#' )
#' s <- c(
#'   0.44, 4.89, 0.92, 1.81, 2.02, 1.55, 3.16, 6.36, 0.66, 2.02, 1.22, 3.96,
#'   7.07, 2.91, 3.38, 2.36, 1.74, 0.06, 5.76, 3.00
#' )
#' ev <- c(1, 0, 0, 1, 0, 1, 0, 1, 1, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 1) #'
#'
#' fakedata <- as.data.frame(cbind(id, u, s, ev))
#' fakedata2ts <- prepare_data(data = fakedata,
#'                             u = "u",
#'                             s_out = "s",
#'                             ev = "ev",
#'                             ds = .5)
#' # Fit a fake model - not optimal smoothing
#' fakemod <- fit2ts(fakedata2ts,
#'   optim_method = "LMMsolver"
#' )
#'
#' # Get hazard
#' get_hazard_2d_LMM(fakemod)
#'
#' # Use a finer grid of points
#' get_hazard_2d_LMM(fakemod,
#'                  plot_grid = list(c(umin = 3, umax = 8.5, du = .1),
#'                                   c(smin = 0, smax = 7.1, ds = .1)))
#'
get_hazard_2d_LMM <- function(fitted_model,
                              plot_grid = NULL,
                              where_slices = NULL,
                              direction = c("u", "s", NULL),
                              tmax = NULL) {
  # ---- Make grid ----
  # check if all information is provided
  if (is.null(plot_grid)) {
    du <- attributes(fitted_model$optimal_model$splRes[[1]]$knots[[1]])$dx
    ds <- attributes(fitted_model$optimal_model$splRes[[1]]$knots[[2]])$dx
    umin <- attributes(fitted_model$optimal_model$splRes[[1]]$knots[[1]])$xmin
    umax <- attributes(fitted_model$optimal_model$splRes[[1]]$knots[[1]])$xmax
    smin <- attributes(fitted_model$optimal_model$splRes[[1]]$knots[[2]])$xmin
    smax <- attributes(fitted_model$optimal_model$splRes[[1]]$knots[[2]])$xmax
    plot_grid <- list(
      c("umin" = umin, "umax" = umax, "du" = du),
      c("smin" = smin, "smax" = smax, "ds" = ds)
    )
  } else {
    if ((length(plot_grid[[1]]) != 3) | (length(plot_grid[[2]]) != 3)) {
      stop("One or more arguments from plot_grid are missing.")
    } else {
      umin <- plot_grid[[1]][1]
      umax <- plot_grid[[1]][2]
      du <- plot_grid[[1]][3]
      smin <- plot_grid[[2]][1]
      smax <- plot_grid[[2]][2]
      ds <- plot_grid[[2]][3]
      if (du <= 0) stop("`du` should be a positive number!")
      if (ds <= 0) stop("`ds` should be a positive number!")
    }
  }
  if (umin < attributes(fitted_model$optimal_model$splRes[[1]]$knots[[1]])$xmin) {
    umin <- attributes(fitted_model$optimal_model$splRes[[1]]$knots[[1]])$xmin
    warning(
      "`umin` is smaller than the lower limit of the domain of Bu. Left boundary adjusted to  =  ",
      attributes(fitted_model$optimal_model$splRes[[1]]$knots[[1]])$xmin
    )
  }
  if (umax > attributes(fitted_model$optimal_model$splRes[[1]]$knots[[1]])$xmax) {
    umax <- attributes(fitted_model$optimal_model$splRes[[1]]$knots[[1]])$xmax
    warning(
      "`umax` is larger than the upper limit of the domain of Bu. Right boundary adjusted to  =  ",
      attributes(fitted_model$optimal_model$splRes[[1]]$knots[[1]])$xmax
    )
  }
  K <- ceiling((umax - umin) / du)
  intu <- seq(umin, umin + K * du, by = du)

  if (smin < attributes(fitted_model$optimal_model$splRes[[1]]$knots[[2]])$xmin) {
    smin <- attributes(fitted_model$optimal_model$splRes[[1]]$knots[[2]])$xmin
    warning(
      "`smin` is smaller than the lower limit of the domain of Bs. Left boundary adjusted to  =  ",
      attributes(fitted_model$optimal_model$splRes[[1]]$knots[[2]])$xmin
    )
  }
  if (smax > attributes(fitted_model$optimal_model$splRes[[1]]$knots[[2]])$xmax) {
    smax <- attributes(fitted_model$optimal_model$splRes[[1]]$knots[[2]])$xmax
    warning(
      "`smax` is larger than the upper limit of the domain of Bs. Right boundary adjusted to  =  ",
      attributes(fitted_model$optimal_model$splRes[[1]]$knots[[2]])$xmax
    )
  }
  K <- ceiling((smax - smin) / ds)
  ints <- seq(smin, smin + K * ds, by = ds)

  new_grid <- list(
    "intu" = intu,
    "umin" = umin,
    "umax" = umax,
    "du" = du,
    "ints" = ints,
    "smin" = smin,
    "smax" = smax,
    "ds" = ds
  )


  # ---- Adjust grid and B-splines if slices are required ----
  if (!is.null(where_slices)) {
    if (is.null(direction)) stop("Direction for cutting slices missing.")
    if (direction == "u") {
      if (min(where_slices) < attributes(fitted_model$optimal_model$splRes[[1]]$knots[[1]])$xmin |
        max(where_slices) > attributes(fitted_model$optimal_model$splRes[[1]]$knots[[1]])$xmax) {
        stop("Desired cutting points outside of range of `B_u`.")
      } else {
        intu <- unique(sort(c(intu, where_slices)))
        new_grid$intu <- intu
      }
    }
    if (direction == "s") {
      if (min(where_slices) < attributes(fitted_model$optimal_model$splRes[[1]]$knots[[2]])$xmin |
        max(where_slices) > attributes(fitted_model$optimal_model$splRes[[1]]$knots[[2]])$xmax) {
        stop("Desired cutting points outside of range of `B_s`.")
      } else {
        news <- unique(sort(c(new_grid$ints, where_slices)))
        new_grid$ints <- news
        grid_us <- expand.grid(u = new_grid$intu, s = new_grid$ints)
        grid_us$t <- with(grid_us, u + s)
        t <- unique(grid_us$t)
        if (!is.null(tmax)) {
          new_grid$intt <- t[t <= tmax]
        } else {
          new_grid$intt <- t
        }
        grid_ts <- expand.grid(t = new_grid$intt, s = new_grid$ints)
        grid_ts$u <- grid_ts$t - grid_ts$s
        newu <- unique(grid_ts$u)
        newu <- newu[newu >= new_grid$umin & newu <= new_grid$umax]
        new_grid$intu <- sort(newu)
      }
    }
  }

  # ---- Get hazard ----
  gridpoints <- expand.grid("u" = new_grid$intu, "s" = new_grid$ints)
  trend2D <- obtainSmoothTrend(fitted_model$optimal_model,
    newdata = gridpoints, includeIntercept = TRUE
  )
  Haz <- matrix(trend2D$ypred, nrow = length(new_grid$intu), ncol = length(new_grid$ints))
  SE_Haz <- matrix(trend2D$se, nrow = length(new_grid$intu), ncol = length(new_grid$ints))
  logHaz <- log(Haz)
  log10Haz <- log10(Haz)
  SE_logHaz <- abs(1 / Haz) * SE_Haz
  SE_log10Haz <- abs(1 / (Haz * log(10))) * SE_Haz


  # ---- Return results in a list ----
  results <- list(
    "new_plot_grid" = new_grid,
    "hazard" = Haz,
    "loghazard" = logHaz,
    "log10hazard" = log10Haz,
    "SE_hazard" = SE_Haz,
    "SE_loghazard" = SE_logHaz,
    "SE_log10hazard" = SE_log10Haz
  )

  return(results)
}
