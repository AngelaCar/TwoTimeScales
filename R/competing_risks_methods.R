# Collection of functions for the competing risks model over two time scales

#' @title Cumulative hazard over two time scales
#'
#' @description Computes the cumulative hazard surface over two time scales
#'      from a fitted model. The function is also called internally from `plot()`
#'      if the user wants to plot the cumulative hazard from a fitted model.
#'
#' @param fitted_model (optional) The output of the function `fit2ts`.
#'   This is an object of class `'haz2ts'` or `'haz2tsLMM'`.
#' @param plot_grid (optional) A list containing the parameters to build a new
#'   finer grid of intervals over `u` and `s` for plotting. This must be of the
#'   form:
#'   * `plot_grid = list(c(umin, umax, du), c(smin, smax, ds))`
#'   where `umin`, `umax` and `smin`, `smax` are the minimum and maximum values
#'   desired for the grid-points over `u` and `s` respectively, and `du`, `ds` are
#'   distances between two adjacent points over `u` and `s` respectively.
#'   Specifying a new denser grid is used to evaluate the B-spline bases used for
#'   estimation on such grid and plot the estimated surfaces with a greater level
#'   of detail.
#'   If not specified, the plotting is done using the same B-splines bases as
#'   for the estimation. The function will check if the parameters for the grid
#'   provided by the user are compatible with those originally used to construct
#'   the B-splines for estimating the model. If not, the grid will be adjusted
#'   accordingly and a warning will be returned.
#' @param cause a character string with a short name for the cause (optional).
#' @param midpoints A Boolean. Default is `FALSE`. If `TRUE`, the estimated
#'        quantities are evaluated at the midpoints of the rectangles
#'        (or parallelograms) of the grids, rather than at each grid-point.
#' @param where_slices A vector of values for the cutting points of the desired
#'   slices of the surface. This option is included mostly for the plotting function.
#'   When using `plot.haz2ts()`, the user selects `which_plot = "cumhaz"` and
#'   `cumhaz_slices = TRUE`, then `where_slices` indicates the location of the
#'   cutting points over the `u` time.
#' @param direction If cross-sectional one-dimensional curves are plotted, this
#'    indicates whether the cutting points are located on the `u` time, or on the
#'    `s` time. For plots of the cumulative hazards, only cutting points over the
#'    `u` time are meaningful.
#' @param tmax The maximum value of `t` that should be plotted.
#'
#' @return A list with the following elements:
#'          * `Haz` a list of estimated hazard and associated SEs
#'           (obtained from the function `get_hazard_2d`);
#'          * `CumHaz` the cumulated hazard estimate over `u` and `s`;
#'          * `cause` (if provided) the short name for the cause.
#'
#' @examples
#' # Create some fake data - the bare minimum
#' id <- 1:20
#' u <- c(5.43, 3.25, 8.15, 5.53, 7.28, 6.61, 5.91, 4.94, 4.25, 3.86, 4.05, 6.86,
#'        4.94, 4.46, 2.14, 7.56, 5.55, 7.60, 6.46, 4.96)
#' s <- c(0.44, 4.89, 0.92, 1.81, 2.02, 1.55, 3.16, 6.36, 0.66, 2.02, 1.22, 3.96,
#'        7.07, 2.91, 3.38, 2.36, 1.74, 0.06, 5.76, 3.00)
#' ev <- c(1, 0, 0, 1, 0, 1, 0, 1, 1, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 1)#'
#'
#' fakedata <- as.data.frame(cbind(id, u, s, ev))
#' fakedata2ts <- prepare_data(u = fakedata$u,
#'                             s_out = fakedata$s,
#'                             ev = fakedata$ev,
#'                             ds = .5)
#' # Fit a fake model - not optimal smoothing
#' fakemod <- fit2ts(fakedata2ts,
#'                   optim_method = "grid_search",
#'                   lrho = list(seq(1 ,1.5 ,.5),
#'                               seq(1 ,1.5 ,.5)))
#'
#' # Obtain the fake cumulated hazard
#' fakecumhaz2ts <- cumhaz2ts(fakemod)
#'
#' @export
#'
cumhaz2ts <- function(fitted_model,
                      plot_grid = NULL,
                      cause = NULL,
                      midpoints = FALSE,
                      where_slices = NULL,
                      direction = c("u", "s", NULL),
                      tmax = NULL) {

  if (inherits(fitted_model, "haz2ts")) {
    Haz <- get_hazard_2d(
      fitted_model = fitted_model,
      plot_grid = plot_grid,
      midpoints = midpoints,
      where_slices = where_slices,
      direction = direction,
      tmax = tmax
    )
  } else {
    if (inherits(fitted_model, "haz2tsLMM")) {
      Haz <- get_hazard_2d_LMM(
        fitted_model = fitted_model,
        plot_grid = plot_grid,
         where_slices = where_slices,
        direction = direction,
        tmax = tmax
      )
    } else {
      stop("'x' must be either a 'haz2ts' object or a 'haz2tsLMM' object")
    }
  }
  # Cumulative Hazards
  ds <- Haz$new_plot_grid$ds
  CumHaz <- t(apply(Haz$hazard  * ds, 1, cumsum))

  res <- list(
    "Haz" = Haz,
    "CumHaz" = CumHaz
  )
  if (!is.null(cause)) {
    res$cause <- cause
  }

  class(res) <- "cumhaz2ts"
  return(res)
}

#' Survival function with two time scales
#'
#' @description
#' Computes the survival matrix, that contains the probability of not
#' experiencing an event (of any cause) by time `s` and fixed entry time `u`.
#' The survival function can be obtained from one fitted model with only one
#' event type, or combining information from several cause-specific hazard
#' in a competing risks model. In the first case, a fitted object of class `'haz2ts'`
#' or `'haz2tsLMM'` can be passed directly as argument to the function. In the
#' competing risks framework, the user should provide a list of cause-specific
#' cumulative hazard matrices. The function is also called internally from `plot()`
#' if the user wants to plot the cumulative hazard from a fitted model.
#'
#'
#' @param cumhaz (optional) a list with all the cause-specific cumulated hazard
#'  matrices (minimum one element needs to be supplied).
#'  If more than one cause-specific cumulated hazard is provided,
#'  then they should all be matrices of the same dimension.
#' @param fitted_model (optional) The output of the function `fit2ts`.
#'  This is an object of class `'haz2ts'` or `'haz2tsLMM'`.
#' @param plot_grid (optional) A list containing the parameters to build a new
#'   finer grid of intervals over `u` and `s` for plotting. This must be of the
#'   form:
#'   * `plot_grid = list(c(umin, umax, du), c(smin, smax, ds))`
#'   where `umin`, `umax` and `smin`, `smax` are the minimum and maximum values
#'   desired for the grid-points over `u` and `s` respectively, and `du`, `ds` are
#'   distances between points over `u` and `s` respectively. Specifying a new
#'   denser grid is used to evaluate the B-spline bases used for estimation on
#'   such grid and plot the estimated surfaces with a greater level of details.
#'   If not specified, the plotting is done using the same B-splines bases as
#'   for the estimation. The function will check if the parameters for the grid
#'   provided by the user are compatible with those originally used to construct
#'   the B-splines for estimating the model. If not, the grid will be adjusted
#'   accordingly and a warning will be returned.
#' @param cause a character string with a short name for the cause (optional).
#' @param midpoints A Boolean. Default is `FALSE`. If `TRUE`, the estimated
#'        quantities are evaluated at the midpoints of the rectangles
#'        (or parallelograms) of the grids, rather than at each grid-point.
#' @param where_slices A vector of values for the cutting points of the desired
#'   slices of the surface. This option is included mostly for the plotting function.
#'   When using `plot.haz2ts()`, the user selects `which_plot = "survival"` and
#'   `surv_slices = TRUE`, then `where_slices` indicates the location of the
#'   cutting points over the `u` time.
#' @param direction If cross-sectional one-dimensional curves are plotted, this
#'    indicates whether the cutting points are located on the `u` time, or on the
#'    `s` time. For plots of the survival function, only cutting points over the
#'    `u` time are meaningful.
#' @param tmax The maximum value of `t` that should be plotted.
#' @return a matrix containing the values of the survival function over `s` and `u`.
#' @examples
#' # Create some fake data - the bare minimum
#' id <- 1:20
#' u <- c(5.43, 3.25, 8.15, 5.53, 7.28, 6.61, 5.91, 4.94, 4.25, 3.86, 4.05, 6.86,
#'        4.94, 4.46, 2.14, 7.56, 5.55, 7.60, 6.46, 4.96)
#' s <- c(0.44, 4.89, 0.92, 1.81, 2.02, 1.55, 3.16, 6.36, 0.66, 2.02, 1.22, 3.96,
#'        7.07, 2.91, 3.38, 2.36, 1.74, 0.06, 5.76, 3.00)
#' ev <- c(1, 0, 0, 1, 0, 1, 0, 1, 1, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 1)
#'
#'
#' fakedata <- as.data.frame(cbind(id, u, s, ev))
#' fakedata2ts <- prepare_data(u = fakedata$u,
#'                             s_out = fakedata$s,
#'                             ev = fakedata$ev,
#'                             ds = .5)
#' # Fit a fake model - not optimal smoothing
#' fakemod <- fit2ts(fakedata2ts,
#'                   optim_method = "grid_search",
#'                   lrho = list(seq(1 , 1.5, .5),
#'                               seq(1 , 1.5, .5)))
#'
#' # Obtain the fake cumulated hazard
#' fakecumhaz2ts <- cumhaz2ts(fakemod)
#' # Fake survival curve
#' fakesurv2ts <- surv2ts(fitted_model = fakemod)
#'
#' @export
#'
surv2ts <- function(cumhaz = NULL,
                    fitted_model = NULL,
                    plot_grid = NULL,
                    cause = NULL,
                    midpoints = FALSE,
                    where_slices = NULL,
                    direction = c("u", "s", NULL),
                    tmax = NULL) {
  if(is.null(cumhaz) & is.null(fitted_model))
    stop("please provide either a list of cumulative hazards (cumhaz) or a fitted model")
  if(!is.null(cumhaz)){
    # check if cumhaz is a list
    if(!is.list(cumhaz)) stop("cumhaz should be a list with at least one element")
    ncauses <- length(cumhaz)
    Surv2ts <- vector("list", length = ncauses)

    for (i in 1:ncauses) {
      Surv2ts[[i]] <- exp(-cumhaz[[i]])
    }

    if (ncauses == 1) {
      names(Surv2ts) <- "Surv2ts"
    } else {
      Surv2ts$Surv2ts <- exp(-(Reduce("+", cumhaz)))
    }
  } else {
    if(!is.null(fitted_model)){
      if (inherits(fitted_model, "haz2ts")) {
        Haz <- get_hazard_2d(
          fitted_model = fitted_model,
          plot_grid = plot_grid,
          midpoints = midpoints,
          where_slices = where_slices,
          direction = direction,
          tmax = tmax
        )
      } else {
        if (inherits(fitted_model, "haz2tsLMM")) {
          Haz <- get_hazard_2d_LMM(
            fitted_model = fitted_model,
            plot_grid = plot_grid,
            where_slices = where_slices,
            direction = direction,
            tmax = tmax
          )
        } else {
          stop("'x' must be either a 'haz2ts' object or a 'haz2tsLMM' object")
        }
      }
      # Cumulative Hazards
      ds <- Haz$new_plot_grid$ds
      CumHaz <- t(apply(Haz$hazard, 1, cumsum) * ds)

      Surv2ts <- list("vector", length = 1) # this is for compatibility
      Surv2ts$Surv2ts <- exp(-CumHaz)
      attr(Surv2ts, "plot_grid") <- Haz$new_plot_grid
    }
  }

  return(Surv2ts)
}

#' Cumulative incidence surface over two time scales
#'
#' @param haz a list of cause-specific hazards
#' @param ds the distance between two consecutive intervals over the `s` time scale.
#'  This has to be equal for all cause-specific hazards
#' @param cause is an optional vector of short names for the causes. It should
#' be of the same length as the number of cause-specific cumulated hazards provided.
#'
#' @return a list with one cumulative incidence matrix for each cause-specific
#'  hazard (named if a vector of short names is passed to `cause`).
#'
#' @examples
#' # Create some fake data - the bare minimum
#' id <- 1:20
#' u <- c(5.43, 3.25, 8.15, 5.53, 7.28, 6.61, 5.91, 4.94, 4.25, 3.86, 4.05, 6.86,
#'        4.94, 4.46, 2.14, 7.56, 5.55, 7.60, 6.46, 4.96)
#' s <- c(0.44, 4.89, 0.92, 1.81, 2.02, 1.55, 3.16, 6.36, 0.66, 2.02, 1.22, 3.96,
#'        7.07, 2.91, 3.38, 2.36, 1.74, 0.06, 5.76, 3.00)
#' ev <- c(1, 0, 0, 1, 0, 1, 0, 1, 1, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 1)#'
#'
#' fakedata <- as.data.frame(cbind(id, u, s, ev))
#' fakedata2ts <- prepare_data(u = fakedata$u,
#'                             s_out = fakedata$s,
#'                             ev = fakedata$ev,
#'                             ds = .5)
#' # Fit a fake model - not optimal smoothing
#' fakemod <- fit2ts(fakedata2ts,
#'                   optim_method = "grid_search",
#'                   lrho = list(seq(1 ,1.5 ,.5),
#'                               seq(1 ,1.5 ,.5)))
#'
#' # Obtain the fake cumulated hazard
#' fakecumhaz2ts <- cumhaz2ts(fakemod)
#' # Fake cumulative incidence function 2ts
#' fakecif2ts <- cuminc2ts(haz = list(fakecumhaz2ts$Haz$hazard),
#'                         ds = .5)
#' @export
#'

# ----- Add option for name of cause
cuminc2ts <- function(haz = list(),
                      ds,
                      cause = NULL) {
  ncauses <- length(haz)
  if (!is.null(cause)) {
    if (length(cause) != ncauses) {
      message("The number of names provided for the causes is not equal to the number of causes provided.")
    }
  }

  # First, calculate the cumulative hazards
  CumHaz <- vector("list", length = ncauses)
  for (i in 1:ncauses){
    CumHaz[[i]] <- t(apply(haz[[i]]  * ds , 1, cumsum))
  }
  # overall survival
  surv <- surv2ts(CumHaz)$Surv2ts
  CIF2ts <- vector("list", length = ncauses)
  for (i in 1:ncauses) {
    CIF2ts[[i]] <- t(apply((haz[[i]] * surv) * ds, 1, cumsum))
    if (!is.null(cause)) {
      names(CIF2ts)[i] <- cause[i]
    }
  }

  class(CIF2ts) <- "cif2ts"
  return(CIF2ts)
}



