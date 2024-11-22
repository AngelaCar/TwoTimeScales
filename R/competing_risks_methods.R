# Collection of functions for the competing risks model over two time scales

#' @title Cumulated hazard over two time scales
#'
#' @description Computes the cumulated hazard surface over two time scales
#'
#' @param fitted_model The output of the function `fit2ts`. This is an object of
#'   class `'haz2ts'`.
#' @param plot_grid (optional) A list containing the parameters to build a new
#'   finer grid of intervals over `u` and `s` for plotting. This must be of the
#'   form:
#'   * `plot_grid = list(c(umin, umax, du), c(smin, smax, ds))`
#'   where `umin`, `umax` and `smin`, `smax` are the minimum and maximum values
#'   desired for the intervals over `u` and `s` respectively, and `du`, `ds` are
#'   distances between intervals over `u` and `s` respectively. Specifying a new
#'   denser grid is used to evaluate the B-spline bases used for estimation on
#'   such grid and plot the estimated surfaces with a greater level of details.
#'   If not specified, the plotting is done using the same B-splines bases as
#'   for the estimation. The function will check if the parameters for the grid
#'   provided by the user are compatible with those originally used to construct
#'   the B-splines for estimating the model. If not, the grid will be adjusted
#'   accordingly and a warning will be returned.
#' @param cause a character string with a short name for the cause (optional).
#'
#' @return A list with the following elements:
#'          * `Haz` the 2d hazard surface (obtained from the function `get_hazard_2d`);
#'          * `CumHaz` the cumulated hazard estimate over `u` and `s`;
#'          * `cause` (if provided) the short name for the cause
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
                      cause = NULL) {
  if (inherits(fitted_model, "haz2ts")) {
    Haz <- get_hazard_2d(
      fitted_model = fitted_model,
      plot_grid = plot_grid
    )
  } else {
    if (inherits(fitted_model, "haz2tsLMM")) {
      Haz <- get_hazard_2d_LMM(
        fitted_model = fitted_model,
        plot_grid = plot_grid
      )
    } else {
      stop("'x' must be either a 'haz2ts' object or a 'haz2tsLMM' object")
    }
  }
  # Cumulative Hazards
  ds <- Haz$new_plot_grid$ds
  CumHaz <- t(apply(Haz$hazard, 1, cumsum) * ds)

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
#' experiencing an event of any cause by time `s` and fixed time at entry `u`,
#' from a list of cause-specific cumulated hazard matrices.
#'
#'
#' @param cumhaz a list with all the cause-specific cumulated hazard matrices
#'  (minimum one element needs to be supplied).
#'  If more than one cause-specific cumulated hazard is provided,
#'  then they should all be matrices of the same dimension.
#' @return a matrix containing the values of the survival function over `s` and `u`.
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
#' # Fake survival curve
#' fakesurv2ts <- surv2ts(cumhaz = list(fakecumhaz2ts$CumHaz))
#'
#' @export
#'
surv2ts <- function(cumhaz = list()) {
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
  return(Surv2ts)
}

#' Cumulative incidence surface over two time scales
#'
#' @param haz a list of cause-specific hazards
#' @param surv the overall survival probability surface over two time scales, obtained
#'  from `surv2ts`
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
#' # Fake survival curve
#' fakesurv2ts <- surv2ts(cumhaz = list(fakecumhaz2ts$CumHaz))
#' # Fake cumulative incidence function 2ts
#' fakecif2ts <- cuminc2ts(haz = list(fakecumhaz2ts$Haz$hazard),
#'                         surv = fakesurv2ts$Surv2ts,
#'                         ds = .5)
#' @export
#'

# ----- Add option for name of cause
cuminc2ts <- function(haz = list(),
                      surv,
                      ds,
                      cause = NULL) {
  ncauses <- length(haz)
  if (!is.null(cause)) {
    if (length(cause) != ncauses) {
      message("The number of names provided for the causes is not equal to the number of causes provided.")
    }
  }
  CIF2ts <- vector("list", length = ncauses)
  for (i in 1:ncauses) {
    CIF2ts[[i]] <- t(apply(haz[[i]] * surv, 1, cumsum) * ds)
    if (!is.null(cause)) {
      names(CIF2ts)[i] <- cause[i]
    }
  }

  class(CIF2ts) <- "cif2ts"
  return(CIF2ts)
}

