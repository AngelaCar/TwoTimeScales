# Collection of functions for the competing risks model over two time scales

#' @title Cumulated hazard over two time scales
#'
#' @description Computes the cumulated hazard surface over two time scales
#' (currently only implemented for objects of class `haz2ts`, not `haz2tsLMM`).
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

#' Survival surface (cause-specific) and overall survival with two time scales
#'
#' @description
#' Computes cause-specific survival estimates over two time scales from a list of
#' cause-specific cumulated hazard matrices. Each matrix of cause-specific survival
#' is calculated as \eqn{\exp(-CumHaz)}. Additionally, if the user provides more than
#' one cause-specific cumulated hazard, the function computes the overall survival
#' matrix, that contains the probability of not experiencing an event of any cause
#' by time `s` and fixed time at entry `u`.
#'
#'
#' @param cumhaz a list with all the cause-specific cumulated hazard matrices
#'  (minimum one element needs to be supplied).
#'  If more than one cause-specific cumulated hazard is provided,
#'  then they should all be matrices of the same dimension.
#' @param cause is an optional vector of short names for the causes. It should
#' be of the same length as the number of cause-specific cumulated hazards provided.
#' @return a list with varying length. If only one cumulated hazard is provided
#'          as input, the function will calculate one two-dimensional survival
#'          surface. The output will be a list of one element.
#'          If K (with K>1) cumulated hazards are supplied, then the output will
#'          be a list with K+1 elements:
#'          * the first K elements will be the K cause-specific survival estimates;
#'            if a vector of short names with the causes is passed to the argument
#'            `cause`, then these will be used to name elements in the output;
#'          * `OverSurv` is the overall two-dimensional survival surface.
#' @export
#'
surv2ts <- function(cumhaz = list(),
                    cause = NULL) {
  ncauses <- length(cumhaz)
  if (!is.null(cause)) {
    if (length(cause) != ncauses) {
      message("The number of names provided for the causes is not equal to the number of causes provided.")
    }
  }
  Surv2ts <- vector("list", length = ncauses)

  for (i in 1:ncauses) {
    Surv2ts[[i]] <- exp(-cumhaz[[i]])
    if (!is.null(cause)) {
      names(Surv2ts)[i] <- cause[i]
    }
  }

  if (ncauses > 1) {
    Surv2ts$OverSurv <- exp(-(Reduce("+", cumhaz)))
  }
  return(Surv2ts)
}

#' Cumulative incidence surface over two time scales
#'
#' @param haz a list of cause-specific hazards
#' @param oversurv the overall survival probability surface over two time scales, obtained
#'  from `surv2ts` (the last element in the list of results)
#' @param ds the distance between two consecutive intervals over the `s` time scale.
#'  This has to be equal for all cause-specific hazards
#' @param cause is an optional vector of short names for the causes. It should
#' be of the same length as the number of cause-specific cumulated hazards provided.
#'
#' @return a list with one cumulative incidence matrix for each cause-specific
#'  hazard (named if a vector of short names is passed to `cause`).
#' @export
#'

# ----- Add option for name of cause
cuminc2ts <- function(haz = list(),
                      oversurv,
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
    CIF2ts[[i]] <- t(apply(haz[[i]] * oversurv, 1, cumsum) * ds)
    if (!is.null(cause)) {
      names(CIF2ts)[i] <- cause[i]
    }
  }

  class(CIF2ts) <- "cif2ts"
  return(CIF2ts)
}

