#' Bin data on two time scales
#'
#' @description `exposures_events_2d()` computes individual or aggregated
#'   matrices of exposure times and event counts starting from individual
#'   records of time at entry in the process (measured over the first time
#'   scale), duration at entry in the process (measured over the second time
#'   scale), duration at exit from the process (measured over the second time
#'   scale), and event's indicator.
#'
#' @details The fixed-time variable `u` and the second time scale `s`
#'   are divided into \eqn{nu} and \eqn{ns} intervals, respectively. The extremes of these
#'   intervals are provided as input to the function. First, the fixed-time at
#'   entry is located in one of the nu bins that cover the whole range of
#'   `u`. Then, the time-at-risk for each individual is split according to
#'   the \eqn{ns} bins that span the whole range of values for `s`, and an event
#'   indicator is placed in the bin where the exit time is located. This is done
#'   by calling the function `exposure_events_1d`. If individual matrices of
#'   exposure and events are required, then the function returns two arrays of
#'   dimension \eqn{nu} by \eqn{ns} by \eqn{n}. If aggregated results are preferred, the
#'   individual contributions are summed in each bin to provide a matrix of
#'   total exposure time and a matrix of total event counts, both of dimensions
#'   \eqn{nu} by \eqn{ns}. See also [prepare_data()] to conveniently prepare individual data
#'   for the analysis with one, or two time scales.
#'
#' @inheritParams exposures_events_1d
#' @param u A vector of fixed times at entry in the process, measured over the
#'   first time scale.
#' @param bins_list is a list with the following (necessary) elements
#' (usually prepared by [make_bins()]):
#'   * `bins_u` a vector of extreme values for the bins over the `u` axis
#'   * `bins_s` a vector of extreme values for the bins over the `s` axis
#' @param individual A Boolean. Default is `FALSE`: if `FALSE` computes the matrices
#'   R and Y collectively for all observations; if `TRUE` computes the matrices
#'   R and Y separately for each individual record.
#'
#' @return A list with the following elements:
#' * `R` an array of exposure times: if `individual == TRUE`,
#'   then `R` is an array of dimension \eqn{nu} by \eqn{ns} by \eqn{n},
#'   otherwise is an array of dimension \eqn{nu} by \eqn{ns}
#' * `Y`an array of event counts: if `individual == TRUE`,
#'   then `Y` is an array of
#'   dimension \eqn{nu} by \eqn{ns} by \eqn{n}, otherwise is an array of
#'   dimension \eqn{nu} by \eqn{ns}
#'
#' @examples
#' # ---- Bin colon cancer data by time at randomization and time since recurrence ----
#' # First create vectors of bins (using function `make_bins()`)
#' bins <- make_bins(u = reccolon2ts$timer, s_out = reccolon2ts$timesr,
#' du = 30, ds = 30)
#' # Now bin data (note: the s_in argument is omitted because data are not left truncated)
#' bindata2d <- exposures_events_2d(u = reccolon2ts$timer,
#' s_out = reccolon2ts$timesr, ev = reccolon2ts$status, bins = bins)
#'
#' @author Angela Carollo \email{carollo@@demogr.mpg.de}
#'
#' @export
#'
exposures_events_2d <- function(u,
                             s_in = NULL,
                             s_out,
                             ev,
                             bins_list,
                             individual = FALSE) {

# ---- Checks on arguments ----
  if (length(bins_list) == 1) {
    stop("`bins_list` must be a list with 2 elements.")
  }

  if(is.null(bins_list$bins_u)) stop("`bins_list$bins_u` missing.")

  if (!is.null(s_in) & (length(s_in) != length(s_out))) {
    stop("`s_in`and `s_out` must have the same length.")
  }

  if (length(s_out) != length(u)) {
    stop("`s_out` and `u` must have the same length.")
  }

  if (length(s_out) != length(ev)) {
    stop("`s_out` and `ev` must have the same length.")
  }

  if (any(ev > 1)) {
    stop("only 0 or 1 allowed in vector `ev`.")
  }


  # ---- Preliminary operations ----
  n <- length(s_out)

  if (is.null(s_in)) {
    s_in <- rep(0, n)
    message("`s_in = NULL`. I will use `s_in = 0` for all observations.")
  }

  nu <- bins_list$nu
  ns <- bins_list$ns

  liu <- bins_list$bins_u[1:nu]

  # ---- Calculation of individual exposure and event counts ----
  R <- Y <- array(NA, dim = c(nu, ns, n))

  for (ind in 1:n) {
    R_ind <- Y_ind <- matrix(0, nu, ns)
    ry_ind <- exposures_events_1d(s_in[ind], s_out[ind], ev[ind], bins = bins_list$bins_s)
    where_u <- max(which(liu <= u[ind]))
    R_ind[where_u, ] <- ry_ind$r
    Y_ind[where_u, ] <- ry_ind$y
    R[, , ind] <- R_ind
    Y[, , ind] <- Y_ind
  }

  # ---- Final results ----
  if (individual == FALSE) {
    R <- apply(R, 1:2, sum)
    Y <- apply(Y, 1:2, sum)
  }

  return(list(R = R, Y = Y))
}
