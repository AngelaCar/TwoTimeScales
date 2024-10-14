#' Bin data on one time scale
#'
#' @description `exposure_events_1d()` computes aggregated measures of exposure
#'   times and event counts starting from individual records of time at entry,
#'   time at exit and event's indicator, over one time scale (`s`).
#'
#' @details The time scale `s` is divided into bins of equal size, which are
#'   provided as input to the function. Then, the time-at-risk for each
#'   individual is split according to these bins, and an event indicator is
#'   placed in the bin where the exit time is located. Finally, the individual
#'   contributions are summed in each bin to provide a vector of total exposure
#'   time and total event counts. See also [prepare_data()] to
#'   conveniently prepare individual data for the analysis with one, or two time
#'   scales.
#'
#' @param s_in A vector of (possibly left truncated) times at entry. If this is
#'   not provided by the user, the function will consider a value of 0 for all
#'   observations.
#' @param s_out A vector of times at event or censoring.
#' @param ev A vector of events' indicators (possible values 0/1).
#' @param bins A vector of interval breaks for discretization (see also [make_bins()]).
#'
#' @return A list with the following elements:
#'   * `R`  A matrix of dimension \eqn{n} by \eqn{ns} containing the exposure times for each
#'     individual separately.
#'   * `r`  A vector of exposure times.
#'   * `Y`  A matrix of dimension \eqn{n} by \eqn{ns} containing the event counts for each
#'     individual separately
#'   * `y`  A vector of event counts.
#'
#'   If the length of the input vectors do not match, an error message is
#'   returned.
#'
#' @examples
#' # ---- Bin colon cancer data by time since recurrence ----
#' # First create vector of bins
#' bins1ts <- make_bins(s_in = reccolon2ts$entrys, s_out = reccolon2ts$timesr, ds = 30)
#' bindata <- exposures_events_1d(s_in = reccolon2ts$entrys,
#' s_out = reccolon2ts$timesr, ev = reccolon2ts$status, bins = bins1ts$bins_s)
#'
#' @author Angela Carollo \email{carollo@@demogr.mpg.de} and Paul Eilers \email{p.eilers@@erasmus.nl}
#'
#' @export
#'
#'
exposures_events_1d <- function(s_in = NULL, s_out, ev,
                                bins) {
  # ---- Checks on inputs ----
  if (!is.null(s_in) & (length(s_in) != length(s_out))) {
    stop("`s_in`and `s_out` must be of the same length.\n")
  }

  if (length(s_out) != length(ev)) {
    stop("`s_out` and `ev` must be of the same length.\n")
  }

  if (length(bins) == 1) {
    stop("`bins` must be a vector of length > 1.\n")
  }

  if (any(ev > 1)) {
    stop("only 0 or 1 allowed in vector `ev`.")
  }

  # ---- Preliminary operations ----
  n <- length(s_out)

  if (is.null(s_in)) {
    s_in <- rep(0, n)
  }
  lint <- bins[-length(bins)]
  rint <- bins[-1]

  ds <- lint[2] - lint[1]

  # ---- Compute exposure and event counts ----

  # Exposure
  R1 <- ds * (outer(s_out, rint, ">="))
  R2 <- 1 * (outer(s_out, rint, "<")) * (outer(s_out, lint, ">")) * outer(s_out, lint, "-")
  R3 <- ds * (outer(s_in, rint, ">="))
  R4 <- 1 * (outer(s_in, rint, "<")) * (outer(s_in, lint, ">")) * (outer(s_in, lint, "-"))
  R <- (R1+R2) - (R3+R4)
  r <- colSums(R)

  # Events
  Y <- 1 * (outer(s_out, rint, "<=")) * (outer(s_out, lint, ">")) * outer(ev, rep(1, length(bins) - 1))
  y <- colSums(Y)

  # ---- Return list ----
  return(list(R = R, r = r, Y = Y, y = y))
}
