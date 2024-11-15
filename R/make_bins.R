#' Construct bins over one or more time axes
#'
#' @description `make_bins()` constructs the bins over the time axes and saves the extremes
#' of the bins in a vector.
#'
#' @details  It allows construction of bins over the time scales `t` and
#'   `s` and/or over the fixed-time axis `u`. The time scale
#'    `s` is always required. See also [prepare_data()] to conveniently
#'     prepare individual data for the analysis with one, or two time scales.
#'
#' @param t_in (optional) A vector of entry times on the time scale `t`.
#' @param t_out (optional) A vector of exit times on the time scale `t`.
#' @param u (optional) A vector of fixed-times at entry in the process.
#' @param s_in (optional) A vector of entry times on the time scale `s`.
#' @param s_out A vector of exit times on the time scale `s`.
#' @param min_t (optional) A minimum value for the bins over `t`.
#'   If `NULL`, the minimum of `t_in` will be used.
#' @param min_u (optional) A minimum value for the bins over `u`.
#'   If `NULL`, the minimum of `u` will be used.
#' @param min_s (optional) A minimum value for the bins over `s`.
#'   If `NULL`, the minimum of `s_in` will be used.
#' @param max_t (optional) A maximum value for the bins over `t`.
#'   If `NULL`, the maximum of `t_out` will be used.
#' @param max_u (optional) A maximum value for the bins over `u`.
#'   If `NULL`, the maximum of `u` will be used.
#' @param max_s (optional) A maximum value for the bins over `s`.
#'   If `NULL`, the maximum of `s_out` will be used.
#' @param dt (optional) A scalar giving the length of the intervals on the `t` time scale.
#' @param du (optional) A scalar giving the length of the intervals on the `u` axis.
#' @param ds A scalar giving the length of the intervals on the `s` time scale.
#'
#' @return A list with the following elements:
#' * `bins_t` if `t_out` is provided, this is a vector of bins extremes for the time scale `t`
#' * `midt` if `t_out` is provided, this is a vector with the midpoints of the bins over `t`
#' * `nt` if `t_out` is provided, this is the number of bins over `t`
#' * `bins_u` if `u` is provided, this is a vector of bins extremes for `u` axis
#' * `midu` if `u` is provided, this is a vector with the midpoints of the bins over `u`
#' * `nu` if `u` is provided, this is the number of bins over `u`
#' * `bins_s` is a vector of bins extremes for the time scale `s`
#' * `mids` is a vector with the midpoints of the bins over `s`
#' * `ns` is the number of bins over `s`
#'
#' @export

make_bins <- function(t_in = NULL, t_out = NULL,
                      u = NULL,
                      s_in = NULL, s_out,
                      min_t = NULL, max_t = NULL,
                      min_u = NULL, max_u = NULL,
                      min_s = NULL, max_s = NULL,
                      dt = NULL, du = NULL, ds) {
  # ---- Make bins ----
  # bins for t (if 2D bins over t and s)
  if (!is.null(t_out)) {
    if (is.null(t_in)) {
      t_in <- t_out - s_out
      message("`t_in` not provided. I will use `t_in = t_out - s_in`.")
    }
    if (is.null(dt)) dt <- ds
    if (is.null(min_t)) min_t <- min(t_in)
    if (is.null(max_t)) max_t <- max(t_out)
    if (!is.null(max_t) & (max_t < max(t_out))) {
      message("`max_t < max(t_out)`; will use `max(t_out)`.")
      max_t <- max(t_out)
    }
    if (!is.null(min_t) & (min_t > min(t_in))) {
      message("`min_t > min(t_out)`; will use `min(t_in)`.")
      min_t <- min(t_in)
    }
    K <- ceiling((max_t - min_t) / dt)
    bins_t <- seq(min_t, min_t + K * dt, by = dt)
  }

  # bins for u (if 2D bins over u and s)
  if (!is.null(u)) {
    if (is.null(du)) du <- ds
    if (is.null(min_u)) min_u <- min(u)
    if (is.null(max_u)) max_u <- max(u)
    if (!is.null(max_u) & (max_u < max(u))) {
      message("`max_u < max(u)`; will use `max(u)`.")
      max_u <- max(u)
    }
    if (!is.null(min_u) & (min_u > min(u))) {
      message("`min_u > min(u)`; will use `min(u)`.")
      min_u <- min(u)
    }
    K <- ceiling((max_u - min_u) / du)
    bins_u <- seq(min_u, min_u + K * du, by = du)
  }
  # bins for s
  if (is.null(s_in)) {
    s_in <- rep(0, length(s_out))
    message("`s_in = NULL`. I will use `s_in = 0` for all observations.")
  }
  if (is.null(min_s)) min_s <- min(s_in)
  if (is.null(max_s)) max_s <- max(s_out)
  if (!is.null(max_s) & (max_s < max(s_out))) {
    message("`max_s < max(s_out)`; will use `max(s_out)`.")
    max_s <- max(s_out)
  }
  if (!is.null(min_s) & (min_s > min(s_in))) {
    message("`min_s > min(s_out)`; will use `min(s_in)`.")
    min_s <- min(s_in)
  }
  K <- ceiling((max_s - min_s) / ds)
  bins_s <- seq(min_s, min_s + K * ds, by = ds)


  # ---- Return results ----
  bins <- list()

  if (!is.null(t_out)) {
    bins$bins_t <- bins_t
    bins$midt <- bins_t[-1] - dt / 2
    bins$nt <- length(bins$midt)
  }

  if (!is.null(u)) {
    bins$bins_u <- bins_u
    bins$midu <- bins_u[-1] - du / 2
    bins$nu <- length(bins$midu)
  }

  bins$bins_s <- bins_s
  bins$mids <- bins_s[-1] - ds / 2
  bins$ns <- length(bins$mids)

  return(bins)
}
