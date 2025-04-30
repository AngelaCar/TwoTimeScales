#' Point-wise prediction hazard 2 time scale
#'
#' @param fitted_model An object of class `'haz2ts'` fitted via `fit2ts()`.
#' @param u The value(s) of `u` where prediction is required
#' @param s The value(s) of `s` where prediction is required
#' @param ds (optional) The distance between two consecutive points on the `s` axis.
#'           If not provided, an optimal minimum value will be chosen automatically and
#'           a warning is returned.
#'
#' @return A data.frame with one row and 6 variable: the values of `u` and `s`
#'         for which predictions of `hazard`, `se_hazard`, the cumulative hazard
#'         `cumhaz` and the `survival` probability are obtained
#' @export
#'
#' @examples
#' id <- 1:20
#' u <- c(
#'   5.43, 3.25, 8.15, 5.53, 7.28, 6.61, 5.91, 4.94, 4.25, 3.86, 4.05, 6.86,
#'   4.94, 4.46, 2.14, 7.56, 5.55, 7.60, 6.46, 4.96
#' )
#' s <- c(
#'   0.44, 4.89, 0.92, 1.81, 2.02, 1.55, 3.16, 6.36, 0.66, 2.02, 1.22, 3.96,
#'   7.07, 2.91, 3.38, 2.36, 1.74, 0.06, 5.76, 3.00
#' )
#' ev <- c(1, 0, 0, 1, 0, 1, 0, 1, 1, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 1)
#' fakedata <- as.data.frame(cbind(id, u, s, ev))
#' fakedata2ts <- prepare_data(
#'   u = fakedata$u,
#'   s_out = fakedata$s,
#'   ev = fakedata$ev,
#'   ds = .5
#' )
#' # Fit a fake model - not optimal smoothing
#' fakemod <- fit2ts(fakedata2ts,
#'   optim_method = "grid_search",
#'   lrho = list(
#'     seq(1, 1.5, .5),
#'     seq(1, 1.5, .5)
#'   )
#' )
#' predict_haz2ts_pointwise(fakemod, u = 5, s = 4.44)
predict_haz2ts_pointwise <- function(fitted_model,
                                     u, s, ds = NULL) {
  Bbases <- fitted_model$optimal_model$Bbases
  midu <- attributes(Bbases$Bu)$x
  umin <- attributes(Bbases$Bu)$xl
  umax <- attributes(Bbases$Bu)$xr
  smin <- attributes(Bbases$Bs)$xl
  smax <- attributes(Bbases$Bs)$xr

  # if ds is not provided, find the smallest `ds` given `s`
  if (is.null(ds)) {
    ds <- round(min(s - as.integer(s)), 1)
    if(ds == 0) ds <- .1
    message("chosen interval: ds = ", ds)

  }

  # remake intervals over `s`
  K <- ceiling((smax - smin) / ds)
  ints <- seq(smin, smin + K * ds, by = ds)
  mids <- ints[-1] - ds / 2

  new_grid <- list(
    "intu" = midu,
    "umin" = umin,
    "umax" = umax,
    "ints" = mids,
    "smin" = smin,
    "smax" = smax
  )

  # if u is already = to one of the interval points, do nothing, otherwise
  if (!(u %in% new_grid$intu)) {
    if (u < new_grid$umin | u > new_grid$umax) {
      stop("New `u` outside of range of `B_u`.")
    } else {
      newu <- unique(sort(c(new_grid$intu, u)))
      new_grid$intu <- newu
    }
  }
  Bu <- JOPS::bbase(new_grid$intu,
    xl = umin, xr = umax,
    nseg = attributes(Bbases$Bu)$nseg,
    bdeg = attributes(Bbases$Bu)$bdeg
  )

  # Don't need to do the same for `s` because we have made the intervals small
  # enough
  Bs <- JOPS::bbase(new_grid$ints,
    nseg = attributes(Bbases$Bs)$nseg,
    xl = smin, xr = smax,
    bdeg = attributes(Bbases$Bs)$bdeg
  )

  # Now calculate hazard and only select the value corresponding to the new u,s pair

  # ---- Calculate (baseline) hazard ----
  Eta <- Bu %*% fitted_model$optimal_model$Alpha %*% t(Bs)
  Haz <- exp(Eta)

  # ---- Calculate Standard Errors for the log-hazard ----
  cu <- ncol(Bu)
  cs <- ncol(Bs)
  # Calculate the row tensor products
  TBu <- Rtens(Bu)
  TBs <- Rtens(Bs)

  Cov_Alpha <- array(fitted_model$optimal_model$Cov_Alpha, c(cu, cs, cu, cs))
  Cov_Alpha <- aperm(Cov_Alpha, c(1, 3, 2, 4))
  Cov_Alpha <- matrix(Cov_Alpha, c(cu^2, cs^2))
  Dim <- c(nrow(TBu), nrow(TBs))
  Var_Eta <- matrix(TBu %*% Cov_Alpha %*% t(TBs), Dim)
  SE_Eta <- sqrt(Var_Eta)

  # ---- Calculate Standard Errors for the hazard ----
  SE_Haz <- Haz * SE_Eta

  # ---- Cumulative hazard ----
  CumHaz <- t(apply(Haz * ds, 1, cumsum))

  # select only value for u,s
  grid_us <- expand.grid(u = new_grid$intu, s = new_grid$ints)
  grid_us$hazard <- c(Haz)
  grid_us$se_hazard <- c(SE_Haz)
  grid_us$cumhaz <- c(CumHaz)
  grid_us$survival <- exp(-grid_us$cumhaz)

  which_ints <- findInterval(s, ints)
  selection <- grid_us[grid_us$u == u & grid_us$s == new_grid$ints[which_ints], ]
  return(selection)
}
