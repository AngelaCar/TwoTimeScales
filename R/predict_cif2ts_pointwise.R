#' Point-wise prediction of cumulative incidence over 2 time scale
#'
#' @param fitted_models a list with cause-specific hazard models
#' @param u The value(s) of `u` where prediction is required
#' @param s The value(s) of `s` where prediction is required
#' @param ds (optional) The distance between two consecutive points on the `s` axis.
#'           If not provided, an optimal minimum value will be chosen automatically and
#'           a warning is returned.
#' @return A data.frame with one row containing:
#'         the values of `u` and `s`for which predictions of the overall survival
#'         (`surv`) probability, and the values of the cumulative incidence functions,
#'         one for each cause, are obtained.
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
#'
#' fakedata <- as.data.frame(cbind(id, u, s, ev))
#' # cause 1
#' fakedata2ts1 <- prepare_data(
#'   u = fakedata$u,
#'   s_out = fakedata$s,
#'   ev = fakedata$ev,
#'   min_u = 2, min_s = 0,
#'   ds = .5
#' )
#' # Fit a fake model - not optimal smoothing for cause type 1
#' fakemod1 <- fit2ts(fakedata2ts1,
#'   optim_method = "grid_search",
#'   lrho = list(
#'     seq(1, 1.5, .5),
#'     seq(1, 1.5, .5)
#'   )
#' )
#' # cause 2
#' fakedata2ts2 <- prepare_data(
#'   u = fakedata$u,
#'   s_out = fakedata$s,
#'   ev = 1 - (fakedata$ev),
#'   min_u = 2, min_s = 0,
#'   ds = .5
#' )
#' # Fit a fake model - not optimal smoothing for cause 2
#' fakemod2 <- fit2ts(fakedata2ts2,
#'   optim_method = "grid_search",
#'   lrho = list(
#'     seq(1, 1.5, .5),
#'     seq(1, 1.5, .5)
#'   )
#' )
#' predict_cif2ts_pointwise(
#'   fitted_models = list(fakemod1, fakemod2),
#'   u = 5.2, s = 4.4
#' )
predict_cif2ts_pointwise <- function(fitted_models = list(),
                                     u, s, ds = NULL) {
  # the Bsplines will be the same for each cause specific models, so we just need to
  # do the following step once

  Bbases <- fitted_models[[1]]$optimal_model$Bbases
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
    "midu" = midu,
    "umin" = umin,
    "umax" = umax,
    "mids" = mids,
    "smin" = smin,
    "smax" = smax
  )

  # if u is already = to one of the interval points, do nothing, otherwise
  if (!(u %in% new_grid$midu)) {
    if (u < new_grid$umin | u > new_grid$umax) {
      stop("New `u` outside of range of `B_u`.")
    } else {
      newu <- unique(sort(c(new_grid$midu, u)))
      new_grid$midu <- newu
    }
  }
  Bu <- JOPS::bbase(new_grid$midu,
    xl = umin, xr = umax,
    nseg = attributes(Bbases$Bu)$nseg,
    bdeg = attributes(Bbases$Bu)$bdeg
  )

  Bs <- JOPS::bbase(new_grid$mids,
    nseg = attributes(Bbases$Bs)$nseg,
    xl = smin, xr = smax,
    bdeg = attributes(Bbases$Bs)$bdeg
  )

  # select only value for u,s
  which_ints <- findInterval(s, ints)
  grid_us <- expand.grid(u = new_grid$midu, s = new_grid$mids)
  n <- nrow(grid_us)
  L <- length(fitted_models) # n. of competing causes
  J <- dim(Bu)[1]
  K <- dim(Bs)[1]

  csh_a <- array(NA, dim = c(J, K, L))
  csh <- matrix(0, n, L)
  colnames(csh) <- paste0("csh_", 1:L)
  cscumh_a <- array(NA, dim = c(J, K, L))
  cif_a <- array(NA, dim = c(J, K, L))
  cif <- matrix(0, n, L)
  colnames(cif) <- paste0("cif_", 1:L)

  for (l in 1:L) {
    # ---- Calculate (baseline) hazard ----
    Eta <- Bu %*% fitted_models[[l]]$optimal_model$Alpha %*% t(Bs)
    csh_a[, , l] <- exp(Eta)
    csh[, l] <- c(csh_a[, , l])

    # ---- Cumulative hazard ----
    cscumh_a[, , l] <- t(apply(csh_a[, , l] * ds, 1, cumsum))
  }
  surv <- exp(-(apply(cscumh_a, 1:2, sum)))
  grid_us$surv <- c(surv)

  for (l in 1:L) {
    cif_a[, , l] <- t(apply((csh_a[, , l] * surv) * ds, 1, cumsum))
    cif[, l] <- c(cif_a[, , l])
  }

  grid_us <- cbind(grid_us, cif)
  selection <- grid_us[grid_us$u == u & grid_us$s == new_grid$mids[which_ints], ]

  return(selection)
}
