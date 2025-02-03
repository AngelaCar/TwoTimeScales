#' Point-wise prediction hazard 2 time scale
#'
#' @param fitted_model An object of class `'haz2ts'` fitted via `fit2ts()`.
#' @param u The value(s) of `u` where prediction is required
#' @param s The value(s) of `s` where prediction is required
#'
#' @return A data.frame with one row and 6 variable: the values of `u` and `s`
#'         for which predictions of `hazard`, `se_hazard`, the cumulative hazard
#'         `cumhaz` and the `survival` probability are obtained
#' @export
#'
#' @examples
#' id <- 1:20
#' u <- c(5.43, 3.25, 8.15, 5.53, 7.28, 6.61, 5.91, 4.94, 4.25, 3.86, 4.05, 6.86,
#'        4.94, 4.46, 2.14, 7.56, 5.55, 7.60, 6.46, 4.96)
#' s <- c(0.44, 4.89, 0.92, 1.81, 2.02, 1.55, 3.16, 6.36, 0.66, 2.02, 1.22, 3.96,
#'        7.07, 2.91, 3.38, 2.36, 1.74, 0.06, 5.76, 3.00)
#' ev <- c(1, 0, 0, 1, 0, 1, 0, 1, 1, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 1)
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
#' predict_haz2ts_pointwise(fakemod, u = 5, s = 4.44)

predict_haz2ts_pointwise <- function(fitted_model,
                                     u, s){
  Bbases <- fitted_model$optimal_model$Bbases
  midu <- attributes(Bbases$Bu)$x
  mids <- attributes(Bbases$Bs)$x
  du <- midu[2] - midu[1]
  ds <- mids[2] - mids[1]
  intu <- midu + du / 2
  intu <- c(intu[1] - du, intu)
  umin <- min(intu)
  umax <- max(intu)
  ints <- mids + ds / 2
  ints <- c(ints[1] - ds, ints)
  smin <- min(ints)
  smax <- max(ints)

  new_grid <- list(
    "intu" = intu,
    "umin" = umin,
    "umax" = umax,
    "ints" = ints,
    "smin" = smin,
    "smax" = smax
  )

  # if u is already = to one of the interval points, do nothing, otherwise
  if(!(u %in% new_grid$intu)){
    if(u < new_grid$umin | u > new_grid$umax){
      stop ("New `u` outside of range of `B_u`.")
    } else {
      newu <- unique(sort(c(new_grid$intu, u)))
      new_grid$intu <- newu
    }
  }
  Bu <- JOPS::bbase(new_grid$intu, nseg = attributes(Bbases$Bu)$nseg,
                      bdeg = attributes(Bbases$Bu)$bdeg)

  # do the same for s
  if(!(s %in% new_grid$ints)){
    if(s < new_grid$smin | s > new_grid$smax){
      stop ("New `s` outside of range of `B_s`.")
    } else {
      news <- unique(sort(c(new_grid$ints, s)))
      new_grid$ints <- news
    }
  }
  Bs <- JOPS::bbase(new_grid$ints, nseg = attributes(Bbases$Bs)$nseg,
                    bdeg = attributes(Bbases$Bs)$bdeg)

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
  #B <- kronecker(Bs, Bu)

  Cov_Alpha <- array(fitted_model$optimal_model$Cov_Alpha, c(cu, cs, cu, cs))
  Cov_Alpha <- aperm(Cov_Alpha, c(1, 3, 2, 4))
  Cov_Alpha <- matrix(Cov_Alpha, c(cu^2, cs^2))
  Dim <- c(nrow(TBu), nrow(TBs))
  Var_Eta <- matrix(TBu %*% Cov_Alpha %*% t(TBs), Dim)
  SE_Eta <- sqrt(Var_Eta)

  # ---- Calculate Standard Errors for the hazard ----
  SE_Haz <- Haz * SE_Eta

  # ---- Cumulative hazard ----
  ds <- diff(new_grid$ints)
  ds <- c(ds, max(ds))
  CumHaz <- t(apply(Haz * ds, 1, cumsum) )

  # select only value for u,s
  grid_us <- expand.grid(u = new_grid$intu, s = new_grid$ints)
  grid_us$hazard <- c(Haz)
  grid_us$se_hazard <- c(SE_Haz)
  grid_us$cumhaz <- c(CumHaz)
  grid_us$survival <- exp(-grid_us$cumhaz)

  selection <- grid_us[grid_us$u == u & grid_us$s == s,]
  return(selection)
}
