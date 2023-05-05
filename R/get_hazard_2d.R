#' Get estimated (log-)hazard surface with 2 time scales
#'
#' @description `get_hazard_2d()` takes as input the results of a model
#' estimated by `fit2ts` and it returns the estimated smooth log-hazard
#' and the smooth hazard together with their standard errors.
#'
#' If the model includes covariates, then only the baseline (log-)hazard is returned.
#' It is possible to provide values that define a new grid for plotting.
#' If not specified, the plotting is done using the same B-splines basis as for
#' the estimation. The function will check if the parameters for the grid
#' provided by the user are compatible with those originally used to construct
#' the B-splines for estimating the model. If not, the grid will be adjusted
#' accordingly and a warning will be returned.
#'
#' @inheritParams plot_haz2ts
#'
#' @return A list with the following elements:
#'   * `new_plot_grid` A list of new specifications of the grid for plotting.
#'   * `hazard` A matrix containing the estimated hazard.
#'   * `loghazard` A matrix containing the estimated log-hazard.
#'   * `SE_hazard` A matrix containing the estimated SEs for the hazard
#'   * `SE_loghazard` A matrix containing the estimated SEs for the log-hazard
#' @export
#'
get_hazard_2d <- function(fitted_model, plot_grid = NULL,
                          where_slices = NULL,
                          direction = c("u", "s", NULL)) {
  Bbases <- fitted_model$optimal_model$Bbases
  direction <- match.arg(direction)

  if (!is.null(plot_grid)) {
    # check if all information are provided
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

      if (umin < attributes(Bbases$Bu)$xl) {
        umin <- attributes(Bbases$Bu)$xl
        warning("`umin` is smaller than the lower limit of the domain of Bu. Left boundary adjusted to  =  ", attributes(Bbases$Bu)$xl)
      }
      if (umax > attributes(Bbases$Bu)$xr) {
        umax <- attributes(Bbases$Bu)$xr
        warning("`umax` is larger than the upper limit of the domain of Bu. Right boundary adjusted to  =  ", attributes(Bbases$Bu)$xr)
      }
      K <- ceiling((umax - umin) / du)
      intu <- seq(umin, umin + K * du, by = du)

      if (smin < attributes(Bbases$Bs)$xl) {
        smin <- attributes(Bbases$Bs)$xl
        warning("`smin` is smaller than the lower limit of the domain of Bs. Left boundary adjusted to  =  ", attributes(Bbases$Bs)$xl)
      }
      if (smax > attributes(Bbases$Bs)$xr) {
        smax <- attributes(Bbases$Bs)$xr
        warning("`smax` is larger than the upper limit of the domain of Bs. Right boundary adjusted to  =  ", attributes(Bbases$Bs)$xr)
      }
      K <- ceiling((smax - smin) / ds)
      ints <- seq(smin, smin + K * ds, by = ds)

      # Evaluate old bases in new grid of points
      Bu <- JOPS::bbase(intu, nseg = attributes(Bbases$Bu)$nseg, bdeg = attributes(Bbases$Bu)$bdeg)
      Bs <- JOPS::bbase(ints, nseg = attributes(Bbases$Bs)$nseg, bdeg = attributes(Bbases$Bs)$bdeg)
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
    }
  } else {
    Bu <- Bbases$Bu
    Bs <- Bbases$Bs
    new_grid <- NULL
  }

  if(!is.null(where_slices)){
    if(is.null(direction)) stop("Direction for cutting slices missing.")
    if(direction == "u" ){
      if(min(where_slices) < attributes(Bbases$Bu)$xl | max(where_slices) > attributes(Bbases$Bu)$xr){
        stop ("Desired cutting points outside of range of `B_u`.")
      } else {
        Bu <- JOPS::bbase(where_slices, nseg = attributes(Bbases$Bu)$nseg, bdeg = attributes(Bbases$Bu)$bdeg)
      }
    }
    if(direction == "s"){
      if(min(where_slices) < attributes(Bbases$Bs)$xl | max(where_slices) > attributes(Bbases$Bs)$xr){
        stop ("Desired cutting points outside of range of `B_s`.")
      } else {
        Bs <- JOPS::bbase(where_slices, nseg = attributes(Bbases$Bs)$nseg, bdeg = attributes(Bbases$Bs)$bdeg)
      }
    }

  }

  # ---- Calculate (baseline) hazard ----
  Eta <- Bu %*% fitted_model$optimal_model$Alpha %*% t(Bs)
  Haz <- exp(Eta)

  # ---- Calculate Standard Errors for the log-hazard ----
  cu <- ncol(Bu)
  cs <- ncol(Bs)
  # Calculate the row tensor products
  TBu <- Rtens(Bu)
  TBs <- Rtens(Bs)
  B <- kronecker(Bs, Bu)

  Cov_Alpha <- array(fitted_model$optimal_model$Cov_Alpha, c(cu, cs, cu, cs))
  Cov_Alpha <- aperm(Cov_Alpha, c(1, 3, 2, 4))
  Cov_Alpha <- matrix(Cov_Alpha, c(cu^2, cs^2))
  Dim <- c(nrow(TBu), nrow(TBs))
  Var_Eta <- matrix(TBu %*% Cov_Alpha %*% t(TBs), Dim)
  SE_Eta <- sqrt(Var_Eta)

  # ---- Calculate Standard Errors for the hazard ----
  SE_Haz <- Haz * SE_Eta

  # ---- Return results in a list ----
  results <- list(
    "new_plot_grid" = new_grid,
    "hazard" = Haz,
    "loghazard" = Eta,
    "SE_hazard" = SE_Haz,
    "SE_loghazard" = SE_Eta
  )

  return(results)
}
