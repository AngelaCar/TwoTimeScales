#' Get estimated (log-)hazard surface with 2 time scales
#'
#' @description `get_hazard_2d()` takes as input the results of a model
#' estimated by `fit2ts()`, `fitpgam()` or `fitvcm()` and it returns the
#' estimated smooth log-hazard and the smooth hazard together with their standard errors.
#'
#' It is possible to provide values that define a new grid for evaluation of the
#' estimated hazard.
#' If not specified, the hazard is evaluated on the same grid used for the
#' binning of the data, and therefore the estimation of the model.
#' The function will check if the parameters for the new grid
#' provided by the user are compatible with those originally used to construct
#' the B-splines for estimating the model. If not, the grid will be adjusted
#' accordingly and a warning will be returned.
#'
#' @inheritParams plot.haz2ts
#' @param fitted_model is an object of class `"haz2ts"`, `"haz2tsPGAM"`, or `"haz2tsVCM"`,
#'  the output of the function `fit2ts()`, `fitpgam()` or `fitvcm()`, respectively.
#' @param tmax The maximum value of `t` that should be plotted.
#' @param midpoints A Boolean. Default is `FALSE`. If `TRUE`, the estimated
#'        quantities are evaluated at the midpoints of the rectangles
#'        (or parallelograms) of the grids, rather than at each grid-point.
#'
#' @return A list with the following elements:
#'   * `new_plot_grid` A list of parameters that specify the new grid, of the form
#'      list("intu", "umin", "umax", "du", "ints", "smin", "smax", "ds")
#'   * `nBu` The B-spline basis for `u`, evaluated over the new grid.
#'   * `nBs` The B-spline basis for `s`, evaluated over the new grid.
#'   * `hazard` A matrix containing the estimated hazard values.
#'   * `loghazard` A matrix containing the estimated log-hazard values.
#'   * `log10hazard`A matrix containing the estimated log10-hazard values.
#'   * `SE_hazard` A matrix containing the estimated SEs for the hazard.
#'   * `SE_loghazard` A matrix containing the estimated SEs for the log-hazard.
#'   * `SE_log10haz` A matrix containing the estimated SEs for the log10-hazard.
#' @export
#'
#' @examples
#' # Create some fake data - the bare minimum
#' id <- 1:20
#' u <- c(
#'   5.43, 3.25, 8.15, 5.53, 7.28, 6.61, 5.91, 4.94, 4.25, 3.86, 4.05, 6.86,
#'   4.94, 4.46, 2.14, 7.56, 5.55, 7.60, 6.46, 4.96
#' )
#' s <- c(
#'   0.44, 4.89, 0.92, 1.81, 2.02, 1.55, 3.16, 6.36, 0.66, 2.02, 1.22, 3.96,
#'   7.07, 2.91, 3.38, 2.36, 1.74, 0.06, 5.76, 3.00
#' )
#' ev <- c(1, 0, 0, 1, 0, 1, 0, 1, 1, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 1) #'
#'
#' fakedata <- as.data.frame(cbind(id, u, s, ev))
#' fakedata2ts <- prepare_data(
#'   data = fakedata,
#'   u = "u",
#'   s_out = "s",
#'   ev = "ev",
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
#' # Obtain 2d hazard
#' get_hazard_2d(fakemod)
#'
#' get_hazard_2d(fakemod,
#'   plot_grid = list(
#'     c(umin = 3, umax = 8.5, du = .1),
#'     c(smin = 0, smax = 7.1, ds = .1)
#'   )
#' )
get_hazard_2d <- function(fitted_model,
                          plot_grid = NULL,
                          where_slices = NULL,
                          direction = c(NULL, "u", "s"),
                          tmax = NULL,
                          midpoints = FALSE) {
  Bbases <- fitted_model$optimal_model$Bbases
  direction <- match.arg(direction)

  # ---- Make grid ----
  model_specifications <- list(
    "smin" = attributes(Bbases$Bs)$xl,
    "smax" = attributes(Bbases$Bs)$xr
  )
  if(inherits(fitted_model, "haz2tsVCM")){
    model_specifications$umin = min(fitted_model$ufitting)
    model_specifications$umax = max(fitted_model$ufitting)
  } else {
    model_specifications$umin = attributes(Bbases$Bu)$xl
    model_specifications$umax = attributes(Bbases$Bu)$xr
  }

  if (!is.null(plot_grid)) {
    # check if all information is provided
    if ((length(plot_grid[[1]]) != 3) | (length(plot_grid[[2]]) != 3)) {
      stop("One or more arguments from plot_grid are missing.")
    } else {
      new_grid <- make_grid(plot_grid = plot_grid,
                            class_fitmodel = class(fitted_model),
                            model_specifications = model_specifications,
                            where_slices = where_slices,
                            direction = direction,
                            tmax = tmax,
                            midpoints = midpoints)
    }

  } else { # no new grid is specified
    if (inherits(fitted_model, "haz2tsVCM")) {
      new_grid <- list(
        "intu" = fitted_model$ufitting,
        "umin" = min(fitted_model$ufitting),
        "umax" = max(fitted_model$ufitting)
      )
    } else {
      Bu <- Bbases$Bu
      new_grid <- list(
        "intu" = unique(attributes(Bu)$x),
        "umin" = attributes(Bu)$xl,
        "umax" = attributes(Bu)$xr
      )
    }
      Bs <- Bbases$Bs
      new_grid$ints <- unique(attributes(Bs)$x)
      new_grid$smin <- attributes(Bs)$xl
      new_grid$smax <- attributes(Bs)$xr
      new_grid$du <- new_grid$intu[2] - new_grid$intu[1]
      new_grid$ds <- new_grid$ints[2] - new_grid$ints[1]
      if (midpoints) {
        midu <- new_grid$intu[-1] - new_grid$du / 2
        mids <- new_grid$ints[-1] - new_grid$ds / 2
        new_grid$intu <- midu
        new_grid$ints <- mids
      }

      # If slices are required the grid needs to be adjusted
      if (!is.null(where_slices)) {
        if (is.null(direction)) stop("Direction for cutting slices missing.")
        if (direction == "u") {
          if (min(where_slices) < model_specifications$umin | max(where_slices) > model_specifications$umax) {
            stop("Desired cutting points outside of range of `B_u`.")
          } else {
            newu <- unique(sort(c(new_grid$intu, where_slices)))
            new_grid$intu <- newu
          }
        }
        if (direction == "s") {
          if (min(where_slices) < model_specifications$smin | max(where_slices) > model_specifications$smax) {
            stop("Desired cutting points outside of range of `B_s`.")
          } else {
            news <- unique(sort(c(new_grid$ints, where_slices)))
            new_grid$ints <- news
            grid_us <- expand.grid(u = new_grid$intu, s = new_grid$ints)
            grid_us$t <- with(grid_us, u + s)
            t <- unique(grid_us$t)
            if (!is.null(tmax)) {
              new_grid$intt <- t[t <= tmax]
            } else {
              new_grid$intt <- t
            }
            grid_ts <- expand.grid(t = new_grid$intt, s = new_grid$ints)
            grid_ts$u <- grid_ts$t - grid_ts$s
            newu <- unique(grid_ts$u)
            newu <- newu[newu >= new_grid$umin & newu <= new_grid$umax]
            new_grid$intu <- sort(newu)
          }
        }
      }
    }
 # ---- Re-evalutate B-splines in new grid ----
  if (inherits(fitted_model, "haz2ts")) {
    # Evaluate old bases in new grid of points
    Bu <- JOPS::bbase(new_grid$intu,
                      nseg = attributes(Bbases$Bu)$nseg,
                      bdeg = attributes(Bbases$Bu)$bdeg,
                      xl = attributes(Bbases$Bu)$xl,
                      xr = attributes(Bbases$Bu)$xr
    )
    Bs <- JOPS::bbase(new_grid$ints,
                      nseg = attributes(Bbases$Bs)$nseg,
                      bdeg = attributes(Bbases$Bs)$bdeg,
                      xl = attributes(Bbases$Bs)$xl,
                      xr = attributes(Bbases$Bs)$xr
    )
  } else {
    if(inherits(fitted_model, "haz2tsPGAM")){
      new_grid_long <- expand.grid(new_grid$intu, new_grid$ints)
      names(new_grid_long) <- c("intu", "ints")
      Bu <- JOPS::bbase(new_grid_long$intu,
                        nseg = attributes(Bbases$Bu)$nseg,
                        bdeg = attributes(Bbases$Bu)$bdeg,
                        xl = attributes(Bbases$Bu)$xl,
                        xr = attributes(Bbases$Bu)$xr
      )
      Bs <- JOPS::bbase(new_grid_long$ints,
                        nseg = attributes(Bbases$Bs)$nseg,
                        bdeg = attributes(Bbases$Bs)$bdeg,
                        xl = attributes(Bbases$Bs)$xl,
                        xr = attributes(Bbases$Bs)$xr
      )
      B <- cbind(Bu, Bs)
    } else {
      if(inherits(fitted_model, "haz2tsVCM")){
        Bs <- JOPS::bbase(new_grid$ints,
                          nseg = attributes(Bbases$Bs)$nseg,
                          bdeg = attributes(Bbases$Bs)$bdeg,
                          xl = attributes(Bbases$Bs)$xl,
                          xr = attributes(Bbases$Bs)$xr
        )
      }
    }
  }


  # ---- Obtain Estimates ----
  # ---- If fully interactive model 2D ----
  if (inherits(fitted_model, "haz2ts")) {
    # ---- Calculate (baseline) hazard ----
    Eta <- Bu %*% fitted_model$optimal_model$Alpha %*% t(Bs)
    Haz <- exp(Eta)
    Log10Haz <- log10(Haz)

    # ---- Calculate Standard Errors for the log-hazard ----
    cu <- ncol(Bu)
    cs <- ncol(Bs)
    # Calculate the row tensor products
    TBu <- Rtens(Bu)
    TBs <- Rtens(Bs)
    # B <- kronecker(Bs, Bu)

    Cov_Alpha <- array(fitted_model$optimal_model$Cov_Alpha, c(cu, cs, cu, cs))
    Cov_Alpha <- aperm(Cov_Alpha, c(1, 3, 2, 4))
    Cov_Alpha <- matrix(Cov_Alpha, c(cu^2, cs^2))
    Dim <- c(nrow(TBu), nrow(TBs))
    Var_Eta <- matrix(TBu %*% Cov_Alpha %*% t(TBs), Dim)
    SE_Eta <- sqrt(Var_Eta)

    # ---- Calculate Standard Errors for the hazard ----
    SE_Haz <- Haz * SE_Eta

    # ---- Calculates Standard Errors for log10(hazard) ----
    const <- log(10)
    SE_log10haz <- abs(1 / (Haz * const)) * SE_Haz
  }


  # --- If VCM model ----
  if (inherits(fitted_model, "haz2tsVCM")) {
    eta <- Bs %*% fitted_model$optimal_model$Alpha[, 1]
    se_eta <- Bs %*% fitted_model$optimal_model$SE_Alpha[, 1]
    vc <- Bs %*% fitted_model$optimal_model$Alpha[, 2]
    se_vc <- Bs %*% fitted_model$optimal_model$SE_Alpha[, 2]

    Eta <- matrix(rep(eta, length(new_grid$intu)),
      byrow = T, nrow = length(new_grid$intu),
      ncol = length(new_grid$ints)
    ) +
      outer(new_grid$intu, as.vector(vc))

    Haz <- exp(Eta)
    Log10Haz <- log10(Haz)

    SE_Eta <-
      matrix(rep(se_eta, length(new_grid$intu)),
                                byrow = T, nrow = length(new_grid$intu),
                                ncol = length(new_grid$ints)
    ) +
      outer(new_grid$intu, as.vector(se_vc))

    # ---- Calculate Standard Errors for the hazard ----
    SE_Haz <- Haz * SE_Eta

    # ---- Calculates Standard Errors for log10(hazard) ----
    const <- log(10)
    SE_log10haz <- abs(1 / (Haz * const)) * SE_Haz

    # No Bu for VCM
    Bu <- NA
  }

  # ---- If P-GAM model ----
  if (inherits(fitted_model, "haz2tsPGAM")) {
    eta <- B %*% fitted_model$optimal_model$alpha
    Eta <- matrix(eta, length(new_grid$intu), length(new_grid$ints))
    Haz <- exp(Eta)
    Log10Haz <- log10(Haz)

    # ---- Calculate Standard Errors for the log-hazard ----
    SE_eta <- B %*% fitted_model$optimal_model$SE_alpha
    SE_Eta <- matrix(SE_eta, length(new_grid$intu), length(new_grid$ints))

    # ---- Calculate Standard Errors for the hazard ----
    SE_Haz <- Haz * SE_Eta

    # ---- Calculates Standard Errors for log10(hazard) ----
    const <- log(10)
    SE_log10haz <- abs(1 / (Haz * const)) * SE_Haz
  }


  # ---- Return results in a list ----
  results <- list(
    "new_plot_grid" = new_grid,
    "nBu" = Bu,
    "nBs" = Bs,
    "hazard" = Haz,
    "loghazard" = Eta,
    "log10hazard" = Log10Haz,
    "SE_hazard" = SE_Haz,
    "SE_loghazard" = SE_Eta,
    "SE_log10hazard" = SE_log10haz
  )

  return(results)
}
