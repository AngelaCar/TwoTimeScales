#' Plot method for a haz2ts object.
#'
#' @description  `plot.haz2ts()` is the plot method for objects of class `haz2ts`.
#'  It produces several kinds of plots of the fitted model with two
#'   time scales (see [fit2ts()]), either in the original (t,s) plane, while respecting the
#'   constraint imposed by the relation of the two time scales, or in the
#'   transformed (u,s) plane.
#'
#' @param x The output of the function `fit2ts`. This is an object of
#'   class `"haz2ts"`.
#' @param which_plot The type of plot required. Can be one of `"hazard"`
#'   (default), `"covariates"`, `"SE"`, `"slices"`, `"survival"` or `"cumhaz"`
#'   (see details section).
#' @param plot_grid (optional) A list containing the parameters to build a new
#'   finer grid of intervals over u and s for plotting. This must be of the
#'   form: `plot_grid = list(c(umin, umax, du), c(smin, smax, ds))`, where
#'   `umin`, `umax` and `smin`, `smax` are the minimum and maximum values
#'   desired for the intervals over `u` and `s` respectively, and `du`, `ds` are
#'   distances between intervals over `u` and `s` respectively. Specifying a new
#'   denser grid is used to evaluate the B-spline bases used for estimation on
#'   such grid and plot the estimated surfaces with a greater level of details.
#'   If not specified, the plotting is done using the same B-splines bases as
#'   for the estimation. The function will check if the parameters for the grid
#'   provided by the user are compatible with those originally used to construct
#'   the B-splines for estimating the model. If not, the grid will be adjusted
#'   accordingly and a warning will be returned.
#' @param where_slices A vector of values for the cutting points of the desired
#'   slices of the surface. If `which_plot == "slices"`, please provide this
#'   argument. Please also provide this argument in case `which_plot = "survival`
#'   or `which_plot = "cumhaz` and `surv_slices = TRUE` or `cumhaz_slices = TRUE`,
#'   respectively.
#' @param direction If `which_plot == "slices"`, indicates the direction for
#'   cutting the surface. If `u`, then the surface will be cut at the selected
#'   values of `u` (indicated by `where_slices`), hence obtaining one-dimensional
#'   curves over `s`. If `s`, then the surface will be cut at the selected values
#'   of `s` (indicated by `where_slices`), hence obtaining one-dimensional curves
#'   over `u`.
#' @param plot_options A list with all possible options for any of the plots:
#'   * `loghazard` A Boolean. Default is `FALSE`. If `FALSE` the function
#'     returns a plot of the hazard surface, if `TRUE` the function returns
#'     a plot of the log-hazard surface.
#'   * `log10hazard` A Boolean. Default is `FALSE`. If `TRUE`,
#'     then a log_10 hazard surface is plotted.
#'   * `cut_extrapolated` A Boolean. Default is `TRUE`. Cuts away the
#'     extrapolated area of the (log-)hazard surface before plotting.
#'   * `rectangular_grid` A Boolean. Default is `FALSE`. If `TRUE`, a
#'     rectangular grid is used for plotting also in the (t,s)-plane as opposed
#'     to the grid of parallelograms used as default in the (t,s)-plane.
#'   * `original` A Boolean. Default is `TRUE`. Plot the (log-)hazard (and/or
#'     the SEs) in the (t,s)-plane. If `FALSE`, the (log-)hazard (and/or the SEs)
#'     will be plotted in the (u,s)-plane.
#'   * `tmax` The maximum value of `t` that should be plotted.
#'   * `surv_slices` A Boolean. Default is `FALSE`. If `TRUE` and
#'     `which_plot == "survival"`, plot survival curves over the time `s` for
#'     selected values of `u`, that are cross-sections of the 2D survival surface.
#'   * `cumhaz_slices` A Boolean. Default is `FALSE`. If `TRUE` and
#'     `which_plot == "cumhaz"`, plot cumulative hazards curves over the time `s` for
#'     selected values of `u`, that are cross-sections of the 2D cumulative hazard surface.
#'   * `midpoints` A Boolean. Default is `FALSE`. If `TRUE`, the estimated quantities
#'     (hazard, survival, etc.) will be evaluated in the mid-points of the bins
#'     rather than at the extremes. Set to `TRUE` if plotting estimated number of
#'     events.
#'   * `col_palette` A function defining the color palette. The default palette
#'     is `viridis::rev(plasma())`. Specifying the color palette as a function
#'     allows for greater flexibility than passing the palette as a vector.
#'
#'   * `n_shades` The number of color shades to plot, default is 50.
#'   * `breaks` The vector of breaks for the color legend. If `n_shades` is provided,
#'     this should be of length `n_shades + 1`.
#'   * `show_legend` A Boolean. Default is `TRUE`. If `FALSE` no legend will be
#'     plotted, useful for multi-panel figures with common legend. Works only
#'     for plots on rectangular grid (i.e. transformed (u,s) plane)
#'   * `main` The title of the plot.
#'   * `xlab` The label of the first time axis (plotted on the x axis).
#'   * `ylab` The label of the second time axis (plotted on the y axis).
#'   * `xlim` A vector with two elements defining the limits of the time scale
#'      on the x axis.
#'   * `ylim` A vector with two elements defining the limits of the time scale
#'      on the y axis.
#'   * `contour_lines` A Boolean. Default is `FALSE`. If `TRUE` white contour
#'     lines are added to the surfaces.
#'   * `contour_col` The color for the contour lines. Default is `white`.
#'   * `contour_cex` The magnification to be used for the contour lines.
#'      Default is `.8`.
#'   * `contour_nlev` The number of contour levels desired. Default is 10.
#'   * `cex_main` The magnification to be used for the main title, default is 1.2 .
#'   * `cex_lab` The magnification to be used for the axis labels, default is 1 .
#'   * `HR` A Boolean. If `TRUE` the HRs with their CIs will be plotted.
#'      Default is `FALSE` (plot the `beta` with their CIs).
#'   * `symmetric_CI` A Boolean. Default is `TRUE`. If a plot of the HRs is
#'     required (`HR == TRUE`), then plot symmetrical Confidence Intervals,
#'     based on the SEs for the HRs calculated by delta method.
#'     If `FALSE`, then CIs are obtained by exponentiating the CIs for the betas.
#'   * `confidence` The level of confidence for the CIs. Default is .95 (alpha
#'     = 0.05).
#'   * `col_beta` The color for the plot of the covariates' effects.
#'   * `pch` The symbol for plotting the point estimates.
#'   * `lwd` The line width.
#' @param \dots Further arguments to image.plot or image
#'
#' @details The vignette "visualization" presents and discusses all the different
#'  plotting options for the fitted model over two time scales.
#'  In most of the cases, the user will want to visualize the hazard surface over
#'  the two time scales. This can be plotted on the hazard scale, the log-hazard
#'  scale or the log10-hazard scale, by switching to `TRUE` the corresponding
#'  argument in `plot_options`.
#'  The survival and cumulative hazard functions can be plotted as two-dimensional
#'  surfaces over `u` and `s` or `t` and `s`. However, it is also very informative
#'  to plot them as one-dimensional curves over `s` (cross-sections or slices).
#'  This is done by selecting `which_plot = "survival"` and `surv_slices = TRUE`
#'  in `plot_options`. Additionally, a vector of values for the cutting points
#'  over the `u`-axis should be passed to the argument `where_slices`, together
#'  with setting `direction = u`.
#'  Similar plot is obtained for the cumulative hazard by selecting `which_plot = "cumhaz"`,
#'  `cumhaz_slices = TRUE`, see examples section.
#'  Please, notice that for the survival function and the cumulative hazard, only
#'  cross-sections of the surface for selected values of `u` (over the `s` time)
#'  can be plotted.
#'
#' @return A plot of the fitted model.
#'
#' @importFrom stats qnorm
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
#'
#' # plot the hazard surface
#' plot(fakemod)
#'
#' # plot the survival function as one-dimension curves over `s`
#' plot(fakemod,
#'   which_plot = "survival",
#'   direction = "u",
#'   where_slices = c(4, 6, 8),
#'   plot_options = list(
#'     surv_slices = TRUE
#'   )
#' )
#'
#' # Plot cross-sections of the hazard over `s` for selected values of `u`
#'
#' plot(fakemod,
#'   which_plot = "slices",
#'   where_slices = c(4, 6, 8),
#'   direction = "u",
#'   plot_options = list(
#'     main = "Cross-sections of the hazard",
#'     xlab = "Time",
#'     ylab = "Hazard"
#'   )
#' )
#'
#' @export

plot.haz2ts <- function(x,
                        plot_grid = NULL,
                        which_plot = c(
                          "hazard", "covariates", "SE", "slices",
                          "survival", "cumhaz"
                        ),
                        where_slices = NULL,
                        direction = c(NULL, "u", "s"),
                        plot_options = list(),
                        ...) {
  if (!inherits(x, "haz2ts")) stop("'x' must be a 'haz2ts' object")

  which_plot <- match.arg(which_plot)

  if (which_plot == "covariates" & is.null(x$optimal_model$beta)) {
    stop("Covariates plot required but x does not have covariates' parameters.")
  }

  if ((!is.null(where_slices)) & is.null(direction)) {
    stop("Please provide a direction for the cutting points - see argument 'direction'")
  }


  u <- s <- NULL
  # ---- Options for plotting ----
  opts <- list(
    loghazard = FALSE,
    log10hazard = FALSE,
    cut_extrapolated = TRUE,
    rectangular_grid = TRUE,
    original = FALSE,
    tmax = NULL,
    midpoints = FALSE,
    surv_slices = FALSE,
    cumhaz_slices = FALSE,
    col_palette = NULL,
    n_shades = 50,
    breaks = NULL,
    show_legend = TRUE,
    main = NULL,
    xlab = NULL,
    ylab = NULL,
    xlim = NULL,
    ylim = NULL,
    contour_lines = TRUE,
    contour_col = NULL,
    contour_cex = .8,
    contour_nlev = 10,
    cex_main = 1.2,
    cex_lab = 1,
    HR = FALSE,
    symmetric_CI = TRUE,
    confidence = .95,
    col_beta = "blue",
    pch = 20,
    lwd = 2,
    lty = 1
  )

  Nopts <- names(opts)
  namesPO <- names(plot_options)

  opts[namesPO] <- plot_options
  if (length(namesPO[!namesPO %in% Nopts]) > 0) {
    warning("Undefined entries in `plot_options`. Default settings are used.\n")
    warning(
      "Undefined keyword(s): ",
      paste(namesPO[!namesPO %in% Nopts], collapse = ", ")
    )
  }

  # ---- Plot covariates (then exit) ----
  if (which_plot == "covariates") {
    plt <- covariates_plot(x,
      confidence_lev = opts$confidence,
      plot_options = list(
        HR = opts$HR,
        symmetric_CI = opts$symmetric_CI,
        main = opts$main,
        ylab = opts$ylab,
        col_beta = opts$col_beta,
        pch = opts$pch,
        cex_main = opts$cex_main,
        cex_lab = opts$cex_lab
      )
    )
    return(invisible(plt))
  }

  # ---- Plot (log-) hazard or SEs surface ----
  # Recreate grid from estimation B-splines
  if (is.null(plot_grid)) {
    Bbases <- x$optimal_model$Bbases
    midu <- attributes(Bbases$Bu)$x
    mids <- attributes(Bbases$Bs)$x
    du <- midu[2] - midu[1]
    ds <- mids[2] - mids[1]
    intu <- midu + du / 2
    intu <- c(intu[1] - du / 2, intu)
    umin <- attributes(Bbases$Bu)$xl
    umax <- attributes(Bbases$Bu)$xr
    smin <- attributes(Bbases$Bs)$xl
    smax <- attributes(Bbases$Bs)$xr
    ints <- mids + ds / 2
    ints <- c(ints[1] - ds / 2, ints)

    plot_grid <- list(
      c("umin" = umin, "umax" = umax, "du" = du),
      c("smin" = smin, "smax" = smax, "ds" = ds)
    )
  }

  # ---- Get (baseline) (log-)hazard and hazard ratios if needed ----
  if (which_plot %in% c("hazard", "slices", "SE")) {
    hazard_SE <- get_hazard_2d(
      fitted_model = x,
      plot_grid = plot_grid,
      where_slices = where_slices,
      direction = direction,
      tmax = opts$tmax,
      midpoints = opts$midpoints
    )
    new_grid <- hazard_SE$new_plot_grid

    if (which_plot %in% c("hazard", "slices")) {
      if (opts$loghazard == TRUE) {
        to_plot <- hazard_SE$loghazard
      } else {
        if (opts$log10hazard) {
          to_plot <- hazard_SE$log10hazard
        } else {
          to_plot <- hazard_SE$hazard
        }
      }
    }

    if (which_plot == "SE") {
      if (opts$loghazard == TRUE) {
        to_plot <- hazard_SE$SE_loghazard
      } else {
        if (opts$log10hazard) {
          to_plot <- hazard_SE$SE_log10hazard
        } else {
          to_plot <- hazard_SE$SE_hazard
        }
      }
    }
  }
  if (which_plot == "survival") {
    surv <- surv2ts(
      fitted_model = x,
      plot_grid = plot_grid,
      midpoints = opts$midpoints,
      where_slices = where_slices,
      direction = direction,
      tmax = opts$tmax
    )
    new_grid <- attr(surv, "plot_grid")
    to_plot <- surv$Surv2ts
  }

  if (which_plot == "cumhaz") {
    CumHaz <- cumhaz2ts(
      fitted_model = x,
      plot_grid = plot_grid,
      midpoints = opts$midpoints,
      where_slices = where_slices,
      direction = direction,
      tmax = opts$tmax
    )
    to_plot <- CumHaz$CumHaz
    new_grid <- CumHaz$Haz$new_plot_grid
  }

  if (is.null(opts$tmax)) {
    opts$tmax <- new_grid$umax + new_grid$smax
  }

  # ---- Cut extrapolated hazard ----
  if (opts$cut_extrapolated) {
    cut <- matrix(NA, nrow(to_plot), ncol(to_plot))
    for (row in 1:nrow(to_plot)) {
      for (col in 1:ncol(to_plot)) {
        cut[row, col] <- ifelse((new_grid$ints[col] + new_grid$intu[row] > opts$tmax + new_grid$du) &
          (((new_grid$smax - new_grid$ints[col]) / (new_grid$umax - new_grid$intu[row]) >= -1) |
            ((new_grid$smax - new_grid$ints[col]) / (new_grid$umax - new_grid$intu[row]) == -Inf)),
        NA, 1
        )
      }
    }

    to_plot <- to_plot * cut

    # adjust legend breaks to match cut surface
    if (which_plot %in% c("hazard", "SE", "cumhaz")) {
      if (is.null(opts$breaks)) {
        K <- (max(to_plot, na.rm = T) - min(to_plot, na.rm = T)) / (opts$n_shades + 1)
        opts$breaks <- seq(min(to_plot, na.rm = T),
          min(to_plot, na.rm = T) + K * (opts$n_shades + 1),
          length = (opts$n_shades + 1)
        )
      }
    }
  }


  # ---- If surface plot, create grid of parallelograms corners for plotting ----
  if (opts$original) {
    if (!opts$rectangular_grid) {
      inty <- expand.grid(new_grid$intu, new_grid$ints)
      intx <- utot(inty$Var1, inty$Var2)

      X1 <- matrix(intx$t,
        nrow = length(new_grid$intu),
        ncol = length(new_grid$ints)
      )
      X2 <- matrix(intx$s,
        nrow = length(new_grid$intu),
        ncol = length(new_grid$ints)
      )
    } else {
      #  ------ transform to (t,s)-plane --------
      grid_us <- expand.grid(u = new_grid$intu, s = new_grid$ints)
      grid_us$t <- grid_us$u + grid_us$s
      grid_us$to_plot <- as.vector(to_plot)

      t <- unique(grid_us$t)

      intt <- t[t <= opts$tmax]

      grid_ts <- expand.grid(t = intt, s = new_grid$ints)
      plotgrid_ts <- merge(grid_ts, grid_us, all.x = TRUE)

      to_plot_v <- plotgrid_ts$to_plot
      dim(to_plot_v) <- c(length(new_grid$ints), length(intt))
      to_plot <- t(to_plot_v)

      X1 <- sort(intt)
      X2 <- sort(new_grid$ints)
    }
  } else {
    X1 <- new_grid$intu
    X2 <- new_grid$ints
  }

  # ---- If slices, organize on grid and select only values where slices are ----
  if (which_plot == "slices" | opts$surv_slices | opts$cumhaz_slices) {
    if (direction == "s") {
      grid_us <- expand.grid(u = new_grid$intu, s = new_grid$ints)
      grid_us$t <- grid_us$u + grid_us$s
      grid_us$to_plot <- as.vector(to_plot)
      onlyslic <- subset(grid_us, s %in% where_slices)
      grid_ts <- expand.grid(s = where_slices, t = unique(sort(onlyslic$t)))
      final_grid <- merge(grid_ts, onlyslic, all.x = T)
      to_plot_v <- final_grid$to_plot
      dim(to_plot_v) <- c(length(unique(sort(final_grid$t))), length(where_slices))
      to_plot <- to_plot_v
      X1 <- unique(sort(final_grid$t))
    } else {
      if (is.null(where_slices)) stop("Please provide location for the cut-points over `u` (where_slices)")
      grid_us <- expand.grid(u = new_grid$intu, s = new_grid$ints)
      grid_us$to_plot <- c(to_plot)
      onlyslic <- subset(grid_us, u %in% where_slices)
      to_plot_v <- onlyslic$to_plot
      dim(to_plot_v) <- c(length(where_slices), length(new_grid$ints))
      to_plot <- to_plot_v
    }
  }

  # ---- Plot (log-)hazard ----
  if (which_plot %in% c("hazard", "survival", "cumhaz") &
    !(opts$surv_slices) & !(opts$cumhaz_slices)) {
    plt <- imageplot_2ts(
      x = X1, y = X2, z = to_plot,
      plot_options = list(
        loghazard = opts$loghazard,
        log10hazard = opts$log10hazard,
        original = opts$original,
        rectangular_grid = opts$rectangular_grid,
        col_palette = opts$col_palette,
        n_shades = opts$n_shades,
        breaks = opts$breaks,
        show_legend = opts$show_legend,
        tmax = opts$tmax,
        main = opts$main,
        xlab = opts$xlab,
        ylab = opts$ylab,
        xlim = opts$xlim,
        ylim = opts$ylim,
        contour_lines = opts$contour_lines,
        contour_col = opts$contour_col,
        contour_cex = opts$contour_cex,
        contour_nlev = opts$contour_nlev,
        cex_main = opts$cex_main,
        cex_lab = opts$cex_lab
      )
    )
    return(invisible(plt))
  }

  # ---- Plot SEs ----
  if (which_plot == "SE") {
    plt <- imageplot_SE(
      x = X1, y = X2, z = to_plot,
      plot_options = list(
        loghazard = opts$loghazard,
        log10hazard = opts$log10hazard,
        original = opts$original,
        rectangular_grid = opts$rectangular_grid,
        col_palette = opts$col_palette,
        n_shades = opts$n_shades,
        breaks = opts$breaks,
        show_legend = opts$show_legend,
        tmax = opts$tmax,
        main = opts$main,
        xlab = opts$xlab,
        ylab = opts$ylab,
        xlim = opts$xlim,
        ylim = opts$ylim,
        contour_lines = opts$contour_lines,
        contour_col = opts$contour_col,
        contour_cex = opts$contour_cex,
        contour_nlev = opts$contour_nlev,
        cex_main = opts$cex_main,
        cex_lab = opts$cex_lab
      )
    )
    return(invisible(plt))
  }

  # ---- Plot slices ----
  if (which_plot == "slices") {
    if (direction == "u") x <- new_grid$ints else x <- X1
    plt <- plot_slices(
      x = x,
      y = to_plot,
      direction = direction,
      plot_options = list(
        loghazard = opts$loghazard,
        log10hazard = opts$log10hazard,
        col_palette = opts$col_palette,
        n_shades = length(where_slices),
        main = opts$main,
        xlab = opts$xlab,
        ylab = opts$ylab,
        xlim = opts$xlim,
        ylim = opts$ylim,
        cex_main = opts$cex_main,
        cex_lab = opts$cex_lab,
        lwd = opts$lwd
      )
    )
  }

  # ---- Plot survival / cumulative hazard slices ----
  if ((which_plot == "survival" & opts$surv_slices) |
    (which_plot == "cumhaz" & opts$cumhaz_slices)) {
    x <- new_grid$ints
    plt <- plot_slices(
      x = x,
      y = to_plot,
      direction = direction,
      plot_options = list(
        loghazard = opts$loghazard,
        log10hazard = opts$log10hazard,
        col_palette = opts$col_palette,
        n_shades = length(where_slices),
        main = opts$main,
        xlab = opts$xlab,
        ylab = opts$ylab,
        xlim = opts$xlim,
        ylim = opts$ylim,
        cex_main = opts$cex_main,
        cex_lab = opts$cex_lab,
        lwd = opts$lwd
      )
    )
  }
}
