#' Plot method for a haz2tsLMM object.
#'
#' @description  `plot.haz2tsLMM()` is the plot method for objects of class `haz2tsLMM`.
#'  It produces plots of the fitted model with two time scales (see [fit2ts()]),
#'  fitted via LMMsolver. The two-dimensional plots are limited to the transformed
#'  plane, that is only plots over `u` and `s` axes are produced.
#' @param x The output of the function `fit2ts`. This is an object of
#'   class `"haz2tsLMM"`.
#' @param which_plot The type of plot required. Can be one of `"hazard"`
#'   (default), `"slices"`, `"covariates"` or `"SE"`.
#' @param plot_grid A list containing the parameters to build a new
#'   finer grid of intervals over `u` and `s` for plotting. This must be of the
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
#'   accordingly and a warning will be returned. While for objects of class 'haz2ts'
#'   this is an optional input, we strongly recommend to provide the plotting
#'   grid for object of class 'haz2tsLMM', given that evaluation of the B-splines
#'   works a bit differently in LMMsolver.
#' @param where_slices A vector of values for the cutting points of the desired
#'   slices of the surface. If `which_plot == "slices"`, please provide this
#'   argument.
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
#'     then a log10 hazard surface is plotted.
#'   * `cut_extrapolated` A Boolean. Default is `TRUE`. Cuts away the
#'     extrapolated area of the (log-)hazard surface before plotting.
#'   * `tmax` The maximum value of `t` that should be plotted.
#'   * `col_palette` A function defining the color palette. The default palette
#'     is `viridis::rev(plasma())`. Specifying the color palette as a function
#'     allows for greater flexibility than passing the palette as a vector.
#'     We provide an example on how to create a function from any color palette
#'     below.
#'   * `n_shades` The number of color shades to plot, default is 50.
#'   * `breaks` The vector of breaks for the color legend. If `n_shades` is provided,
#'     this should be of length `n_shades + 1`.
#'   * `show_legend` A Boolean. Default is `TRUE`. If `FALSE` no legend will be
#'     plotted, useful for multi-panel figures with common legend. Works only
#'     for plots on rectangular grid (i.e. transformed (u,s) plane).
#'   * `main` The title of the plot.
#'   * `xlab` The label of the first time axis (plotted on the x axis).
#'   * `ylab` The label of the second time axis (plotted on the y axis).
#'   * `xlim` A vector with two elements defining the limits of the time scale
#'      on the x axis.
#'   * `ylim` A vector with two elements defining the limits of the time scale
#'      on the y axis.
#'   * `contour_lines` A Boolean. Default is `FALSE`. If `TRUE` contour
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
#' @details The function `obtainSmoothTrend` from the R-package `LMMsolver` is
#'          used here. We refer the interested readers to https://biometris.github.io/LMMsolver/
#'          for more details on `LMMsolver` and its usage.
#'
#' @return A plot of the fitted model.
#'
#' @importFrom stats qnorm
#'
#' @export
#'
#' @examples
#' # Create a color pallete function from a RColorBrewer palette, using the function
#' # colorRampPalette from grDevices.
#' \dontrun{
#' mypal <- function(n) {
#'   colorRampPalette(RColorBrewer::brewer.pal(9, "YlGnBu"))(n)
#' }
#' # if mod_haz is a fitted model of class haz2ts, the following code will
#' # produce a cross-sections plot of the hazard over `s` for selected values
#' # of `u`, with the palette specified above
#'
#' plot(mod_haz,
#'   which_plot = "slices",
#'   where_slices = c(30, 60, 90, 180, 365, 1000, 2000),
#'   direction = "u",
#'   plot_options = list(
#'     col_palette = mypal,
#'     main = "Cross-sections of the hazard",
#'     xlab = "Time since recurrence",
#'     ylab = "Hazard"
#'   )
#' )
#' }
#'
plot.haz2tsLMM <- function(x,
                           plot_grid = NULL,
                           which_plot = c("hazard", "covariates", "SE", "slices"),
                           where_slices = NULL,
                           direction = c(NULL, "u", "s"),
                           plot_options = list(),
                           ...) {
  if (!inherits(x, "haz2tsLMM")) stop("'x' must be a 'haz2tsLMM' object")

  which_plot <- match.arg(which_plot)

  if (which_plot == "covariates" & x$covariates == "no") {
    stop("Covariates plot required but x does not have covariates' parameters.")
  }

  u <- s <- NULL

  # ---- Options for plotting ----
  opts <- list(
    loghazard = FALSE,
    log10hazard = FALSE,
    rectangular_grid = TRUE,
    cut_extrapolated = TRUE,
    original = FALSE,
    tmax = NULL,
    col_palette = NULL,
    n_shades = NULL,
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
    lwd = 2
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
  if (is.null(plot_grid)) {
    du <- attributes(x$optimal_model$splRes[[1]]$knots[[1]])$dx
    ds <- attributes(x$optimal_model$splRes[[1]]$knots[[2]])$dx
    umin <- attributes(x$optimal_model$splRes[[1]]$knots[[1]])$xmin
    umax <- attributes(x$optimal_model$splRes[[1]]$knots[[1]])$xmax
    smin <- attributes(x$optimal_model$splRes[[1]]$knots[[2]])$xmin
    smax <- attributes(x$optimal_model$splRes[[1]]$knots[[2]])$xmax
    plot_grid <- list(
      c("umin" = umin, "umax" = umax, "du" = du),
      c("smin" = smin, "smax" = smax, "ds" = ds)
    )
  }

  # ---- Get (baseline) (log-)hazard and hazard ratios if needed ----
  if (which_plot != "covariates") {

    hazard_SE <- get_hazard_2d_LMM(fitted_model = x,
                                   where_slices = where_slices,
                                   direction = direction,
                                   plot_grid = plot_grid,
                                   tmax = opts$tmax
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
        to_plot <- hazard_SE$SE_log10hazard
      } else {
        if (opts$log10hazard) {
          to_plot <- hazard_SE$SE_log10Haz
        } else {
          to_plot <- hazard_SE$SE_hazard
        }
      }
    }
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
                                  (new_grid$smax - new_grid$ints[col]) / (new_grid$umax - new_grid$intu[row]) >= -1, NA, 1)
      }
    }

    to_plot <- to_plot * cut
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
  if (which_plot == "slices") {
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
      grid_us <- expand.grid(u = new_grid$intu, s = new_grid$ints)
      grid_us$to_plot <- c(to_plot)
      onlyslic <- subset(grid_us, unique(u) %in% where_slices)
      to_plot_v <- onlyslic$to_plot
      dim(to_plot_v) <- c(length(where_slices), length(new_grid$ints))
      to_plot <- to_plot_v
    }
  }
  # ---- Plot (log-)hazard ----
  if (which_plot == "hazard") {
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
      ),
      ...
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
      ),
      ...
    )
    return(invisible(plt))
  }

  # ---- Plot slices ----
  if (which_plot == "slices") {
    if (direction == "u") x <- new_grid$ints else x <- X1
    plt <- plot_slices(
      x = x,
      y = to_plot_v,
      direction = direction,
      plot_options = list(
        loghazard = opts$loghazard,
        log10hazard = opts$log10hazard,
        col_palette = opts$col_palette,
        n_shades = opts$n_shades,
        main = opts$main,
        xlab = opts$xlab,
        ylab = opts$ylab,
        xlim = opts$xlim,
        ylim = opts$ylim,
        cex_main = opts$cex_main,
        cex_lab = opts$cex_lab,
        lwd = opts$lwd
      ),
      ...
    )
  }
}
