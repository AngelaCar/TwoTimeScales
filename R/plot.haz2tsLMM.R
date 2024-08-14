#' Plot method for a haz2ts object.
#'
#' @description  `plot.haz2tsLMM()` is the plot method for objects of class `haz2tsLMM`.
#'  It produces plots of the fitted model with two time scales (see [fit2ts()]),
#'  fitted via LMMsolver. It only produces a plot over the `u` and `s` axes.
#' @param x The output of the function `fit2ts`. This is an object of
#'   class `"haz2tsLMM"`.
#' @param which_plot The type of plot required. Can be one of `"hazard"`
#'   (default), `"covariates"` or `"SE"`.
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
#' @param plot_options A list with all possible options for any of the plots:
#'   * `loghazard` A Boolean. Default is `FALSE`. If `FALSE` the function
#'     returns a plot of the hazard surface, if `TRUE` the function returns
#'     a plot of the log-hazard surface.
#'   * `log10hazard` A Boolean. Default is `FALSE`. If `TRUE`,
#'     then a log_10 hazard surface is plotted.
#'   * `cut_extrapolated` A Boolean. Default is `TRUE`. Cuts away the
#'     extrapolated area of the (log-)hazard surface before plotting.
#'   * `tmax` The maximum value of `t` that should be plotted.
#'   * `col_palette` A function defining the color palette. The default palette
#'     is `viridis::rev(plasma())`.
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
#'
#' @return A plot of the fitted model.
#'
#' @importFrom stats qnorm
#'
#' @export

plot.haz2tsLMM <- function(x,
                        plot_grid = NULL,
                        which_plot = c("hazard", "covariates", "SE"),
                        plot_options = list(),
                        ...) {
  if (!inherits(x, "haz2tsLMM")) stop("'x' must be a 'haz2tsLMM' object")

  which_plot <- match.arg(which_plot)

  if (which_plot == "covariates" & x$covariates == "no") {
    stop("Covariates plot required but x does not have covariates' parameters.")
  }

  # ---- Options for plotting ----
  opts <- list(
    loghazard = FALSE,
    log10hazard = FALSE,
    cut_extrapolated = TRUE,
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

  # ---- Get (baseline) (log-)hazard and hazard ratios if needed ----
  if (which_plot != "covariates") { # the only case in which we don't need to call
    # Recreate grid from estimation B-splines
    if (is.null(plot_grid)) {
      intu <- x$optimal_model$splRes[[1]]$knots[[1]]
      umin <- attributes(intu)$xmin
      umax <- attributes(intu)$xmax
      intu <- intu[intu >= umin & intu <= umax]
      ints <- x$optimal_model$splRes[[1]]$knots[[2]]
      smin <- attributes(ints)$xmin
      smax <- attributes(ints)$xmax
      ints <- ints[ints >= smin & ints <= smax]
      du <- intu[2] - intu[1]
      ds <- ints[2] - ints[1]
      new_grid <- expand.grid(intu, ints)
      names(new_grid) <- c("u", "s")
    } else {
      umin <- plot_grid[[1]][1]
      umax <- plot_grid[[1]][2]
      du <- plot_grid[[1]][3]
      smin <- plot_grid[[2]][1]
      smax <- plot_grid[[2]][2]
      ds <- plot_grid[[2]][3]
      if (umin < attributes(x$optimal_model$splRes[[1]]$knots[[1]])$xmin) {
        umin <- attributes(x$optimal_model$splRes[[1]]$knots[[1]])$xmin
        warning("`umin` is smaller than the lower limit of the domain of Bu. Left boundary adjusted to  =  ", attributes(x$optimal_model$splRes[[1]]$knots[[1]])$xmin)
      }
      if (umax > attributes(x$optimal_model$splRes[[1]]$knots[[1]])$xmax) {
        umax <- attributes(x$optimal_model$splRes[[1]]$knots[[1]])$xmax
        warning("`umax` is larger than the upper limit of the domain of Bu. Right boundary adjusted to  =  ", attributes(x$optimal_model$splRes[[1]]$knots[[1]])$xmax)
      }
      K <- ceiling((umax - umin) / du)
      intu <- seq(umin, umin + K * du, by = du)

      if (smin < attributes(x$optimal_model$splRes[[1]]$knots[[2]])$xmin) {
        smin <- attributes(x$optimal_model$splRes[[1]]$knots[[2]])$xmin
        warning("`smin` is smaller than the lower limit of the domain of Bs. Left boundary adjusted to  =  ", attributes(x$optimal_model$splRes[[1]]$knots[[2]])$xmin)
      }
      if (smax > attributes(x$optimal_model$splRes[[1]]$knots[[2]])$xmax) {
        smax <- attributes(x$optimal_model$splRes[[1]]$knots[[2]])$xmax
        warning("`smax` is larger than the upper limit of the domain of Bs. Right boundary adjusted to  =  ", attributes(x$optimal_model$splRes[[1]]$knots[[2]])$xmax)
      }
      K <- ceiling((smax - smin) / ds)
      ints <- seq(smin, smin + K * ds, by = ds)

      new_grid <- expand.grid(intu, ints)
      names(new_grid) <- c("u", "s")
    }

    # get_hazard_2d
    trend2D <- obtainSmoothTrend(x$optimal_model,
                                 newdata = new_grid,  includeIntercept = TRUE)
    Haz <- matrix(trend2D$ypred, nrow = length(intu), ncol = length(ints))
    SE_Haz <- matrix(trend2D$se, nrow = length(intu), ncol = length(ints))
    logHaz <- log(Haz)
    SE_logHaz <- abs(1/Haz) * SE_Haz
    SE_log10Haz <- abs(1/(Haz * log(10))) * SE_Haz

    if (which_plot %in% c("hazard")) {
      if (opts$loghazard == TRUE) {
        to_plot <- log(Haz)
      } else {
        if(opts$log10hazard){
          to_plot <- log10(Haz)
        } else {
          to_plot <- Haz
        }
      }
    }

    if (which_plot == "SE") {
      if (opts$loghazard == TRUE) {
        to_plot <- SE_logHaz
      } else {
        if(opts$log10hazard){
          to_plot <- SE_log10Haz
        } else {
          to_plot <- SE_Haz
        }
      }
    }
  }

  if (is.null(opts$tmax)) {
    opts$tmax <- umax + smax
  }

  # ---- Cut extrapolated hazard ----
  if (opts$cut_extrapolated) {
    cut <- matrix(NA, nrow(to_plot), ncol(to_plot))
    for (row in 1:nrow(to_plot)) {
      for (col in 1:ncol(to_plot)) {
        cut[row, col] <- ifelse((ints[col] + intu[row] > opts$tmax + du) &
                                  (smax - ints[col]) / (umax - intu[row]) >= -1, NA, 1)
        }
    }

    to_plot <- to_plot * cut
  }


  # ---- If surface plot, create grid of parallelograms corners for plotting ----

    X1 <- intu
    X2 <- ints

  # ---- Plot (log-)hazard ----
  if (which_plot == "hazard") {
    plt <- imageplot_2ts(
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
      )
    )
    return(invisible(plt))
  }

}
