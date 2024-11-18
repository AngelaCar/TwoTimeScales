#' Image Plot of Standard Errors for the 2ts hazard
#'
#' `imageplot_SE()` plots an image of the SEs of the two time scales hazard with
#'   contour lines.
#'
#'
#' @param x The coordinates for the x-axis. This is a vector of intervals
#'   over the `u` axis (default), a matrix with the corner points of
#'   the parallelograms over the `t` time scale, or a vector of intervals for the `t`
#'   time scale.
#' @param y The coordinates for the y-axis. This is a vector of intervals
#'   over the `s` time scale (default), or a matrix with the corner points of
#'   the parallelograms over the `s` time scale.
#' @param z The values of the surface to plot, organized in a matrix with
#'   dimensions compatible with those of `x` and `y`. These can be the SEs for
#'   the hazard, the SEs for the log-hazard or the SEs for the log10-hazard.
#' @param plot_options A list of options for the plot:
#'  * `loghazard` A Boolean. Default is `FALSE`. If `FALSE` the function
#'     returns a plot of the standard errors of the hazard surface,
#'     if `TRUE` the function returns a plot of the standard errors of the
#'     log-hazard surface.
#'  * `log10hazard` A Boolean. Default is `FALSE`. If `TRUE`, the function
#'     returns a plot of the standard errors of the log10-hazard surface
#'  * `original` A Boolean. Default is `TRUE`. Plot the (log-)hazard in the
#'    (t,s)-plane. If `FALSE`, the (log-)hazard will be plotted in the (u,s)-plane.
#'  * `rectangular_grid` A Boolean. Default is `FALSE`. If `TRUE`, a
#'     rectangular grid is used for plotting also in the (t,s)-plane as opposed
#'     to the grid of parallelograms used as default in the (t,s)-plane.
#'  * `col_palette` A function defining the color palette. The default palette
#'    is `rev(colorspace::sequential_hcl(n = 50, "Red-Purple"))`.
#'  * `n_shades` The number of color shades to plot, default is 50.
#'  * `breaks` The vector of breaks for the color legend. If `n_shades` is provided,
#'    this should be of length `n_shades + 1`. Otherwise, `n_shades` will be
#'    recalculated accordingly.
#'   * `show_legend` A Boolean. Default is `TRUE`. If `FALSE` no legend will be
#'     plotted, useful for multi-panel figures with common legend. Works only
#'     for plots on rectangular grid!
#'  * `tmax` The maximum value of `t` that should be plotted.
#'  * `main` The title of the plot.
#'  * `xlab` The label of the first time axis (plotted on the x axis).
#'  * `ylab` The label of the second time axis (plotted on the y axis).
#'  * `xlim` A vector with two elements defining the limits of the time scale
#'     on the x axis.
#'  * `ylim` A vector with two elements defining the limits of the time scale
#'    on the y axis.
#'  * `cex_main` The magnification to be used for the main title, default is `1.2`.
#'  * `cex_lab` The magnification to be used for the axis labels, default is `1`.
#'  * `contour_lines` A Boolean. Default is `FALSE`. If `TRUE` white contour
#'    lines are added to the surfaces.
#'  * `contour_col` The color for the contour lines. Default is `white`.
#'  * `contour_cex` The magnification to be used for the contour lines.
#'    Default is `.8`.
#'  * `contour_nlev` The number of contour levels. Default is `10`.
#' @param \dots Further arguments to image.plot or image
#'
#' @return An image plot of the SEs for the (log-) hazard surface.
#'
#' @importFrom fields image.plot
#' @importFrom graphics axis box contour image
#' @importFrom colorspace sequential_hcl
#'
#' @export
#'

imageplot_SE <- function(x, y, z,
                         plot_options = list(),
                         ...) {
  # ---- Set options for plotting ----
  opts <- list(
    loghazard = FALSE,
    log10hazard = FALSE,
    original = FALSE,
    rectangular_grid = TRUE,
    col_palette = NULL,
    n_shades = NULL,
    breaks = NULL,
    show_legend = TRUE,
    tmax = NULL,
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
    cex_lab = 1
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

  # ---- Set breaks and color palette ----
  if (is.null(opts$n_shades)) {opts$n_shades <- 50}
  if (is.null(opts$breaks)) {
    K <- (max(z, na.rm = T)-min(z, na.rm = T))/(opts$n_shades + 1)
    opts$breaks <- seq(min(z, na.rm = T), min(z, na.rm = T) + K*(opts$n_shades + 1),
                       length = (opts$n_shades + 1))
  } else {
    opts$n_shades <- length(opts$breaks) - 1
  }
  if (is.null(opts$col_palette)) {
    col_palette <- rev(colorspace::sequential_hcl(n = opts$n_shades, "Red-Purple"))
  } else {
    col_palette <- opts$col_palette(n = opts$n_shades)
  }
  if(is.null(opts$contour_col)) opts$contour_col <- "grey"

  # ---- Title and labels ----
  if (is.null(opts$main)) {
    if (opts$loghazard) {
      opts$main <- "SEs log-hazard"
    } else {
      if (opts$log10hazard) {
        opts$main <- "Ses log10-hazard"
      } else {
        opts$main <- "SEs hazard"
      }
    }
  }

  if (is.null(opts$xlab)) opts$xlab <- ifelse(opts$original, "t", "u")
  if (is.null(opts$ylab)) opts$ylab <- "s"

  # ---- Axes limits ----
  if (is.null(opts$xlim)) {
    if (is.null(opts$tmax)) {
      opts$xlim <- c(min(unique(x)), max(unique(x)))
    } else {
      if(!opts$original){
        opts$xlim <- c(min(unique(x)), max(unique(x)))
      } else {
        opts$xlim <- c(min(unique(x)), opts$tmax)
      }
    }
  }
  if(is.null(opts$ylim)) opts$ylim <- c(min(unique(y)), max(unique(y)))

  # ---- No contour if original and not rectangular_grid
  if(opts$original & opts$rectangular_grid) opts$contour_lines <- FALSE

  # ---- Plot ----
  plt <- {
    if(opts$show_legend){
      image.plot(
      x, y,
      z,
      xlim = opts$xlim,
      ylim = opts$ylim,
      col = col_palette,
      breaks = opts$breaks,
      midpoint = !(opts$rectangular_grid),
      main = opts$main,
      cex.main = opts$cex_main,
      xlab = opts$xlab,
      ylab = opts$ylab,
      cex.lab = opts$cex_lab,
      ...
    )} else {
      if(! opts$rectangular_grid) stop("Cannot plot without a legend using a non rectangular grid...")
      image(
        x, y,
        z,
        xlim = opts$xlim,
        ylim = opts$ylim,
        col = col_palette,
        breaks = opts$breaks,
        main = opts$main,
        cex.main = opts$cex_main,
        xlab = opts$xlab,
        ylab = opts$ylab,
        cex.lab = opts$cex_lab,
        ...
      )}
    if (opts$contour_lines) {
      contour(x, y,
              z,
              col = opts$contour_col,
              add = TRUE,
              labcex = opts$contour_cex,
              nlevels = opts$contour_nlev
      )
    }
    box()
  }

  return(plt)
}