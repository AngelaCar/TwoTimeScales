#' Plot slices of the (log-) hazard
#'
#' @description `plot_slices()` plots slices of the (log-)hazard with two time
#'   scales, at selected values of one of the two time dimensions.
#'
#'
#' @param x A vector of values for the x-axis. This is a vector of values over
#'   the axis opposite to the one where the sliced are cut.
#' @param y A matrix of (log-)hazard values.
#' @param direction Either `"u"` or `"s"`.
#' @param plot_options A list of options for the plot:
#'  * `loghazard` A Boolean. Default is `FALSE`. If `FALSE` the function
#'     returns a plot of the hazard surface, if `TRUE` the function returns
#'     a plot of the log-hazard surface.
#'  * `col_palette` A function defining the color palette. The default palette
#'    is `grDevices::gray.colors()`.
#'  * `main` The title of the plot.
#'  * `xlab` The label of the first time axis (plotted on the x axis).
#'  * `ylab` The label of the second time axis (plotted on the y axis).
#'  * `xlim` A vector with two elements defining the limits of the time scale
#'     on the x axis.
#'  * `ylim` A vector with two elements defining the limits of the time scale
#'    on the y axis.
#'  * `cex_main` The magnification to be used for the main title, default is `1.2`.
#'  * `cex_lab` The magnification to be used for the axis labels, default is `1`.
#'  * `lwd` The line's width. Default is `2`.
#'
#' @return A plot of the slices of the hazard cut at selected points.
#'
#' @importFrom grDevices gray.colors
#' @importFrom graphics matplot
#' @export
#'
plot_slices <- function(x, y,
                        direction,
                        plot_options = list()) {
  if (direction == "s") {
    to_plot <- y
  } else {
    to_plot <- t(y)
  }

  # ---- Set options for plotting ----
  opts <- list(
    loghazard = FALSE,
    col_palette = NULL,
    main = NULL,
    xlab = NULL,
    ylab = NULL,
    xlim = NULL,
    ylim = NULL,
    cex_main = 1.2,
    cex_lab = 1,
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

  # ---- Set color palette ----
  if (is.null(opts$col_palette)) {
    n_shades <- ncol(to_plot)
    col_palette <- grDevices::gray.colors(n_shades)
  } else {
    col_palette <- opts$col_palette(n_shades)
  }

  # ---- Title and labels ----
  if (is.null(opts$main)) opts$main <- ifelse(opts$loghazard, "log-hazard", "hazard")
  if (is.null(opts$xlab)) opts$xlab <- ifelse(direction == "s", "t", "s")
  if (is.null(opts$ylab)) opts$ylab <- ifelse(opts$loghazard, "log-hazard", "hazard")

  # ---- Axes limits ----
  if (is.null(opts$xlim)) opts$xlim <- c(min(unique(x)), max(unique(x)))
  if (is.null(opts$ylim)) opts$ylim <- c(min(to_plot, na.rm = T), max(to_plot, na.rm = T))

  # ---- Plot ----
  matplot(x, to_plot,
    type = "l",
    lwd = opts$lwd,
    lty = 1,
    col = col_palette,
    main = opts$main,
    cex.main = opts$cex_main,
    xlab = opts$xlab,
    ylab = opts$ylab,
    cex.lab = opts$cex_lab,
    xlim = opts$xlim,
    ylim = opts$ylim
  )
}
