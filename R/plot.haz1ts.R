#' Plot method for a haz1ts object.
#'
#' @description  `plot.haz1ts()` is a plot method for objects of class `haz1ts`.
#  Plot the estimated smooth (log-) hazard with one time scale.
#'
#' @param x The output of the function `fit1ts`.
#' @param plot_grid (optional) A named vector containing the parameters to build a new
#'   grid of intervals over `s` for plotting the estimated hazard on a finer
#'   grid. This must be of the form: `plot_grid = c(smin, smax, ds)`,
#'   where `smin`, `smax` are the minimum and maximum values desired for the
#'   intervals over `s`, and `ds` is the distance between intervals over `s`. If
#'   not specified, the plotting is done using the same B-splines basis as for
#'   the estimation. The function will check if the parameters for the grid
#'   provided by the user are compatible with those originally used to construct
#'   the B-splines for estimating the model. If not, the grid will be adjusted
#'   accordingly and a warning will be returned.
#' @param which_plot The type of plot required. Can be one of `"hazard"`
#'   (default) or `"covariates"`.
#' @param plot_options A list with all possible options for any of the plots:
#'   * `loghazard` A Boolean. Default is `FALSE`. If `FALSE` the function
#'     returns a plot of the hazard curve, if `TRUE` the function returns
#'     a plot of the log-hazard curve.
#'   * `log10hazard` A Boolean. Default is `FALSE`. If `TRUE` it returns a plot
#'      of the log10-hazard curve.
#'   * `col` The color of the curve plotted. Default is `"black"`.
#'   * `add_CI` A Boolean. If `TRUE` (default) the confidence bands will be added.
#'   * `col_CI` The color for the confidence bands. The default is the same color
#'     of the curve, with a 50% transparancy level.
#'   * `main` The title of the plot.
#'   * `xlab` The label of the time axis (plotted on the x axis).
#'   * `ylab` The label of the y-axis (hazard, log-hazard or log10-hazard).
#'   * `xlim` A vector with two elements defining the limits of the time scale
#'      on the x axis.
#'   * `ylim` A vector with two elements defining the limits of function plotted
#'      on the y axis (hazard, log-hazard or log10-hazard).
#'   * `xmin` The minimum value on the x-axis.
#'   * `ymin` The minimum value on the y-axis.
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
#' @param \dots Further arguments to plot.
#'
#' @return A plot of the type required.
#'
#' @importFrom grDevices adjustcolor
#' @importFrom graphics abline lines par points polygon segments
#' @importFrom stats qnorm
#'
#' @export
#'
#' @examples
#'## preparing data - no covariates
#' dt1ts <- prepare_data(data = reccolon2ts,
#'                       s_in = "entrys",
#'                       s_out = "timesr",
#'                       events = "status",
#'                       ds = 180)
#'
#' ## fitting the model with fit1ts() - default options
#'
#' mod1 <- fit1ts(dt1ts)
#'
#' plot(mod1)
#'


plot.haz1ts <- function(x,
                        which_plot = c("hazard", "covariates"),
                        plot_grid = NULL,
                        plot_options = list(),
                        ...) {

  if (!inherits(x, "haz1ts")) stop("'x' must be a 'haz1ts' object")

  which_plot <- match.arg(which_plot)

  if (which_plot == "covariates" & is.null(x$optimal_model$beta))
    stop("Covariates plot required but x does not have covariates' parameters.")

  # ---- Options for plotting ----
  opts <- list(
    loghazard = FALSE,
    log10hazard = FALSE,
    col = "black",
    add_CI = TRUE,
    col_CI = NULL,
    main = NULL,
    xlab = "s",
    ylab = NULL,
    xlim = NULL,
    ylim = NULL,
    cex_main = 1.2,
    cex_lab = 1,
    HR = FALSE,
    symmetric_CI = TRUE,
    confidence = .95,
    col_beta = "blue",
    pch = 20
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
  if(which_plot == "covariates"){
    plt <- covariates_plot(x,
                           confidence_lev = opts$confidence,
                           plot_options = list(
                             HR = opts$HR,
                             symmetric_CI = opts$symmetric_CI,
                             main = opts$main,
                             ylab = opts$ylab,
                             ylim = opts$ylim,
                             col_beta = opts$col_beta,
                             pch = opts$pch,
                             cex_main = opts$cex_main,
                             cex_lab = opts$cex_lab
                           ))
    return(invisible(plt))

  }

  # ---- Plot (log-)hazard curve ----
  # Recreate grid from estimation B-splines
  if (is.null(plot_grid)) { # recreate intervals from B-splines bases
    Bbases <- x$optimal_model$Bbases
    mids <- attributes(Bbases$Bs)$x
    ds <- mids[2] - mids[1]
    ints <- mids + ds / 2
    ints <- c(ints[1] - ds, ints)
    smin <- attributes(Bbases$Bs)$xl
    smax <- attributes(Bbases$Bs)$xr
    plot_grid <- c("smin" = smin, "smax" = smax, "ds" = ds)
  }

  # ---- Z-value level for confidence intervals ----
  cialp <- (1-opts$confidence)/2
  zval <- abs(qnorm(cialp))

  # ---- Get (baseline) (log-)hazard and hazard ratios if needed ----
  if(which_plot == "hazard"){
    hazard_SE <- get_hazard_1d(x, plot_grid)
    new_grid <- hazard_SE$new_plot_grid
    if(is.null(opts$xlim)) opts$xlim <- c(new_grid$smin, new_grid$smax)

    if(opts$loghazard){
      to_plot <- hazard_SE$loghazard
      lci <- hazard_SE$loghazard - zval * hazard_SE$SE_loghazard
      uci <- hazard_SE$loghazard + zval * hazard_SE$SE_loghazard

      if(is.null(opts$main)) opts$main <- "log-hazard"
      if(is.null(opts$ylab)) opts$ylab <- "log-hazard"
    } else {
      if(opts$log10hazard) {
        to_plot <- hazard_SE$log10hazard
        lci <- hazard_SE$log10hazard - zval * hazard_SE$SE_log10hazard
        uci <- hazard_SE$log10hazard + zval * hazard_SE$SE_log10hazard

        if(is.null(opts$main)) opts$main <- "log10-hazard"
        if(is.null(opts$ylab)) opts$ylab <- "log10-hazard"
      } else {
      to_plot <- hazard_SE$hazard
      lci <- exp(hazard_SE$loghazard - zval * hazard_SE$SE_loghazard)
      uci <- exp(hazard_SE$loghazard + zval * hazard_SE$SE_loghazard)
      if(is.null(opts$main)) opts$main <- "hazard"
      if(is.null(opts$ylab)) opts$ylab <- "hazard"
        }
      }

    if(is.null(opts$col_CI)) col_CI <- adjustcolor(opts$col, alpha.f = 0.5)

    if(is.null(opts$ylim)) {
      if(opts$add_CI){
        opts$ylim <- c(min(lci, na.rm = T), max(uci, na.rm = T))
      } else {
        opts$ylim <- c(min(to_plot), max(to_plot))
      }
    }
  }


  # ---- Plot (log-) hazard ----
  if(which_plot == "hazard"){
    plot(new_grid$ints,
         to_plot,
         type = "l",
         #lwd = 2,
         xlim = opts$xlim,
         ylim = opts$ylim,
         col = opts$col,
         main = opts$main,
         cex.main = opts$cex_main,
         xlab = opts$xlab,
         ylab = opts$ylab,
         cex.lab = opts$cex_lab,
         axes = F,
         ...
    )
    axis(1, pos = min(opts$ylim))
    axis(2, pos = new_grid$smin)
    if (opts$add_CI) {
      lines(new_grid$ints, lci,
            lwd = 1,
            col = col_CI
      )
      lines(new_grid$ints, uci,
            lwd = 1,
            col = col_CI
      )
      polygon(c(rev(new_grid$ints), new_grid$ints), c(rev(lci), uci),
              col = col_CI, border = NA
      )

    }
  }
}

