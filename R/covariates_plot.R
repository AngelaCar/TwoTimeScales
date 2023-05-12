#' Plot of the covariates' effects
#'
#' @description `covariates_plot()` produces a plot of the covariates' effects
#'   with confidence intervals, or of the Hazard Ratios with confidence intervals.
#'
#' @param fitted_model A list returned by the function `fit2ts` or `fit1ts`.
#' @param confidence_lev The level of confidence for the CIs. Default is .95 (alpha
#'   = 0.05).
#' @param plot_options A list of options for the plot:
#'  * `HR` A Boolean. If `TRUE` the HRs with their CIs will be plotted.
#'    Default is `FALSE` (plot the `beta` with their CIs).
#'  * `symmetric_ci` A Boolean. Default is `TRUE`. If a plot of the HRs is
#'   required (`HR == TRUE`), then plot symmetrical Confidence Intervals,
#'   based on the SEs for the HRs calculated by delta method.
#'   If `FALSE`, then CIs are obtained by exponentiating the CIs for the betas.
#'  * `main` The title of the plot.
#'  * `ylab` The label of the second time axis (plotted on the y axis).
#'  * `ylim` A vector with two elements defining the limits of the time scale
#'      on the y axis.
#'  * `col_beta` The color for the plot of the covariates' effects.
#'  * `pch` The symbol for plotting the point estimates.
#'  * `cex_main` The magnification factor for the main of the plot.
#' @return A plot of the covariates' effects.
#' @importFrom stats qnorm
#' @importFrom graphics abline axis
#' @export
#'

covariates_plot <- function(fitted_model,
                            confidence_lev = .95,
                            plot_options = list()){

  # ---- Get estimates ----
  HR_SE <- get_hr(fitted_model)

  # ---- Set plotting options ----
  opts <- list(
    HR = FALSE,
    symmetric_CI = FALSE,
    main = NULL,
    ylab = NULL,
    ylim = NULL,
    col_beta = "blue",
    pch = 20,
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

  if(is.null(opts$main_cov)) {
    opts$main_cov <- ifelse(opts$HR, "hazard ratios", "betas")
  }
  if(is.null(opts$ylab)) {
    opts$ylab <- ifelse(opts$HR, "hazard ratios", "betas")
  }

  # ---- Z-value level for confidence intervals ----
  cialp <- (1-confidence_lev)/2
  zval <- abs(qnorm(cialp))

  # ---- what to plot ----
  if(opts$HR & opts$symmetric_CI){
    y <- HR_SE$HR
    se_y <- HR_SE$SE_HR
  } else {
    y <- HR_SE$beta
    se_y <- HR_SE$SE_beta
  }

  namesb <- attributes(fitted_model$optimal_model$beta)$dimnames[[1]]

  if(opts$HR & !(opts$symmetric_CI)){
    LCIs <- exp(y - zval * se_y)
    UCIs <- exp(y + zval * se_y)
    y <- exp(y)
  } else {
    LCIs <- y - zval * se_y
    UCIs <- y + zval * se_y
  }

  if(is.null(opts$ylim)) {
      opts$ylim <- c(min(LCIs, na.rm = T)+.005, max(UCIs, na.rm = T)+.005)
  }
  # ---- Plot ----
  p <- length(HR_SE$beta)
  x <- 1:p

  plt <- {
    plot(
    x, y,
    xlab = "",
    ylab = opts$ylab,
    xaxt = "n",
    xlim = c(0, p + 1),
    ylim = opts$ylim,
    pch = opts$pch,
    col = opts$col_beta,
    main = opts$main,
    cex.main = opts$cex_main,
    cex.lab = opts$cex_lab
  )
  abline(h = 0, col = "grey", lty = 2)
  segments(x, LCIs, x, UCIs,
           col = opts$col_beta
  )
  segments((x - .15), LCIs, (x + .15), LCIs,
           col = opts$col_beta
  )
  segments((x - .15), UCIs, (x + .15), UCIs,
           col = opts$col_beta
  )
  axis(1,
       at = x,
       tick = F,
       labels = namesb,
       las = 2
  )
  }
  return(plt)
}
