#' Plot of the covariates' effects
#'
#' @description `covariates_plot()` produces a plot of the covariates' effects (\eqn{\hat\beta})
#'   with confidence intervals, or of the Hazard Ratios (\eqn{\exp(\hat\beta)}) with confidence intervals.
#'
#' @param fitted_model A list returned by the function `fit2ts` or `fit1ts`.
#' @param confidence_lev The level of confidence for the CIs. Default is 0.95 (\eqn{\alpha
#'   = 0.05}).
#' @param plot_options A list of options for the plot:
#'  * `HR` A Boolean. If `TRUE` the HRs with their CIs will be plotted.
#'    Default is `FALSE` (plot the `beta` with their CIs).
#'  * `symmetric_ci` A Boolean. Default is `TRUE`. If a plot of the HRs is
#'   required (`HR == TRUE`), then plot symmetrical Confidence Intervals,
#'   based on the SEs for the HRs calculated by delta method.
#'   If `FALSE`, then CIs are obtained by exponentiating the CIs for the betas.
#'  * `main` The title of the plot.
#'  * `ylab` The label for the y-axis.
#'  * `ylim` A vector with two elements defining the limits for the y-axis.
#'  * `col_beta` The color for the plot of the covariates' effects.
#'  * `pch` The symbol for plotting the point estimates.
#'  * `cex_main` The magnification factor for the main of the plot.
#' @param \dots further arguments passed to plot()
#'
#' @return A plot of the covariates' effects. The different covariates are plotted
#'          on the x-axis, and on the y-axis the effects on the coefficient- or
#'          on the HR-scale are plotted. The main estimate is represented by a
#'          point and the CIs are added as vertical bars.
#' @importFrom stats qnorm coef
#' @importFrom graphics abline axis
#'
#' @examples
#' # Create some fake data - the bare minimum
#' id <- 1:20
#' u <- c(5.43, 3.25, 8.15, 5.53, 7.28, 6.61, 5.91, 4.94, 4.25, 3.86, 4.05, 6.86,
#'        4.94, 4.46, 2.14, 7.56, 5.55, 7.60, 6.46, 4.96)
#' s <- c(0.44, 4.89, 0.92, 1.81, 2.02, 1.55, 3.16, 6.36, 0.66, 2.02, 1.22, 3.96,
#'        7.07, 2.91, 3.38, 2.36, 1.74, 0.06, 5.76, 3.00)
#' ev <- c(1, 0, 0, 1, 0, 1, 0, 1, 1, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 1)
#' x1 <- c(0, 0, 0, 1, 0, 1, 1, 1, 1, 0, 0, 0, 1, 0, 1, 1, 0, 0, 1, 0)
#'
#' fakedata <- as.data.frame(cbind(id, u, s, ev, x1))
#' covs <- subset(fakedata, select = c("x1"))
#' fakedata2ts <- prepare_data(u = fakedata$u,
#'                             s_out = fakedata$s,
#'                             ev = fakedata$ev,
#'                             ds = .5,
#'                             individual = TRUE,
#'                             covs = covs)
#' # Fit a fake model - not optimal smoothing
#' fakemod <- fit2ts(fakedata2ts,
#'                   optim_method = "grid_search",
#'                   lrho = list(seq(1 ,1.5 ,.5),
#'                               seq(1 ,1.5 ,.5)))
#' # Covariates plot with default options
#' covariates_plot(fakemod)
#'
#' # Plot the hazard ratios instead
#' covariates_plot(fakemod,
#'                 plot_options = list(
#'                 HR = TRUE))
#'
#' # Change confidence level
#' covariates_plot(fakemod,
#'                 confidence_lev = .99)
#' @export
#'

covariates_plot <- function(fitted_model,
                            confidence_lev = .95,
                            plot_options = list(), ...){

  # ---- Get estimates ----
  if(inherits(fitted_model, "haz2ts") | inherits(fitted_model, "haz1ts")){
    HR_SE <- get_hr(fitted_model)
    namesb <- attributes(fitted_model$optimal_model$beta)$dimnames[[1]]

  } else {
    namesCov <- fitted_model$optimal_model$term.labels.f[-c(1,length(fitted_model$optimal_model$term.labels.f))]
    coefLMM <- coef(fitted_model$optimal_model, se = T)[namesCov]
    coeftab <- matrix(0, length(namesCov), 3)
    for(ind in 1:length(namesCov)){
      coeftab[ind,2] <- round(as.numeric(coefLMM[[ind]]$value), 4)
      coeftab[ind,3] <- as.numeric(coefLMM[[ind]]$se)
    }
    coeftab <- as.data.frame(coeftab)
    #coeftab[,1] <- namesCov
    colnames(coeftab) <- c("coef", "beta", "SE_beta")
    HR_SE <- list(
      "beta" = coeftab$beta,
      "SE_beta" = coeftab$SE_beta
    )
    HR_SE$HR <- exp(HR_SE$beta)
    HR_SE$SE_HR <- HR_SE$HR * HR_SE$SE_beta
    namesb <- namesCov
  }
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
    cex.lab = opts$cex_lab,
    ...
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
