#' Fit a log-additive model over two time scales
#'
#' @description
#' Fits a log-additive model over two time scales. If covariates are included in
#' the model, `fitpgam()` will fit a proportional hazards model, where the baseline
#' hazard is a surface over the two time scales, where the effect of the time scales
#' is additive on the log-scale. The model is described in Carollo et al. (2025).
#'
#'
#' @param data2ts (optional) an object of class `"data2ts"` created by the function
#'   `prepare_data()`. Proving this input is the easiest way to use the function
#'   `fitpgam()`. However, the user can also provide the input data together with
#'   a list of bins, as explained by the following parameters' descriptions.
#'
#' @inheritParams grid_search_2d
#' @param bins a list with the specification for the bins. This is created by
#'   the function `prepare_data()`. If a list prepared externally from such function
#'   if provided, it should contain the following elements:
#'     * `bins_u` A vector of bins extremes for the time scale `u`.
#'     * `midu` A vector with the midpoints of the bins over `u`.
#'     * `nu` The number of bins over `u`.
#'     * `bins_s` A vector of bins extremes for the time scale `s`.
#'     * `mids` A vector with the midpoints of the bins over `s`.
#'     * `ns` The number of bins over `s`.
#' @param Bbases_spec A list with the specification for the B-splines basis
#'   with the following elements:
#'   * `bdeg` The degree of the B-splines basis. Default is 3 (for cubic B-splines).
#'   * `nseg_u` The number of segments for the B-splines over `u`. Default is 10.
#'   * `min_u` (optional) The lower limit of the domain of `Bu`.
#'     Default is `min(bins_u)`.
#'   * `max_u` (optional) The upper limit of the domain of `Bu`.
#'     Default is `max(bins_u)`.
#'   * `nseg_s` The number of segments for the B-splines over `s`. Default is 10.
#'   * `min_s` (optional) The lower limit of the domain of `Bs`.
#'     Default is `min(bins_s)`.
#'   * `max_s` (optional) The upper limit of the domain of `Bs`.
#'     Default is `max(bins_s)`.
#' @param pord The order of the penalty. Default is 2.
#' @param optim_method The method to be used for optimization:
#'   `"ucminf"` (default) for the numerical optimization of the AIC (or BIC),
#'    `"grid_search"` for a grid search of the minimum AIC (or BIC)
#'     over a grid of \eqn{\log_{10}(\rho_u)} and \eqn{\log_{10}(\rho_s)} values.
#' @param lrho A vector of two elements if `optim_method == "ucminf"`.
#'   Default is `c(0,0)`. A list of two vectors of values for \eqn{\log_{10}(\rho_u)}
#'   and \eqn{\log_{10}(\rho_s)} if `optim_method == "grid_search"`. In the latter case,
#'   if a list with two vectors is not provided, a default sequence of
#'   values is used for both \eqn{\log_{10}(\rho_u)} and \eqn{\log_{10}(\rho_s)}.
#'
#' @return An object of class `haz2tsPGAM`, that is a list with the following elements:
#'   * `optimal_model` A list with :
#'     * `alpha` The vector of estimated P-splines coefficients of length
#'       \eqn{c_u + c_s}.
#'     * `SE_alpha` The vector of estimated standard errors of the `alpha` coefficients,
#'       of dimension \eqn{c_u + c_s}.
#'     * `beta` (if covariates) The vector of length \eqn{p} of estimated covariates coefficients.
#'     * `se_beta` The vector of length \eqn{p} of estimated Standard Errors for the `beta`
#'        coefficients.
#'     * `eta0` (if covariates) The vector of values of the baseline linear predictor (log-hazard).
#'     * `eta` (if covariates) The vector of values of the baseline linear predictor
#'       (log-hazard) of length \eqn{n_u * n_s}.
#'     * `H` The hat-matrix.
#'     * `deviance` The deviance.
#'     * `ed` The (total) effective dimension of the model.
#'     * `aic` The value of the AIC.
#'     * `bic` The value of the BIC.
#'     * `Bbases` a list with the B-spline bases `Bu` and `Bs`
#'   * `optimal_logrho` A vector with the optimal values of \eqn{\log_{10}(\rho_u)} and
#'     \eqn{\log_{10}(\rho_s)}.
#'   * `P_optimal` The optimal penalty matrix P.
#'   * `EDu` The effective dimensions along the `u` axis.
#'   * `EDs` The effective dimensions along the `s` axis.
#'   * `AIC` (if `par_gridsearch$return_aic == TRUE`) The matrix of AIC values.
#'   * `BIC` (if `par_gridsearch$return_bic == TRUE`) The matrix of BIC values.
#'
#' @details The name P-GAM follows from the way this class of models is named in
#' Eilers and Marx (2021) Practical Smoothing. The Joys of P-splines.
#' In Carollo et al. (2025) this model is referred to as log-additive.
#'
#'   @references Carollo, A., Putter, H., Eilers, P. H. C., & Gampe, J. (2025).
#'   Analysis of Time-to-Event Data With Two Time Scales.
#'   An Application to Transitions out of Cohabitation.
#'   Sociological Methods & Research, 0(0). \doi{10.1177/00491241251374193}
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
#' fitpgam(fakedata2ts,
#'   optim_method = "grid_search",
#'   lrho = list(seq(1, 1.5, .5), seq(1, 1.5, .5))
#' )
fitpgam <- function(data2ts = NULL,
                    Y = NULL, R = NULL, Z = NULL,
                    bins = NULL,
                    Bbases_spec = list(),
                    pord = 2,
                    ridge = 1e-6,
                    lrho = c(0, 0),
                    optim_method = c("ucminf", "grid_search"),
                    optim_criterion = c("aic", "bic"),
                    control_algorithm = list(),
                    par_gridsearch = list()) {
  optim_method <- match.arg(optim_method)
  optim_criterion <- match.arg(optim_criterion)

  # ---- Check arguments ----
  if (is.null(data2ts)) {
    if (is.null(Y) | is.null(R)) {
      stop("Please provide either an object created by the function `prepare_data`,
           or an array of event counts `Y` and an array of exposure times `R`.")
    }
    if (is.null(bins)) stop("`bins_list` is missing. Please provide a list with the bins specifications.")
  } else {
    Y <- data2ts$bindata$Y
    R <- data2ts$bindata$R
    Z <- data2ts$bindata$Z
    bins <- data2ts$bins
  }

  # ---- Controls for iterative process ----
  con <- list(
    maxiter = 20,
    conv_crit = 1e-5,
    verbose = FALSE,
    monitor_ev = FALSE,
    xtol = 1e-5
  )
  Ncon <- names(con)
  namesCon <- names(control_algorithm)

  con[namesCon] <- control_algorithm
  if (length(namesCon[!namesCon %in% Ncon]) > 0) {
    warning("Undefined entries in control! Default settings are used.\n")
    warning(
      "Undefined keyword(s): ",
      paste(namesCon[!namesCon %in% Ncon], collapse = ", ")
    )
  }

  # ---- Controls for grid search, if selected as optimization method ----
  if (optim_method == "grid_search") {
    if (!is.list(lrho)) {
      lru <- seq(-3, 3, by = .25)
      lrs <- lru
      message(
        "Grid search method selected, but no list with grid values over `log_10(rho_u)` and `log_10(rho_s)` is provided.
              I will use a default grid: \n `log_10(rho_u)`: ",
        lru, "\n `log_10(rho_s)`: ", lrs
      )
    } else {
      lru <- lrho[[1]]
      lrs <- lrho[[2]]
    }
    gsp <- list(
      plot_aic = FALSE,
      plot_bic = FALSE,
      return_aic = TRUE,
      return_bic = TRUE,
      col = grey.colors(n = 10),
      plot_contour = FALSE,
      mark_optimal = FALSE,
      main_aic = "AIC grid",
      main_bic = "BIC grid"
    )
    Ngsp <- names(gsp)
    namesGsp <- names(par_gridsearch)

    gsp[namesGsp] <- par_gridsearch
    if (length(namesGsp[!namesGsp %in% Ngsp]) > 0) {
      warning("Undefined entries in parameters for the grid search. Default settings are used.\n")
      warning(
        "Undefined keyword(s): ",
        paste(namesGsp[!namesGsp %in% Ngsp], collapse = ", ")
      )
    }
  }

  # ---- Checks for dimensions of input data ----
  check_inputs(
    R = R,
    Y = Y,
    Z = Z
  )

  # ---- Reformat data for pgam model ----
  r <- c(R)
  y <- c(Y)
  midu <- data2ts$bins$midu
  mids <- data2ts$bins$mids
  grid <- expand.grid(midu, mids)
  data.grid <- cbind(r, y, grid)
  names(data.grid) <- c("r", "y", "midu", "mids")

  nu <- length(midu)
  ns <- length(mids)

  # ---- Bbases specification ----
  Bbases <- list(
    bdeg = 3,
    nseg_u = 10,
    nseg_s = 10,
    min_u = min(bins$bins_u),
    max_u = max(bins$bins_u),
    min_s = min(bins$bins_s),
    max_s = max(bins$bins_s)
  )

  NBbases <- names(Bbases)
  namesBbases <- names(Bbases_spec)

  Bbases[namesBbases] <- Bbases_spec
  if (length(namesBbases[!namesBbases %in% NBbases]) > 0) {
    warning("Undefined entries in Bbases_spec. Default settings are used.\n")
    warning(
      "Undefined keyword(s): ",
      paste(namesBbases[!namesBbases %in% NBbases], collapse = ", ")
    )
  }

  # ---- Construct B-splines ----
  Bu <- JOPS::bbase(
    x = data.grid$midu, nseg = Bbases$nseg_u,
    xl = Bbases$min_u, xr = Bbases$max_u, bdeg = Bbases$bdeg
  )
  nbu <- ncol(Bu)
  Bs <- JOPS::bbase(
    x = data.grid$mids, nseg = Bbases$nseg_s,
    xl = Bbases$min_s, xr = Bbases$max_s, bdeg = Bbases$bdeg
  )
  nbs <- ncol(Bs)

  B <- cbind(Bu, Bs)
  nb <- ncol(B)

  # ---- Penalty parameters ----
  Du <- diff(diag(nbu), diff = pord)
  Ds <- diff(diag(nbs), diff = pord)
  DutDu <- t(Du) %*% Du
  DstDs <- t(Ds) %*% Ds

  # ---- Optimization ----

  if (optim_method == "grid_search") {
    optimal_model <- grid_search_pgam(
      lru = lru, lrs = lrs,
      r = r, y = y, Z = Z,
      B = B, nb = nb, nbu = nbu, nbs = nbs,
      DutDu = DutDu, DstDs = DstDs,
      ridge = ridge,
      optim_criterion = optim_criterion,
      control_algorithm = con,
      par_gridsearch = gsp
    )
    results <- optimal_model
  }
  if (optim_method == "ucminf") {
    optimal_model <- fitpgammodel_ucminf(
      y = y, r = r,
      Z = Z,
      optim_criterion = optim_criterion,
      lrho = lrho,
      B = B, nb = nb, nbu = nbu, nbs = nbs,
      DutDu = DutDu, DstDs = DstDs,
      ridge = ridge,
      control_algorithm = con
    )
    results <- optimal_model
  }

  # re-write the B-splines argument
  results$optimal_model$Bbases <- list("Bu" = Bu, "Bs" = Bs)
  return(results)
}
