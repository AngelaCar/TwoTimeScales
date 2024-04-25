#' Fit a smooth hazard model with two time scales
#'
#' @description `fit2ts()` fits a smooth hazard model with two time scales.
#'
#'   Two methods are implemented for the search of the optimal smoothing
#'   parameters (and therefore optimal model): a numerical optimization of the
#'   AIC or BIC of the model and a search for the minimum AIC or BIC of the
#'   model over a grid of `log_10` values for the smoothing parameters.
#'   Construction of the B-splines bases and of the penalty matrix is
#'   incorporated within the function. If a matrix of covariates is provided,
#'   the function will estimate a model with covariates.
#'
#' @param data2ts (optional) an object of class created by the function
#'   `prepare_data()`. Proving this input is the easiest way to use the function
#'   `fit2ts`. However, the user can also provide the input data together with
#'   a list of bins, as explained by the following parameters' descriptions.
#' @inheritParams grid_search_2d
#' @param bins a list with the specification for the bins. This is created by
#'   the function `prepare_data`. If a list prepared externally from such function
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
#'   `"ucminf"` (default) for the numerical optimization of the AIC (or BIC)
#'   and `"grid_search"` for a grid search of the minimum AIC (or BIC)
#'   over a grid of `log_10(rho_u)` and `log_10(rho_s)` values.
#' @param lrho A vector of two elements if `optim_method == "ucminf"`.
#'   Default is `c(0,0)`. A list of two vectors of values for `log_10(rho_u)`
#'   and `log_10(rho_s)` if `optim_method == "grid_search"`. In the latter case,
#'   if a list with two vectors is not provided, a default sequence of
#'   values is used for both `log_10(rho_u)` and `log_10(rho_s)`.
#'
#' @return An object of class `h2tsfit` with the following elements:
#'   * `optimal_model` A list with :
#'     * `Alpha` The matrix of estimated P-splines coefficients of dimension
#'       cu by cs.
#'     * `Cov_alpha` The variance-covariance matrix of the `Alpha` coefficients,
#'       of dimension cucs by cucs.
#'     * `beta` The vector of length p of estimated covariates coefficients
#'        (if model with covariates).
#'     * `Cov_beta` The variance-covariance matrix of the `beta` coefficients,
#'       of dimension p by p (if model with covariates).
#'     * `SE_beta` The vector of length p of estimated Standard Errors for the `beta`
#'       coefficients (if model with covariates)..
#'     * `Eta` or `Eta0` The matrix of values of the (baseline) linear predictor
#'       (log-hazard) of dimension nu by ns.
#'     * `H` The hat-matrix.
#'     * `deviance` The deviance.
#'     * `ed` The effective dimension of the model.
#'     * `aic` The value of the AIC.
#'     * `bic` The value of the BIC.
#'     * `Bbases` a list with the B-spline bases `Bu` and `Bs`
#'   * `optimal_logrho` A vector with the optimal values of `log10(rho_u)` and
#'     `log10(rho_s)`.
#'   * `P_optimal` The optimal penalty matrix P.
#'   * `AIC` (if `par_gridsearch$return_aic == TRUE`) The matrix of AIC values.
#'   * `BIC` (if `par_gridsearch$return_bic == TRUE`) The matrix of BIC values.
#'
#' @import JOPS
#' @importFrom grDevices grey.colors
#' @export
fit2ts <- function(data2ts = NULL,
                   Y = NULL, R = NULL, Z = NULL,
                   bins = NULL,
                   Bbases_spec = list(),
                   pord = 2,
                   optim_method = c("ucminf", "grid_search"),
                   optim_criterion = c("aic", "bic"),
                   lrho = c(0, 0),
                   Wprior = NULL,
                   ridge = 0,
                   control_algorithm = list(),
                   par_gridsearch = list()) {
  # ---- Check arguments ----
  optim_method <- match.arg(optim_method)
  optim_criterion <- match.arg(optim_criterion)

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
    monitor_ev = FALSE
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
    R = R, Y = Y,
    Z = Z,
    Wprior = Wprior
  )

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
  Bu <- JOPS::bbase(x = bins$midu, nseg = Bbases$nseg_u,
                    xl = Bbases$min_u, xr = Bbases$max_u, bdeg = Bbases$bdeg)
  nbu <- ncol(Bu)
  Bs <- JOPS::bbase(x = bins$mids, nseg = Bbases$nseg_s,
                    xl = Bbases$min_s, xr = Bbases$max_s, bdeg = Bbases$bdeg)
  nbs <- ncol(Bs)

  # ---- Penalty parameters ----
  Du <- diff(diag(nbu), diff = pord)
  Iu <- diag(nbu)
  Ds <- diff(diag(nbs), diff = pord)
  Is <- diag(nbs)

  # ---- Optimization ----

  if (optim_method == "grid_search") {
    optimal_model <- grid_search_2d(
      lru = lru, lrs = lrs,
      R = R, Y = Y, Z = Z,
      Bu = Bu, Bs = Bs,
      Iu = Iu, Is = Is,
      Du = Du, Ds = Ds,
      Wprior = Wprior,
      ridge = ridge,
      optim_criterion = optim_criterion,
      control_algorithm = con,
      par_gridsearch = gsp
    )
  }
  if (optim_method == "ucminf") {
    optimal_model <- fit2tsmodel_ucminf(
      Y, R,
      Z = Z,
      optim_criterion = optim_criterion,
      lrho = lrho,
      Bu = Bu, Bs = Bs,
      Iu = Iu, Is = Is,
      Du = Du, Ds = Ds,
      Wprior = Wprior,
      ridge = ridge,
      control_algorithm = con
    )
  }

  # ---- Save results in list and return list ----
  results <- optimal_model
  #class(results) <- "haz2ts"

  return(results)
}
