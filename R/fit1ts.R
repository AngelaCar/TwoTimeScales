#' Fit a smooth hazard model with one time scale
#'
#' @description `fit1ts()` fits a smooth hazard model with one time scale.
#'
#'   Two methods are implemented for the search of the optimal smoothing
#'   parameter (and therefore optimal model): a numerical optimization of the
#'   AIC or BIC of the model and a search for the minimum AIC or BIC of the
#'   model over a grid of `log_10` values for the smoothing parameter.
#'   Construction of the B-splines basis and of the penalty matrix is
#'   incorporated within the function. If a matrix of covariates is provided,
#'   the function will estimate a model with covariates.
#'
#' @param data1ts (optional) an object created by the function
#'   `prepare_data()`. Providing this input is the easiest way to use the function
#'   `fit1ts`. However, the user can also provide the input data together with
#'   a list of bins, as explained by the following parameters' descriptions.
#' @inheritParams grid_search_1d
#' @param bins a list with the specification for the bins. This is created by
#'   the function `prepare_data`. Alternatively, a list with the following elements
#'   can be provided:
#'     * `bins_s` is a vector of bins extremes for the time scale `s`.
#'     * `mids` is a vector with the midpoints of the bins over `s`.
#'     * `ns` is the number of bins over `s`.
#' @param Bbases_spec A list with the specification for the B-splines basis
#'   with the following elements:
#'   * `bdeg` The degree of the B-splines basis. Default is 3 (for cubic B-splines).
#'   * `nseg_s` The number of segments for the B-splines over `s`. Default is 10.
#'   * `min_s` (optional) The lower limit of the domain of `Bs`.
#'     Default is `min(bins_s)`.
#'   * `max_s` (optional) The upper limit of the domain of `Bs`.
#'     Default is `max(bins_s)`.
#' @param pord The order of the penalty. Default is 2.
#' @param ridge A ridge penalty parameter: default is 0.
#' @param optim_method The method to be used for optimization:
#'   `"ucminf"` (default) for the numerical optimization of the AIC (or BIC)
#'   and `"grid_search"` for a grid search of the minimum AIC (or BIC)
#'   over a grid of `log_10(rho_s)` values.
#' @param lrho A number if `optim_method == "ucminf"`, default is 0.
#'   A vector of values for `log_10(rho_s)` if `optim_method == "grid_search"`.
#'   In the latter case, if a vector is not provided, a default sequence of
#'   values is used for `log_10(rho_s)`.
#'
#' @return An object of class `haz1ts` with the following elements:
#'   * `optimal_model` A list with:
#'     * `alpha` The vector of estimated P-splines coefficients of length cs.
#'     * `SE_alpha` The vector of estimated Standard Errors for the `alpha` coefficients,
#'        of length cs.
#'     * `beta` The vector of estimated covariates coefficients of length p
#'       (if model with covariates).
#'     * `se_beta` The vector of estimated Standard Errors for the
#'       `beta` coefficients of length p (if model with covariates).
#'     * `eta` or `eta0`. The vector of values of the (baseline) linear predictor
#'       (log-hazard).
#'     * `H` The hat-matrix.
#'     * `Cov` The full variance-covariance matrix.
#'     * `deviance` The deviance.
#'     * `ed` The effective dimension of the model.
#'     * `aic` The value of the AIC.
#'     * `bic` The value of the BIC.
#'     * `Bbases` a list with the B-spline basis `Bs` (this is a list for
#'       compatibility with functions in 2d).
#'   * `optimal_logrho` The optimal value of `log10(rho_s)`.
#'   * `P_optimal` The optimal penalty matrix P.
#'   * `AIC` (if `par_gridsearch$return_aic == TRUE`) The vector of AIC values.
#'   * `BIC` (if `par_gridsearch$return_bic == TRUE`) The vector of BIC values.
#'
#' @import JOPS
#' @export
#'
#' @examples
#' \dontrun{
#' ## preparing data - no covariates
#' dt1ts <- prepare_data(s_in = reccolon2ts$entrys,
#'                         s_out = reccolon2ts$timesr,
#'                         events = reccolon2ts$status,
#'                         ds = 30)
#'
#' ## fitting the model with fit1ts() - default options
#'
#' mod1 <- fit1ts(dt1ts)
#'
#' ## preparing the data - covariates
#'
#' covs <- subset(reccolon2ts, select = c("rx", "node4", "sex"))
#' dt1ts_cov <- prepare_data(s_in = reccolon2ts$entrys,
#'                           s_out = reccolon2ts$timesr,
#'                           events = reccolon2ts$status,
#'                           ds = 30,
#'                           individual = TRUE,
#'                           covs = covs)
#'
#' ## fitting the model with fit1ts() - default options
#'
#' mod2 <- fit1ts(d)
#'
#' ##
#' }
#'
#'
#'
fit1ts <- function(data1ts = NULL,
                   y = NULL, r = NULL,
                   Z = NULL,
                   bins = NULL,
                   Bbases_spec = list(),
                   Wprior = NULL,
                   pord = 2,
                   optim_method = c("ucminf", "grid_search"),
                   optim_criterion = c("aic", "bic"),
                   lrho = 0,
                   ridge = 0,
                   control_algorithm = list(),
                   par_gridsearch = list()) {

  # ---- Check all arguments ----
  optim_method <- match.arg(optim_method)
  optim_criterion <- match.arg(optim_criterion)

  if (is.null(data1ts)) {
    if (is.null(y) | is.null(r)) {
      stop("Please provide either an object created by the function `prepare_data`, or an array of event counts y and an array of exposure r.")
    }
    if (is.null(bins)) stop("Please provide a list with the bins specifications.")
  } else {
    y <- data1ts$bindata[[2]]
    r <- data1ts$bindata[[1]]
    Z <- data1ts$bindata$Z
    bins <- data1ts$bins
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
    if (length(lrho) == 1) {
      lrho <- seq(-3, 3, by = .25)
      message("Grid search method selected, but no vector with grid values over `log_10(rho_s)` is provided.
              I will use a default grid, ", lrho)
    }
    gsp <- list(
      plot_aic = FALSE,
      plot_bic = FALSE,
      return_aic = TRUE,
      return_bic = TRUE,
      mark_optimal = TRUE,
      main_aic = "AIC grid",
      main_bic = "BIC grid"
    )
    Ngsp <- names(gsp)
    namesGsp <- names(par_gridsearch)

    gsp[namesGsp] <- par_gridsearch
    if (length(namesGsp[!namesGsp %in% Ngsp]) > 0) {
      warning("Undefined entries in parameters for the grid search! Default settings are used.\n")
      warning(
        "Undefined keyword(s): ",
        paste(namesGsp[!namesGsp %in% Ngsp], collapse = ", ")
      )
    }
  }

  # ---- Bbases specification ----
  Bbases <- list(
    bdeg = 3,
    nseg_s = 10,
    min_s = min(bins$bins_s),
    max_s = max(bins$bins_s)
  )

  NBbases <- names(Bbases)
  namesBbases <- names(Bbases_spec)

  Bbases[namesBbases] <- Bbases_spec
  if (length(namesBbases[!namesBbases %in% NBbases]) > 0) {
    warning("undefined entries in Bbases_spec! Default settings are used.\n")
    warning(
      "undefined keyword(s): ",
      paste(namesBbases[!namesBbases %in% NBbases], collapse = ", ")
    )
  }

  # ---- Construct B-splines ----
  Bs <- JOPS::bbase(x = bins$mids, nseg = Bbases$nseg_s, xl = Bbases$min_s, xr = Bbases$max_s, bdeg = Bbases$bdeg)
  nbs <- ncol(Bs)

  # ---- Penalty parameters ----
  Ds <- diff(diag(nbs), diff = pord)

  # ---- Optimization ----

  if (optim_method == "grid_search") {
    optimal_model <- grid_search_1d(
      lrho = lrho,
      r = r, y = y,
      Z = Z,
      Bs = Bs,
      Ds = Ds,
      Wprior = Wprior,
      optim_criterion = optim_criterion,
      control_algorithm = con,
      par_gridsearch = gsp
    )
  }

  if (optim_method == "ucminf") {
    optimal_model <- fit1tsmodel_ucminf(
      lrho = lrho,
      y = y, r = r,
      Z = Z,
      Bs = Bs,
      Ds = Ds,
      Wprior = Wprior,
      optim_criterion = optim_criterion,
      control_algorithm = con
    )
  }

  # ---- Save results in list and return list ----
  results <- optimal_model
  class(results) <- "haz1ts"


  return(results)
}
