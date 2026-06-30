#' Fit a Varying Coefficient Model (VCM) over two time scales
#'
#' @description
#' `fitvcm()` fits a varying coefficient model over two time scales with P-splines.
#'            The model is described in Carollo et al. (2025).
#'
#'
#' @param data2ts (optional) an object of class `"data2ts"` created by the function
#'   `prepare_data()`. Proving this input is the easiest way to use the function
#'   `fitvcm()`. However, the user can also provide the input data together with
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
#'   * `nseg_s` The number of segments for the B-splines over `s`. Default is 10.
#'   * `min_s` (optional) The lower limit of the domain of `Bs`.
#'     Default is `min(bins_s)`.
#'   * `max_s` (optional) The upper limit of the domain of `Bs`.
#'     Default is `max(bins_s)`.
#' @param pord The order of the penalty. Default is 2.
#' @param kappa A ridge penalty.
#' @param lsmpar The starting values for the two smoothing parameters.
#'               A vector of two elements, default is `c(0,0)`.
#' @param control_algorithm A list with optional values for the parameters of
#'   the iterative processes:
#'   * `maxiter` The maximum number of iteration for the IWSL algorithm.
#'     Default is 20.
#'   * `conv_crit` The convergence criteria, expressed as difference between
#'     estimates at iteration i and i+1. Default is `1e-5`.
#'   * `verbose` A Boolean. Default is `FALSE`. If `TRUE` monitors the iteration
#'     process.
#'   * `monitor_ev` A Boolean. Default is `FALSE`. If `TRUE` monitors the
#'     evaluation of the model over the `log_10(rho_s)` values.
#'
#' @return An object of class `haz2tsVCM`, that is a list with the following elements:
#'   * `optimal_model` A list with :
#'     * `Alpha` The matrix of estimated P-splines coefficients of dimension
#'       \eqn{c_s} by 2.
#'     * `Cov_alpha` The variance-covariance matrix of the `Alpha` coefficients,
#'       of dimension \eqn{c_uc_s} by \eqn{c_uc_s}.
#'     * `Eta` The matrix of values of the baseline linear predictor
#'       (log-hazard) of dimension \eqn{n_u} by \eqn{n_s}.
#'     * `H` The hat-matrix.
#'     * `deviance` The deviance.
#'     * `ed` The effective dimension of the model.
#'     * `aic` The value of the AIC.
#'     * `bic` The value of the BIC.
#'     * `Bbases` a list with the B-spline bases `Bu` and `Bs`
#'   * `optimal_logsmpar` A vector with the optimal values of \eqn{\log_{10}(\theta)} and \eqn{\log_{10}(\phi)}.
#'   * `P_optimal` The optimal penalty matrix P.
#'   * `AIC` (if `par_gridsearch$return_aic == TRUE`) The matrix of AIC values.
#'   * `BIC` (if `par_gridsearch$return_bic == TRUE`) The matrix of BIC values.
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
#' fitvcm(fakedata2ts)
#'
fitvcm <- function(data2ts = NULL,
                   Y = NULL, R = NULL, Z = NULL,
                   bins = NULL,
                   Bbases_spec = list(),
                   pord = 2,
                   kappa = 1e-10,
                   lsmpar = c(0, 0),
                   control_algorithm = list()) {
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
    maxiter_vcm = 15,
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

  # ---- Checks for dimensions of input data ----
  check_inputs(
    R = R, Y = Y,
    Z = Z
  )

  # ---- Reformat data for vcm model ----
  y <- as.vector(Y)
  expose <- as.vector(R)
  ufitting <- rep(bins$midu, ncol(Y))
  sfitting <- rep(bins$mids, each = nrow(Y))

  m <- length(y)

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
    warning("Undefined entries in Bbases_spec. Default settings are used.\n")
    warning(
      "Undefined keyword(s): ",
      paste(namesBbases[!namesBbases %in% NBbases], collapse = ", ")
    )
  }

  # ---- Construct B-splines ----
  Bs <- JOPS::bbase(
    x = sfitting, nseg = Bbases$nseg_s,
    xl = Bbases$min_s, xr = Bbases$max_s,
    bdeg = Bbases$bdeg
  )
  nbs <- ncol(Bs)
  Ds <- diff(diag(nbs), diff = pord)
  P <- t(Ds) %*% Ds

  # Build VCM basis and block diag penalty matrix
  C <- cbind(Bs, diag(ufitting) %*% Bs)
  K <- diag(2 * nbs) * kappa


  # Initialize
  eta <- log((y + 0.1) / (expose + 0.2))
  a <- 0
  smpar <- 10^lsmpar

  # Iterate for lambdas (Schall algorithm)
  for (it2 in 1:con$maxiter) {
    Q <- kronecker(diag(smpar), P) + K

    # Fit VCM
    a1 <- a
    for (it in 1:con$maxiter_vcm) {
      mu <- expose * exp(eta)
      z <- y - mu + mu * eta
      W <- diag(c(mu))
      S <- t(C) %*% W %*% C
      anew <- solve(S + Q, t(C) %*% z)
      da <- max(abs(anew - a))
      if (da < con$conv_crit) break
      a <- anew
      eta <- C %*% a
    }

    # Compute effective dimensions and sums of squares of coefficients, aic, bic
    y[y == 0] <- 10^-4
    mu_c <- mu
    mu_c[mu_c == 0] <- 10^-4
    dev <- 2 * sum(y * log(y / mu_c))
    G <- solve(S + Q, S)
    ssa <- eds <- rep(0, 2)
    g <- diag(G)
    ed <- sum(g)
    VarCov <- solve(S + Q)
    for (k in 1:2) {
      r <- (k - 1) * nbs + (1:nbs)
      ssa[k] <- sum(a[r]^2)
      eds[k] <- sum(g[r])
    }
    aic <- dev + 2 * ed
    n_obs <- sum(expose > 0)
    bic <- dev + ed * log(n_obs)

    sp <- smpar
    smpar <- eds / ssa
    dsp <- max(abs(log10(sp) - log10(smpar)))
    if (con$monitor_ev) cat(it2, dsp, log10(smpar), "\n")
  }

  # ---- Compute final estimates ----
  Alpha <- matrix(a, nbs, 2) # coefficients for log-baseline [, 1] and
  # for varying coefficient [, 2]
  SE_Alpha <- matrix(sqrt(diag(VarCov)), nbs, 2) # ses of the coefficients for
  # log-baseline [, 1] and for varying coefficient [, 2]

  s <- bins$mids
  new_Bs <- JOPS::bbase(
    x = s, nseg = Bbases$nseg_s,
    xl = Bbases$min_s, xr = Bbases$max_s,
    bdeg = Bbases$bdeg
  )
  u <- bins$midu
  eta <- new_Bs %*% Alpha[, 1]
  vc <- new_Bs %*% Alpha[, 2]

  Eta <- matrix(rep(eta, length(u)),
    byrow = T,
    nrow = length(u),
    ncol = length(s)
  ) +
    outer(u, as.vector(vc))

  # ---- Return results ----
  optimal_model <- list(
    Alpha = Alpha,
    SE_Alpha = SE_Alpha,
    VarCov = VarCov,
    deviance = dev,
    Eta = Eta,
    ed = ed,
    aic = aic,
    bic = bic,
    Bbases = Bbases
  )

  optimal_model$Bbases$Bs <- Bs
  results <- list(
    "optimal_model" = optimal_model,
    "optimal_logsmpar" = log10(smpar),
    "P_optim" = P,
    "nevents" = sum(Y),
    "ns" = dim(Y)[1],
    "cs" = nbs,
    "ufitting" = u
  )
  class(results) <- "haz2tsVCM"
  return(results)
}
