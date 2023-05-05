#' Iterative Weighted Least Squares algorithm for 1ts model
#'
#' @description `iwls_1d()` fits the 1ts model with IWLS algorithm.
#'
#' @param r A vector of exposure times of length ns.
#' @param y A vector of event counts of length ns.
#' @param Bs A matrix of B-splines for the `s` time scale of dimension ns by cs.
#' @param P The penalty matrix of dimension cs by cs.
#' @param control_algorithm is a list with optional values for the parameters of iterative processes:
#' @param control_algorithm A list with optional values for the parameters of
#'   the iterative processes:
#'   * `maxiter` The maximum number of iteration for the IWSL algorithm.
#'     Default is 20.
#'   * `conv_crit` The convergence criteria, expressed as difference between
#'     estimates at iteration i and i+1. Default is `1e-5`.
#'   * `verbose` A Boolean. Default is `FALSE`. If `TRUE` monitors the iteration
#'     process.
#'
#' @return A list with the following elements:
#' * `alpha` The vector of estimated P-splines coefficients of length cs.
#' * `SE_alpha` The vector of estimated Standard Errors for the `alpha` coefficients,
#'   of length cs.
#' * `H` The hat-matrix.
#' * `Cov` The full variance-covariance matrix.
#' * `deviance` The deviance.
#' * `ed` The effective dimension of the model.
#' * `aic` The value of the AIC.
#' * `bic` The value of the BIC.
#' * `Bbases` a list with the B-spline basis `Bs` (this is a list for
#'    compatibility with functions in 2d).
#'
#' @export

iwls_1d <- function(r, y,
                    Bs, P,
                    control_algorithm = list(maxiter = 20,
                                             conv_crit = 1e-5,
                                             verbose = FALSE)) {

  # ---- Preliminary operations ----

  # Controls iterative process
  maxiter <- control_algorithm$maxiter
  conv_crit <- control_algorithm$conv_crit
  verbose <- control_algorithm$verbose

  # Initialize the algorithm
  m  <- ncol(Bs)
  a0 <- log(sum(y) / sum(r))
  alpha <- rep(a0, m)

  # ---- IWLS algorithm ----

  for (iter in 1:maxiter) {                     # Start of IWLS algorithm
    eta <- Bs %*% alpha
    mu <- r * exp(eta)
    w <- c(mu)
    z <- eta * mu + (y - mu)
    BWB <- t(Bs) %*% (w * Bs)
    BWBpP <- BWB + P
    Bz <- t(Bs) %*% z
    alpha.old <- alpha
    alpha <- solve(BWBpP, Bz)

    # Check for convergence
    delta <- max(abs(alpha - alpha.old))

    # Monitor
    if (verbose) cat(iter, delta, '\n')
    if (delta > conv_crit & iter == maxiter)
      warning("Max number of iterations ", iter, " reached but the algorithm did not converge.")
    if (delta <= conv_crit) break         # End of IWLS algorithm
  }

  # ---- Optimal model ----

  # Compute optimal quantities
  eta <- Bs %*% alpha
  mu <- r * exp(eta)
  w <- c(mu)

  # Calculate AIC
  y[y == 0] <- 10 ^ -4
  mu_c <- mu
  mu_c[mu_c==0] <- 10 ^ -4
  dev <- 2 * sum(y * log(y / mu_c))
  H <- solve(BWBpP, BWB)
  ed <- sum(diag(H))
  aic <- dev + 2 * ed
  n_obs <- length(y)
  bic <- dev + ed * log(n_obs)

  # Calculate Variance-covariance matrix and SEs
  Cov <- solve(BWBpP)
  se <- sqrt(diag(Cov))

  # ---- Return results ----
  Bbases <- list("Bs" = Bs)

  results <- list(
    alpha = alpha,
    SE_alpha = se,
    eta = eta,
    H = H,
    deviance = dev,
    ed = ed,
    aic = aic,
    bic = bic,
    Bbases = Bbases)
  return(results)
}
