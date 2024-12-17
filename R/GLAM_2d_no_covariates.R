#' Fit the 2d GLAM without covariates
#'
#' @description `GLAM_2d_no_covariates()` fits a GLAM for the hazard with two time
#'   scales, without covariates.
#'
#' @param R A matrix of exposure times of dimension nu by ns.
#' @param Y A matrix of event counts of dimension nu by ns.
#' @inheritParams GLAM_2d_covariates
#'
#' @return A list with the following elements:
#'   * `Alpha` The matrix of estimated P-splines coefficients of dimension
#'     cu by cs.
#'   * `Cov_alpha` The variance-covariance matrix of the `Alpha` coefficients,
#'     of dimension cucs by cucs.
#'   * `Eta0` The matrix of values of the baseline linear predictor (log-hazard)
#'     of dimension nu by ns.
#'   * `H` The hat-matrix.
#'   * `deviance` The deviance.
#'   * `ed` The effective dimension of the model.
#'   * `aic` The value of the AIC.
#'   * `bic` The value of the BIC.
#'   * `Bbases` a list with the B-spline bases `Bu` and `Bs`.
#'

GLAM_2d_no_covariates <- function(R, Y,
                                  Bu, Bs,
                                  Wprior = NULL,
                                  P,
                                  ridge = 0,
                                  control_algorithm = list(
                                    maxiter = 20,
                                    conv_crit = 1e-5,
                                    verbose = FALSE
                                  )) {
  # ---- Preparatory steps ----
  # Dimensions (they have been already checked for consistency)
  nu <- nrow(R)
  ns <- ncol(R)
  cu <- ncol(Bu)
  cs <- ncol(Bs)

  # Controls iterative process
  maxiter <- control_algorithm$maxiter
  conv_crit <- control_algorithm$conv_crit
  verbose <- control_algorithm$verbose

  # Penalty matrix with ridge on diagonal
  P_ridge <- P + ridge * diag(nrow(P))

  # If Wprior is null, initialize a matrix of 1 of the same dimensions as R
  if (is.null(Wprior)) {
    Wprior <- matrix(1, nu, ns)
  }

  # Calculate the row tensor products
  TBu <- Rtens(Bu)
  TBs <- Rtens(Bs)

  # ---- IWLS algorithm ----
  # Initialize the algorithm
  alpha <- rep(log((sum(Y) / sum(R))), cu * cs)

  for (iter in 1:maxiter) { # Start of IWLS algorithm
    Alpha <- matrix(alpha, cu, cs)
    Eta <- RHt(Bs, RHt(Bu, Alpha))
    Mu <- R * exp(Eta)
    W <- Mu * Wprior

    # Left side equation
    BtWB <- t(TBu) %*% W %*% TBs
    dim(BtWB) <- c(cu, cu, cs, cs)
    BtWB <- aperm(BtWB, c(1, 3, 2, 4))
    dim(BtWB) <- c(cu * cs, cu * cs)
    BtWBpP <- BtWB + P_ridge
    Inv <- solve(BtWBpP)

    # Right side equation
    Z <- Eta * Mu + (Y - Mu)
    RS <- t(Bu) %*% (Wprior * Z) %*% Bs
    rs <- RS
    dim(rs) <- c(cu * cs, 1)

    # Update alpha
    old.alpha <- alpha
    alpha <- Inv %*% rs

    # Check for convergence
    delta <- max(abs(alpha - old.alpha))
    # Monitor
    if (verbose) cat(iter, delta, "\n")
    if (delta > conv_crit & iter == maxiter) {
      warning("Max number of iterations ", iter, " reached but the algorithm did not converge.")
    }
    if (delta <= conv_crit) break
  } # End of IWLS algorithm

  # ---- Compute optimal quantities after convergence ----
  Alpha <- matrix(alpha, cu, cs)
  Eta <- RHt(Bs, RHt(Bu, Alpha))
  Mu <- R * exp(Eta)

  # Compute Effective Dimension
  BtWB <- t(TBu) %*% (Mu * Wprior) %*% TBs
  dim(BtWB) <- c(cu, cu, cs, cs)
  BtWB <- aperm(BtWB, c(1, 3, 2, 4))
  dim(BtWB) <- c(cu * cs, cu * cs)
  L <- BtWB + P_ridge
  G <- solve(L, BtWB)
  ed <- sum(diag(G))

  # Compute hat-matrix H
  GG <- solve(L)
  Cov <- GG
  dim(GG) <- c(cu, cs, cu, cs)
  GG <- aperm(GG, c(1, 3, 2, 4))
  dim(GG) <- c(cu * cu, cs * cs)
  H <- Mu * (TBu %*% GG %*% t(TBs))

  # Compute AIC
  Y_c <- Y
  Y_c[Y_c == 0] <- 1e-7
  Mu_c <- Mu
  Mu_c[Mu_c == 0] <- 1e-7
  dev <- 2 * sum(Y_c * log(Y_c / Mu_c))
  aic <- dev + 2 * ed
  n_obs <- sum(R > 0)
  bic <- dev + ed * log(n_obs)

  # Variance-covariance matrix and SEs
  SE <- sqrt(diag(Cov))
  SE_Alpha <- matrix(SE, cu, cs)

  Bbases <- list("Bu" = Bu, "Bs" = Bs)

  # ---- Save results in list ----
  # results
  results <- list(
    Alpha = Alpha,
    Cov_Alpha = Cov,
    deviance = dev,
    Eta = Eta,
    H = H,
    ed = ed,
    aic = aic,
    bic = bic,
    Bbases = Bbases,
    nevents = sum(Y)
  )
  return(results)
}
