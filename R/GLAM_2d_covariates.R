#' Fit the 2d GLAM with covariates
#'
#' @description `GLAM_2d_covariates()` fits a GLAM for the hazard with two time
#'   scales, with covariates.
#'
#' @param R A 3d-array of dimensions nu by ns by n containing exposure times.
#' @param Y A 3d-array of dimensions nu by ns by n containing events' indicators.
#' @param Z (optional) A regression matrix of covariates values of dimensions
#'   n by p.
#' @param Bu A matrix of B-splines for the `u` time scale of dimension nu by cu.
#' @param Bs A matrix of B-splines for the `s` time scale of dimension ns by cs.
#' @param P The penalty matrix of dimension cucs by cucs.
#' @param Wprior An optional matrix of a-priori weights.
#' @param ridge A ridge penalty parameter: default is 0.
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
#'   * `Alpha` The matrix of estimated P-splines coefficients of dimension
#'     cu by cs.
#'   * `Cov_alpha` The variance-covariance matrix of the `Alpha` coefficients,
#'     of dimension cucs by cucs.
#'   * `beta` The vector of length p of estimated covariates coefficients.
#'   * `Cov_beta` The variance-covariance matrix of the `beta` coefficients,
#'     of dimension p by p.
#'   * `SE_beta` The vector of length p of estimated Standard Errors for the `beta`
#'      coefficients.
#'   * `Eta0` The matrix of values of the baseline linear predictor (log-hazard)
#'     of dimension nu by ns.
#'   * `H` The hat-matrix.
#'   * `deviance` The deviance.
#'   * `ed` The effective dimension of the model.
#'   * `aic` The value of the AIC.
#'   * `bic` The value of the BIC.
#'   * `Bbases` a list with the B-spline bases `Bu` and `Bs`.
#'
#' @export
#'
GLAM_2d_covariates <- function(R, Y,
                               Bu, Bs, Z,
                               Wprior = NULL,
                               P,
                               ridge = 0,
                               control_algorithm = list(
                                 maxiter = 20,
                                 conv_crit = 1e-5,
                                 verbose = FALSE
                               )) {
  # ---- Preparatory steps ----

  # Controls iterative process
  maxiter <- control_algorithm$maxiter
  conv_crit <- control_algorithm$conv_crit
  verbose <- control_algorithm$verbose

  # Dimensions
  nu <- dim(R)[1]
  ns <- dim(R)[2]
  cu <- ncol(Bu)
  cs <- ncol(Bs)
  n <- dim(R)[3]
  p <- ncol(Z)

  # Weights
  # if Wprior is null, initialize an array of 1 of dimensions nu, ns, n
  if (is.null(Wprior)) {
    Weight <- array(1, dim = c(nu, ns, n))
  } else { # or we transform the matrix Wprior in a tri-dimensional array
    Weight <- array(Wprior, dim = c(nu, ns, n))
  }

  # Penalty matrix with ridge penalty on diagonal
  P_ridge <- P + ridge * spam::as.spam(diag(nrow(P)))

  # Vectors of 1
  one.n <- matrix(1, n, 1)
  one.p <- matrix(1, p, 1)
  one.nu <- matrix(1, nu, 1)
  one.ns <- matrix(1, ns, 1)

  # Row tensor products and kronecker product
  TBu <- Rtens(Bu)
  TBs <- Rtens(Bs)
  B.kron <- kronecker(Bs, Bu)

  # ---- GLAM algorithm ----
  # Initialize the parameters vector
  theta <- rep(NA, length = (cu * cs + p))
  theta[1:(cu * cs)] <- log((sum(Y) / sum(R))) # B-splines coefficients
  theta[(cu * cs + 1):(cu * cs + p)] <- 0  # Covariates coefficients

  for (iter in 1:maxiter) { # Start of IWLS algorithm
    Alpha <- array(theta[1:(cu * cs)], dim = c(cu, cs, 1))
    beta <- array(theta[(cu * cs + 1):(cu * cs + p)], dim = c(1, 1, p))
    Base <- RHt(one.n, RHt(Bs, RHt(Bu, Alpha)))
    Risk <- RHt(Z, RHt(one.ns, RHt(one.nu, beta)))
    Eta <- Base + Risk
    Mu <- R * exp(Eta)
    W <- Mu * Weight
    vx <- apply(W, 3, sum)
    U <- matrix(1, nu, ns)
    U <- apply(W, c(1, 2), sum)

    # Compute inner product and inverse
    # The blocks A, B, Bt and C are the blocks of the inner product
    # Block A
    A <- t(TBu) %*% U %*% TBs
    dim(A) <- c(cu, cu, cs, cs)
    A <- aperm(A, c(1, 3, 2, 4))
    dim(A) <- c(cu * cs, cu * cs)
    ApP <- A + P_ridge

    # BlockB and BlockB transposed
    B <- matrix(0, nu * ns, p)
    for (i in 1:n) B <- B + outer(c(W[, , i]), c(Z[i, ]))
    B <- t(B.kron) %*% B
    Bt <- t(B)

    # Block C
    C <- t(Z) %*% (vx * Z)
    Cinv <- solve(C)

    # Inversion
    CinvBt <- solve(C, Bt)
    BCinv <- B %*% Cinv

    XX <- solve(ApP - (B %*% CinvBt))
    Inv11 <- XX
    Inv12 <- - XX %*% BCinv
    Inv21 <- - CinvBt %*% XX
    Inv22 <- Cinv + CinvBt %*% XX %*% BCinv

    Inv <- cbind(rbind(Inv11, Inv21), rbind(Inv12, Inv22))

    # Right side
    ZZ <- Eta * Mu + (Y - Mu)
    ZZ <- ZZ * Weight
    rs11 <- array(RHt(t(one.n), RHt(t(Bs), RHt(t(Bu), ZZ))),
      dim = cu * cs * 1
    )
    rs12 <- array(RHt(t(Z), RHt(t(one.ns), RHt(t(one.nu), ZZ))),
      dim = 1 * 1 * p
    )
    rs <- c(rs11, rs12)

    # Update theta
    old.theta <- theta
    theta <- Inv %*% rs

    # Check for convergence
    delta <- max(abs(theta - old.theta))
    # Monitor
    if (verbose) cat(iter, delta, '\n')
    if (delta > conv_crit & iter == maxiter)
      warning("Max number of iterations ", iter, " reached but the algorithm did not converge.")
    if (delta <= conv_crit) break
  } # End of algorithm

  # ---- Compute optimal quantities after convergence ----
  Alpha <- array(theta[1:(cu * cs)], dim = c(cu, cs, 1))
  beta <- array(theta[(cu * cs + 1):(cu * cs + p)],
                dim = c(1, 1, p))
  Base <- RHt(one.n, RHt(Bs, RHt(Bu, Alpha)))
  Risk <- RHt(Z, RHt(one.ns, RHt(one.nu, beta)))
  Eta <- Base + Risk
  Mu <- R * exp(Eta)
  W <- Mu * Weight

  # Compute Effective Dimension
  # Block A
  vx <- apply(W, 3, sum)
  U <- matrix(1, nu, ns)
  U <- apply(W, c(1, 2), sum)
  A <- t(TBu) %*% U %*% TBs
  dim(A) <- c(cu, cu, cs, cs)
  A <- aperm(A, c(1, 3, 2, 4))
  dim(A) <- c(cu * cs, cu * cs)

  # BlockB and BlockB transposed
  B <- matrix(0, nu * ns, p)
  for (i in 1:n) B <- B + outer(c(W[, , i]), c(Z[i, ]))
  B <- t(B.kron) %*% B
  Bt <- t(B)

  # Block C
  C <- t(Z) %*% (vx * Z)

  # Compute hat-matrix H
  Right <- cbind(rbind(A, Bt), rbind(B, C))
  Left <- cbind(rbind(A + P_ridge, Bt), rbind(B, C))
  H <- solve(Left, Right)
  ed <- sum(diag(H))

  # Calculate AIC
  Y_c <- Y
  Y_c[Y_c == 0] <- 1e-7
  Mu_c <- Mu
  Mu_c[Mu_c == 0] <- 1e-7
  dev <- 2 * sum(Y_c * log(Y_c / Mu_c))
  aic <- dev + 2 * ed
  #n_obs <- prod(dim(Y))
  n_obs <- sum(R > 0)
  bic <- dev + ed * log(n_obs)

  # Variance-covariance matrix and SEs
  Cov <- solve(Left)
  SE <- sqrt(diag(Cov))
  Cov_Alpha <- Cov[1:(cu * cs), 1:(cu * cs)]
  Cov_beta <- Cov[(cu * cs + 1):(cu * cs + p), (cu * cs + 1):(cu * cs + p)]

  # ---- Save results in list ----
  Alpha <- matrix(Alpha, nrow = cu, ncol = cs)
  beta <- array(theta[(cu * cs + 1):(cu * cs + p)],
                dim = c(p))
  names(beta) <- attributes(Z)$dimnames[[2]]
  SE_Alpha <- array(SE[1:(cu * cs)], dim = c(cu, cs))
  SE_beta <- array(SE[(cu * cs + 1):(cu * cs + p)], dim = c(p))

  Bbases <- list("Bu" = Bu, "Bs" = Bs)
  # results
  results <- list(
    Alpha = Alpha,
    Cov_Alpha = Cov_Alpha,
    beta = beta,
    Cov_beta = Cov_beta,
    SE_beta = SE_beta,
    Eta0 = Base,
    H = H,
    deviance = dev,
    ed = ed,
    aic = aic,
    bic = bic,
    Bbases = Bbases,
    nevents = sum(Y)
  )
  return(results)
}
