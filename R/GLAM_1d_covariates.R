#' Fit the 1d GLAM with covariates
#'
#' @description `GLAM_1d_covariates()` fits a GLAM for the hazard with one time
#'   scale, with covariates.
#'
#' @param R A 2d-array of dimensions ns by n containing exposure times.
#' @param Y A 2d-array of dimensions ns by n containing events' indicators.
#' @param Bs A matrix of B-splines for the `s` time scale of dimension ns by cs.
#' @param Z A regression matrix of covariates values of dimensions n by p.
#' @param Wprior An optional vector of length ns of a-priori weights.
#' @param P The penalty matrix of dimension cs by cs.
#' @param control_algorithm A list with optional values for the parameters of
#'   iterative processes:
#'     *`maxiter` The maximum number of iterations for the IWSL algorithm,
#'     default is 20 .
#'     * `conv_crit` The convergence criteria, expressed as difference between
#'     estimates at iteration i and i+1, default is `1e-5`.
#'     * `verbose` A Boolean. Default is `FALSE`. If `TRUE`, monitor the
#'     iteration process.
#'
#' @return A list with the following elements:
#' * `alpha` The vector of estimated P-splines coefficients of length cs.
#' * `SE_alpha` The vector of estimated Standard Errors for the `alpha` coefficients,
#'   of length cs.
#' * `beta` The vector of length p of estimated covariates coefficients.
#' * `se_beta` The vector of length p of estimated Standard Errors for the `beta`
#'    coefficients.
#' * `eta0` The vector of values of the baseline linear predictor (log-hazard).
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
#'
GLAM_1d_covariates <- function(R, Y,
                               Bs, Z = Z,
                               Wprior = NULL,
                               P,
                               control_algorithm = list(maxiter = 20,
                                                        conv_crit = 1e-5,
                                                        verbose = FALSE)) {
  # ---- Preparatory steps ----

  # Controls iterative process
  maxiter <- control_algorithm$maxiter
  conv_crit <- control_algorithm$conv_crit
  verbose <- control_algorithm$verbose

  # Dimensions
  ns <- dim(R)[1]
  cs <- ncol(Bs)
  n <- dim(R)[2]
  p <- ncol(Z)

  # Vectors of 1
  one.n <- matrix(1, n, 1)
  one.p <- matrix(1, p, 1)
  one.ns <- matrix(1, ns, 1)

  # Weights
  # if Wprior is null, initialize an array of 1 of dimensions ns, n
  if (is.null(Wprior)) {
    Weight <- array(1, dim = c(ns, n))
  } else { # or we transform the vector Wprior in a two-dimensional array
    Weight <- array(Wprior, dim = c(ns, n))
  }

  # Row tensor product
  TBs <- Rtens(Bs)
  B1 <- Rtens(one.ns, Bs)
  B2 <- Rtens(Z, one.n)

  # ---- GLAM algorithm ----
  # Initialize the parameters vector
  theta <- rep(NA, length = (cs + p))
  theta[1:cs] <- log((sum(Y) / sum(R))) # B-splines coefficients
  theta[(cs + 1):(cs + p)] <- 0  # Covariates coefficients

  for (iter in 1:maxiter) { # Start of IWLS algorithm
    alpha <- array(theta[1:cs], dim = c(cs, 1))
    beta <- array(theta[(cs + 1):(cs + p)], dim = c(1, p))
    Base <- RHt(one.n, RHt(Bs, alpha))
    Risk <- RHt(Z, RHt(one.ns, beta))
    Eta <- Base + Risk
    Mu <- R * exp(Eta)
    W <- Mu * Weight
    vx <- apply(W, 2, sum)

    # Compute inner product and inverse
    # The blocks A, B, Bt and C are the blocks of the inner product
    # Block A
    A <- t(TBs) %*% W %*% one.n
    dim(A) <- c(cs, cs)
    ApP <- A + P

    # BlockB and BlockB transposed
    B <- t(B1) %*% W %*% B2
    dim(B) <- c(cs, p)
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
    rs11 <- array(RHt(t(one.n), RHt(t(Bs), ZZ)),
                  dim = cs * 1
    )
    rs12 <- array(RHt(t(Z), RHt(t(one.ns), ZZ)),
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
  alpha <- array(theta[1:cs], dim = c(cs, 1))
  beta <- array(theta[(cs + 1):(cs + p)],
                dim = c(1, p))
  Base <- RHt(one.n, RHt(Bs, alpha))
  Risk <- RHt(Z, RHt(one.ns, beta))
  Eta <- Base + Risk
  Mu <- R * exp(Eta)
  W <- Mu * Weight
  eta0 <- RHt(Bs, alpha)

  # Compute Effective Dimension
  # Block A
  A <- t(TBs) %*% W %*% one.n
  dim(A) <- c(cs, cs)
  ApP <- A + P

  # BlockB and BlockB transposed
  B <- t(B1) %*% W %*% B2
  dim(B) <- c(cs, p)
  Bt <- t(B)

  # Block C
  C <- t(Z) %*% (vx * Z)

  # Compute hat-matrix H
  Right <- cbind(rbind(A, Bt), rbind(B, C))
  Left <- cbind(rbind(ApP, Bt), rbind(B, C))
  H <- solve(Left, Right)
  ed <- sum(diag(H))

  # Calculate AIC
  Y_c <- Y
  Y_c[Y_c == 0] <- 1e-7
  Mu_c <- Mu
  Mu_c[Mu_c == 0] <- 1e-7
  dev <- 2 * sum(Y_c * log(Y_c / Mu_c))
  aic <- dev + 2 * ed
  #n_obs <- prod(dim(Y)) #
  n_obs <- sum(R > 0)
  bic <- dev + ed * log(n_obs)

  # Variance-covariance matrix and SEs
  Cov <- solve(Left)
  SE <- sqrt(diag(Cov))
  Cov_alpha <- Cov[1:cs, 1:cs]

  # ---- Save results in list ----
  alpha <- array(alpha, dim = cs)
  beta <- array(beta, dim = p)
  names(beta) <- attributes(Z)$dimnames[[2]]
  SE_alpha <- array(SE[1:cs], dim = cs)
  SE_beta <- array(SE[(cs + 1):(cs + p)], dim = p)

  Bbases <- list("Bs" = Bs)

  # results
  results <- list(
    alpha = alpha,
    SE_alpha = SE_alpha,
    beta = beta,
    SE_beta = SE_beta,
    Cov = Cov,
    deviance = dev,
    Eta = Eta,
    H = H,
    aic = aic,
    bic = bic,
    ed = ed,
    Bbases = Bbases
  )
  return(results)
}
