#' Summary function for object of class 'haz2tsLMM'
#'
#' @param x an object of class 'haz2tsLMM' returned by the function [fit2ts()]
#' @param ... further arguments
#'
#' @importFrom stats coef
#' @return a printed summary of the fitted model, including optimal smoothing
#'          paramters, the effective dimension ED and the AIC/BIC. For model
#'          with covariates, a regression table is also returned.
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
#' fakedata2ts <- prepare_data(data = fakedata,
#'                             u = "u",
#'                             s_out = "s",
#'                             ev = "ev",
#'                             ds = .5)
#' # Fit a fake model - not optimal smoothing
#' fakemod <- fit2ts(fakedata2ts,
#'   optim_method = "LMMsolver"
#' )
#' summary(fakemod)
#' @export
#'
haz2tsLMM_summary <- function(x, ...) {
  if (!inherits(x, "haz2tsLMM")) stop("'x' must be a 'haz2tsLMM' object")

  nevents <- x$nevents
  nu <- x$nu
  ns <- x$ns
  cu <- x$cu
  cs <- x$cs
  AIC <- x$AIC_BIC$AIC
  BIC <- x$AIC_BIC$BIC
  ED <- x$AIC_BIC$ED

  # extract smoothing parameters and transform them to be compatible with those
  # obtained from GLAM
  sc_lambda <- as.numeric(x$optimal_model$theta[c(1, 2)] / x$optimal_model$theta[3])
  dx1 <- attr(x$optimal_model$splRes[[1]]$knots[[1]], which = "dx")
  dx2 <- attr(x$optimal_model$splRes[[1]]$knots[[2]], which = "dx")

  # transform back to standard penalties:
  pord <- x$optimal_model$splRes[[1]]$pord
  rho <- rep(NA, 2)
  rho[1] <- (1 / dx1)^(2 * pord - 1) * sc_lambda[1]
  rho[2] <- (1 / dx2)^(2 * pord - 1) * sc_lambda[2]
  log10rhou <- log10(rho[1])
  log10rhos <- log10(rho[2])

  res <- list(
    nevents = nevents, nu = nu, ns = ns, cu = cu, cs = cs, AIC = AIC, BIC = BIC, ED = ED,
    rho = rho, log10rhou = log10rhou, log10rhos = log10rhos
  )

  if (x$covariates == "yes") {
    namesCov <- x$optimal_model$term.labels.f[-c(1, length(x$optimal_model$term.labels.f))]
    coefLMM <- coef(x$optimal_model, se = T)[namesCov]
    coeftab <- matrix(0, length(namesCov), 3)
    for (ind in 1:length(namesCov)) {
      coeftab[ind, 2] <- round(as.numeric(coefLMM[[ind]]$value), 4)
      coeftab[ind, 3] <- as.numeric(coefLMM[[ind]]$se)
    }
    coeftab <- as.data.frame(coeftab)
    coeftab[, 1] <- namesCov
    colnames(coeftab) <- c("coef", "beta", "SE_beta")
    coeftab$HR <- exp(coeftab[, "beta"])
    coeftab$SE_HR <- coeftab$HR * coeftab$SE_beta
    coeftab$lowerci <- coeftab$HR - 1.96 * coeftab$SE_HR
    coeftab$upperci <- coeftab$HR + 1.96 * coeftab$SE_HR
    res$coeftab <- coeftab
  }

  cat("Number of events = ", res$nevents, "\n")
  cat("Model specifications:\n")
  cat("  nu = ", res$nu, "\n")
  cat("  ns = ", res$ns, "\n")
  cat("  cu = ", res$cu, "\n")
  cat("  cs = ", res$cs, "\n")
  cat("\nOptimal smoothing: \n")
  cat("  log10(rho_u) = ", res$log10rhou, "\n")
  cat("  log10(rho_s) = ", res$log10rhos, "\n")
  cat("  rho_u = ", res$rho[1], "\n")
  cat("  rho_s = ", res$rho[2], "\n")
  cat("\n")
  if(is.null(res$coeftab)) cat("Model with no covariates") else {
    colnames(res$coeftab) <- c("", "beta", "se(beta)", "exp(beta)",
                               "se(exp(beta))", "lower .95", "upper.95")
    print(res$coeftab[, colnames(res$coeftab) %in% c("", "beta", "se(beta)", "exp(beta)",
                                                    "lower .95", "upper.95")] )
  }
  cat("\n\n")
  cat("Model diagnostics: \n")
  cat("  AIC = ", res$AIC, "\n")
  cat("  BIC = ", res$BIC, "\n")
  cat("  ED = ",  res$ED)


  return(invisible())
}
