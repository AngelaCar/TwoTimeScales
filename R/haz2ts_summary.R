#' Summary function for object of class 'haz2ts'
#'
#' @param x an object of class 'haz2ts' returned by the function [fit2ts()]
#' @param ... further arguments
#'
#' @return a printed summary of the fitted model
#' @export
#'
haz2ts_summary <- function(x,...){
  if (!inherits(x, "haz2ts")) stop("'x' must be a 'haz2ts' object")

  nevents <- x$optimal_model$nevents
  nu <- nrow(x$optimal_model$Bbases$Bu)
  ns <- nrow(x$optimal_model$Bbases$Bs)
  cu <- ncol(x$optimal_model$Bbases$Bu)
  cs <- ncol(x$optimal_model$Bbases$Bs)
  AIC <- x$optimal_model$aic
  BIC <- x$optimal_model$bic
  ED <- x$optimal_model$ed
  log10rhou <- x$optimal_logrho[1]
  log10rhos <- x$optimal_logrho[2]
  rhou <- 10^log10rhou
  rhos <- 10^log10rhos

  if(!is.null(x$optimal_model$beta)){ # model with covariates
    hr <- get_hr(x)
  }

  res <- list(nevents=nevents,nu=nu,ns=ns,cu=cu,cs=cs,AIC=AIC,BIC=BIC,ED=ED,
              log10rhou=log10rhou,log10rhos=log10rhos,rhou=rhou,rhos=rhos)
  if(!is.null(x$optimal_model$beta)){
    coeftab <- cbind(coef= hr$beta, secoef=hr$SE_beta,
                     expcoef=hr$HR, seexpcoef=hr$SE_HR)
    res$coeftab <- coeftab
  }

  cat("Model specifications:\n")
  cat("  number of events = ", res$nevents, "\n")
  cat("  nu = ", res$nu, "\n")
  cat("  ns = ", res$ns, "\n")
  cat("  cu = ", res$cu, "\n")
  cat("  cs = ", res$cs, "\n")
  cat("\nOptimal smoothing: \n")
  cat("  log10(rho_u) = ", res$log10rhou, "\n")
  cat("  log10(rho_s) = ", res$log10rhos, "\n")
  cat("  rho_u = ", res$rhou, "\n")
  cat("  rho_s = ", res$rhos, "\n")
  cat("\n")
  if(is.null(res$coeftab)) cat("Model with no covariates") else {
    colnames(res$coeftab) <- c("coef", "se(coef)", "exp(coef)", "se(exp(coef))")
    print(res$coeftab)
  }
  cat("\n\n")
  cat("Model fit: \n")
  cat("  AIC = ", res$AIC, "\n")
  cat("  BIC = ", res$BIC, "\n")

  return(invisible())
}
