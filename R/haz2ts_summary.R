#' Summary function for object of class 'haz2ts'
#'
#' @param x an object of class 'haz2ts' returned by the function [fit2ts()]
#' @param ... further arguments
#'
#' @return a printed summary of the fitted model
#'
#' @examples
#' # Create some fake data - the bare minimum
#' id <- 1:20
#' u <- c(5.43, 3.25, 8.15, 5.53, 7.28, 6.61, 5.91, 4.94, 4.25, 3.86, 4.05, 6.86,
#'        4.94, 4.46, 2.14, 7.56, 5.55, 7.60, 6.46, 4.96)
#' s <- c(0.44, 4.89, 0.92, 1.81, 2.02, 1.55, 3.16, 6.36, 0.66, 2.02, 1.22, 3.96,
#'        7.07, 2.91, 3.38, 2.36, 1.74, 0.06, 5.76, 3.00)
#' ev <- c(1, 0, 0, 1, 0, 1, 0, 1, 1, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 1)#'
#'
#' fakedata <- as.data.frame(cbind(id, u, s, ev))
#' fakedata2ts <- prepare_data(data = fakedata,
#'                             u = "u",
#'                             s_out = "s",
#'                             ev = "ev",
#'                             ds = .5)
#' # Fit a fake model - not optimal smoothing
#' fakemod <- fit2ts(fakedata2ts,
#'                   optim_method = "grid_search",
#'                   lrho = list(seq(1 ,1.5 ,.5),
#'                               seq(1 ,1.5 ,.5)))
#' summary(fakemod)
#'
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
                     expcoef=hr$HR,
                     lowerci=hr$HR - 1.96*hr$SE_HR,
                     upperci=hr$HR + 1.96*hr$SE_HR)
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
  cat("  rho_u = ", res$rhou, "\n")
  cat("  rho_s = ", res$rhos, "\n")
  cat("\n")
  if(is.null(res$coeftab)) cat("Model with no covariates") else {
    colnames(res$coeftab) <- c("beta", "se(beta)", "exp(beta)", "lower .95", "upper.95")
    print(res$coeftab)
  }
  cat("\n\n")
  cat("Model diagnostics: \n")
  cat("  AIC = ", res$AIC, "\n")
  cat("  BIC = ", res$BIC, "\n")
  cat("  ED = ",  res$ED)

  return(invisible())
}
