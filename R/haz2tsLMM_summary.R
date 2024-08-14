#' Summary function for object of class 'haz2tsLMM'
#'
#' @param x an object of class 'haz2tsLMM' returned by the function [fit2ts()]
#' @param ... further arguments
#'
#' @return a printed summary of the fitted model
#' @export
#'
haz2tsLMM_summary <- function(x,...){
  if (!inherits(x, "haz2tsLMM")) stop("'x' must be a 'haz2tsLMM' object")

  nevents <- x$nevents
  nu <- x$nu
  ns <- x$ns
  cu <- x$cu
  cs <- x$cs
  AIC <- x$AIC_BIC$AIC
  BIC <- x$AIC_BIC$BIC
  ED <- x$AIC_BIC$ED
  #log10rhou <- x$optimal_model$
  #log10rhos <- x$optimal_logrho[2]
  #rhou <- 10^log10rhou
  #rhos <- 10^log10rhos

  # if(!is.null(x$optimal_model$beta)){ # model with covariates
  #   hr <- get_hr(x)
  # }

  res <- list(nevents=nevents,nu=nu,ns=ns,cu=cu,cs=cs,AIC=AIC,BIC=BIC,ED=ED)
  if(x$covariates == "yes"){
    namesCov <- x$optimal_model$term.labels.f[-c(1,length(x$optimal_model$term.labels.f))]
    coefLMM <- coef(x$optimal_model, se = T)[namesCov]
    coeftab <- matrix(0, length(namesCov), 3)
      for(ind in 1:length(namesCov)){
        coeftab[ind,2] <- round(as.numeric(coefLMM[[ind]]$value), 4)
        coeftab[ind,3] <- as.numeric(coefLMM[[ind]]$se)
      }
    coeftab <- as.data.frame(coeftab)
    coeftab[,1] <- namesCov
    colnames(coeftab) <- c("coef", "beta", "SE_beta")
    coeftab$HR <- exp(coeftab[,'beta'])
    coeftab$SE_HR <- coeftab$HR * coeftab$SE_beta

    res$coeftab <- coeftab
  }

  cat("Model specifications:\n")
  cat("  number of events = ", res$nevents, "\n")
  cat("  nu = ", res$nu, "\n")
  cat("  ns = ", res$ns, "\n")
  cat("  cu = ", res$cu, "\n")
  cat("  cs = ", res$cs, "\n")
  # cat("\nOptimal smoothing: \n")
  # cat("  log10(rho_u) = ", res$log10rhou, "\n")
  # cat("  log10(rho_s) = ", res$log10rhos, "\n")
  # cat("  rho_u = ", res$rhou, "\n")
  # cat("  rho_s = ", res$rhos, "\n")
  cat("\n")
  if(is.null(res$coeftab)) cat("Model with no covariates") else {
    colnames(res$coeftab) <- c("coef", "beta", "se(beta)", "exp(beta)", "se(exp(beta))")
    print(res$coeftab)
  }
  cat("\n\n")
  cat("Model fit: \n")
  cat("  AIC = ", res$AIC, "\n")
  cat("  BIC = ", res$BIC, "\n")

  return(invisible())
}
