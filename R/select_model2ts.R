#' Selects best model among a group of fitted models
#'
#' @description
#' `select_model2ts()` takes as input a list of hazard models with two time scales,
#'                     and returns a table with the best fitting model indicated.
#'                     The selection is based on minimization of the AIC or BIC,
#'                     so the models do not need to be nested. It can be used, for
#'                     example, to identify the best model among the log-additive
#'                     model, the varying coefficient model, and the full 2D model.
#'
#' @param model_list is a list of fitted objects, of class `haz2ts`, `hat2tsLMM`,
#'                   `haz2tsPGAM`, or `haz2tsVCM`.
#' @param sel_criteria to determine the best model - it can be `'aic'` or `'bic'`.
#'
#' @return A data.frame with as many rows as the number of models being compared,
#'         and the following columns:
#'         * Model: the model name (see details and examples)
#'         * Type: the class of the model
#'         * AIC
#'         * BIC
#'         * Best: indicates which of the model is the best fitting one with
#'                 respect to the criteria indicated
#' @details In the model list, it is possible to provide a character name for
#'          each model. In such case, these names will also be reported in the
#'          results' table. An example is provided (See examples).
#'
#' @export
#'
#' @examples
#' id <- 1:20
#' u <- c(
#'   5.43, 3.25, 8.15, 5.53, 7.28, 6.61, 5.91, 4.94, 4.25, 3.86, 4.05, 6.86,
#'   4.94, 4.46, 2.14, 7.56, 5.55, 7.60, 6.46, 4.96
#' )
#' s <- c(
#'   0.44, 4.89, 0.92, 1.81, 2.02, 1.55, 3.16, 6.36, 0.66, 2.02, 1.22, 3.96,
#'   7.07, 2.91, 3.38, 2.36, 1.74, 0.06, 5.76, 3.00
#' )
#' ev <- c(1, 0, 0, 1, 0, 1, 0, 1, 1, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 1)
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
#' full2d <- fit2ts(fakedata2ts,
#'   optim_method = "grid_search",
#'   lrho = list(seq(1, 1.5, .5), seq(1, 1.5, .5))
#' )
#' pgam <- fitpgam(fakedata2ts,
#'   optim_method = "grid_search",
#'   optim_criterion = "aic",
#'   lrho = list(seq(1, 1.5, .5), seq(1, 1.5, .5))
#' )
#' vcm <- fitvcm(fakedata2ts)
#'
#' select_model2ts(model_list = list(
#'   "full interaction" = full2d,
#'   "additive" = pgam,
#'   "varying coeff" = vcm
#' ))
#'

#'
select_model2ts <- function(model_list,
                            sel_criteria = c("aic", "bic")) {
  # ---- Check arguments ----
  sel_criteria <- match.arg(sel_criteria)

  if (!is.list(model_list) | length(model_list) == 1) {
    stop("`model_list` should be a list of at list two fitted models")
  }

  # ---- Extract model components ----
  # length of the list = number of models to compare
  nmod <- length(model_list)

  models <- data.frame(matrix(NA, nrow = nmod, ncol = 5))
  names(models) <- c("Model", "Type", "AIC", "BIC", "Best")

  fitmod_names <- names(model_list)
  if (is.null(fitmod_names)) {
    fitmod_names <- paste0("model_", 1:nmod)
  }
  models$Model <- fitmod_names

  for (i in 1:nmod) {
    models$Type[i] <- class(model_list[[i]])
    models$AIC[i] <- model_list[[i]]$optimal_model$aic
    models$BIC[i] <- model_list[[i]]$optimal_model$bic
  }

  if (sel_criteria == "aic") {
    models$Best[which(models$AIC == min(models$AIC))] <- "*"
  } else {
    models$Best[which(models$BIC == min(models$BIC))] <- "*"
  }

  return(models)
}
