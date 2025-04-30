#' Predict overall survival and cumulative incidence functions
#'
#' @param models A list with cause-specific hazard models of class `'haz2ts'`,
#'               fitted via `fit2ts()`. The models have to be fitted on the same
#'               grid of points over the time axes.
#' @param newdata A dataframe with columns cointaing the values of the variable `u`
#'                and the time scale `s` for which predictions are to be obtained.
#' @param u The name of the variable in `newdata` containing values for the variable `u`.
#' @param s The name of the variable in `newdata` containing values for the variable `s`.
#' @param ds (optional) The distance between two consecutive points on the `s` axis.
#'           If not provided, an optimal minimum value will be chosen automatically and
#'           a warning is returned.
#'
#' @return A dataframe cointaing the values of `u` and `s` in `newdata`,
#'        the predicted survival probability and the values of the cumulative incidence
#'        functions (cif), for each combination of `u` and `s`. There will be as many values of
#'        the cif as the number of cause-specific hazard models.
#'
#' @details
#' Predictions of cumulated quantities can be provided only within intervals of values on the `s` time scale.
#'
#' @export
#'
#' @examples
#' id <- 1:20
#' u <- c(5.43, 3.25, 8.15, 5.53, 7.28, 6.61, 5.91, 4.94, 4.25, 3.86, 4.05, 6.86,
#'        4.94, 4.46, 2.14, 7.56, 5.55, 7.60, 6.46, 4.96)
#' s <- c(0.44, 4.89, 0.92, 1.81, 2.02, 1.55, 3.16, 6.36, 0.66, 2.02, 1.22, 3.96,
#'        7.07, 2.91, 3.38, 2.36, 1.74, 0.06, 5.76, 3.00)
#' ev <- c(1, 0, 0, 1, 0, 1, 0, 1, 1, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 1)#'
#' fakedata <- as.data.frame(cbind(id, u, s, ev))
#' # cause 1
#' fakedata2ts1 <- prepare_data(u = fakedata$u,
#'                             s_out = fakedata$s,
#'                             ev = fakedata$ev,
#'                             ds = .5)
#' # Fit a fake model - not optimal smoothing for cause type 1
#' fakemod1 <- fit2ts(fakedata2ts1,
#'                   optim_method = "grid_search",
#'                   lrho = list(seq(1 ,1.5 ,.5),
#'                               seq(1 ,1.5 ,.5)))
#' # cause 2
#' fakedata2ts2 <- prepare_data(u = fakedata$u,
#'                               s_out = fakedata$s,
#'                               ev = 1-(fakedata$ev),
#'                               ds = .5)
#' # Fit a fake model - not optimal smoothing for cause 2
#' fakemod2 <- fit2ts(fakedata2ts2,
#'                     optim_method = "grid_search",
#'                     lrho = list(seq(1 ,1.5 ,.5),
#'                                 seq(1 ,1.5 ,.5)))
#'
#' newdata <- as.data.frame(cbind("u" =  c(2.5,3.4,6), "s" = c(.2,.5,1.3)))
#'
#' predict_comprisks2ts(models = list(fakemod1, fakemod2),
#'                      newdata,
#'                      u = "u", s = "s")

predict_comprisks2ts <- function(models = list(),
                                 newdata,
                                 u, s, ds = NULL) {
  L <- length(models) # n. of competing causes
  n <- nrow(newdata)
  predicted <- newdata
  predicted$survival <- NULL
  names_cif <- paste0("cif_", 1:L)
  cifs <- matrix(0, 0, L)

  if(is.null(ds)){
    ds <- round(min((newdata[, which(names(newdata) == s)] -
                       as.integer(newdata[, which(names(newdata) == s)]))),1)
    if(ds == 0) ds <- .1 # make a minimal choice if ds is 0
    message("chosen interval: ds = ", ds)

  }

  for (i in 1:n) {
    pred <- predict_cif2ts_pointwise(
        fitted_models = models,
        u = newdata[i, which(names(newdata) == u)],
        s = newdata[i, which(names(newdata) == s)],
        ds = ds
      )
    predicted[i, "survival"] <- pred$surv
    cifs <- rbind(cifs, subset(pred, select = names_cif))
  }
  predicted <- cbind(predicted, cifs)
  return(predicted)
}
