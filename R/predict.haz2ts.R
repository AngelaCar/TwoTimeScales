#' Prediction method for objects of class `'haz2ts'`
#'
#' @param x an object of class `'haz2ts'`, the output of the function `fit2ts()`.
#' @param newdata (optional) A dataframe with columns cointaing the values of the variable `u`
#'                and the time scale `s` for which predictions are to be obtained.
#' @param originaldata (optional) The original dataset. Provide it to obtain individual predictions
#'                for each observation in the data.
#' @param u The name of the variable in `newdata`, or in `originaldata` containing
#'          values for the variable `u`.
#' @param s The name of the variable in `newdata`, or in `originaldata` containing
#'          values for the variable `s`. Note that over the `s` axis predictions
#'          are provided only within intervals of values, as it is necessary to
#'          approximate cumulated quantities.
#' @param id (optional) The name of the variable in `newdata`, or in `originaldata` containing
#'          the identification of each observation. It is not required for predictions on a
#'          new dataset.
#' @param ds (optional) The distance between two consecutive points on the `s` axis.
#'           If not provided, an optimal minimum value will be chosen automatically and
#'           a warning is returned.
#'
#' @return A dataframe. This can be the original dataframe (`originaldata`), where only the
#'         variables `id`, `u` and `s` are selected, or the new data frame (`newdata`),
#'         together with the predicted values for the hazard `hazard` and its
#'         standard errors `se_hazard`, the cumulative hazard `cumhazard` and the
#'         survival probability `survival`.
#'
#' @details
#' Predictions of cumulated quantities can be provided only within intervals of values on the `s` time scale.
#' @export
#'
#' @examples
#' # Create the same fake data as in other examples
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
#' fakedata <- as.data.frame(cbind(id, u, s, ev))
#' fakedata2ts <- prepare_data(
#'   u = fakedata$u,
#'   s_out = fakedata$s,
#'   ev = fakedata$ev,
#'   ds = .5
#' )
#' # Fit a fake model - not optimal smoothing
#' fakemod <- fit2ts(fakedata2ts,
#'   optim_method = "grid_search",
#'   lrho = list(
#'     seq(1, 1.5, .5),
#'     seq(1, 1.5, .5)
#'   )
#' )
#' # Create a new dataset for prediction
#' newdata <- as.data.frame(cbind("u" = c(2.5, 3.4, 6), "s" = c(.2, .5, 1.3)))
#'
#' # First - predict on original data
#' predict(object = fakemod,
#'   originaldata = fakedata, u = "u", s = "s", id = "id"
#' )
#'
#' # Now - predict on new dataset
#' predict(object = fakemod,
#'   newdata = newdata, u = "u", s = "s"
#' )
#'
predict.haz2ts <- function(x,
                           newdata = NULL,
                           originaldata = NULL,
                           u, s, id = NULL, ds = NULL) {
  if (!inherits(x, "haz2ts")) stop("'x' must be a 'haz2ts' object")

  if (!is.null(originaldata)) {
    n <- nrow(originaldata)
    var <- c(id, u, s)
    predicted <- subset(originaldata, select = var)
    predicted$hazard <- NULL
    predicted$cumhazard <- NULL
    predicted$se_hazard <- NULL
    predicted$survival <- NULL

    if(is.null(ds)){
      ds <- round(min((originaldata[, which(names(originaldata) == s)] -
                         as.integer(originaldata[, which(names(originaldata) == s)]))),1)
      if(ds == 0) ds <- .1 # make a minimal choice if ds is 0
    }

    for (i in 1:n) {
      pred <- predict_haz2ts_pointwise(
        fitted_model = x,
        u = originaldata[i, which(names(originaldata) == u)],
        s = originaldata[i, which(names(originaldata) == s)],
        ds = ds
      )
      predicted[i, "hazard"] <- pred$hazard
      predicted[i, "se_hazard"] <- pred$se_hazard
      predicted[i, "cumhazard"] <- pred$cumhaz
      predicted[i, "survival"] <- pred$survival
    }
  } else {
    if (!is.null(newdata)) {
      n <- nrow(newdata)
      predicted <- newdata
      predicted$hazard <- NULL
      predicted$cumhazard <- NULL
      predicted$se_hazard <- NULL
      predicted$survival <- NULL

      if(is.null(ds)){
        ds <- round(min((newdata[, which(names(newdata) == s)] -
                           as.integer(newdata[, which(names(newdata) == s)]))),1)
        if(ds == 0) ds <- .1 # make a minimal choice if ds is 0
        message("chosen interval: ds = ", ds)
      }

      for (i in 1:n) {
        pred <- predict_haz2ts_pointwise(
          fitted_model = x,
          u = newdata[i, which(names(newdata) == u)],
          s = newdata[i, which(names(newdata) == s)],
          ds = ds
        )
        predicted[i, "hazard"] <- pred$hazard
        predicted[i, "se_hazard"] <- pred$se_hazard
        predicted[i, "cumhazard"] <- pred$cumhaz
        predicted[i, "survival"] <- pred$survival
      }
    }
  }
  return(predicted)
}
