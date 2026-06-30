# Method dispatch for different object types
# Summary methods for `'haz2ts'`, `'haz2tsPGAM'`, `'haz2tsVCM'` and `'haz2tsLMM'` objects

##
#' @export
summary.haz2ts <- function(object,...) {
  haz2ts_summary(x = object,...)
}

##
#' @export
summary.haz2tsPGAM <- function(object,...) {
  haz2tsPGAM_summary(x = object,...)
}

##
#' @export
summary.haz2tsVCM <- function(object,...) {
  haz2tsVCM_summary(x = object,...)
}

##
#' @export
summary.haz2tsLMM <- function(object,...){
  haz2tsLMM_summary(x = object, ...)
}

# Prediction method for `'haz2ts'` objects
#' @export
predict.haz2ts <- function(object,...){
  predict_haz2ts(x = object, ...)
}
