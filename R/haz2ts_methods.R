# Summary methods for `'haz2ts'` and `'haz2tsLMM'` objects
#' @export
summary.haz2ts <- function(object,...) {
  haz2ts_summary(x = object,...)
}

##
#' @export
summary.haz2tsLMM <- function(object,...){
  haz2tsLMM_summary(x = object, ...)
}
