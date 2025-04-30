#' Transforms the time scale u back to t
#'
#' @description This function transforms values of u back to values of t.
#' It is used by the function `plot_haz2ts`.
#'
#' @noRd
#' @keywords internal

utot <- function(u, s){
  t <- u + s
  s <- s
  return(as.data.frame(cbind(t, s)))
}
