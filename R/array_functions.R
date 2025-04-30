#' Array functions for the computation of the GLAM model
#'
#' Defines a family of functions for the computation of the GLAM model
#'
#' @param X a matrix
#' @param Y (optional) a matrix
#' @param A an array
#'
#' @noRd
#' @keywords internal
#'

# Row tensor of a matrix X and a matrix Y
Rtens <- function(X, Y=NULL) {
  if(is.null(Y)) Y=X
  one.y <- matrix(1, 1, ncol(Y))
  one.x <- matrix(1, 1, ncol(X))
  kronecker(X, one.y) * kronecker(one.x, Y)
}

# H-transform of an array A by a matrix X
Ht <- function(X, A) {
  d <- dim(A)
  M <- matrix(A, nrow = d[1])
  XM <- X %*% M
  array(XM, c(nrow(XM), d[-1]))
}

# Rotation of an array A
Rotate <- function(A) {
  d <- 1:length(dim(A))
  d1 <- c(d[-1], d[1])
  aperm(A, d1)
}

# Rotated H-transform of an array A by a matrix X
RHt <- function(X, A) Rotate(Ht(X, A))
