#' Check inputs for the GLAM model
#'
#' @description `check_inputs()` checks the inputs provided by the user for the
#'   estimation of the GLAM model.
#'
#'   It checks the dimensions of the data matrices and if provided, the matrix
#'   containing the covariate values. If the user provides a matrix of weights,
#'   the dimension of this matrix is also checked. This is an internal
#'   function called by the function `fit2ts`.
#'
#' @param R A matrix of exposure times of dimension nu by ns, if `Z` is not
#'   provided, or a three-dimensional array of dimensions nu by ns by n, if `Z` is provided.
#' @param Y A matrix of events of dimension nu by ns, if `Z` is not provided, or
#'   a three-dimensional array of dimensions nu by ns by n, if `Z` is provided.
#' @param Z An optional regression matrix of covariates.
#'   If provided, it must have n rows.
#' @param Wprior An optional matrix of a-priori weights provided by the user.
#'   This must be of dimension nu by ns.
#'
#' @noRd
#'
check_inputs <- function(R, Y,
                         Z = NULL,
                         Wprior = NULL) {
  # check if R and Y have the same dimensions:
  if (any(dim(R) != dim(Y))) {
    stop("`R` and `Y` must be of the same dimension.")
  }

  # if !is.na(Z) then R and Y should be array of dimensions nu,ns,n and
  # Z should have number of rows = n

  if (!is.null(Z)) {
    if (length(dim(R)) < 3) {
      stop("`Z` is provided but `R` is not a three-dimensional array.")
    } else {
      if (nrow(Z) != dim(R)[3]) {
        stop("The number of rows of Z should be equal to the third dimension of R.")
      }
    }
  }

  # check on weights dimensions

  if (!is.null(Wprior)) {
    if (any(dim(Wprior) != dim(R)[c(1, 2)]))
      stop("The weight matrix Wprior must be of the same dimension as R and Y.")
  }
}
