# Collection of functions for the competing risks model over two time scales

#' @title Cumulated hazard over two time scales
#'
#' @description Computes the cumulated hazard surface over two time scales
#' (currently only implemented for objects of class `haz2ts`, not `haz2tsLMM`).
#'
#' @inheritParams plot.haz2ts
#' @param cause a character string with a short name for the cause (optional).
#'
#' @return a list with the hazard estimates (obtained from function `get_hazard_2d`),
#'  the cumulated hazard estimate and (if provided) the name for the cause
#'  (for cause-specific hazards).
#' @export
#'
cum_hazard2ts <- function(fitted_model,
                          plot_grid = NULL,
                          cause = NULL){

  Haz <- get_hazard_2d(fitted_model = fitted_model,
                       plot_grid = plot_grid)
  # Cumulative Hazards
  ds <- Haz$new_plot_grid$ds
  CumHaz <- t(apply(Haz$hazard, 1, cumsum) * ds)

  res <- list(
    "Haz" = Haz,
    "CumHaz" = CumHaz
  )
  if(!is.null(cause)){
    res$cause <- cause
  }

  class(res) <- "cumhaz2ts"
  return(res)
}

#' Survival functions (cause-specific) and overall survival with two time scales
#'
#' @description
#' Computes cause-specific survival functions over two time scales from a list of
#' cause-specific cumulated hazard matrices. Each matrix of cause-specific survival
#' is calculated as exp(-CumHaz). Additionally, if the user provides more than
#' one cause-specific cumulated hazard, the function computes the overall survival
#' matrix, that contains the probability of not experiencing an event of any cause
#' by time `s` and fixed time at entry `u`.
#'
#'
#' @param cumhaz a list with all the cause-specific cumulated hazard matrices
#'  (minimum 1). If more than one cause-specific cumulated hazard is provided,
#'  then they should all be matrices of the same dimension.
#' @param cause is an optional vector of short names for the causes. It should
#' be of the same length as the number of cause-specific cumulated hazards provided.
#' @return a list with cause-specific survival probability matrices and
#' one overall survival probability matrix (if at least two cause-specific
#' cumulated hazards are provided)
#' @export
#'
surv2ts <- function(cumhaz = list(),
                    cause = NULL){

  ncauses <- length(cumhaz)
  if(!is.null(cause)){
    if(length(cause) != ncauses)
      message("The number of names provided for the causes is not equal to the number of causes provided.")
  }
  Surv2ts <- vector("list", length = ncauses)

  for(i in 1:ncauses){
    Surv2ts[[i]] <- exp(-cumhaz[[i]])
    names(Surv2ts[[i]]) <- cause[i]
  }

  if(ncauses > 1){
    Surv2ts$OverSurv <- exp(-(Reduce("+", cumhaz)))
  }
  return(Surv2ts)
}

#' Cumulative incidence function over two time scales
#'
#' @param haz a list of cause-specific hazards
#' @param oversurv the overall survival function over two time scales, obtained
#'  from `surv2ts` (the last element in the list of results)
#' @param ds the distance between two consecutive intervals over the `s` time scale.
#'  This has to be equal for all cause-specific hazards
#' @param cause is an optional vector of short names for the causes. It should
#' be of the same length as the number of cause-specific cumulated hazards provided.
#'
#' @return a list with one cumulative incidence matrix for each cause-specific
#'  hazard
#' @export
#'

# ----- Add option for name of cause
cuminc2ts <- function(haz = list(),
                      oversurv,
                      ds,
                      cause = NULL){
  ncauses <- length(haz)
  if(!is.null(cause)){
    if(length(cause) != ncauses)
      message("The number of names provided for the causes is not equal to the number of causes provided.")
  }
  CIF2ts <- vector("list", length = ncauses)
  for(i in 1:ncauses){
    CIF2ts[[i]] <- t(apply(haz[[i]] * oversurv, 1, cumsum) * ds)
    names(CIF2ts[[i]]) <- cause[i]
  }

  return(CIF2ts)
}

