#' Process data to fit model with LMMsolver
#'
#' @inheritParams fit2ts
#'
#' @return A dataset in long form to fit the model with LMMsolver
#'
#' @importFrom reshape2 melt
#' @export
#'

prepare_data_LMMsolver <- function(Y = Y, R = R, Z = NULL, bins = bins){

  if(is.null(Z)){
    datalong = reshape2::melt(R)
    colnames(datalong) = c('u_ind', 's_ind', 'r')
    datalong$y = reshape2::melt(Y)$value
    datalong$u <- bins$midu[datalong$u_ind]
    datalong$s <- bins$mids[datalong$s_ind]

    # Remove bins that are not at risk
    datalong = datalong[datalong$r > 0, ]
    attr(datalong, "bininfo") <- bins
  } else {
    datalong = melt(R)
    colnames(datalong) = c('u_ind', 's_ind', 'id', 'r')
    datalong$y = reshape2::melt(Y)$value
    datalong$u <- bins$midu[datalong$u_ind]
    datalong$s <- bins$mids[datalong$s_ind]

    # Remove bins that are not at risk
    datalong = datalong[datalong$r > 0, ]
    datalong = cbind(datalong, Z[datalong$id, ])
    attr(datalong, "bininfo") <- bins
  }
  return(datalong)
}
