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
    # if covariates, check names first
    # if mathematical symbols in names, str2lang(x) will return an error
    # so they have to be substitute with something else
      namesZ <- colnames(Z)
      if(length(which(grepl("+", namesZ, fixed = T))) != 0){
        namesZ[which(grepl("+", namesZ, fixed = T))] <- sub('+','&',
                                                            x = namesZ[which(grepl("+", namesZ, fixed = T))] ,
                                                            fixed = T)
      }
      if(length(which(grepl("-", namesZ, fixed = T))) != 0){
        namesZ[which(grepl("-", namesZ, fixed = T))] <- sub('-','_',
                                                            x = namesZ[which(grepl("-", namesZ, fixed = T))] ,
                                                            fixed = T)
      }
      if(length(which(grepl("<", namesZ, fixed = T))) != 0){
        namesZ[which(grepl("<", namesZ, fixed = T))] <- sub('<','less',
                                                            x = namesZ[which(grepl("<", namesZ, fixed = T))] ,
                                                            fixed = T)
      }
      if(length(which(grepl(">", namesZ, fixed = T))) != 0){
        namesZ[which(grepl(">", namesZ, fixed = T))] <- sub('>','more',
                                                            x = namesZ[which(grepl(">", namesZ, fixed = T))] ,
                                                            fixed = T)
      }
      colnames(Z) <- namesZ
      datalong = cbind(datalong, Z[datalong$id, ])
      attr(datalong, "cov_names") <- namesZ
  }
    attr(datalong, "bininfo") <- bins

  return(datalong)
}
