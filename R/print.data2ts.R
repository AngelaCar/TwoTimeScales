#' Print method for a `data2ts` object
#'
#' @description Print method for an object of class `data2ts`
#'
#' @param x of class `data2ts`, as prepared by \code{\link{prepare_data}}
#' @param \dots Further arguments to print
#'
#' @return No return value
#'
#' @examples
#' # Bin the colon cancer data over s (time since recurrence)
#'
#' dt1ts <- prepare_data(data = reccolon2ts,
#'                       s_in = "entrys",
#'                       s_out = "timesr",
#'                       events = "status",
#'                       ds = 180)
#'
#' print(dt1ts)
#'
#' # Bin the colon cancer data over u (time at recurrence) and s (time since recurrence)
#' dt2ts <- prepare_data(data = reccolon2ts,
#'                       u = "timer",
#'                       s_in = "entrys",
#'                       s_out = "timesr",
#'                       events = "status",
#'                       ds = 180)
#' print(dt2ts)
#'
#' @author Angela Carollo \email{carollo@@demogr.mpg.de}
#'
#' @importFrom  utils str
#'
#' @export


print.data2ts <- function(x,...)
{
  if (!inherits(x, "data2ts"))
    stop("'x' must be an 'data2ts' object")
  cat("An object of class 'data2ts'\n\nData:\n")
  print(str(x, max.level = 1))
  cat("\nRange covered by the bins: \n")
  if(length(x$bins)>3)
    print(lapply(x$bins[c(1,4)], range))
  else
    print(lapply(x$bins[1], range))
  cat("\nNumber of bins: \n")
  if(length(x$bins)>3)
    print(x$bins[c(3,6)])
  else
    print(x$bins[c(3)])
  cat("\nOverview of the binned data:\n")
  cat(paste0("Total exposure time: ", lapply(x$bindata[1], sum)))
  cat("\n")
  cat(paste0("Total number of events: ", lapply(x$bindata[2], sum)))
  # if covariates
  if(!is.null(x$bindata$Z)){
    cat("\nCovariates:\n")
    print(colnames(x$bindata$Z))
  }
  return(invisible())
}

