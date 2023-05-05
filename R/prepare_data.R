#' Prepare raw data by binning them in 1d or 2d
#'
#' @description `prepare_data()` prepares the raw individual time-to-event data
#' for hazard estimation in 1d or 2d.
#'
#' Given the raw data, this function first constructs the bins over one or two
#'   axes and then computes the aggregated (or individual)
#'   vectors or matrices of exposure times and events indicators.
#'
#' @inheritParams make_bins
#' @param events A vector of event's indicators (possible values 0/1).
#' @param individual A Boolean. Default is `FALSE`: if `FALSE` computes the matrices `R` and `Y`
#'   collectively for all observations; if `TRUE` computes the matrices `R` and `Y` separately for each individual record.
#' @param covs A data.frame with the variables to be used as covariates.
#'   The function will create dummy variables for any factor variable passed as argument in `covs`.
#'   If a variable of class character is passed as argument, it will be converted to factor.
#'
#' @return A list with the following elements:
#' * `bins` a list:
#'   * `bins_t` if `t_out` is provided, this is a vector of bins extremes for the time scale `t`.
#'   * `mid_t` if `t_out` is provided, this is a vector with the midpoints of the bins over `t`.
#'   * `nt` if `t_out` is provided, this is the number of bins over `t`.
#'   * `bins_u` if `u` is provided, this is a vector of bins extremes for `u` axis.
#'   * `midu` if `u` is provided, this is a vector with the midpoints of the bins over `u`.
#'   * `nu` if `u` is provided, this is the number of bins over `u`.
#'   * `bins_s` is a vector of bins extremes for the time scale `s`.
#'   * `mids` is a vector with the midpoints of the bins over `s`.
#'   * `ns` is the number of bins over `s`.
#' * `bindata`:
#'   * `r` or `R` an array of exposure times: if binning the data over
#'      one time scale only this is a vector.
#'      If binning the data over two time scales and if `individual == TRUE`
#'      then `R` is an array of dimension nu by ns by n, otherwise it is an
#'      array of dimension nu by ns
#'   * `y` or `Y` an array of event counts: if binning the data over one time
#'      scale only this is a vector.
#'      If binning the data over two time scales and if `individual == TRUE`
#'      then `Y` is an array of dimension nu by ns by n, otherwise it is an
#'      array of dimension nu by ns
#'   * `Z` A matrix of covariates' values to be used in the model,
#'      of dimension n by p
#'
#' @examples
#'
#' # Bin data over s = time since recurrence only, with intervals of length 30 days
#' # aggregated data (no covariates)
#' binned_data <- prepare_data(s_out = reccolon2ts$timesr, events = reccolon2ts$status, ds = 30)
#' # Visualize vector of event counts
#' print(binned_data$bindata$y)
#' # Visualize midpoints of the bins
#' print(binned_data$bins$mids)
#' # Visualize number of bins
#' print(binned_data$bins$ns)
#'
#'
#' # Now using ds = .3 and the same variable measured in years
#' binned_data <- prepare_data(s_out = reccolon2ts$timesr_y, events = reccolon2ts$status, ds = .3)
#' # Visualize vector of exposure timess
#' print(binned_data$bindata$r)
#'
#'
#' # Bin data over u = time at recurrence and s = time since recurrence, measured in days
#' # aggregated data (no covariates)
#' # Note that if we do not provide du this is taken to be equal to ds
#' binned_data <- prepare_data(
#'   u = reccolon2ts$timer, s_out = reccolon2ts$timesr,
#'   events = reccolon2ts$status, ds = 30
#' )
#'
#' # Visualize matrix of event counts
#' print(binned_data$bindata$Y)
#'
#' # Visualize midpoints of bins over u
#' print(binned_data$bins$midu)
#'
#'
#' # Bin data over u = time at recurrence and s = time since recurrence, measured in day
#' # individual-level data required
#' # we provide two covariates: nodes (numerical) and rx (factor)
#' covs <- subset(reccolon2ts, select = c("nodes", "rx"))
#' binned_data <- prepare_data(
#'   u = reccolon2ts$timer, s_out = reccolon2ts$timesr,
#'   events = reccolon2ts$status, ds = 30, individual = TRUE, covs = covs
#' )
#'
#' # Visualize structure of binned data
#' print(str(binned_data$bindata))
#'
#' @export

prepare_data <- function(t_in = NULL, t_out = NULL,
                         u = NULL,
                         s_in = NULL, s_out,
                         events,
                         min_t = NULL, max_t = NULL,
                         min_u = NULL, max_u = NULL,
                         min_s = NULL, max_s = NULL,
                         dt = NULL, du = NULL, ds,
                         individual = FALSE,
                         covs = NULL) {

  # The check on all inputs is done by the individual functions called by
  # prepare_data(), and will be here omitted.

  # ---- Create bins ----
  bins <- make_bins(
    t_in = t_in, t_out = t_out,
    u = u,
    s_in = s_in, s_out = s_out,
    min_t = min_t, max_t = max_t,
    min_u = min_u, max_u = max_u,
    min_s = min_s, max_s = max_s,
    dt = dt, du = du, ds = ds
  )

  # ---- Exposure and events calculation ----
  n <- length(s_out)

  # Only s provided
  if (is.null(t_out) & is.null(u)) {

    ry <- exposures_events_1d(
      s_in = s_in,
      s_out = s_out,
      ev = events,
      bins = bins$bins_s
    )
    if(individual){
      bindata <- list(
        R = t(ry$R),
        Y = t(ry$Y)
      )
    } else {
      bindata <- list(
        r = ry$r,
        y = ry$y
      )
    }
  }

  # t and s provided
  if (!is.null(t_out)) {
    if (individual) {
      ns <- bins$ns
      nt <- bins$nu

      R <- Y <- array(NA, dim = c(nt, ns, n))

      for (ind in 1:n) {
        RY_ind <- exposures_events_Lexis(
          t_in = t_in[ind],
          t_out = t_out[ind],
          s_in = s_in[ind],
          s_out = s_out[ind],
          ev = events[ind],
          bins_list = bins
        )
        R[, , ind] <- RY_ind$R
        Y[, , ind] <- RY_ind$Y
      }
      bindata <- list(
        R = R,
        Y = Y
      )
    } else {
      bindata <- exposures_events_Lexis(
        t_in = t_in,
        t_out = t_out,
        s_in = s_in,
        s_out = s_out,
        ev = events,
        bins_list = bins
      )
    }
  }

  # u and s provided
  if (!is.null(u)) {
    bindata <- exposures_events_2d(
      u = u,
      s_in = s_in,
      s_out = s_out,
      ev = events,
      individual = individual,
      bins_list = bins
    )
  }

  # ---- Prepare covariates' matrix ----
  # Returns warning is individual == FALSE

  if (!is.null(covs)) {

    if (individual == FALSE) {
      warning("Covariates provided but `individual == FALSE`.")
      message("I will nevertheless prepare the regression matrix `Z`.")
    }

    Z <- matrix(NA,nrow=nrow(covs),ncol=0)
    namesZ <- vector(mode="character", length=0L)
    namesc <- attributes(covs)$names
    # how many columns?
    for (i in 1:ncol(covs)) {
      covtemp <- covs[, i]
      if (is.numeric(covtemp)) {
        Z <- cbind(Z, covtemp)
        namesZ <- c(namesZ, namesc[i])
      } # if numeric copy the variable
      else {
        if (is.factor(covtemp)) {
          nlev <- nlevels(covtemp)
          ndum <- nlev - 1
          lev <- levels(covtemp)
          for (j in 2:(nlev)) {
            z <- as.numeric(covtemp == lev[j])
            Z <- cbind(Z, z)
          }
          namesZ <- c(namesZ, lev[2:nlev])
        } else {
          if (is.character(covtemp)) {
            covtemp <- factor(covtemp,
                              levels = sort(unique(covtemp)))
            nlev <- nlevels(covtemp)
            ndum <- nlev - 1
            lev <- levels(covtemp)
            for (j in 2:(nlev)) {
              z <- as.numeric(covtemp == lev[j])
              Z <- cbind(Z, z)
            }
            namesZ <- c(namesZ, lev[2:nlev])
          }
        }
      }
      colnames(Z) <- namesZ
    }

    bindata$Z <- Z
  }

  if (individual == TRUE & is.null(covs)) {
    warning("Individual-level data requested, but covariates are not provided.\n
            I assume that the regression matrix `Z` is prepared externally.")
  }


  # ---- Return results ----
  data.binned <- list(
    "bins" = bins,
    "bindata" = bindata
  )

  class(data.binned) <- "data2ts"
  return(data.binned)
}
