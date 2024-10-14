#' Bin data on the Lexis diagram
#'
#' @description `exposures_events_Lexis()` computes aggregated matrices of
#' exposure times and event counts over two time
#' scales, on the Lexis diagram.
#'
#' The time scales are `t` and `s`. This function uses functions from
#' the package `popEpi` and from the package `Epi`, and code shared by Bendix Carstensen
#' on the website bendixcarstensen.com. See also [prepare_data()]
#' to conveniently prepare individual data for the analysis with one,
#' or two time scales.
#'
#' @inheritParams make_bins
#' @param ev A vector of event's indicators (possible values 0/1).
#' @param bins_list A list with the following (necessary) elements:
#' * `bins_t` a vector of extreme values for the bins over the `t` axis.
#' * `nt` the number of bins over `t`.
#' * `bins_s` a vector of extreme values for the bins over the `t` axis.
#' * `ns` the number of bins over `s`.
#' @return A list with the following elements:
#' * `R` an array of exposure times of dimension \eqn{nt} by \eqn{ns}
#' * `Y` an array of event counts of dimension \eqn{nt} by \eqn{ns}
#'
#' @author Angela Carollo \email{carollo@@demogr.mpg.de}
#' @references Carstensen B, Plummer M, Laara E, Hills M (2022).
#' Epi: A Package for Statistical Analysis in Epidemiology.
#' R package version 2.47.1, https://CRAN.R-project.org/package=Epi.
#'
#' Miettinen J, Rantanen M, Seppa K (2023).
#' popEpi: Functions for Epidemiological Analysis using Population Data.
#' R package version 0.4.11, https://cran.r-project.org/package=popEpi.
#'
#'
#' @export

exposures_events_Lexis <-
  function(t_in = NULL,
           t_out,
           s_in = NULL,
           s_out,
           ev,
           bins_list) {


    # ---- Checks on arguments ----
    if (length(bins_list) == 1) {
      stop("`bins_list` must be a list with 2 elements.")
    }

    if(is.null(bins_list$bins_t)) stop("`bins_list$bins_t` missing.")

    if (!is.null(s_in) & (length(s_in) != length(s_out))) {
      stop("`s_in`and `s_out` must have the same length.")
    }

    if (!is.null(t_in) & (length(t_in) != length(t_out))) {
      stop("`t_in`and `t_out` must have the same length.")
    }

    if (length(s_out) != length(t_out)) {
      stop("`s_out` and `t_out` must have the same length.")
    }

    if (length(s_out) != length(ev)) {
      stop("`s_out` and `ev` must have the same length.")
    }

    if (any(ev > 1)) {
      stop("only 0 or 1 allowed in vector `ev`.")
    }

    if(is.null(t_in)){
      t_in <- t_out - s_out
      message("`t_in` not provided. I will use `t_in = t_out - s_in`.")
    }

    if (is.null(s_in)) {
      s_in <- rep(0, length(s_out))
      message("`s_in = NULL`. I will use `s_in = 0` for all observations.")
    }

    # ---- Compute R and Y ---
    # Transform to Lexis object
    LexData <- Epi::Lexis(
      entry = list(t = t_in, s = s_in),
      exit = list(t = t_out, s = s_out),
      exit.status = ev
    )
    LexDataSplit <-
      popEpi::splitMulti(LexData,
        breaks = list(
          t = bins_list$bins_t,
          s = bins_list$bins_s
        )
      )
    LexDataSplit$which.it <-
      cut(
        LexDataSplit$t,
        breaks = bins_list$bins_t,
        include.lowest = T,
        right = T
      )
    LexDataSplit$which.is <-
      cut(
        LexDataSplit$s,
        breaks = bins_list$bins_s,
        include.lowest = T,
        right = T
      )

    # Tabulate events and person-years
    Tabulation <-
      with(LexDataSplit, xtabs(cbind(
        Y = lex.Xst,
        R = lex.dur
      ) ~
        which.it +
        which.is))
    Y <- Tabulation[, , "Y"]
    R <- Tabulation[, , "R"]

    Y <- matrix(Y,
      nrow = bins_list$nt,
      ncol = bins_list$ns
    )
    R <- matrix(R,
      nrow = bins_list$nt,
      ncol = bins_list$ns
    )

    return(list(R = R, Y = Y))
  }
