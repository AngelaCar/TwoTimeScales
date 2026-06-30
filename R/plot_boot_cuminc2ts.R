#' Image Plot of Bootstrap Uncertainty for Cumulative Incidence Functions
#'
#' `plot_boot_cuminc2ts()` plots an image of the pointwise bootstrap
#'   uncertainty measures (lower bound, upper bound, or standard errors)
#'   for the cumulative incidence function of a single cause, over two
#'   time scales.
#'
#' @param boot_object The object returned by \code{boot_cuminc2ts()}, a list
#'   with elements \code{results} (one per cause) and \code{grid}.
#' @param which_cause A character string or integer indicating which cause
#'   to plot. If a character string, it must match one of the names in
#'   \code{boot_object$results}. If an integer, it is used as an index.
#' @param which_quantity A character string specifying which uncertainty
#'   measure to plot. One of \code{"lower"}, \code{"upper"}, or \code{"se"}.
#'   Default is \code{"se"}.
#' @param plot_options A list of options passed to \code{imageplot_SE()}.
#'   See \code{\link{imageplot_SE}} for the full list of options.
#' @param \dots Further arguments passed to \code{imageplot_SE()}.
#'
#' @return An image plot of the selected bootstrap uncertainty measure for
#'   the cumulative incidence function of the selected cause.
#'
#' @seealso \code{\link{boot_cuminc2ts}}, \code{\link{imageplot_SE}}
#'
#' @examples
#' # --- Fake data -----------------------------------------------------------
#' set.seed(1234)
#' n <- 30
#' fakedata <- data.frame(
#'   id     = 1:n,
#'   u      = round(runif(n, min = 24, max = 58), 2),
#'   s_out  = round(runif(n, min = 0.5, max = 10), 2),
#'   cause1 = c(rep(1, 8), rep(0, 22)),
#'   cause2 = c(rep(0, 8), rep(1, 7), rep(0, 15))
#' )
#'
#' # --- Bootstrap -----------------------------------------------------------
#' boot_cif <- boot_cuminc2ts(
#'   data     = fakedata,
#'   causes   = c("cause1", "cause2"),
#'   cause_names = c("cause1", "cause2"),
#'   prepare_data_args = list(
#'     u      = "u",
#'     s_out  = "s_out",
#'     min_u  = 24, max_u = 58,
#'     min_s  = 0,  max_s = 10,
#'     du     = 1,  ds    = .5
#'   ),
#'   fit2ts_args = list(
#'     Bbases_spec = list(
#'       bdeg   = 3,
#'       nseg_u = 7, min_u = 24, max_u = 58,
#'       nseg_s = 3, min_s = 0, max_s = 10
#'     ),
#'     optim_criterion = "bic"
#'   ),
#'   cumhaz2ts_args = list(
#'     plot_grid = list(
#'       c(umin = 24, umax = 58, du = .5),
#'       c(smin = 0,  smax = 10, ds = .2)
#'     )
#'   ),
#'   ds         = .2,
#'   nboot          = 10,
#'   seed       = 1234,
#'   conf_level = 0.95,
#'   parallel   = FALSE
#' )
#'
#' # --- Plot standard errors for cause 1 ------------------------------------
#' plot_boot_cuminc2ts(
#'   boot_object   = boot_cif,
#'   which_cause   = "cause1",
#'   which_quantity = "se",
#'   plot_options  = list(
#'     main = "Bootstrap SE - Cause 1",
#'     xlab = "Age",
#'     ylab = "Time since entry"
#'   )
#' )
#'
#' @export
plot_boot_cuminc2ts <- function(boot_object,
                                which_cause,
                                which_quantity = c("se", "lower", "upper"),
                                plot_options = list(),
                                ...) {

  # ---- Input checks --------------------------------------------------------
  if (!is.list(boot_object) ||
      !all(c("results_cif", "grid") %in% names(boot_object))) {
    stop(paste("`boot_object` must be the object returned by",
               "`boot_cuminc2ts()`, with elements `results_cif` and `grid`."))
  }

  if (is.null(boot_object$grid)) {
    stop("Grid coordinates are missing from `boot_object$grid`.")
  }

  which_quantity <- match.arg(which_quantity)

  # Resolve which_cause to a name
  cause_names <- names(boot_object$results_cif)

  if (is.numeric(which_cause)) {
    if (which_cause < 1 || which_cause > length(cause_names)) {
      stop(sprintf(
        "`which_cause` index %d is out of range. There are %d causes.",
        which_cause, length(cause_names)
      ))
    }
    which_cause <- cause_names[which_cause]
  } else if (is.character(which_cause)) {
    if (!which_cause %in% cause_names) {
      stop(sprintf(
        "`which_cause` '%s' not found. Available causes are: %s.",
        which_cause, paste(cause_names, collapse = ", ")
      ))
    }
  } else {
    stop("`which_cause` must be a character string or an integer index.")
  }

  # ---- Extract components --------------------------------------------------
  cause_result <- boot_object$results_cif[[which_cause]]
  z            <- cause_result[[which_quantity]]
  x            <- boot_object$grid$u
  y            <- boot_object$grid$s

  if (is.null(z)) {
    stop(sprintf(
      "Quantity '%s' is missing for cause '%s'.",
      which_quantity, which_cause
    ))
  }

  # ---- Default plot title if not provided ----------------------------------
  if (is.null(plot_options$main)) {
    quantity_label <- switch(which_quantity,
                             "se"    = "Bootstrap SE",
                             "lower" = "Lower confidence bound",
                             "upper" = "Upper confidence bound"
    )
    plot_options$main <- sprintf("%s - %s", quantity_label, which_cause)
  }

  # ---- Plot ----------------------------------------------------------------
  imageplot_SE(
    x            = x,
    y            = y,
    z            = z,
    plot_options = plot_options,
    ...
  )
}
