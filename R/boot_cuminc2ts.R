#' Bootstrap confidence intervals for cumulative incidence functions with two time scales
#'
#' @description
#' `boot_cuminc2ts()` is a wrapper function to perform non-parametric bootstrap
#' in order to obtain uncertainty measures for the cumulative incidence functions
#' (and the overall survival function) over two time scales. It allows parallel
#' computing.
#'
#' @param data A data.frame containing the original individual-level data.
#' @param causes A character vector of column names in \code{data}, one for each
#'   competing cause (binary 0/1 indicators).
#' @param cause_names A character vector of names for the causes, passed to
#'   \code{cumhaz2ts()} and \code{cuminc2ts()}. Defaults to the values in
#'   \code{causes}.
#' @param prepare_data_args A named list of arguments passed to
#'   \code{prepare_data()}, excluding \code{data} and \code{events}.
#' @param fit2ts_args A named list of arguments passed to \code{fit2ts()},
#'   excluding \code{data2ts}.
#' @param cumhaz2ts_args A named list of arguments passed to
#'   \code{cumhaz2ts()}, excluding \code{fitted_model} and \code{cause}.
#' @param ds The bin width in the \code{s} direction, passed to
#'   \code{cuminc2ts()}.
#' @param nboot Integer. Number of bootstrap replicates. Default is 200.
#' @param seed Integer or NULL. Random seed for reproducibility. Default is
#'   NULL.
#' @param conf_level Numeric. Confidence level for the intervals. Default is
#'   0.95.
#' @param parallel Logical. Whether to use parallel computation. Default is
#'   FALSE.
#' @param ncpus Integer. Number of CPU cores to use when \code{parallel =
#'   TRUE}. Default is 2.
#'
#' @return A list of two lists: the first list contains the results of the bootstrap
#'  with one element per cause. Each element is itself a list
#'   containing:
#'   \describe{
#'     \item{\code{lower}}{Matrix of lower confidence bounds.}
#'     \item{\code{upper}}{Matrix of upper confidence bounds.}
#'     \item{\code{se}}{Matrix of pointwise standard errors.}
#'     \item{\code{boot_replicates}}{List of \code{nboot} matrices of bootstrap
#'       CIF estimates.}
#'   }
#'   The second list is the grid of points in correspondence of which the estimates
#'   are obtained. This is returned to facilitate plotting of these quantities.
#'
#' @details
#' It may happen that, as a consequence of the resampling, the range of values
#' for both `u` and  `s` differ across the various bootstrap samples. To avoid
#' problems it is necessary to bin the data, fit the model, and compute the
#' cumulative incidence functions on the same range of values across all the samples.
#' To do so, simply specify `min_u`, `max_u`, `min_s`, and `max_s` wherever required.
#' This is shown in the example below.
#'
#' @importFrom stats quantile sd
#' @importFrom utils setTxtProgressBar txtProgressBar
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
#' \donttest{
#' # --- prepare_data for each cause -----------------------------------------
#' rec2ts <- prepare_data(
#'   data   = fakedata,
#'   u      = "u",
#'   s_out  = "s_out",
#'   events = "cause1",
#'   min_u  = 24, max_u = 58,
#'   min_s  = 0,  max_s = 10,
#'   du     = 1,  ds    = .5
#' )
#'
#' death2ts <- prepare_data(
#'   data   = fakedata,
#'   u      = "u",
#'   s_out  = "s_out",
#'   events = "cause2",
#'   min_u  = 24, max_u = 58,
#'   min_s  = 0,  max_s = 10,
#'   du     = 1,  ds    = .5
#' )
#'
#' # --- Fit cause-specific hazard models ------------------------------------
#' mod_cause1 <- fit2ts(
#'   rec2ts,
#'   Bbases_spec = list(
#'     bdeg   = 3,
#'     nseg_u = 7, min_u = 24, max_u = 58,
#'     nseg_s = 3, min_s = 0,  max_s = 10
#'   ),
#'   optim_criterion = "bic"
#' )
#'
#' mod_cause2 <- fit2ts(
#'   death2ts,
#'   Bbases_spec = list(
#'     bdeg   = 3,
#'     nseg_u = 7, min_u = 24, max_u = 58,
#'     nseg_s = 3, min_s = 0,  max_s = 10
#'   ),
#'   optim_criterion = "bic"
#' )
#'}
#' # --- Bootstrap confidence intervals --------------------------------------
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
#'       nseg_s = 3, min_s = 0,  max_s = 10
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
#'   nboot      = 10,
#'   seed       = 1234,
#'   conf_level = 0.95,
#'   parallel   = FALSE
#' )
#'
#' # --- Inspect output ------------------------------------------------------
#' # Names of causes
#' names(boot_cif)
#'
#' # Dimensions of the lower confidence bound matrix for cause 1
#' dim(boot_cif$results_cif[["cause1"]]$lower)
#'
#' # Pointwise standard errors for cause 2 (first few rows and columns)
#' boot_cif$results_cif[["cause2"]]$se[1:5, 1:5]
#'
#' @export
#'
boot_cuminc2ts <- function(data,
                           causes,
                           cause_names = NULL,
                           prepare_data_args = list(),
                           fit2ts_args = list(),
                           cumhaz2ts_args = list(),
                           ds,
                           nboot = 200,
                           seed = NULL,
                           conf_level = 0.95,
                           parallel = FALSE,
                           ncpus = 2) {

  # ---- Input checks --------------------------------------------------------
  if (!is.data.frame(data)) {
    stop("`data` must be a data.frame.")
  }
  if (!is.character(causes) || length(causes) < 2) {
    stop("`causes` must be a character vector with at least two cause names.")
  }
  if (!all(causes %in% names(data))) {
    stop("All elements of `causes` must be column names in `data`.")
  }
  if (!is.list(prepare_data_args)) {
    stop("`prepare_data_args` must be a named list.")
  }
  if (!is.list(fit2ts_args)) {
    stop("`fit2ts_args` must be a named list.")
  }
  if (!is.list(cumhaz2ts_args)) {
    stop("`cumhaz2ts_args` must be a named list.")
  }
  if (!is.numeric(nboot) || nboot < 1) {
    stop("`nboot` must be a positive integer.")
  }
  if (!is.numeric(conf_level) || conf_level <= 0 || conf_level >= 1) {
    stop("`conf_level` must be a number strictly between 0 and 1.")
  }
  if (!is.logical(parallel)) {
    stop("`parallel` must be logical.")
  }
  if (parallel && (!is.numeric(ncpus) || ncpus < 1)) {
    stop("`ncpus` must be a positive integer when `parallel = TRUE`.")
  }

  # ---- Setup ---------------------------------------------------------------
  if (is.null(cause_names)) cause_names <- causes
  if (!is.null(seed)) set.seed(seed)

  n <- nrow(data)
  alpha <- 1 - conf_level
  n_causes <- length(causes)

  # ---- Define single bootstrap iteration -----------------------------------
  one_boot <- function(b) {

    # 1. Resample individuals
    idx <- sample(seq_len(n), size = n, replace = TRUE)
    boot_data <- data[idx, ]

    # 2. & 3. prepare_data + fit2ts for each cause
    haz_list <- vector("list", n_causes)

    for (k in seq_len(n_causes)) {

      # prepare_data for cause k
      pd_args_k <- c(list(data = boot_data, events = causes[k]),
                     prepare_data_args)
      data2ts_k <- tryCatch(
        do.call(prepare_data, pd_args_k),
        error = function(e) NULL
      )
      if (is.null(data2ts_k)) return(NULL)

      # fit2ts for cause k
      fit_args_k <- c(list(data2ts = data2ts_k), fit2ts_args)
      mod_k <- tryCatch(
        do.call(fit2ts, fit_args_k),
        error = function(e) NULL
      )
      if (is.null(mod_k)) return(NULL)

      # 4. cumhaz2ts for cause k
      ch_args_k <- c(list(fitted_model = mod_k, cause = cause_names[k]),
                     cumhaz2ts_args)
      H_k <- tryCatch(
        do.call(cumhaz2ts, ch_args_k),
        error = function(e) NULL
      )
      if (is.null(H_k)) return(NULL)

      haz_list[[k]] <- H_k$Haz$hazard
    }

    Surv2ts <- exp(-(Reduce("+", haz_list)))


      # 6. cuminc2ts
    cif_b <- tryCatch(
      cuminc2ts(haz = haz_list, ds = ds, cause = cause_names),
      error = function(e) NULL
    )
    cif_b$Surv2ts <- Surv2ts
    return(cif_b)
  }

  # ---- Run bootstrap -------------------------------------------------------
  if (parallel) {
    cl <- parallel::makeCluster(ncpus)
    on.exit(parallel::stopCluster(cl), add = TRUE)

    # Export everything the worker nodes need
    parallel::clusterExport(
      cl,
      varlist = c("data", "causes", "cause_names", "prepare_data_args",
                  "fit2ts_args", "cumhaz2ts_args", "ds", "n", "n_causes",
                  "prepare_data", "fit2ts", "cumhaz2ts", "cuminc2ts",
                  "one_boot"),
      envir = environment()
    )
    parallel::clusterSetRNGStream(cl, seed)

    message(sprintf("Running %d bootstrap replicates on %d cores...", nboot, ncpus))
    boot_results <- parallel::parLapply(cl, seq_len(nboot), one_boot)

  } else {

    pb <- txtProgressBar(min = 0, max = nboot, style = 3)
    boot_results <- vector("list", nboot)

    for (b in seq_len(nboot)) {
      boot_results[[b]] <- one_boot(b)
      setTxtProgressBar(pb, b)
    }
    close(pb)
  }

  # ---- Handle failed replicates --------------------------------------------
  failed <- vapply(boot_results, is.null, logical(1))
  n_failed <- sum(failed)

  if (n_failed > 0) {
    warning(sprintf(
      "%d out of %d bootstrap replicates failed and were discarded.",
      n_failed, nboot
    ))
  }

  boot_results <- boot_results[!failed]
  n_success <- length(boot_results)

  if (n_success < 2) {
    stop("Fewer than 2 bootstrap replicates succeeded. Cannot compute intervals.")
  }

  # ---- Compute summaries per cause -----------------------------------------
  results_cif <- vector("list", n_causes)
  names(results_cif) <- cause_names

  # Extract grid coordinates from the first successful replicate
  # These are fixed across all replicates since plot_grid is fixed

  first_haz <- tryCatch(
    {
      pd_args <- c(list(data = data, events = causes[1]), prepare_data_args)
      data2ts <- do.call(prepare_data, pd_args)
      fit_args <- c(list(data2ts = data2ts), fit2ts_args)
      mod <- do.call(fit2ts, fit_args)
      ch_args <- c(list(fitted_model = mod, cause = cause_names[1]),
                   cumhaz2ts_args)
      do.call(cumhaz2ts, ch_args)
    },
    error = function(e) NULL
  )

  grid <- if (!is.null(first_haz)) {
    list(
      u = first_haz$plot_grid$intu,
      s = first_haz$plot_grid$ints
    )
  } else {
    warning("Could not extract grid coordinates from fitted model.")
    NULL
  }

  # boot_results is a list of nboot cuminc2ts objects (lists of matrices)
  # We need to reorganize by cause first, then compute quantiles and SE

  for (k in seq_len(n_causes)) {

    # Extract the k-th CIF matrix from each successful replicate
    # cuminc2ts returns a list of matrices, one per cause
    boot_mats <- lapply(boot_results, function(cif_b) cif_b[[k]])

    # Stack into a 3D array: rows x cols x B
    cif_array <- tryCatch(
      simplify2array(boot_mats),
      error = function(e) {
        stop(sprintf(
          "Bootstrap CIF matrices for cause '%s' have inconsistent dimensions.",
          cause_names[k]
        ))
      }
    )


    # Pointwise quantiles and SE
    lower_mat <- apply(cif_array, c(1, 2), quantile,
                       probs = alpha / 2, na.rm = TRUE)
    upper_mat <- apply(cif_array, c(1, 2), quantile,
                       probs = 1 - alpha / 2, na.rm = TRUE)
    se_mat    <- apply(cif_array, c(1, 2), sd, na.rm = TRUE)

    results_cif[[k]] <- list(
      lower           = lower_mat,
      upper           = upper_mat,
      se              = se_mat,
      boot_replicates = boot_mats
    )
  }

  surv_mats <- lapply(boot_results, function(cif_b) cif_b$Surv2ts)
  # Stack into a 3D array: rows x cols x B
  surv_array <- tryCatch(
    simplify2array(surv_mats),
    error = function(e) {
      stop(
        "Bootstrap overall survival matrix has inconsistent dimensions."
      )
    }
  )

  surv_se_mat <- apply(surv_array, c(1, 2), sd, na.rm = TRUE)
  surv_lower_mat <- apply(surv_array, c(1, 2), quantile,
                          probs = alpha / 2, na.rm = TRUE)
  surv_upper_mat <- apply(surv_array, c(1, 2), quantile,
                          probs = 1 - alpha / 2, na.rm = TRUE)
  results_surv = list(
    lower      = surv_lower_mat,
    upper      = surv_upper_mat,
    se         = surv_se_mat,
    boot_replicates = surv_mats
  )

  message(sprintf(
    "Done. %d successful replicates out of %d attempted.",
    n_success, nboot
  ))

  return(list(
    results_cif = results_cif,
    results_surv = results_surv,
    grid    = grid
  ))
}
