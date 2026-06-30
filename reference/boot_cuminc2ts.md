# Bootstrap confidence intervals for cumulative incidence functions with two time scales

`boot_cuminc2ts()` is a wrapper function to perform non-parametric
bootstrap in order to obtain uncertainty measures for the cumulative
incidence functions (and the overall survival function) over two time
scales. It allows parallel computing.

## Usage

``` r
boot_cuminc2ts(
  data,
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
  ncpus = 2
)
```

## Arguments

- data:

  A data.frame containing the original individual-level data.

- causes:

  A character vector of column names in `data`, one for each competing
  cause (binary 0/1 indicators).

- cause_names:

  A character vector of names for the causes, passed to
  [`cumhaz2ts()`](https://angelacar.github.io/TwoTimeScales/reference/cumhaz2ts.md)
  and
  [`cuminc2ts()`](https://angelacar.github.io/TwoTimeScales/reference/cuminc2ts.md).
  Defaults to the values in `causes`.

- prepare_data_args:

  A named list of arguments passed to
  [`prepare_data()`](https://angelacar.github.io/TwoTimeScales/reference/prepare_data.md),
  excluding `data` and `events`.

- fit2ts_args:

  A named list of arguments passed to
  [`fit2ts()`](https://angelacar.github.io/TwoTimeScales/reference/fit2ts.md),
  excluding `data2ts`.

- cumhaz2ts_args:

  A named list of arguments passed to
  [`cumhaz2ts()`](https://angelacar.github.io/TwoTimeScales/reference/cumhaz2ts.md),
  excluding `fitted_model` and `cause`.

- ds:

  The bin width in the `s` direction, passed to
  [`cuminc2ts()`](https://angelacar.github.io/TwoTimeScales/reference/cuminc2ts.md).

- nboot:

  Integer. Number of bootstrap replicates. Default is 200.

- seed:

  Integer or NULL. Random seed for reproducibility. Default is NULL.

- conf_level:

  Numeric. Confidence level for the intervals. Default is 0.95.

- parallel:

  Logical. Whether to use parallel computation. Default is FALSE.

- ncpus:

  Integer. Number of CPU cores to use when `parallel = TRUE`. Default is
  2.

## Value

A list of two lists: the first list contains the results of the
bootstrap with one element per cause. Each element is itself a list
containing:

- `lower`:

  Matrix of lower confidence bounds.

- `upper`:

  Matrix of upper confidence bounds.

- `se`:

  Matrix of pointwise standard errors.

- `boot_replicates`:

  List of `nboot` matrices of bootstrap CIF estimates.

The second list is the grid of points in correspondence of which the
estimates are obtained. This is returned to facilitate plotting of these
quantities.

## Details

It may happen that, as a consequence of the resampling, the range of
values for both `u` and `s` differ across the various bootstrap samples.
To avoid problems it is necessary to bin the data, fit the model, and
compute the cumulative incidence functions on the same range of values
across all the samples. To do so, simply specify `min_u`, `max_u`,
`min_s`, and `max_s` wherever required. This is shown in the example
below.

## Examples

``` r
# --- Fake data -----------------------------------------------------------
set.seed(1234)
n <- 30
fakedata <- data.frame(
  id     = 1:n,
  u      = round(runif(n, min = 24, max = 58), 2),
  s_out  = round(runif(n, min = 0.5, max = 10), 2),
  cause1 = c(rep(1, 8), rep(0, 22)),
  cause2 = c(rep(0, 8), rep(1, 7), rep(0, 15))
)

# \donttest{
# --- prepare_data for each cause -----------------------------------------
rec2ts <- prepare_data(
  data   = fakedata,
  u      = "u",
  s_out  = "s_out",
  events = "cause1",
  min_u  = 24, max_u = 58,
  min_s  = 0,  max_s = 10,
  du     = 1,  ds    = .5
)
#> `s_in = NULL`. I will use `s_in = 0` for all observations.
#> `s_in = NULL`. I will use `s_in = 0` for all observations.

death2ts <- prepare_data(
  data   = fakedata,
  u      = "u",
  s_out  = "s_out",
  events = "cause2",
  min_u  = 24, max_u = 58,
  min_s  = 0,  max_s = 10,
  du     = 1,  ds    = .5
)
#> `s_in = NULL`. I will use `s_in = 0` for all observations.
#> `s_in = NULL`. I will use `s_in = 0` for all observations.

# --- Fit cause-specific hazard models ------------------------------------
mod_cause1 <- fit2ts(
  rec2ts,
  Bbases_spec = list(
    bdeg   = 3,
    nseg_u = 7, min_u = 24, max_u = 58,
    nseg_s = 3, min_s = 0,  max_s = 10
  ),
  optim_criterion = "bic"
)

mod_cause2 <- fit2ts(
  death2ts,
  Bbases_spec = list(
    bdeg   = 3,
    nseg_u = 7, min_u = 24, max_u = 58,
    nseg_s = 3, min_s = 0,  max_s = 10
  ),
  optim_criterion = "bic"
)
# }
# --- Bootstrap confidence intervals --------------------------------------
boot_cif <- boot_cuminc2ts(
  data     = fakedata,
  causes   = c("cause1", "cause2"),
  cause_names = c("cause1", "cause2"),
  prepare_data_args = list(
    u      = "u",
    s_out  = "s_out",
    min_u  = 24, max_u = 58,
    min_s  = 0,  max_s = 10,
    du     = 1,  ds    = .5
  ),
  fit2ts_args = list(
    Bbases_spec = list(
      bdeg   = 3,
      nseg_u = 7, min_u = 24, max_u = 58,
      nseg_s = 3, min_s = 0,  max_s = 10
    ),
    optim_criterion = "bic"
  ),
  cumhaz2ts_args = list(
    plot_grid = list(
      c(umin = 24, umax = 58, du = .5),
      c(smin = 0,  smax = 10, ds = .2)
    )
  ),
  ds         = .2,
  nboot      = 10,
  seed       = 1234,
  conf_level = 0.95,
  parallel   = FALSE
)
#>   |                                                                              |                                                                      |   0%
#> `s_in = NULL`. I will use `s_in = 0` for all observations.
#> `s_in = NULL`. I will use `s_in = 0` for all observations.
#> `s_in = NULL`. I will use `s_in = 0` for all observations.
#> `s_in = NULL`. I will use `s_in = 0` for all observations.
#>   |                                                                              |=======                                                               |  10%
#> `s_in = NULL`. I will use `s_in = 0` for all observations.
#> `s_in = NULL`. I will use `s_in = 0` for all observations.
#> `s_in = NULL`. I will use `s_in = 0` for all observations.
#> `s_in = NULL`. I will use `s_in = 0` for all observations.
#>   |                                                                              |==============                                                        |  20%
#> `s_in = NULL`. I will use `s_in = 0` for all observations.
#> `s_in = NULL`. I will use `s_in = 0` for all observations.
#> Warning: Max number of iterations 20 reached but the algorithm did not converge.
#> Warning: Max number of iterations 20 reached but the algorithm did not converge.
#> Warning: Max number of iterations 20 reached but the algorithm did not converge.
#> `s_in = NULL`. I will use `s_in = 0` for all observations.
#> `s_in = NULL`. I will use `s_in = 0` for all observations.
#>   |                                                                              |=====================                                                 |  30%
#> `s_in = NULL`. I will use `s_in = 0` for all observations.
#> `s_in = NULL`. I will use `s_in = 0` for all observations.
#> `s_in = NULL`. I will use `s_in = 0` for all observations.
#> `s_in = NULL`. I will use `s_in = 0` for all observations.
#>   |                                                                              |============================                                          |  40%
#> `s_in = NULL`. I will use `s_in = 0` for all observations.
#> `s_in = NULL`. I will use `s_in = 0` for all observations.
#> `s_in = NULL`. I will use `s_in = 0` for all observations.
#> `s_in = NULL`. I will use `s_in = 0` for all observations.
#>   |                                                                              |===================================                                   |  50%
#> `s_in = NULL`. I will use `s_in = 0` for all observations.
#> `s_in = NULL`. I will use `s_in = 0` for all observations.
#> `s_in = NULL`. I will use `s_in = 0` for all observations.
#> `s_in = NULL`. I will use `s_in = 0` for all observations.
#>   |                                                                              |==========================================                            |  60%
#> `s_in = NULL`. I will use `s_in = 0` for all observations.
#> `s_in = NULL`. I will use `s_in = 0` for all observations.
#> `s_in = NULL`. I will use `s_in = 0` for all observations.
#> `s_in = NULL`. I will use `s_in = 0` for all observations.
#>   |                                                                              |=================================================                     |  70%
#> `s_in = NULL`. I will use `s_in = 0` for all observations.
#> `s_in = NULL`. I will use `s_in = 0` for all observations.
#> `s_in = NULL`. I will use `s_in = 0` for all observations.
#> `s_in = NULL`. I will use `s_in = 0` for all observations.
#>   |                                                                              |========================================================              |  80%
#> `s_in = NULL`. I will use `s_in = 0` for all observations.
#> `s_in = NULL`. I will use `s_in = 0` for all observations.
#> `s_in = NULL`. I will use `s_in = 0` for all observations.
#> `s_in = NULL`. I will use `s_in = 0` for all observations.
#>   |                                                                              |===============================================================       |  90%
#> `s_in = NULL`. I will use `s_in = 0` for all observations.
#> `s_in = NULL`. I will use `s_in = 0` for all observations.
#> `s_in = NULL`. I will use `s_in = 0` for all observations.
#> `s_in = NULL`. I will use `s_in = 0` for all observations.
#>   |                                                                              |======================================================================| 100%
#> Warning: 2 out of 10 bootstrap replicates failed and were discarded.
#> `s_in = NULL`. I will use `s_in = 0` for all observations.
#> `s_in = NULL`. I will use `s_in = 0` for all observations.
#> Done. 8 successful replicates out of 10 attempted.

# --- Inspect output ------------------------------------------------------
# Names of causes
names(boot_cif)
#> [1] "results_cif"  "results_surv" "grid"        

# Dimensions of the lower confidence bound matrix for cause 1
dim(boot_cif$results_cif[["cause1"]]$lower)
#> [1] 69 51

# Pointwise standard errors for cause 2 (first few rows and columns)
boot_cif$results_cif[["cause2"]]$se[1:5, 1:5]
#>             [,1]        [,2]        [,3]        [,4]        [,5]
#> [1,] 0.002609207 0.005180806 0.007717291 0.010221142 0.012694819
#> [2,] 0.002444048 0.004861595 0.007254613 0.009625085 0.011975006
#> [3,] 0.002287397 0.004557973 0.006813294 0.009054959 0.011284593
#> [4,] 0.002138211 0.004267947 0.006390467 0.008507070 0.010619098
#> [5,] 0.001995467 0.003989552 0.005983281 0.007977729 0.009974023
```
