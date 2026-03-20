# Grid search for the optimal 1ts model

`grid_search_1d()` performs a grid search for the minimum AIC or BIC of
the one time scale model.

It finds the optimal values of `log_10(rho_s)` and returns the estimated
optimal model.

## Usage

``` r
grid_search_1d(
  r,
  y,
  Z = NULL,
  lrho,
  Bs,
  Ds,
  Wprior = NULL,
  optim_criterion = c("aic", "bic"),
  control_algorithm = list(maxiter = 20, conv_crit = 1e-05, verbose = FALSE, monitor_ev =
    FALSE),
  par_gridsearch = list(plot_aic = FALSE, plot_bic = FALSE, return_aic = TRUE, return_bic
    = TRUE, mark_optimal = TRUE, main_aic = "AIC grid", main_bic = "BIC grid")
)
```

## Arguments

- r:

  A vector of exposure times of length ns, or an array of dimension ns
  by n.

- y:

  A vector of event counts of length ns, or an array of dimension ns by
  n.

- Z:

  (optional) A regression matrix of covariates of dimension n by p.

- lrho:

  A vector of `log_10(rho_s)` values.

- Bs:

  A matrix of B-splines for the time scale `s`.

- Ds:

  The difference matrix of the penalty.

- Wprior:

  An optional vector of a-priori weights.

- optim_criterion:

  The criterion to be used for optimization: `"aic"` (default) or
  `"bic"`. BIC penalizes model complexity more strongly than AIC, so
  that its usage is recommended when a smoother fit is preferable (see
  also Camarda, 2012).

- control_algorithm:

  A list with optional values for the parameters of the iterative
  processes:

  - `maxiter` The maximum number of iteration for the IWSL algorithm.
    Default is 20.

  - `conv_crit` The convergence criteria, expressed as difference
    between estimates at iteration i and i+1. Default is `1e-5`.

  - `verbose` A Boolean. Default is `FALSE`. If `TRUE` monitors the
    iteration process.

  - `monitor_ev` A Boolean. Default is `FALSE`. If `TRUE` monitors the
    evaluation of the model over the `log_10(rho_s)` values.

- par_gridsearch:

  A list of parameters for the grid_search:

  - `plot_aic` A Boolean. Default is `FALSE`. If `TRUE`, plot the AIC
    values over the grid of `log_10(rhos)` values.

  - `plot_bic` A Boolean. Default is `FALSE`. If `TRUE`, plot the BIC
    values over the grid of `log_10(rhos)` values.

  - `return_aic` A Boolean. Default is `TRUE`. Return the AIC values.

  - `return_bic` A Boolean. Default is `TRUE`. Return the BIC values.

  - `mark_optimal` A Boolean. Default is `TRUE`. If the plot of the AIC
    or BIC values is returned, marks the optimal `log_10(rho_s)` in the
    plot.

  - `main_aic` The title of the AIC plot. Default is `"AIC grid"`.

  - `main_bic` The title of the BIC plot. Default is `"BIC grid"`.

## Value

An list with the following elements:

- `optimal_model` A list containing the results of the optimal model.

- `optimal_logrho` The optimal value of `log10(rho_s)`.

- `P_optimal` The optimal penalty matrix P.

- `AIC` (if `par_gridsearch$return_aic == TRUE`) The vector of AIC
  values.

- `BIC` (if `par_gridsearch$return_bic == TRUE`) The vector of BIC
  values.

## References

Camarda, C. G. (2012). "MortalitySmooth: An R Package for Smoothing
Poisson Counts with P-Splines." Journal of Statistical Software, 50(1),
1–24. https://doi.org/10.18637/jss.v050.i01
