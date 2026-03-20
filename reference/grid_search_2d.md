# Grid search for the optimal 2ts model

`grid_search_2d()` performs a grid search for the minimum AIC or BIC of
the two time scales model.

It finds the optimal values of `log_10(rho_u)` and `log_10(rho_s)` and
returns the estimated optimal model.

## Usage

``` r
grid_search_2d(
  lru,
  lrs,
  R,
  Y,
  Bu,
  Bs,
  Z = NULL,
  Iu,
  Is,
  Du,
  Ds,
  Wprior = NULL,
  ridge = 0,
  optim_criterion = c("aic", "bic"),
  control_algorithm = list(maxiter = 20, conv_crit = 1e-05, verbose = FALSE, monitor_ev =
    FALSE),
  par_gridsearch = list(plot_aic = FALSE, plot_bic = FALSE, return_aic = TRUE, return_bic
    = TRUE, col = grey.colors(n = 10), plot_contour = FALSE, mark_optimal = TRUE,
    main_aic = "AIC grid", main_bic = "BIC grid")
)
```

## Arguments

- lru:

  A vector of `log_10(rho_u)` values.

- lrs:

  A vector of `log_10(rho_s)` values.

- R:

  A matrix (or 3d-array) of exposure times of dimension nu by ns (or nu
  by ns by n).

- Y:

  A matrix (or 3d-array) of event counts of dimension nu by ns (or nu by
  ns by n).

- Bu:

  A matrix of B-splines for the `u` time scale of dimension nu by cu.

- Bs:

  A matrix of B-splines for the `s` time scale of dimension ns by cs.

- Z:

  (optional) A regression matrix of covariates values of dimensions n by
  p.

- Iu:

  An identity matrix of dimension nbu by nbu.

- Is:

  An identity matrix of dimension nbs by nbs.

- Du:

  The difference matrix over `u`.

- Ds:

  The difference matrix over `s`.

- Wprior:

  An optional matrix of a-priori weights.

- ridge:

  A ridge penalty parameter: default is 0. This is useful when, in some
  cases the algorithm shows convergence problems. In this case, set to a
  small number, for example `1e-4`.

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
    values over the grid of `log_10(rho_u)` and `log_10(rho_s)` values.

  - `plot_bic` A Boolean. Default is `FALSE`. If `TRUE`, plot the BIC
    values over the grid of `log_10(rho_u)` and `log_10(rho_s)` values.

  - `return_aic` A Boolean. Default is `TRUE`. Return the AIC values.

  - `return_bic` A Boolean. Default is `TRUE`. Return the BIC values.

  - `col` The color palette to be used for the AIC/BIC plot. Default is
    `grDevices::gray.colors(n=10)`.

  - `plot_contour` A Boolean. Default is `TRUE`. Adds white contour
    lines to the AIC/BIC plot.

  - `mark_optimal` A Boolean. Default is `TRUE`. If the plot of the AIC
    or BIC values is returned, marks the optimal combination of
    `log_10(rho_u)` and `log_10(rho_s)` in the plot.

  - `main_aic` The title of the AIC plot. Default is `"AIC grid"`.

  - `main_bic` The title of the BIC plot. Default is `"BIC grid"`.

## Value

An object of class `"haz2ts"`, that is a list with the following
elements:

- `optimal_model` A list containing the results of the optimal model.

- `optimal_logrho` The optimal couple of `log_10(rho_u)` and
  `log_10(rho_s)` values.

- `P_optimal` The optimal penalty matrix P.

- `AIC` (if `par_gridsearch$return_aic == TRUE`) The vector of AIC
  values.

- `BIC` (if `par_gridsearch$return_bic == TRUE`) The vector of BIC
  values.

## References

Camarda, C. G. (2012). "MortalitySmooth: An R Package for Smoothing
Poisson Counts with P-Splines." Journal of Statistical Software, 50(1),
1–24. https://doi.org/10.18637/jss.v050.i01
