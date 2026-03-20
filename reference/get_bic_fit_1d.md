# Return the BIC of 1ts model

`get_bic_fit_1d()` fits the 1ts model with or without individual level
covariates and it returns the BIC of the model. See also
[`fit1tsmodel_ucminf()`](https://angelacar.github.io/TwoTimeScales/reference/fit1tsmodel_ucminf.md)
and
[`fit1ts()`](https://angelacar.github.io/TwoTimeScales/reference/fit1ts.md).

## Usage

``` r
get_bic_fit_1d(
  lrho,
  r,
  y,
  Z = NULL,
  Bs,
  Ds,
  Wprior = NULL,
  control_algorithm = list(maxiter = 20, conv_crit = 1e-05, verbose = FALSE, monitor_ev =
    FALSE)
)
```

## Arguments

- lrho:

  A starting value for \\\log\_{10}(\varrho_s)\\. Default is 0.

- r:

  A vector of exposure times of length ns, or an array of dimension ns
  by n.

- y:

  A vector of event counts of length ns, or an array of dimension ns by
  n.

- Z:

  (optional) A regression matrix of covariates of dimension n by p.

- Bs:

  A matrix of B-splines for the time scale `s`.

- Ds:

  The difference matrix of the penalty.

- Wprior:

  An optional vector of a-priori weights.

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

## Value

the `bic` value of the fitted model.
