# Return the AIC of P-GAM model

`get_aic_fit_pgam()` fits the log-additive model with two time scales
with or without individual level covariates and it returns the AIC of
the model. See also
[`fitpgammodel_ucminf()`](https://angelacar.github.io/TwoTimeScales/reference/fitpgammodel_ucminf.md)
and
[`fitpgam()`](https://angelacar.github.io/TwoTimeScales/reference/fitpgam.md).

## Usage

``` r
get_aic_fit_pgam(
  lrho,
  r,
  y,
  Z = NULL,
  B,
  nb,
  nbu,
  nbs,
  DutDu,
  DstDs,
  ridge = 1e-06,
  control_algorithm = list(maxiter = 20, conv_crit = 1e-05, verbose = FALSE, monitor_ev =
    FALSE)
)
```

## Arguments

- lrho:

  A vector of two elements, the initial values for
  \\\log\_{10}(\varrho_u)\\ and \\\log\_{10}(\varrho_s)\\.

- r:

  A vector of exposure times of length `nu` by `ns`.

- y:

  A vector of event counts of length `nu` by `ns`.

- Z:

  (optional) A regression matrix of covariates values of dimensions n by
  p.

- B:

  A matrix of B-splines.

- nb:

  The number of B-splines.

- nbu:

  The number of B-splines for the `u` dimension.

- nbs:

  The number of B-splines for the `s` dimension.

- DutDu:

  \\D_u^TD_u\\ to multiply \\\varrho_u\\.

- DstDs:

  \\D_s^TD_s\\ to multiply \\\varrho_s\\.

- ridge:

  A ridge penalty parameter: default is 0. This is useful when, in some
  cases the algorithm shows convergence problems. In this case, set to a
  small number, for example `1e-6`.

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

The `aic` value of the fitted model.
