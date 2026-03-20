# Return the BIC of 2ts model

`get_bic_fit_2d()` fits the 2ts model with or without individual level
covariates and it returns the BIC of the model. See also
[`fit2tsmodel_ucminf()`](https://angelacar.github.io/TwoTimeScales/reference/fit2tsmodel_ucminf.md)
and
[`fit2ts()`](https://angelacar.github.io/TwoTimeScales/reference/fit2ts.md).

## Usage

``` r
get_bic_fit_2d(
  lrho,
  R,
  Y,
  Z = NULL,
  Bu,
  Bs,
  Iu,
  Is,
  Du,
  Ds,
  Wprior = NULL,
  ridge = 0,
  control_algorithm = list(maxiter = 20, conv_crit = 1e-05, verbose = FALSE, monitor_ev =
    FALSE)
)
```

## Arguments

- lrho:

  A vector of two elements, the initial values for
  \\\log\_{10}(\varrho_u)\\ and \\\log\_{10}(\varrho_s)\\.

- R:

  A matrix (or 3d-array) of exposure times of dimension nu by ns (or nu
  by ns by n).

- Y:

  A matrix (or 3d-array) of event counts of dimension nu by ns (or nu by
  ns by n).

- Z:

  (optional) A regression matrix of covariates values of dimensions n by
  p.

- Bu:

  A matrix of B-splines for the `u` time scale of dimension nu by cu.

- Bs:

  A matrix of B-splines for the `s` time scale of dimension ns by cs.

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

The `bic` value of the fitted model.
