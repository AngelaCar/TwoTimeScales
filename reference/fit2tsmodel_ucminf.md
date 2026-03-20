# Numerical optimization of the 2ts model

`fit2tsmodel_ucminf()` performs a numerical optimization of the AIC or
BIC of the two time scales model.

It finds the optimal values of \\\log\_{10}(\rho_u)\\ and
\\\log\_{10}(\rho_s)\\ and returns the estimated optimal model. See also
[`ucminf::ucminf()`](https://rdrr.io/pkg/ucminf/man/ucminf.html).

## Usage

``` r
fit2tsmodel_ucminf(
  Y,
  R,
  Z = NULL,
  optim_criterion = c("aic", "bic"),
  lrho = c(0, 0),
  Bu,
  Bs,
  Iu,
  Is,
  Du,
  Ds,
  Wprior = NULL,
  ridge = 0,
  control_algorithm = list()
)
```

## Arguments

- Y:

  A matrix (or 3d-array) of event counts of dimension nu by ns (or nu by
  ns by n).

- R:

  A matrix (or 3d-array) of exposure times of dimension nu by ns (or nu
  by ns by n).

- Z:

  (optional) A regression matrix of covariates values of dimensions n by
  p.

- optim_criterion:

  The criterion to be used for optimization: `"aic"` (default) or
  `"bic"`.

- lrho:

  A vector of two elements, the initial values for
  \\\log\_{10}(\varrho_u)\\ and \\\log\_{10}(\varrho_s)\\.

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

An object of class `haz2ts` with the following elements:

- `optimal_model` A list containing the results of the optimal model.

- `optimal_logrho` A vector with the optimal values of `log10(rho_u)`
  and `log10(rho_s)`.

- `P_optimal` The optimal penalty matrix P.

## References

Nielsen H, Mortensen S (2024). *ucminf: General-Purpose Unconstrained
Non-Linear Optimization*. R package version 1.2.2,
<https://CRAN.R-project.org/package=ucminf>
