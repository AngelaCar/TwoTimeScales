# Iterative Weighted Least Squares algorithm for 1ts model

`iwls_1d()` fits the 1ts model with IWLS algorithm.

## Usage

``` r
iwls_1d(
  r,
  y,
  Bs,
  P,
  Wprior = NULL,
  control_algorithm = list(maxiter = 20, conv_crit = 1e-05, verbose = FALSE)
)
```

## Arguments

- r:

  A vector of exposure times of length ns.

- y:

  A vector of event counts of length ns.

- Bs:

  A matrix of B-splines for the `s` time scale of dimension ns by cs.

- P:

  The penalty matrix of dimension cs by cs.

- Wprior:

  An optional vector of length ns of a-priori weights.

- control_algorithm:

  A list with optional values for the parameters of the iterative
  processes:

  - `maxiter` The maximum number of iteration for the IWSL algorithm.
    Default is 20.

  - `conv_crit` The convergence criteria, expressed as difference
    between estimates at iteration i and i+1. Default is `1e-5`.

  - `verbose` A Boolean. Default is `FALSE`. If `TRUE` monitors the
    iteration process.

## Value

A list with the following elements:

- `alpha` The vector of estimated P-splines coefficients of length cs.

- `SE_alpha` The vector of estimated Standard Errors for the `alpha`
  coefficients, of length cs.

- `H` The hat-matrix.

- `Cov` The full variance-covariance matrix.

- `deviance` The deviance.

- `ed` The effective dimension of the model.

- `aic` The value of the AIC.

- `bic` The value of the BIC.

- `Bbases` a list with the B-spline basis `Bs` (this is a list for
  compatibility with functions in 2d).
