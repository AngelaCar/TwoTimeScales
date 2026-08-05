# Fit the 1d GLAM with covariates

`GLAM_1d_covariates()` fits a GLAM for the hazard with one time scale,
with covariates.

## Usage

``` r
GLAM_1d_covariates(
  R,
  Y,
  Bs,
  Z = Z,
  ns = NULL,
  n = NULL,
  Wprior = NULL,
  P,
  control_algorithm = list(maxiter = 20, conv_crit = 1e-05, verbose = FALSE)
)
```

## Arguments

- R:

  A 2d-array of dimensions \\n_s\\ by \\n\\ containing exposure times.

- Y:

  A 2d-array of dimensions \\n_s\\ by \\n\\ containing event indicators.

- Bs:

  A matrix of B-splines for the `s` time scale of dimension \\n_s\\ by
  \\c_s\\.

- Z:

  A regression matrix of covariates values of dimensions \\n\\ by \\p\\.

- Wprior:

  An optional vector of length \\n_s\\ of a-priori weights.

- P:

  The penalty matrix of dimension \\c_s\\ by \\c_s\\.

- control_algorithm:

  A list with optional values for the parameters of iterative processes:
  \*`maxiter` The maximum number of iterations for the IWLS algorithm,
  default is 20 . \* `conv_crit` The convergence criteria, expressed as
  difference between estimates at iteration i and i+1, default is
  `1e-5`. \* `verbose` A Boolean. Default is `FALSE`. If `TRUE`, monitor
  the iteration process.

## Value

A list with the following elements:

- `alpha` The vector of estimated P-splines coefficients of length
  \\c_s\\.

- `SE_alpha` The vector of estimated Standard Errors for the `alpha`
  coefficients, of length \\c_s\\.

- `beta` The vector of length \\p\\ of estimated covariates
  coefficients.

- `se_beta` The vector of length \\p\\ of estimated Standard Errors for
  the `beta` coefficients.

- `eta0` The vector of values of the baseline linear predictor
  (log-hazard).

- `H` The hat-matrix.

- `Cov` The full variance-covariance matrix.

- `deviance` The deviance.

- `ed` The effective dimension of the model.

- `aic` The value of the AIC.

- `bic` The value of the BIC.

- `Bbases` a list with the B-spline basis `Bs` (this is a list for
  compatibility with functions in 2d).
