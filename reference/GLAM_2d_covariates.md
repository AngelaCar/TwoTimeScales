# Fit the 2d GLAM with covariates

`GLAM_2d_covariates()` fits a GLAM for the hazard with two time scales,
with covariates.

## Usage

``` r
GLAM_2d_covariates(
  R,
  Y,
  Bu,
  Bs,
  Z,
  Wprior = NULL,
  P,
  ridge = 0,
  control_algorithm = list(maxiter = 20, conv_crit = 1e-05, verbose = FALSE)
)
```

## Arguments

- R:

  A 3d-array of dimensions nu by ns by n containing exposure times.

- Y:

  A 3d-array of dimensions nu by ns by n containing event indicators.

- Bu:

  A matrix of B-splines for the `u` time scale of dimension nu by cu.

- Bs:

  A matrix of B-splines for the `s` time scale of dimension ns by cs.

- Z:

  (optional) A regression matrix of covariates values of dimensions n by
  p.

- Wprior:

  An optional matrix of a-priori weights.

- P:

  The penalty matrix of dimension cucs by cucs.

- ridge:

  A ridge penalty parameter: default is 0.

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

- `Alpha` The matrix of estimated P-splines coefficients of dimension cu
  by cs.

- `Cov_alpha` The variance-covariance matrix of the `Alpha`
  coefficients, of dimension cucs by cucs.

- `beta` The vector of length p of estimated covariates coefficients.

- `Cov_beta` The variance-covariance matrix of the `beta` coefficients,
  of dimension p by p.

- `SE_beta` The vector of length p of estimated Standard Errors for the
  `beta` coefficients.

- `Cov_Alpha_beta` The matrix with the covariance values between the
  `Alpha` and `beta` coefficients.

- `Eta0` The matrix of values of the baseline linear predictor
  (log-hazard) of dimension nu by ns.

- `H` The hat-matrix.

- `deviance` The deviance.

- `ed` The effective dimension of the model.

- `aic` The value of the AIC.

- `bic` The value of the BIC.

- `Bbases` a list with the B-spline bases `Bu` and `Bs`.
