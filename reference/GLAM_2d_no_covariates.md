# Fit the 2d GLAM without covariates

`GLAM_2d_no_covariates()` fits a GLAM for the hazard with two time
scales, without covariates.

## Usage

``` r
GLAM_2d_no_covariates(
  R,
  Y,
  Bu,
  Bs,
  Wprior = NULL,
  P,
  ridge = 0,
  control_algorithm = list(maxiter = 20, conv_crit = 1e-05, verbose = FALSE)
)
```

## Arguments

- R:

  A matrix of exposure times of dimension nu by ns.

- Y:

  A matrix of event counts of dimension nu by ns.

- Bu:

  A matrix of B-splines for the `u` time scale of dimension nu by cu.

- Bs:

  A matrix of B-splines for the `s` time scale of dimension ns by cs.

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

- `Eta0` The matrix of values of the baseline linear predictor
  (log-hazard) of dimension nu by ns.

- `H` The hat-matrix.

- `deviance` The deviance.

- `ed` The effective dimension of the model.

- `aic` The value of the AIC.

- `bic` The value of the BIC.

- `Bbases` a list with the B-spline bases `Bu` and `Bs`.
