# Numerical optimization of the 1ts model

`fit1tsmodel_ucminf()` performs a numerical optimization of the AIC or
BIC of the one time scale model.

It finds the optimal values of \\\log\_{10}(\varrho_s)\\ and returns the
estimated optimal model. See also
[`ucminf::ucminf()`](https://rdrr.io/pkg/ucminf/man/ucminf.html).

## Usage

``` r
fit1tsmodel_ucminf(
  r,
  y,
  Z = NULL,
  lrho = 0,
  Bs,
  Ds,
  Wprior = NULL,
  optim_criterion = c("aic", "bic"),
  control_algorithm = list()
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

  A starting value for \\\log\_{10}(\varrho_s)\\. Default is 0.

- Bs:

  A matrix of B-splines for the time scale `s`.

- Ds:

  The difference matrix of the penalty.

- Wprior:

  An optional vector of a-priori weights.

- optim_criterion:

  The criterion to be used for optimization: `"aic"` (default) or
  `"bic"`.

- control_algorithm:

  A list with optional values for the parameters of the iterative
  processes:

  - `maxiter` The maximum number of iteration for the IWLS algorithm.
    Default is 20.

  - `conv_crit` The convergence criteria, expressed as difference
    between estimates at iteration i and i+1. Default is `1e-5`.

  - `verbose` A Boolean. Default is `FALSE`. If `TRUE` monitors the
    iteration process.

  - `monitor_ev` A Boolean. Default is `FALSE`. If `TRUE` monitors the
    evaluation of the model over the `log_10(rho_s)` values.

  - `xtol` The relative tolerance to stop the algorithm. For details see
    `help(ucminf)`. Here is it set to `1e-5`.

## Value

An object of class `haz1ts` with the following elements:

- `optimal_model` A list containing the results of the optimal model.

- `optimal_logrho` The optimal value of `log10(rho_s)`.

- `P_optimal` The optimal penalty matrix P.

## References

Nielsen H, Mortensen S (2024). *ucminf: General-Purpose Unconstrained
Non-Linear Optimization*. R package version 1.2.2,
<https://CRAN.R-project.org/package=ucminf>
