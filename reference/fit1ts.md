# Fit a smooth hazard model with one time scale

`fit1ts()` fits a smooth hazard model with one time scale.

Three methods are implemented for the search of the optimal smoothing
parameter (and therefore optimal model): a numerical optimization of the
AIC or BIC of the model, a search for the minimum AIC or BIC of the
model over a grid of \\\log\_{10}\\ values for the smoothing parameter
and the estimation using a sparse mixed model representation of
P-splines. Construction of the B-splines basis and of the penalty matrix
is incorporated within the function. If a matrix of covariates is
provided, the function will estimate a model with covariates.

## Usage

``` r
fit1ts(
  data1ts = NULL,
  y = NULL,
  r = NULL,
  Z = NULL,
  bins = NULL,
  Bbases_spec = list(),
  Wprior = NULL,
  pord = 2,
  optim_method = c("ucminf", "grid_search", "LMMsolver"),
  optim_criterion = c("aic", "bic"),
  lrho = 0,
  ridge = 0,
  control_algorithm = list(),
  par_gridsearch = list()
)
```

## Arguments

- data1ts:

  (optional) an object created by the function
  [`prepare_data()`](https://angelacar.github.io/TwoTimeScales/reference/prepare_data.md).
  Providing this input is the easiest way to use the function
  `fit1ts()`. However, the user can also provide the input data together
  with a list of bins, as explained by the following parameters'
  descriptions.

- y:

  A vector of event counts of length ns, or an array of dimension ns by
  n.

- r:

  A vector of exposure times of length ns, or an array of dimension ns
  by n.

- Z:

  (optional) A regression matrix of covariates of dimension n by p.

- bins:

  a list with the specification for the bins. This is created by the
  function
  [`prepare_data()`](https://angelacar.github.io/TwoTimeScales/reference/prepare_data.md).
  Alternatively, a list with the following elements can be provided: \*
  `bins_s` is a vector of intervals for the time scale `s`. \* `mids` is
  a vector with the midpoints of the intervals over `s`. \* `ns` is the
  number of bins over `s`.

- Bbases_spec:

  A list with the specification for the B-splines basis with the
  following elements:

  - `bdeg` The degree of the B-splines basis. Default is 3 (for cubic
    B-splines).

  - `nseg_s` The number of segments for the B-splines over `s`. Default
    is 10.

  - `min_s` (optional) The lower limit of the domain of `Bs`. Default is
    `min(bins_s)`.

  - `max_s` (optional) The upper limit of the domain of `Bs`. Default is
    `max(bins_s)`.

- Wprior:

  An optional vector of a-priori weights.

- pord:

  The order of the penalty. Default is 2.

- optim_method:

  The method to be used for optimization: `"ucminf"` (default) for the
  numerical optimization of the AIC (or BIC), `"grid_search"` for a grid
  search of the minimum AIC (or BIC) over a grid of
  \\\log\_{10}(\varrho_s)\\ values, and `"LMMsolver"` to solve the model
  as sparse linear mixed model using the package LMMsolver.

- optim_criterion:

  The criterion to be used for optimization: `"aic"` (default) or
  `"bic"`. BIC penalizes model complexity more strongly than AIC, so
  that its usage is recommended when a smoother fit is preferable (see
  also Camarda, 2012).

- lrho:

  A number if `optim_method == "ucminf"`, default is 0. A vector of
  values for \\\log\_{10}(\varrho_s)\\ if
  `optim_method == "grid_search"`. In the latter case, if a vector is
  not provided, a default sequence of values is used for
  \\\log\_{10}(\varrho_s)\\ .

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

  - `monitor_ev` A Boolean. Default is `FALSE`. If `TRUE` monitors the
    evaluation of the model over the `log_10(rho_s)` values.

  - `xtol` The relative tolerance to stop the algorithm, only relevant
    for numerical optimization algorithm. For details see
    `help(ucminf)`. Here is it set to `1e-5`.

- par_gridsearch:

  A list of parameters for the grid_search:

  - `plot_aic` A Boolean. Default is `FALSE`. If `TRUE`, plot the AIC
    values over the grid of `log_10(rhos)` values.

  - `plot_bic` A Boolean. Default is `FALSE`. If `TRUE`, plot the BIC
    values over the grid of `log_10(rhos)` values.

  - `return_aic` A Boolean. Default is `TRUE`. Return the AIC values.

  - `return_bic` A Boolean. Default is `TRUE`. Return the BIC values.

  - `mark_optimal` A Boolean. Default is `TRUE`. If the plot of the AIC
    or BIC values is returned, marks the optimal `log_10(rho_s)` in the
    plot.

  - `main_aic` The title of the AIC plot. Default is `"AIC grid"`.

  - `main_bic` The title of the BIC plot. Default is `"BIC grid"`.

## Value

An object of class `haz1ts`, or of class `haz1tsLMM`. For objects of
class `haz1ts` this is

- `optimal_model` A list with:

  - `alpha` The vector of estimated P-splines coefficients of length
    \\c_s\\.

  - `SE_alpha` The vector of estimated Standard Errors for the `alpha`
    coefficients, of length \\c_s\\.

  - `beta` The vector of estimated covariate coefficients of length
    \\p\\ (if model with covariates).

  - `se_beta` The vector of estimated Standard Errors for the `beta`
    coefficients of length \\p\\ (if model with covariates).

  - `eta` or `eta0`. The vector of values of the (baseline) linear
    predictor (log-hazard) of length \\n_s\\.

  - `H` The hat-matrix.

  - `Cov` The full variance-covariance matrix.

  - `deviance` The deviance.

  - `ed` The effective dimension of the model.

  - `aic` The value of the AIC.

  - `bic` The value of the BIC.

  - `Bbases` a list with the B-spline basis `Bs` (this is a list for
    compatibility with functions in 2d).

- `optimal_logrho` The optimal value of `log10(rho_s)`.

- `P_optimal` The optimal penalty matrix P.

- `AIC` (if `par_gridsearch$return_aic == TRUE`) The vector of AIC
  values.

- `BIC` (if `par_gridsearch$return_bic == TRUE`) The vector of BIC
  values.

Objects of class `haz1tsLMM` have a slight different structure. They are
a list with:

- `optimal_model` an object of class `LMMsolve`

- `AIC_BIC` a list with, among other things, the AIC and BIC values and
  the ED of the model

- `n_events` the number of events

- `ns` the number of bins over the s-axis

- `cs` the number of B-splines over the s-axis

- `covariates` an indicator for PH model

## Details

Some functions from the R-package `LMMsolver` are used here. We refer
the interested readers to https://biometris.github.io/LMMsolver/ for
more detail on `LMMsolver` and its usage.

## References

Boer, Martin P. 2023. “Tensor Product P-Splines Using a Sparse Mixed
Model Formulation.” Statistical Modelling 23 (5-6): 465–79.
https://doi.org/10.1177/1471082X231178591.#'

## Examples

``` r
## preparing data - no covariates
dt1ts <- prepare_data(
  data = reccolon2ts,
  s_in = "entrys",
  s_out = "timesr",
  events = "status",
  ds = 180
)

## fitting the model with fit1ts() - default options, that is ucminf optimization

mod1 <- fit1ts(dt1ts)

## fitting with LMMsolver
mod2 <- fit1ts(dt1ts,
  optim_method = "LMMsolver"
)

## preparing the data - covariates

dt1ts_cov <- prepare_data(
  data = reccolon2ts,
  s_in = "entrys",
  s_out = "timesr",
  events = "status",
  ds = 180,
  individual = TRUE,
  covs = c("rx", "node4", "sex")
)

## fitting the model with fit1ts() - grid search over only two log_10(rho_s) values

mod3 <- fit1ts(dt1ts_cov,
  optim_method = "grid_search",
  lrho = c(1, 1.5)
)
```
