# Point-wise prediction of cumulative incidence over 2 time scale

Point-wise prediction of cumulative incidence over 2 time scale

## Usage

``` r
predict_cif2ts_pointwise(fitted_models = list(), u, s, ds = NULL)
```

## Arguments

- fitted_models:

  a list with cause-specific hazard models

- u:

  The value(s) of `u` where prediction is required

- s:

  The value(s) of `s` where prediction is required

- ds:

  (optional) The distance between two consecutive points on the `s`
  axis. If not provided, an optimal minimum value will be chosen
  automatically and a warning is returned.

## Value

A data.frame with one row containing: the values of `u` and `s`for which
predictions of the overall survival (`surv`) probability, and the values
of the cumulative incidence functions, one for each cause, are obtained.

## Examples

``` r
id <- 1:20
u <- c(
  5.43, 3.25, 8.15, 5.53, 7.28, 6.61, 5.91, 4.94, 4.25, 3.86, 4.05, 6.86,
  4.94, 4.46, 2.14, 7.56, 5.55, 7.60, 6.46, 4.96
)
s <- c(
  0.44, 4.89, 0.92, 1.81, 2.02, 1.55, 3.16, 6.36, 0.66, 2.02, 1.22, 3.96,
  7.07, 2.91, 3.38, 2.36, 1.74, 0.06, 5.76, 3.00
)
ev <- c(1, 0, 0, 1, 0, 1, 0, 1, 1, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 1)

fakedata <- as.data.frame(cbind(id, u, s, ev))
# cause 1
fakedata2ts1 <- prepare_data(
  u = fakedata$u,
  s_out = fakedata$s,
  ev = fakedata$ev,
  min_u = 2, min_s = 0,
  ds = .5
)
#> `s_in = NULL`. I will use `s_in = 0` for all observations.
#> `s_in = NULL`. I will use `s_in = 0` for all observations.
# Fit a fake model - not optimal smoothing for cause type 1
fakemod1 <- fit2ts(fakedata2ts1,
  optim_method = "grid_search",
  lrho = list(
    seq(1, 1.5, .5),
    seq(1, 1.5, .5)
  )
)
# cause 2
fakedata2ts2 <- prepare_data(
  u = fakedata$u,
  s_out = fakedata$s,
  ev = 1 - (fakedata$ev),
  min_u = 2, min_s = 0,
  ds = .5
)
#> `s_in = NULL`. I will use `s_in = 0` for all observations.
#> `s_in = NULL`. I will use `s_in = 0` for all observations.
# Fit a fake model - not optimal smoothing for cause 2
fakemod2 <- fit2ts(fakedata2ts2,
  optim_method = "grid_search",
  lrho = list(
    seq(1, 1.5, .5),
    seq(1, 1.5, .5)
  )
)
predict_cif2ts_pointwise(
  fitted_models = list(fakemod1, fakemod2),
  u = 5.2, s = 4.4
)
#> chosen interval: ds = 0.4
#>       u   s      surv     cif_1    cif_2
#> 161 5.2 4.6 0.2133407 0.3975849 0.341374
```
