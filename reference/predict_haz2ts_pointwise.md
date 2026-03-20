# Point-wise prediction hazard 2 time scale

Point-wise prediction hazard 2 time scale

## Usage

``` r
predict_haz2ts_pointwise(fitted_model, upred, spred, zpred = NULL, ds = NULL)
```

## Arguments

- fitted_model:

  An object of class `'haz2ts'` fitted via
  [`fit2ts()`](https://angelacar.github.io/TwoTimeScales/reference/fit2ts.md).

- upred:

  The value(s) of `u` where prediction is required

- spred:

  The value(s) of `s` where prediction is required

- zpred:

  (optional) is a vector of values for the covariates in the model

- ds:

  (optional) The distance between two consecutive points on the `s`
  axis. If not provided, an optimal minimum value will be chosen
  automatically and a warning is returned.

## Value

A data.frame with one row and 6 variable: the values of `u` and `s` for
which predictions of `hazard`, `se_hazard`, the cumulative hazard
`cumhaz` and the `survival` probability are obtained

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
fakedata2ts <- prepare_data(
  u = fakedata$u,
  s_out = fakedata$s,
  ev = fakedata$ev,
  ds = .5
)
#> `s_in = NULL`. I will use `s_in = 0` for all observations.
#> `s_in = NULL`. I will use `s_in = 0` for all observations.
# Fit a fake model - not optimal smoothing
fakemod <- fit2ts(fakedata2ts,
  optim_method = "grid_search",
  lrho = list(
    seq(1, 1.5, .5),
    seq(1, 1.5, .5)
  )
)
predict_haz2ts_pointwise(fakemod, upred = 5, spred = 4.44)
#> chosen interval: ds = 0.4
#>     u   s basehazard se_basehazard    cumhaz  survival
#> 161 5 4.6  0.1632987     0.1021754 0.7466594 0.4739472
```
