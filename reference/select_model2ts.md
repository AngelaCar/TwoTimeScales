# Selects best model among a group of fitted models

`select_model2ts()` takes as input a list of hazard models with two time
scales, and returns a table with the best fitting model indicated. The
selection is based on minimization of the AIC or BIC, so the models do
not need to be nested. It can be used, for example, to identify the best
model among the log-additive model, the varying coefficient model, and
the full 2D model.

## Usage

``` r
select_model2ts(model_list, sel_criteria = c("aic", "bic"))
```

## Arguments

- model_list:

  is a list of fitted objects, of class `haz2ts`, `hat2tsLMM`,
  `haz2tsPGAM`, or `haz2tsVCM`.

- sel_criteria:

  to determine the best model - it can be `'aic'` or `'bic'`.

## Value

A data.frame with as many rows as the number of models being compared,
and the following columns: \* Model: the model name (see details and
examples) \* Type: the class of the model \* AIC \* BIC \* Best:
indicates which of the model is the best fitting one with respect to the
criteria indicated

## Details

In the model list, it is possible to provide a character name for each
model. In such case, these names will also be reported in the results'
table. An example is provided (See examples).

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
  data = fakedata,
  u = "u",
  s_out = "s",
  ev = "ev",
  ds = .5
)
#> `s_in = NULL`. I will use `s_in = 0` for all observations.
#> `s_in = NULL`. I will use `s_in = 0` for all observations.
# Fit a fake model - not optimal smoothing
full2d <- fit2ts(fakedata2ts,
  optim_method = "grid_search",
  lrho = list(seq(1, 1.5, .5), seq(1, 1.5, .5))
)
pgam <- fitpgam(fakedata2ts,
  optim_method = "grid_search",
  optim_criterion = "aic",
  lrho = list(seq(1, 1.5, .5), seq(1, 1.5, .5))
)
vcm <- fitvcm(fakedata2ts)

select_model2ts(model_list = list(
  "full interaction" = full2d,
  "additive" = pgam,
  "varying coeff" = vcm
))
#>              Model       Type      AIC      BIC Best
#> 1 full interaction     haz2ts 42.86643 52.93603 <NA>
#> 2         additive haz2tsPGAM 39.98714 49.07529    *
#> 3    varying coeff  haz2tsVCM 47.23913 66.04327 <NA>

```
