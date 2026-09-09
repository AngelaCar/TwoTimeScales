# Prediction method for objects of class `'haz2ts'`

Prediction method for objects of class `'haz2ts'`

## Usage

``` r
predict_haz2ts(
  x,
  newdata = NULL,
  originaldata = NULL,
  u,
  s,
  z = NULL,
  id = NULL,
  ds = NULL
)
```

## Arguments

- x:

  an object of class `'haz2ts'`, the output of the function
  [`fit2ts()`](https://angelacar.github.io/TwoTimeScales/reference/fit2ts.md).

- newdata:

  (optional) A dataframe with columns cointaing the values of the
  variable `u` and the time scale `s` for which predictions are to be
  obtained.

- originaldata:

  (optional) The original dataset. Provide it to obtain individual
  predictions for each observation in the data.

- u:

  The name of the variable in `newdata`, or in `originaldata` containing
  values for the variable `u`.

- s:

  The name of the variable in `newdata`, or in `originaldata` containing
  values for the variable `s`. Note that over the `s` axis predictions
  are provided only within intervals of values, as it is necessary to
  approximate cumulated quantities.

- z:

  Covariates value

- id:

  (optional) The name of the variable in `newdata`, or in `originaldata`
  containing the identification of each observation. It is not required
  for predictions on a new dataset.

- ds:

  (optional) The distance between two consecutive points on the `s`
  axis. If not provided, an optimal minimum value will be chosen
  automatically and a warning is returned.

## Value

A dataframe. This can be the original dataframe (`originaldata`), where
only the variables `id`, `u` and `s` are selected, or the new data frame
(`newdata`), together with the predicted values for the hazard `hazard`
and its standard errors `se_hazard`, the cumulative hazard `cumhazard`
and the survival probability `survival`.

## Details

Predictions of cumulated quantities can be provided only within
intervals of values on the `s` time scale.

## Examples

``` r
# Create the same fake data as in other examples
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
z1 <- rbinom(n = 20, size = 20, prob = .5)
z2 <- rnorm(n = 20)

fakedata <- as.data.frame(cbind(id, u, s, ev, z1, z2))
fakedata2ts <- prepare_data(
  u = fakedata$u,
  s_out = fakedata$s,
  ev = fakedata$ev,
  ds = .5, individual = TRUE,
  covs = subset(fakedata, select = c("z1", "z2"))
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
# Create a new dataset for prediction
newdata <- as.data.frame(cbind("u" = c(2.5, 3.4, 6),
                               "s" = c(.2, .5, 1.3)))

# First - predict on original data
predict(object = fakemod,
  originaldata = fakedata, u = "u", s = "s", id = "id"
)
#>    id    u    s   hazard se_hazard  cumhazard     survival
#> 1   1 5.43 0.44 2.495792  6.015473  1.2278831 2.929120e-01
#> 2   2 3.25 4.89 6.939030 22.570593 26.1694168 4.312873e-12
#> 3   3 8.15 0.92 1.006956  2.376716  0.9836509 3.739434e-01
#> 4   4 5.53 1.81 2.694720  6.555217  4.7842196 8.360646e-03
#> 5   5 7.28 2.02 1.466431  3.360636  2.9086045 5.455180e-02
#> 6   6 6.61 1.55 1.814980  4.185762  2.7691807 6.271337e-02
#> 7   7 5.91 3.16 2.573711  6.358072  7.4474787 5.829094e-04
#> 8   8 4.94 6.36 4.598209 13.286604 23.1955021 8.439591e-11
#> 9   9 4.25 0.66 3.502920  9.188226  2.3770289 9.282596e-02
#> 10 10 3.86 2.02 4.414879 12.085321  8.3771526 2.300641e-04
#> 11 11 4.05 1.22 3.901952 10.427029  4.7673836 8.502597e-03
#> 12 12 6.86 3.96 1.910181  4.910188  6.8247124 1.086588e-03
#> 13 13 4.94 7.07 4.878881 14.529867 26.5257677 3.019990e-12
#> 14 14 4.46 2.91 4.092417 10.932506 10.8364998 1.966835e-05
#> 15 15 2.14 3.38 7.631547 25.374380 20.6521687 1.073686e-09
#> 16 16 7.56 2.36 1.346345  3.115455  3.0317925 4.822911e-02
#> 17 17 5.55 1.74 2.657399  6.445267  4.4860518 1.126503e-02
#> 18 18 7.60 0.06 1.167410  2.894301  0.1167410 8.898156e-01
#> 19 19 6.46 5.76 2.488046  7.238958 12.1169853 5.465881e-06
#> 20 20 4.96 3.00 3.539978  9.161398  9.7687610 5.721119e-05


# Now - predict on new dataset
predict(object = fakemod,
  newdata = newdata, u = "u", s = "s"
)
#> chosen interval: ds = 0.2
#> Warning: Right boundary adjusted to max(x) = 7.5
#> Warning: Right boundary adjusted to max(x) = 7.5
#> Warning: Right boundary adjusted to max(x) = 7.5
#>     u   s   hazard cumhazard se_hazard   survival
#> 1 2.5 0.2 4.668981  1.843764 14.667175 0.15822080
#> 2 3.4 0.5 4.106394  2.407444 11.646665 0.09004517
#> 3 6.0 1.3 2.210598  2.966389  5.204505 0.05148888

# Now - predict including covariates
predict(object = fakemod,
        originaldata = fakedata,
        u = "u", s = "s", id = "id",
        z = c("z1", "z2")
)
#>    id    u    s z1         z2     hazard  se_hazard  cumhazard  survival
#> 1   1 5.43 0.44 10 -0.3993150 0.09035658 0.06939913 0.04445375 0.9565198
#> 2   2 3.25 4.89 10 -1.4319202 0.19840261 0.30059923 0.74824306 0.4731972
#> 3   3 8.15 0.92  8 -2.0101657 0.04810397 0.07941872 0.04699065 0.9540963
#> 4   4 5.53 1.81  7  0.3383292 0.30406833 0.21214832 0.53984451 0.5828389
#> 5   5 7.28 2.02 11  0.6512870 0.04888114 0.05306512 0.09695369 0.9075980
#> 6   6 6.61 1.55 10  2.4331515 0.12554268 0.12154676 0.19154501 0.8256825
#> 7   7 5.91 3.16  7  1.1913382 0.35293304 0.30699435 1.02127286 0.3601362
#> 8   8 4.94 6.36 12  0.9224428 0.11809294 0.12494138 0.59571568 0.5511680
#> 9   9 4.25 0.66  9  0.6210974 0.22112358 0.13618539 0.15005115 0.8606639
#> 10 10 3.86 2.02 12  0.3563335 0.09962275 0.06813167 0.18903237 0.8277597
#> 11 11 4.05 1.22 10  0.0711464 0.15730200 0.08867782 0.19219073 0.8251495
#> 12 12 6.86 3.96 11  0.7318026 0.06485547 0.08372746 0.23171625 0.7931712
#> 13 13 4.94 7.07 10  0.6197897 0.22296375 0.25287463 1.21222161 0.2975355
#> 14 14 4.46 2.91  9 -0.2425597 0.21205708 0.12212841 0.56151578 0.5703439
#> 15 15 2.14 3.38 12  2.1871192 0.26169019 0.39561389 0.70817486 0.4925423
#> 16 16 7.56 2.36 10 -0.5817275 0.04675201 0.05705196 0.10527943 0.9000730
#> 17 17 5.55 1.74 10  0.7000802 0.12369182 0.05840095 0.20880864 0.8115505
#> 18 18 7.60 0.06  7  1.4921766 0.17148220 0.26250407 0.01714822 0.9829980
#> 19 19 6.46 5.76  9  0.5265534 0.15370158 0.24275372 0.74853901 0.4730572
#> 20 20 4.96 3.00 10  1.0377721 0.17799433 0.08243671 0.49118496 0.6119009
#>    basehazard se_basehazard
#> 1    2.495792      6.015473
#> 2    6.939030     22.570593
#> 3    1.006956      2.376716
#> 4    2.694720      6.555217
#> 5    1.466431      3.360636
#> 6    1.814980      4.185762
#> 7    2.573711      6.358072
#> 8    4.598209     13.286604
#> 9    3.502920      9.188226
#> 10   4.414879     12.085321
#> 11   3.901952     10.427029
#> 12   1.910181      4.910188
#> 13   4.878881     14.529867
#> 14   4.092417     10.932506
#> 15   7.631547     25.374380
#> 16   1.346345      3.115455
#> 17   2.657399      6.445267
#> 18   1.167410      2.894301
#> 19   2.488046      7.238958
#> 20   3.539978      9.161398
# If one wants to predict with only one of the covariates at a different
# value than the baseline, the other one(s) should be fixed at their
# baseline levels

newdata2 <- as.data.frame(cbind("u" = c(2.5, 3.4, 6),
                                "s" = c(.2, .5, 1.3),
                                "z1" = c(1, 2, 3),
                                "z2" = c(0, 0, 0)))
predict(object = fakemod,
        newdata = newdata2,
        u = "u", s = "s", id = "id",
        z = c("z1", "z2")
)
#> chosen interval: ds = 0.2
#> Warning: Right boundary adjusted to max(x) = 7.5
#> Warning: Right boundary adjusted to max(x) = 7.5
#> Warning: Right boundary adjusted to max(x) = 7.5
#>     u   s z1 z2    hazard cumhazard se_hazard  survival basehazard
#> 1 2.5 0.2  1  0 3.3811239  1.335194  9.817081 0.2631072   4.668981
#> 2 3.4 0.5  2  0 2.1534688  1.262508  5.036644 0.2829435   4.106394
#> 3 6.0 1.3  3  0 0.8395117  1.126536  1.336703 0.3241541   2.210598
#>   se_basehazard
#> 1     14.667175
#> 2     11.646665
#> 3      5.204505
```
