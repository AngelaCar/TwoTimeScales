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
#>    id    u    s     hazard se_hazard   cumhazard  survival
#> 1   1 5.43 0.44 0.15064076 0.3503795 0.075690953 0.9271027
#> 2   2 3.25 4.89 0.16024763 0.4399220 0.944039214 0.3890532
#> 3   3 8.15 0.92 0.07158399 0.1788405 0.070027233 0.9323684
#> 4   4 5.53 1.81 0.14173158 0.3182666 0.275991442 0.7588194
#> 5   5 7.28 2.02 0.09375046 0.2156988 0.192703827 0.8247262
#> 6   6 6.61 1.55 0.11058440 0.2512125 0.176579627 0.8381320
#> 7   7 5.91 3.16 0.12581855 0.2815624 0.418567856 0.6579885
#> 8   8 4.94 6.36 0.12657482 0.2922401 0.948071371 0.3874876
#> 9   9 4.25 0.66 0.18954905 0.4464218 0.134577717 0.8740849
#> 10 10 3.86 2.02 0.18498223 0.4285434 0.413319942 0.6614506
#> 11 11 4.05 1.22 0.18953429 0.4418416 0.254480627 0.7753191
#> 12 12 6.86 3.96 0.10467557 0.2599607 0.414870077 0.6604261
#> 13 13 4.94 7.07 0.12292868 0.2906392 1.035213696 0.3551505
#> 14 14 4.46 2.91 0.16133749 0.3645245 0.525305710 0.5913745
#> 15 15 2.14 3.38 0.19781634 0.5505536 0.775177016 0.4606222
#> 16 16 7.56 2.36 0.08827455 0.2072859 0.204614059 0.8149618
#> 17 17 5.55 1.74 0.14158628 0.3184558 0.260611424 0.7705803
#> 18 18 7.60 0.06 0.08131508 0.2155693 0.008131508 0.9919015
#> 19 19 6.46 5.76 0.11057708 0.3023623 0.656795417 0.5185103
#> 20 20 4.96 3.00 0.14896262 0.3318964 0.495850781 0.6090525


# Now - predict on new dataset
predict(object = fakemod,
  newdata = newdata, u = "u", s = "s"
)
#> chosen interval: ds = 0.2
#> Warning: Right boundary adjusted to max(x) = 7.5
#> Warning: Right boundary adjusted to max(x) = 7.5
#> Warning: Right boundary adjusted to max(x) = 7.5
#>     u   s    hazard  cumhazard se_hazard  survival
#> 1 2.5 0.2 0.2456699 0.09908983 0.6601345 0.9056613
#> 2 3.4 0.5 0.2166524 0.13168108 0.5342636 0.8766205
#> 3 6.0 1.3 0.1291541 0.18221844 0.2933662 0.8334193

# Now - predict including covariates
predict(object = fakemod,
        originaldata = fakedata,
        u = "u", s = "s", id = "id",
        z = c("z1", "z2")
)
#>    id    u    s z1          z2     hazard  se_hazard   cumhazard  survival
#> 1   1 5.43 0.44  8 -0.32083991 0.15401676 0.11467254 0.077387261 0.9255314
#> 2   2 3.25 4.89  9  1.47100572 0.09916912 0.17218110 0.584217957 0.5575417
#> 3   3 8.15 0.92  9  1.70432940 0.04154101 0.05747600 0.040637611 0.9601770
#> 4   4 5.53 1.81 12  0.04324404 0.12680413 0.08178189 0.246923482 0.7812005
#> 5   5 7.28 2.02  8 -0.33265732 0.09616415 0.08985181 0.197665157 0.8206446
#> 6   6 6.61 1.55 12 -1.82223542 0.16543108 0.13017410 0.264158052 0.7678522
#> 7   7 5.91 3.16  9  1.41126240 0.07915520 0.06773494 0.263330212 0.7684881
#> 8   8 4.94 6.36  8 -0.83758243 0.14921601 0.15184529 1.117658544 0.3270447
#> 9   9 4.25 0.66  9 -1.12376279 0.23979641 0.17359287 0.170252780 0.8434516
#> 10 10 3.86 2.02  9  3.04376589 0.07421478 0.07234452 0.165823756 0.8471955
#> 11 11 4.05 1.22 10  0.23502131 0.16352948 0.08399861 0.219564943 0.8028680
#> 12 12 6.86 3.96  8 -0.03325861 0.09886754 0.14028658 0.391850599 0.6758051
#> 13 13 4.94 7.07 11 -2.73221952 0.23827465 0.26815605 2.006571559 0.1344488
#> 14 14 4.46 2.91 13 -0.09979059 0.14890993 0.12369124 0.484842274 0.6157943
#> 15 15 2.14 3.38 11  0.97603173 0.13800393 0.22109161 0.540791894 0.5822870
#> 16 16 7.56 2.36 13  0.41386892 0.07072130 0.09456018 0.163926886 0.8488041
#> 17 17 5.55 1.74 14  0.91232216 0.09805859 0.10887597 0.180491972 0.8348594
#> 18 18 7.60 0.06  8  1.98373220 0.04405452 0.07067912 0.004405452 0.9956042
#> 19 19 6.46 5.76 10  1.16910851 0.07375371 0.14469556 0.438075418 0.6452771
#> 20 20 4.96 3.00  9 -0.50873702 0.15907189 0.07567617 0.529501415 0.5888985
#>    basehazard se_basehazard
#> 1  0.15064076     0.3503795
#> 2  0.16024763     0.4399220
#> 3  0.07158399     0.1788405
#> 4  0.14173158     0.3182666
#> 5  0.09375046     0.2156988
#> 6  0.11058440     0.2512125
#> 7  0.12581855     0.2815624
#> 8  0.12657482     0.2922401
#> 9  0.18954905     0.4464218
#> 10 0.18498223     0.4285434
#> 11 0.18953429     0.4418416
#> 12 0.10467557     0.2599607
#> 13 0.12292868     0.2906392
#> 14 0.16133749     0.3645245
#> 15 0.19781634     0.5505536
#> 16 0.08827455     0.2072859
#> 17 0.14158628     0.3184558
#> 18 0.08131508     0.2155693
#> 19 0.11057708     0.3023623
#> 20 0.14896262     0.3318964
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
#>     u   s z1 z2    hazard  cumhazard se_hazard  survival basehazard
#> 1 2.5 0.2  1  0 0.2436439 0.09827264 0.6077892 0.9064018  0.2456699
#> 2 3.4 0.5  2  0 0.2130936 0.12951808 0.4376349 0.8785187  0.2166524
#> 3 6.0 1.3  3  0 0.1259850 0.17774724 0.2055417 0.8371540  0.1291541
#>   se_basehazard
#> 1     0.6601345
#> 2     0.5342636
#> 3     0.2933662
```
