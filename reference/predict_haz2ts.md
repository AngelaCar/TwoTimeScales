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
#>    id    u    s    hazard se_hazard  cumhazard    survival
#> 1   1 5.43 0.44 0.6072060 0.9978781 0.30316460 0.738477536
#> 2   2 3.25 4.89 1.0887850 2.4576861 4.20725546 0.014887171
#> 3   3 8.15 0.92 0.3411957 0.6802598 0.35591421 0.700532715
#> 4   4 5.53 1.81 0.5968252 0.9280079 1.13709740 0.320748678
#> 5   5 7.28 2.02 0.3874782 0.6418565 0.87010487 0.418907616
#> 6   6 6.61 1.55 0.4722751 0.7679877 0.78067430 0.458097012
#> 7   7 5.91 3.16 0.5298979 0.8166548 1.76471520 0.171235547
#> 8   8 4.94 6.36 0.7034043 1.2114875 4.30176843 0.013544585
#> 9   9 4.25 0.66 0.6853213 1.1000959 0.47074009 0.624539879
#> 10 10 3.86 2.02 0.7745896 1.2464374 1.51119197 0.220646817
#> 11 11 4.05 1.22 0.7217774 1.1485244 0.90004271 0.406552297
#> 12 12 6.86 3.96 0.3883457 0.6927739 1.73515579 0.176372721
#> 13 13 4.94 7.07 0.7146268 1.2903305 4.79865338 0.008240837
#> 14 14 4.46 2.91 0.7426744 1.1786959 2.10624500 0.121694070
#> 15 15 2.14 3.38 1.1009874 2.4340152 2.86832878 0.056793762
#> 16 16 7.56 2.36 0.3520825 0.5975389 0.92185570 0.397780195
#> 17 17 5.55 1.74 0.5949336 0.9279192 1.07389717 0.341674355
#> 18 18 7.60 0.06 0.4149519 0.9039435 0.04149519 0.959353946
#> 19 19 6.46 5.76 0.4072842 0.8343537 2.67575332 0.068854939
#> 20 20 4.96 3.00 0.6713441 1.0377240 2.03705732 0.130411907


# Now - predict on new dataset
predict(object = fakemod,
  newdata = newdata, u = "u", s = "s"
)
#> chosen interval: ds = 0.2
#> Warning: Right boundary adjusted to max(x) = 7.5
#> Warning: Right boundary adjusted to max(x) = 7.5
#> Warning: Right boundary adjusted to max(x) = 7.5
#>     u   s    hazard cumhazard se_hazard  survival
#> 1 2.5 0.2 0.6671570 0.2630913 1.2605187 0.7686717
#> 2 3.4 0.5 0.6930127 0.4074481 1.1658347 0.6653460
#> 3 6.0 1.3 0.5464481 0.7740255 0.8777995 0.4611530

# Now - predict including covariates
predict(object = fakemod,
        originaldata = fakedata,
        u = "u", s = "s", id = "id",
        z = c("z1", "z2")
)
#>    id    u    s z1            z2     hazard  se_hazard   cumhazard  survival
#> 1   1 5.43 0.44  8  1.4051088798 0.16457602 0.12947586 0.082169182 0.9211161
#> 2   2 3.25 4.89 14 -0.7954609493 0.11085250 0.19526014 0.428353404 0.6515811
#> 3   3 8.15 0.92 10 -1.5665144652 0.06672547 0.09789766 0.069603877 0.9327632
#> 4   4 5.53 1.81 10 -1.0405791109 0.11671731 0.10325603 0.222374900 0.8006152
#> 5   5 7.28 2.02  9  1.0199337428 0.08920858 0.07773283 0.200323034 0.8184663
#> 6   6 6.61 1.55 10 -0.7020819780 0.09235984 0.07599797 0.152671504 0.8584117
#> 7   7 5.91 3.16  9  0.9733157770 0.12199765 0.08172573 0.406287927 0.6661183
#> 8   8 4.94 6.36  6 -0.0768176526 0.26422743 0.27763623 1.615920118 0.1987078
#> 9   9 4.25 0.66 10  0.8929249245 0.13402392 0.08264367 0.092059644 0.9120507
#> 10 10 3.86 2.02  9 -0.7775030885 0.17833269 0.14554851 0.347919626 0.7061556
#> 11 11 4.05 1.22  5  0.4367971056 0.31918856 0.26412356 0.398022081 0.6716472
#> 12 12 6.86 3.96  9  0.4134439348 0.08940829 0.11258568 0.399482504 0.6706670
#> 13 13 4.94 7.07 10  0.9763417720 0.13975501 0.16422535 0.938442095 0.3912369
#> 14 14 4.46 2.91 12  1.1465004990 0.10479587 0.07863102 0.297203971 0.7428925
#> 15 15 2.14 3.38  9  1.2172716875 0.25347880 0.36920406 0.660371346 0.5166594
#> 16 16 7.56 2.36 10  0.0004800131 0.06885454 0.06812545 0.180281451 0.8350352
#> 17 17 5.55 1.74 11  0.7551250562 0.09882923 0.05654155 0.178393739 0.8366130
#> 18 18 7.60 0.06  9  0.3424035105 0.09553381 0.13084702 0.009553381 0.9904921
#> 19 19 6.46 5.76  9  0.1684728243 0.09376849 0.16076826 0.616034966 0.5400816
#> 20 20 4.96 3.00 10  1.3970665088 0.13129049 0.08795798 0.398374334 0.6714106
#>    basehazard se_basehazard
#> 1   0.6072060     0.9978781
#> 2   1.0887850     2.4576861
#> 3   0.3411957     0.6802598
#> 4   0.5968252     0.9280079
#> 5   0.3874782     0.6418565
#> 6   0.4722751     0.7679877
#> 7   0.5298979     0.8166548
#> 8   0.7034043     1.2114875
#> 9   0.6853213     1.1000959
#> 10  0.7745896     1.2464374
#> 11  0.7217774     1.1485244
#> 12  0.3883457     0.6927739
#> 13  0.7146268     1.2903305
#> 14  0.7426744     1.1786959
#> 15  1.1009874     2.4340152
#> 16  0.3520825     0.5975389
#> 17  0.5949336     0.9279192
#> 18  0.4149519     0.9039435
#> 19  0.4072842     0.8343537
#> 20  0.6713441     1.0377240
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
#> 1 2.5 0.2  1  0 0.5667047 0.2234783 1.0020773 0.7997323  0.6671570
#> 2 3.4 0.5  2  0 0.5000331 0.2939882 0.7006621 0.7452853  0.6930127
#> 3 6.0 1.3  3  0 0.3349156 0.4743967 0.3803199 0.6222603  0.5464481
#>   se_basehazard
#> 1     1.2605187
#> 2     1.1658347
#> 3     0.8777995
```
