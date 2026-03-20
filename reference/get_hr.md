# Get the Hazard Ratios with their Standard Errors

`get_hr()` takes as input the results of a model with covariates
estimated by
[`fit2ts()`](https://angelacar.github.io/TwoTimeScales/reference/fit2ts.md)
or
[`fit1ts()`](https://angelacar.github.io/TwoTimeScales/reference/fit1ts.md)
and returns the estimated hazard ratios together with their standard
errors.

## Usage

``` r
get_hr(fitted_model)
```

## Arguments

- fitted_model:

  A list returned by the function
  [`fit2ts()`](https://angelacar.github.io/TwoTimeScales/reference/fit2ts.md)
  or
  [`fit1ts()`](https://angelacar.github.io/TwoTimeScales/reference/fit1ts.md).

## Value

A list with the following elements:

- `HR` A vector of hazard ratios (calculated as \\\exp(\hat\beta)\\).

- `SE_HR` A vector of Standard Errors for the hazard ratios calculated
  via the delta method.

- `beta` A vector of the estimated \\\hat\beta\\ coefficients.

- `SE_beta` A vector of the Standard Errors for the beta coefficients.

## Examples

``` r
# Create some fake data - the bare minimum
id <- 1:20
u <- c(5.43, 3.25, 8.15, 5.53, 7.28, 6.61, 5.91, 4.94, 4.25, 3.86, 4.05, 6.86,
       4.94, 4.46, 2.14, 7.56, 5.55, 7.60, 6.46, 4.96)
s <- c(0.44, 4.89, 0.92, 1.81, 2.02, 1.55, 3.16, 6.36, 0.66, 2.02, 1.22, 3.96,
       7.07, 2.91, 3.38, 2.36, 1.74, 0.06, 5.76, 3.00)
ev <- c(1, 0, 0, 1, 0, 1, 0, 1, 1, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 1)
x1 <- c(0, 0, 0, 1, 0, 1, 1, 1, 1, 0, 0, 0, 1, 0, 1, 1, 0, 0, 1, 0)

fakedata <- as.data.frame(cbind(id, u, s, ev, x1))
fakedata2ts <- prepare_data(data = fakedata,
                            u = "u",
                            s_out = "s",
                            ev = "ev",
                            ds = .5,
                            individual = TRUE,
                            covs = "x1")
#> `s_in = NULL`. I will use `s_in = 0` for all observations.
#> `s_in = NULL`. I will use `s_in = 0` for all observations.
# Fit a fake model - not optimal smoothing
fakemod <- fit2ts(fakedata2ts,
                  optim_method = "grid_search",
                  lrho = list(seq(1, 1.5, .5),
                              seq(1, 1.5, .5)))
get_hr(fakemod)
#> $beta
#>         x1 
#> -0.3833005 
#> 
#> $SE_beta
#> [1] 0.7593286
#> 
#> $HR
#>       x1 
#> 0.681608 
#> 
#> $SE_HR
#>        x1 
#> 0.5175645 
#> 
```
