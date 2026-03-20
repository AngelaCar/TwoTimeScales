# Summary function for object of class 'haz2tsLMM'

Summary function for object of class 'haz2tsLMM'

## Usage

``` r
haz2tsLMM_summary(x, ...)
```

## Arguments

- x:

  an object of class 'haz2tsLMM' returned by the function
  [`fit2ts()`](https://angelacar.github.io/TwoTimeScales/reference/fit2ts.md)

- ...:

  further arguments

## Value

a printed summary of the fitted model, including optimal smoothing
paramters, the effective dimension ED and the AIC/BIC. For model with
covariates, a regression table is also returned.

## Examples

``` r
# Create some fake data - the bare minimum
id <- 1:20
u <- c(
  5.43, 3.25, 8.15, 5.53, 7.28, 6.61, 5.91, 4.94, 4.25, 3.86, 4.05, 6.86,
  4.94, 4.46, 2.14, 7.56, 5.55, 7.60, 6.46, 4.96
)
s <- c(
  0.44, 4.89, 0.92, 1.81, 2.02, 1.55, 3.16, 6.36, 0.66, 2.02, 1.22, 3.96,
  7.07, 2.91, 3.38, 2.36, 1.74, 0.06, 5.76, 3.00
)
ev <- c(1, 0, 0, 1, 0, 1, 0, 1, 1, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 1) #'

fakedata <- as.data.frame(cbind(id, u, s, ev))
fakedata2ts <- prepare_data(data = fakedata,
                            u = "u",
                            s_out = "s",
                            ev = "ev",
                            ds = .5)
#> `s_in = NULL`. I will use `s_in = 0` for all observations.
#> `s_in = NULL`. I will use `s_in = 0` for all observations.
# Fit a fake model - not optimal smoothing
fakemod <- fit2ts(fakedata2ts,
  optim_method = "LMMsolver"
)
summary(fakemod)
#> Number of events =  8 
#> Model specifications:
#>   nu =  13 
#>   ns =  15 
#>   cu =  13 
#>   cs =  13 
#> 
#> Optimal smoothing: 
#>   log10(rho_u) =  -0.2779936 
#>   log10(rho_s) =  6.832826 
#>   rho_u =  0.5272377 
#>   rho_s =  6804970 
#> 
#> Model with no covariates
#> 
#> Model diagnostics: 
#>   AIC =  41.81362 
#>   BIC =  53.51648 
#>   ED =  4.887498
```
