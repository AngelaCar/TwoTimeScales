# Summary function for object of class 'haz2tsPGAM'

Summary function for object of class 'haz2tsPGAM'

## Usage

``` r
haz2tsPGAM_summary(x, ...)
```

## Arguments

- x:

  an object of class 'haz2tsPGAM' returned by the function
  [`fitpgam()`](https://angelacar.github.io/TwoTimeScales/reference/fitpgam.md)

- ...:

  further arguments

## Value

a printed summary of the fitted model

## Examples

``` r
# Create some fake data - the bare minimum
id <- 1:20
u <- c(5.43, 3.25, 8.15, 5.53, 7.28, 6.61, 5.91, 4.94, 4.25, 3.86, 4.05, 6.86,
       4.94, 4.46, 2.14, 7.56, 5.55, 7.60, 6.46, 4.96)
s <- c(0.44, 4.89, 0.92, 1.81, 2.02, 1.55, 3.16, 6.36, 0.66, 2.02, 1.22, 3.96,
       7.07, 2.91, 3.38, 2.36, 1.74, 0.06, 5.76, 3.00)
ev <- c(1, 0, 0, 1, 0, 1, 0, 1, 1, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 1)#'

fakedata <- as.data.frame(cbind(id, u, s, ev))
fakedata2ts <- prepare_data(data = fakedata,
                            u = "u",
                            s_out = "s",
                            ev = "ev",
                            ds = .5)
#> `s_in = NULL`. I will use `s_in = 0` for all observations.
#> `s_in = NULL`. I will use `s_in = 0` for all observations.
# Fit a fake model - not optimal smoothing
fakemod <- fitpgam(fakedata2ts,
                  optim_method = "grid_search",
                  lrho = list(seq(1, 1.5, .5),
                              seq(1, 1.5, .5)))
summary(fakemod)
#> Number of events =  
#> Model specifications:
#>   nu =  195 
#>   ns =  195 
#>   cu =  13 
#>   cs =  13 
#> 
#> Optimal smoothing: 
#>   log10(rho_u) =  1 
#>   log10(rho_s) =  1.5 
#>   rho_u =  10 
#>   rho_s =  31.62278 
#> 
#> Model with no covariates
#> 
#> Model diagnostics: 
#>   AIC =  39.98714 
#>   BIC =  49.07529 
#>   ED =  3.79551  EDu =  2.01417  EDs =  1.78134
```
