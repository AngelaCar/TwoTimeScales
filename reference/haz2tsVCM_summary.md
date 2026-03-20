# Summary function for object of class 'haz2tsVCM'

Summary function for object of class 'haz2tsVCM'

## Usage

``` r
haz2tsVCM_summary(x, ...)
```

## Arguments

- x:

  an object of class 'haz2tsVCM' returned by the function
  [`fitvcm()`](https://angelacar.github.io/TwoTimeScales/reference/fitvcm.md)

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
fakemod <- fitvcm(fakedata2ts)
summary(fakemod)
#> Number of events =  
#> Model specifications:
#>   ns =  195 
#>   cs =  13 
#> 
#> Optimal smoothing: 
#>   log10(theta) =  -0.7157954 
#>   log10(phi) =  0.8022907 
#>   theta =  0.1923998 
#>   phi =  6.342942 
#> 
#> Model with no covariates
#> 
#> Model diagnostics: 
#>   AIC =  47.23913 
#>   BIC =  66.04327 
#>   ED =  7.853223
```
