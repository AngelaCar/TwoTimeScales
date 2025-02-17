
<!-- README.md is generated from README.Rmd. Please edit that file -->

# TwoTimeScales

<!-- badges: start -->

[![R-CMD-check](https://github.com/AngelaCar/TwoTimeScales/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/AngelaCar/TwoTimeScales/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

`TwoTimeScales` provides a suite of functions for the analysis of time
to event data with two time scales.

The package provides tools to estimate a smooth hazard that varies over
two time scales and also, if covariates are available, to estimate a
proportional hazards model with such a two-dimensional baseline hazard.
The hazard of one event is modelled by means of two-dimensional
P-splines model of Poisson data with offset. Regression model to study
the effect of covariates on the baseline hazard are possible in a
proportional hazards fashion. Functions to plot the (baseline) hazard
are also provided.

Additionally, functions to estimate cumulated surfaces (such as
cumulated hazard, survival function and incidence function) over two
time scales are included.

## Installation

### From CRAN

Just run install.packages(“TwoTimeScales”)

### From GitHub

You can install the development version of TwoTimeScales from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("AngelaCar/TwoTimeScales")
```

Vignettes need to be explicitly requested. To install a version of the
package including the vignettes, please run:

``` r
devtools::install_github("AngelaCar/TwoTimeScales",
                          dependencies = TRUE, build_vignettes = TRUE)
```

Note: This will take more time than installing the version of the
package without vignettes. Alternatively, we suggest to read the
vignettes online here:
[TwoTimeScales](https://angelacar.github.io/TwoTimeScales/)

## How to use TwoTimeScales

A general introduction to the usage of the package is given in:

``` r
vignette("TwoTimeScales")
```

or here
[Introduction](https://angelacar.github.io/TwoTimeScales/articles/TwoTimeScales.html).

Further examples on can be found in the vignettes.

``` r
vignette("onetime")
vignette("twotimes")
vignette("visualization")
```

## Acknowledgments

We are very grateful to Guillermo Briseno Sanchez, José Carlos Andrade
Santacruz and Lisa Rieker for testing the package and their helpful
feedback.
