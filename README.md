
<!-- README.md is generated from README.Rmd. Please edit that file -->

# TwoTimeScales

<!-- badges: start -->
<!-- badges: end -->

`TwoTimeScales` provides a collection of functions for the analysis of
time to event data with two time scales.

The main goal of the analysis is to estimate a smooth two dimensional
hazard curve, while also respecting the constraint imposed by the
relation between the two time scales. Analyses are based on the
P-splines model of Poisson data with offset. Regression model to study
the effect of covariates on the baseline hazard are possible in a
proportional hazards fashion. Functions to plot the (baseline) hazard
are also provided.

## Installation

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
package without vignettes.

## How to use TwoTimeScales

A general introduction to the functionalities of the package is given
in:

``` r
vignette("TwoTimeScales")
```

Further examples on can be found in the vignettes.

``` r
vignette("onetime")
vignette("twotimes")
vignette("visualization")
```

Or read them in the package website.
