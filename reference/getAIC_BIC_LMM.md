# Calculates AIC and BIC from object fitted via LMMsolver

`getAIC_BIC_LMM` is an utility function that takes an object of class
`'LMMsolve'` fitted via
[`fit1ts()`](https://angelacar.github.io/TwoTimeScales/reference/fit1ts.md)
or
[`fit2ts()`](https://angelacar.github.io/TwoTimeScales/reference/fit2ts.md)
and calculates AIC, BIC and ED.

## Usage

``` r
getAIC_BIC_LMM(fit, offset)
```

## Arguments

- fit:

  An object of class `"LMMsolve"`

- offset:

  The vector of exposure times from dataLMM

## Value

A list with: \* `ED` effective dimension of the full model; \* `EDbase`
effective dimension of the baseline hazard only; \* `Dev` deviance; \*
`AIC` the aic; \* `BIC` the bic; \* `n_beta` the number of estimated
covariate parameters (if PH model).
