# Changelog

## TwoTimeScales 1.2.0

CRAN release: 2026-03-25

### New functions

- [`fitpgam()`](https://angelacar.github.io/TwoTimeScales/reference/fitpgam.md)
  and related functions added to fit a log-additive model over two time
  scales.
- [`fitvcm()`](https://angelacar.github.io/TwoTimeScales/reference/fitvcm.md)
  and related functions added to fit a varying-coefficient model over
  two time scales.
- [`select_model2ts()`](https://angelacar.github.io/TwoTimeScales/reference/select_model2ts.md)
  added to compare several fitted models and identify the best-fitting
  one.
- [`make_grid()`](https://angelacar.github.io/TwoTimeScales/reference/make_grid.md)
  facilitates the creation of a new plotting grid.

### Tests

- Added unit tests for
  [`prepare_data()`](https://angelacar.github.io/TwoTimeScales/reference/prepare_data.md),
  [`fit1ts()`](https://angelacar.github.io/TwoTimeScales/reference/fit1ts.md),
  and
  [`fit2ts()`](https://angelacar.github.io/TwoTimeScales/reference/fit2ts.md).

## TwoTimeScales 1.1.1

### Updates

- `predict.haz2ts()` can now predict also including covariates fixed at
  arbitrary values

- [`GLAM_2d_covariates()`](https://angelacar.github.io/TwoTimeScales/reference/GLAM_2d_covariates.md)
  returns covariances between the alpha and beta parameters

- All functions that use ucminf to minimize the AIC/BIC of the model wrt
  the smoothing parameter(s) now have a smaller value for the option
  `xtol`. Additionally, it can also be changed in the control lists.

- Fixed a small typo in `plot_haz1ts()` that did not allow to plot
  confidence intervals in color specified by user.

- Fixed problem with variable names in
  [`prepare_data_LMMsolver()`](https://angelacar.github.io/TwoTimeScales/reference/prepare_data_LMMsolver.md)

- [`print.data2ts()`](https://angelacar.github.io/TwoTimeScales/reference/print.data2ts.md)
  now prints rounded values for total exposure time

## TwoTimeScales 1.1.0

### New functions

- `predict.haz2ts()` method added to objects of class `'haz2ts'`. It
  allows prediction of the hazard, its standard errors, the cumulative
  hazard and the survival probability, from a fitted model of type
  `'haz2ts'`, for arbitrary values of the time scales. It can also be
  used to obtain individual predictions for the original data points.

- `predict_comprisk2ts()` allows prediction for competing risks models.
  It takes as input a list of cause-specific hazard models over two time
  scales, all fitted with
  [`fit2ts()`](https://angelacar.github.io/TwoTimeScales/reference/fit2ts.md),
  on the same grid, and a new data.frame with values of the two time
  scales for which predictions are requested.

### Minor changes

- `DESCRIPTION` has been updated, to correct a small typo.

## TwoTimeScales 1.0.0

CRAN release: 2024-12-23

- Initial CRAN submission.
