# TwoTimeScales 1.1.1

## Updates
* `predict.haz2ts()` can now predict also including covariates fixed at arbitrary 
  values.

* `GLAM_2d_covariates()` returns covariances between the alpha and beta parameters

* All functions that use ucminf to minimize the AIC/BIC of the model wrt the smoothing
  parameter(s) now have a smaller value for the option `xtol`. 
  Additionally, it can also be changed in the control lists.
  
* Fixed a small typo in `plot_haz1ts()` that did not allow to plot confidence intervals in
  color specified by user.
  
* Fixed problem with variable names in `prepare_data_LMMsolver()`

* `print.data2ts()` now prints rounded values for total exposure time

# TwoTimeScales 1.1.0

## New functions
* `predict.haz2ts()` method added to objects of class `'haz2ts'`.
  It allows prediction of the hazard, its standard errors, the cumulative hazard
  and the survival probability, from a fitted model of type `'haz2ts'`, for 
  arbitrary values of the time scales. It can also be used to obtain individual
  predictions for the original data points.
  
* `predict_comprisk2ts()` allows prediction for competing risks models. 
  It takes as input a list of cause-specific hazard models over two time scales,
  all fitted with `fit2ts()`, on the same grid, and a new data.frame with values
  of the two time scales for which predictions are requested.
  
## Minor changes
* `DESCRIPTION` has been updated, to correct a small typo.

# TwoTimeScales 1.0.0

* Initial CRAN submission.
