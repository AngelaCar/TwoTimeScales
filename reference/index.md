# Package index

## Data

- [`reccolon2ts`](https://angelacar.github.io/TwoTimeScales/reference/reccolon2ts.md)
  : Data from the chemotherapy for stace B/C colon cancer study

## Data processing

- [`prepare_data()`](https://angelacar.github.io/TwoTimeScales/reference/prepare_data.md)
  : Prepare raw data by binning them in 1d or 2d

- [`prepare_data_LMMsolver()`](https://angelacar.github.io/TwoTimeScales/reference/prepare_data_LMMsolver.md)
  : Process data to fit model with LMMsolver

- [`make_bins()`](https://angelacar.github.io/TwoTimeScales/reference/make_bins.md)
  : Construct bins over one or more time axes

- [`exposures_events_1d()`](https://angelacar.github.io/TwoTimeScales/reference/exposures_events_1d.md)
  : Bin data on one time scale

- [`exposures_events_2d()`](https://angelacar.github.io/TwoTimeScales/reference/exposures_events_2d.md)
  : Bin data on two time scales

- [`exposures_events_Lexis()`](https://angelacar.github.io/TwoTimeScales/reference/exposures_events_Lexis.md)
  : Bin data on the Lexis diagram

- [`print(`*`<data2ts>`*`)`](https://angelacar.github.io/TwoTimeScales/reference/print.data2ts.md)
  :

  Print method for a `data2ts` object

## One time scale

- [`fit1ts()`](https://angelacar.github.io/TwoTimeScales/reference/fit1ts.md)
  : Fit a smooth hazard model with one time scale
- [`fit1tsmodel_ucminf()`](https://angelacar.github.io/TwoTimeScales/reference/fit1tsmodel_ucminf.md)
  : Numerical optimization of the 1ts model
- [`grid_search_1d()`](https://angelacar.github.io/TwoTimeScales/reference/grid_search_1d.md)
  : Grid search for the optimal 1ts model
- [`get_aic_fit_1d()`](https://angelacar.github.io/TwoTimeScales/reference/get_aic_fit_1d.md)
  : Return the AIC of 1ts model
- [`get_bic_fit_1d()`](https://angelacar.github.io/TwoTimeScales/reference/get_bic_fit_1d.md)
  : Return the BIC of 1ts model
- [`iwls_1d()`](https://angelacar.github.io/TwoTimeScales/reference/iwls_1d.md)
  : Iterative Weighted Least Squares algorithm for 1ts model
- [`GLAM_1d_covariates()`](https://angelacar.github.io/TwoTimeScales/reference/GLAM_1d_covariates.md)
  : Fit the 1d GLAM with covariates
- [`get_hazard_1d()`](https://angelacar.github.io/TwoTimeScales/reference/get_hazard_1d.md)
  : Get estimated (log-)hazard values with 1 time scale
- [`get_hazard_1d_LMM()`](https://angelacar.github.io/TwoTimeScales/reference/get_hazard_1d_LMM.md)
  : Get estimated (log-)hazard values with 1 time scale
- [`get_hr()`](https://angelacar.github.io/TwoTimeScales/reference/get_hr.md)
  : Get the Hazard Ratios with their Standard Errors

## Two time scales

- [`fit2ts()`](https://angelacar.github.io/TwoTimeScales/reference/fit2ts.md)
  : Fit a smooth hazard model with two time scales

- [`fit2tsmodel_ucminf()`](https://angelacar.github.io/TwoTimeScales/reference/fit2tsmodel_ucminf.md)
  : Numerical optimization of the 2ts model

- [`fitpgam()`](https://angelacar.github.io/TwoTimeScales/reference/fitpgam.md)
  : Fit a log-additive model over two time scales

- [`fitvcm()`](https://angelacar.github.io/TwoTimeScales/reference/fitvcm.md)
  : Fit a Varying Coefficient Model (VCM) over two time scales

- [`grid_search_2d()`](https://angelacar.github.io/TwoTimeScales/reference/grid_search_2d.md)
  : Grid search for the optimal 2ts model

- [`get_aic_fit_2d()`](https://angelacar.github.io/TwoTimeScales/reference/get_aic_fit_2d.md)
  : Return the AIC of 2ts model

- [`get_bic_fit_2d()`](https://angelacar.github.io/TwoTimeScales/reference/get_bic_fit_2d.md)
  : Return the BIC of 2ts model

- [`getAIC_BIC_LMM()`](https://angelacar.github.io/TwoTimeScales/reference/getAIC_BIC_LMM.md)
  : Calculates AIC and BIC from object fitted via LMMsolver

- [`GLAM_2d_no_covariates()`](https://angelacar.github.io/TwoTimeScales/reference/GLAM_2d_no_covariates.md)
  : Fit the 2d GLAM without covariates

- [`GLAM_2d_covariates()`](https://angelacar.github.io/TwoTimeScales/reference/GLAM_2d_covariates.md)
  : Fit the 2d GLAM with covariates

- [`get_hazard_2d()`](https://angelacar.github.io/TwoTimeScales/reference/get_hazard_2d.md)
  : Get estimated (log-)hazard surface with 2 time scales

- [`get_hazard_2d_LMM()`](https://angelacar.github.io/TwoTimeScales/reference/get_hazard_2d_LMM.md)
  : Get estimated (log-)hazard surface with 2 time scales

- [`get_hr()`](https://angelacar.github.io/TwoTimeScales/reference/get_hr.md)
  : Get the Hazard Ratios with their Standard Errors

- [`haz2ts_summary()`](https://angelacar.github.io/TwoTimeScales/reference/haz2ts_summary.md)
  : Summary function for object of class 'haz2ts'

- [`haz2tsLMM_summary()`](https://angelacar.github.io/TwoTimeScales/reference/haz2tsLMM_summary.md)
  : Summary function for object of class 'haz2tsLMM'

- [`haz2tsPGAM_summary()`](https://angelacar.github.io/TwoTimeScales/reference/haz2tsPGAM_summary.md)
  : Summary function for object of class 'haz2tsPGAM'

- [`haz2tsVCM_summary()`](https://angelacar.github.io/TwoTimeScales/reference/haz2tsVCM_summary.md)
  : Summary function for object of class 'haz2tsVCM'

- [`make_grid()`](https://angelacar.github.io/TwoTimeScales/reference/make_grid.md)
  : Make a grid of points to evaluate B-splines

- [`predict_haz2ts_pointwise()`](https://angelacar.github.io/TwoTimeScales/reference/predict_haz2ts_pointwise.md)
  : Point-wise prediction hazard 2 time scale

- [`predict_haz2ts()`](https://angelacar.github.io/TwoTimeScales/reference/predict_haz2ts.md)
  :

  Prediction method for objects of class `'haz2ts'`

- [`select_model2ts()`](https://angelacar.github.io/TwoTimeScales/reference/select_model2ts.md)
  : Selects best model among a group of fitted models

## Plotting functions

- [`plot(`*`<haz1ts>`*`)`](https://angelacar.github.io/TwoTimeScales/reference/plot.haz1ts.md)
  : Plot method for a haz1ts object.
- [`plot(`*`<haz1tsLMM>`*`)`](https://angelacar.github.io/TwoTimeScales/reference/plot.haz1tsLMM.md)
  : Plot method for a haz1ts object.
- [`plot(`*`<haz2ts>`*`)`](https://angelacar.github.io/TwoTimeScales/reference/plot.haz2ts.md)
  : Plot method for a haz2ts object.
- [`plot(`*`<haz2tsLMM>`*`)`](https://angelacar.github.io/TwoTimeScales/reference/plot.haz2tsLMM.md)
  : Plot method for a haz2tsLMM object.
- [`plot(`*`<haz2tsPGAM>`*`)`](https://angelacar.github.io/TwoTimeScales/reference/plot.haz2tsPGAM.md)
  : Plot method for a haz2tsPGAM object.
- [`plot(`*`<haz2tsVCM>`*`)`](https://angelacar.github.io/TwoTimeScales/reference/plot.haz2tsVCM.md)
  : Plot method for a haz2tsVCM object.
- [`plot_slices()`](https://angelacar.github.io/TwoTimeScales/reference/plot_slices.md)
  : Plot slices of the (log-) hazard
- [`imageplot_2ts()`](https://angelacar.github.io/TwoTimeScales/reference/imageplot_2ts.md)
  : Image Plot of 2ts hazard
- [`imageplot_SE()`](https://angelacar.github.io/TwoTimeScales/reference/imageplot_SE.md)
  : Image Plot of Standard Errors for the 2ts hazard
- [`covariates_plot()`](https://angelacar.github.io/TwoTimeScales/reference/covariates_plot.md)
  : Plot of the covariates' effects

## Competing risks functions

- [`cumhaz2ts()`](https://angelacar.github.io/TwoTimeScales/reference/cumhaz2ts.md)
  : Cumulative hazard over two time scales
- [`surv2ts()`](https://angelacar.github.io/TwoTimeScales/reference/surv2ts.md)
  : Survival function with two time scales
- [`cuminc2ts()`](https://angelacar.github.io/TwoTimeScales/reference/cuminc2ts.md)
  : Cumulative incidence surface over two time scales
- [`predict_cif2ts_pointwise()`](https://angelacar.github.io/TwoTimeScales/reference/predict_cif2ts_pointwise.md)
  : Point-wise prediction of cumulative incidence over 2 time scale
- [`predict_comprisks2ts()`](https://angelacar.github.io/TwoTimeScales/reference/predict_comprisks2ts.md)
  : Predict overall survival and cumulative incidence functions
