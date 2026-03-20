# Plot method for a haz1ts object.

`plot.haz1tsLMM()` is a plot method for objects of class `haz1tsLMM`.

## Usage

``` r
# S3 method for class 'haz1tsLMM'
plot(
  x,
  which_plot = c("hazard", "covariates"),
  plot_grid = NULL,
  plot_options = list(),
  ...
)
```

## Arguments

- x:

  The output of the function `fit1ts`.

- which_plot:

  The type of plot required. Can be one of `"hazard"` (default) or
  `"covariates"`.

- plot_grid:

  (optional) A named vector containing the parameters to build a new
  grid of intervals over `s` for plotting the estimated hazard on a
  finer grid. This must be of the form: `plot_grid = c(smin, smax, ds)`,
  where `smin`, `smax` are the minimum and maximum values desired for
  the intervals over `s`, and `ds` is the distance between intervals
  over `s`. If not specified, the plotting is done using the same
  B-splines basis as for the estimation. The function will check if the
  parameters for the grid provided by the user are compatible with those
  originally used to construct the B-splines for estimating the model.
  If not, the grid will be adjusted accordingly and a warning will be
  returned.

- plot_options:

  A list with all possible options for any of the plots:

  - `loghazard` A Boolean. Default is `FALSE`. If `FALSE` the function
    returns a plot of the hazard curve, if `TRUE` the function returns a
    plot of the log-hazard curve.

  - `log10hazard` A Boolean. Default is `FALSE`. If `TRUE` it returns a
    plot of the log10-hazard curve.

  - `col` The color of the curve plotted. Default is `"black"`.

  - `add_CI` A Boolean. If `TRUE` (default) the confidence bands will be
    added.

  - `col_CI` The color for the confidence bands. The default is the same
    color of the curve, with a 50% transparency level.

  - `main` The title of the plot.

  - `xlab` The label of the time axis (plotted on the x axis).

  - `ylab` The label of the y-axis (hazard, log-hazard or log10-hazard).

  - `xlim` A vector with two elements defining the limits of the time
    scale on the x axis.

  - `ylim` A vector with two elements defining the limits of function
    plotted on the y axis (hazard, log-hazard or log10-hazard).

  - `xmin` The minimum value on the x-axis.

  - `ymin` The minimum value on the y-axis.

  - `cex_main` The magnification to be used for the main title, default
    is 1.2 .

  - `cex_lab` The magnification to be used for the axis labels, default
    is 1 .

  - `HR` A Boolean. If `TRUE` the HRs with their CIs will be plotted.
    Default is `FALSE` (plot the `beta` with their CIs).

  - `symmetric_CI` A Boolean. Default is `TRUE`. If a plot of the HRs is
    required (`HR == TRUE`), then plot symmetrical Confidence Intervals,
    based on the SEs for the HRs calculated by delta method. If `FALSE`,
    then CIs are obtained by exponentiating the CIs for the betas.

  - `confidence` The level of confidence for the CIs. Default is .95
    (alpha = 0.05).

  - `col_beta` The color for the plot of the covariates' effects.

  - `pch` The symbol for plotting the point estimates.

  - `lwd` The line width.

  - `lty` The line type.

- ...:

  Further arguments to plot.

## Value

A plot of the type required.

## Details

The function `obtainSmoothTrend` from the R-package `LMMsolver` is used
here. We refer the interested readers to
https://biometris.github.io/LMMsolver/ for more details on `LMMsolver`
and its usage.

## Examples

``` r
## preparing data - no covariates
dt1ts <- prepare_data(data = reccolon2ts,
                      s_in = "entrys",
                      s_out = "timesr",
                      events = "status",
                      ds = 180)

## fitting the model with fit1ts() - default options

mod1 <- fit1ts(dt1ts,
optim_method = "LMMsolver")
plot(mod1)
```
