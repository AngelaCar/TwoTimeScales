# Plot method for a haz2tsVCM object.

`plot.haz2tsVCM()` is the plot method for objects of class `haz2tsVCM`.
It produces several kinds of plots of the fitted varying coefficients
model with two time scales (see
[`fitvcm()`](https://angelacar.github.io/TwoTimeScales/reference/fitvcm.md)),
either on the original (t,s) plane, while respecting the constraint
imposed by the relation of the two time scales, or on the transformed
(u,s) plane.

## Usage

``` r
# S3 method for class 'haz2tsVCM'
plot(
  x,
  plot_grid = NULL,
  which_plot = c("hazard", "SE", "slices", "survival", "cumhaz"),
  where_slices = NULL,
  direction = c(NULL, "u", "s"),
  plot_options = list(),
  ...
)
```

## Arguments

- x:

  The output of the function
  [`fitvcm()`](https://angelacar.github.io/TwoTimeScales/reference/fitvcm.md).
  This is an object of class `"haz2tsVCM"`.

- plot_grid:

  (optional) A list containing the parameters to build a new finer grid
  of intervals over u and s for plotting. This must be of the form:
  `plot_grid = list(c(umin, umax, du), c(smin, smax, ds))`, where
  `umin`, `umax` and `smin`, `smax` are the minimum and maximum values
  desired for the intervals over `u` and `s` respectively, and `du`,
  `ds` are distances between intervals over `u` and `s` respectively.
  Specifying a new denser grid is used to evaluate the B-spline bases
  used for estimation on such grid and plot the estimated surfaces with
  a greater level of details. If not specified, the plotting is done
  using the same B-splines bases as for the estimation. The function
  will check if the parameters for the grid provided by the user are
  compatible with those originally used to construct the B-splines for
  estimating the model. If not, the grid will be adjusted accordingly
  and a warning will be returned.

- which_plot:

  The type of plot required. Can be one of `"hazard"` (default), `"SE"`,
  `"slices"`, `"survival"` or `"cumhaz"` (see details section).

- where_slices:

  A vector of values for the cutting points of the desired slices of the
  surface. If `which_plot == "slices"`, please provide this argument.
  Please also provide this argument in case `which_plot = "survival` or
  `which_plot = "cumhaz` and `surv_slices = TRUE` or
  `cumhaz_slices = TRUE`, respectively.

- direction:

  If `which_plot == "slices"`, indicates the direction for cutting the
  surface. If `u`, then the surface will be cut at the selected values
  of `u` (indicated by `where_slices`), hence obtaining one-dimensional
  curves over `s`. If `s`, then the surface will be cut at the selected
  values of `s` (indicated by `where_slices`), hence obtaining
  one-dimensional curves over `u`.

- plot_options:

  A list with all possible options for any of the plots:

  - `loghazard` A Boolean. Default is `FALSE`. If `FALSE` the function
    returns a plot of the hazard surface, if `TRUE` the function returns
    a plot of the log-hazard surface.

  - `log10hazard` A Boolean. Default is `FALSE`. If `TRUE`, then a
    log_10 hazard surface is plotted.

  - `cut_extrapolated` A Boolean. Default is `TRUE`. Cuts away the
    extrapolated area of the (log-)hazard surface before plotting.

  - `rectangular_grid` A Boolean. Default is `FALSE`. If `TRUE`, a
    rectangular grid is used for plotting also in the (t,s)-plane as
    opposed to the grid of parallelograms used as default in the
    (t,s)-plane.

  - `original` A Boolean. Default is `TRUE`. Plot the (log-)hazard
    (and/or the SEs) in the (t,s)-plane. If `FALSE`, the (log-)hazard
    (and/or the SEs) will be plotted in the (u,s)-plane.

  - `tmax` The maximum value of `t` that should be plotted.

  - `surv_slices` A Boolean. Default is `FALSE`. If `TRUE` and
    `which_plot == "survival"`, plot survival curves over the time `s`
    for selected values of `u`, that are cross-sections of the 2D
    survival surface.

  - `cumhaz_slices` A Boolean. Default is `FALSE`. If `TRUE` and
    `which_plot == "cumhaz"`, plot cumulative hazards curves over the
    time `s` for selected values of `u`, that are cross-sections of the
    2D cumulative hazard surface.

  - `midpoints` A Boolean. Default is `FALSE`. If `TRUE`, the estimated
    quantities (hazard, survival, etc.) will be evaluated in the
    mid-points of the bins rather than at the extremes. Set to `TRUE` if
    plotting estimated number of events.

  - `col_palette` A function defining the color palette. The default
    palette is `viridis::rev(plasma())`. Specifying the color palette as
    a function allows for greater flexibility than passing the palette
    as a vector.

  - `n_shades` The number of color shades to plot, default is 50.

  - `breaks` The vector of breaks for the color legend. If `n_shades` is
    provided, this should be of length `n_shades + 1`.

  - `show_legend` A Boolean. Default is `TRUE`. If `FALSE` no legend
    will be plotted, useful for multi-panel figures with common legend.
    Works only for plots on rectangular grid (i.e. transformed (u,s)
    plane)

  - `main` The title of the plot.

  - `xlab` The label of the first time axis (plotted on the x axis).

  - `ylab` The label of the second time axis (plotted on the y axis).

  - `xlim` A vector with two elements defining the limits of the time
    scale on the x axis.

  - `ylim` A vector with two elements defining the limits of the time
    scale on the y axis.

  - `contour_lines` A Boolean. Default is `FALSE`. If `TRUE` white
    contour lines are added to the surfaces.

  - `contour_col` The color for the contour lines. Default is `white`.

  - `contour_cex` The magnification to be used for the contour lines.
    Default is `.8`.

  - `contour_nlev` The number of contour levels desired. Default is 10.

  - `cex_main` The magnification to be used for the main title, default
    is 1.2 .

  - `cex_lab` The magnification to be used for the axis labels, default
    is 1 .

  - `lwd` The line width.

  - `lty` The line type.

- ...:

  Further arguments to image.plot or image

## Value

A plot of the fitted model.

## Details

The vignette "visualization" presents and discusses all the different
plotting options for the fitted model over two time scales. In most of
the cases, the user will want to visualize the hazard surface over the
two time scales. This can be plotted on the hazard scale, the log-hazard
scale or the log10-hazard scale, by switching to `TRUE` the
corresponding argument in `plot_options`. The survival and cumulative
hazard functions can be plotted as two-dimensional surfaces over `u` and
`s` or `t` and `s`. However, it is also very informative to plot them as
one-dimensional curves over `s` (cross-sections or slices). This is done
by selecting `which_plot = "survival"` and `surv_slices = TRUE` in
`plot_options`. Additionally, a vector of values for the cutting points
over the `u`-axis should be passed to the argument `where_slices`,
together with setting `direction = u`. Similar plot is obtained for the
cumulative hazard by selecting `which_plot = "cumhaz"`,
`cumhaz_slices = TRUE`, see examples section. Please, notice that for
the survival function and the cumulative hazard, only cross-sections of
the surface for selected values of `u` (over the `s` time) can be
plotted.

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
fakedata2ts <- prepare_data(
  u = fakedata$u,
  s_out = fakedata$s,
  ev = fakedata$ev,
  ds = .5
)
#> `s_in = NULL`. I will use `s_in = 0` for all observations.
#> `s_in = NULL`. I will use `s_in = 0` for all observations.
# Fit a fake model - not optimal smoothing
fakemod <- fitvcm(fakedata2ts)

# plot the hazard surface
plot(fakemod)


# plot the survival function as one-dimension curves over `s`
plot(fakemod,
  which_plot = "survival",
  direction = "u",
  where_slices = c(4, 6, 8),
  plot_options = list(
    surv_slices = TRUE
  )
)


# Plot cross-sections of the hazard over `s` for selected values of `u`

plot(fakemod,
  which_plot = "slices",
  where_slices = c(4, 6, 8),
  direction = "u",
  plot_options = list(
    main = "Cross-sections of the hazard",
    xlab = "Time",
    ylab = "Hazard"
  )
)

```
