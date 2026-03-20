# Cumulative hazard over two time scales

Computes the cumulative hazard surface over two time scales from a
fitted model. The function is also called internally from
[`plot()`](https://rdrr.io/r/graphics/plot.default.html) if the user
wants to plot the cumulative hazard from a fitted model.

## Usage

``` r
cumhaz2ts(
  fitted_model,
  plot_grid = NULL,
  cause = NULL,
  midpoints = FALSE,
  where_slices = NULL,
  direction = c("u", "s", NULL),
  tmax = NULL
)
```

## Arguments

- fitted_model:

  (optional) An object of class `'haz2ts'`, `'haz2tsLMM'`,
  `'haz2tsPGAM'`, or `'haz2tsVCM'`.

- plot_grid:

  (optional) A list containing the parameters to build a new finer grid
  of intervals over `u` and `s` for plotting. This must be of the form:

  - `plot_grid = list(c(umin, umax, du), c(smin, smax, ds))` where
    `umin`, `umax` and `smin`, `smax` are the minimum and maximum values
    desired for the grid-points over `u` and `s` respectively, and `du`,
    `ds` are distances between two adjacent points over `u` and `s`
    respectively. Specifying a new denser grid is used to evaluate the
    B-spline bases used for estimation on such grid and plot the
    estimated surfaces with a greater level of detail. If not specified,
    the plotting is done using the same B-splines bases as for the
    estimation. The function will check if the parameters for the grid
    provided by the user are compatible with those originally used to
    construct the B-splines for estimating the model. If not, the grid
    will be adjusted accordingly and a warning will be returned.

- cause:

  a character string with a short name for the cause (optional).

- midpoints:

  A Boolean. Default is `FALSE`. If `TRUE`, the estimated quantities are
  evaluated at the midpoints of the rectangles (or parallelograms) of
  the grids, rather than at each grid-point.

- where_slices:

  A vector of values for the cutting points of the desired slices of the
  surface. This option is included mostly for the plotting function.
  When using
  [`plot.haz2ts()`](https://angelacar.github.io/TwoTimeScales/reference/plot.haz2ts.md),
  the user selects `which_plot = "cumhaz"` and `cumhaz_slices = TRUE`,
  then `where_slices` indicates the location of the cutting points over
  the `u` time.

- direction:

  If cross-sectional one-dimensional curves are plotted, this indicates
  whether the cutting points are located on the `u` time, or on the `s`
  time. For plots of the cumulative hazards, only cutting points over
  the `u` time are meaningful.

- tmax:

  The maximum value of `t` that should be plotted.

## Value

A list with the following elements: \* `Haz` a list of estimated hazard
and associated SEs (obtained from the function
[`get_hazard_2d()`](https://angelacar.github.io/TwoTimeScales/reference/get_hazard_2d.md));
\* `CumHaz` the cumulated hazard estimate over `u` and `s`; \* `cause`
(if provided) the short name for the cause.

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
fakedata2ts <- prepare_data(u = fakedata$u,
                            s_out = fakedata$s,
                            ev = fakedata$ev,
                            ds = .5)
#> `s_in = NULL`. I will use `s_in = 0` for all observations.
#> `s_in = NULL`. I will use `s_in = 0` for all observations.
# Fit a fake model - not optimal smoothing
fakemod <- fit2ts(fakedata2ts,
                  optim_method = "grid_search",
                  lrho = list(seq(1 ,1.5 ,.5),
                              seq(1 ,1.5 ,.5)))

# Obtain the fake cumulated hazard
fakecumhaz2ts <- cumhaz2ts(fakemod)
```
