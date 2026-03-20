# Survival function with two time scales

Computes the survival matrix containing the probability of not
experiencing an event (of any cause) by time `s` and fixed entry time
`u`. The survival function can be obtained from one fitted model with
only one event type, or combining information from several
cause-specific hazard in a competing risks model. In the first case, a
fitted object of class `'haz2ts'`, `'haz2tsLMM'`, `'haz2tsPGAM'`, or
`'haz2tsVCM'` can be passed directly as argument to the function. In the
competing risks framework, the user should provide a list of
cause-specific cumulative hazard matrices. The function is also called
internally from [`plot()`](https://rdrr.io/r/graphics/plot.default.html)
if the user wants to plot the cumulative hazard from a fitted model.

## Usage

``` r
surv2ts(
  cumhaz = NULL,
  fitted_model = NULL,
  plot_grid = NULL,
  cause = NULL,
  midpoints = FALSE,
  where_slices = NULL,
  direction = c("u", "s", NULL),
  tmax = NULL
)
```

## Arguments

- cumhaz:

  (optional) a list with all the cause-specific cumulated hazard
  matrices (minimum one element needs to be supplied). If more than one
  cause-specific cumulated hazard is provided, then they should all be
  matrices of the same dimension.

- fitted_model:

  (optional) An object of class `'haz2ts'`, `'haz2tsLMM'`,
  `'haz2tsPGAM'`, or `'haz2tsVCM'`.

- plot_grid:

  (optional) A list containing the parameters to build a new finer grid
  of intervals over `u` and `s` for plotting. This must be of the form:

  - `plot_grid = list(c(umin, umax, du), c(smin, smax, ds))` where
    `umin`, `umax` and `smin`, `smax` are the minimum and maximum values
    desired for the grid-points over `u` and `s` respectively, and `du`,
    `ds` are distances between points over `u` and `s` respectively.
    Specifying a new denser grid is used to evaluate the B-spline bases
    used for estimation on such grid and plot the estimated surfaces
    with a greater level of details. If not specified, the plotting is
    done using the same B-splines bases as for the estimation. The
    function will check if the parameters for the grid provided by the
    user are compatible with those originally used to construct the
    B-splines for estimating the model. If not, the grid will be
    adjusted accordingly and a warning will be returned.

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
  the user selects `which_plot = "survival"` and `surv_slices = TRUE`,
  then `where_slices` indicates the location of the cutting points over
  the `u` time.

- direction:

  If cross-sectional one-dimensional curves are plotted, this indicates
  whether the cutting points are located on the `u` time, or on the `s`
  time. For plots of the survival function, only cutting points over the
  `u` time are meaningful.

- tmax:

  The maximum value of `t` that should be plotted.

## Value

a matrix containing the values of the survival function over `s` and
`u`.

## Examples

``` r
# Create some fake data - the bare minimum
id <- 1:20
u <- c(5.43, 3.25, 8.15, 5.53, 7.28, 6.61, 5.91, 4.94, 4.25, 3.86, 4.05, 6.86,
       4.94, 4.46, 2.14, 7.56, 5.55, 7.60, 6.46, 4.96)
s <- c(0.44, 4.89, 0.92, 1.81, 2.02, 1.55, 3.16, 6.36, 0.66, 2.02, 1.22, 3.96,
       7.07, 2.91, 3.38, 2.36, 1.74, 0.06, 5.76, 3.00)
ev <- c(1, 0, 0, 1, 0, 1, 0, 1, 1, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 1)


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
                  lrho = list(seq(1 , 1.5, .5),
                              seq(1 , 1.5, .5)))

# Obtain the fake cumulated hazard
fakecumhaz2ts <- cumhaz2ts(fakemod)
# Fake survival curve
fakesurv2ts <- surv2ts(fitted_model = fakemod)
```
