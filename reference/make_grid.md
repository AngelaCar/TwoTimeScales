# Make a grid of points to evaluate B-splines

Make a grid of points to evaluate B-splines

## Usage

``` r
make_grid(
  plot_grid,
  model_specifications,
  class_fitmodel,
  where_slices = NULL,
  direction = c(NULL, "u", "s"),
  tmax = NULL,
  midpoints = FALSE
)
```

## Arguments

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

- model_specifications:

  A list with `'umin'`, `'umax'`, `'smin'`, and `'smax'` used for
  fitting.

- class_fitmodel:

  The class of the fitted model. Can be `'haz2ts'`, `'haz2tsPGAM'`, or
  `'haz2tsVCM'`.

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

- tmax:

  The maximum value of `t` that should be plotted.

- midpoints:

  A Boolean. Default is `FALSE`. If `TRUE`, the estimated quantities are
  evaluated at the midpoints of the rectangles (or parallelograms) of
  the grids, rather than at each grid-point.

## Value

A list of specification for the grid

## Examples

``` r
make_grid(plot_grid = list(
    c(umin = 3, umax = 8.5, du = .1),
    c(smin = 0, smax = 7.1, ds = .1)
  ),
  model_specification = list(umin = 2, umax = 8.5, smin = 0, smax = 7.5),
  class_fitmodel = "haz2ts")
#> $intu
#>  [1] 3.0 3.1 3.2 3.3 3.4 3.5 3.6 3.7 3.8 3.9 4.0 4.1 4.2 4.3 4.4 4.5 4.6 4.7 4.8
#> [20] 4.9 5.0 5.1 5.2 5.3 5.4 5.5 5.6 5.7 5.8 5.9 6.0 6.1 6.2 6.3 6.4 6.5 6.6 6.7
#> [39] 6.8 6.9 7.0 7.1 7.2 7.3 7.4 7.5 7.6 7.7 7.8 7.9 8.0 8.1 8.2 8.3 8.4 8.5
#> 
#> $umin
#> umin 
#>    3 
#> 
#> $umax
#> umax 
#>  8.5 
#> 
#> $du
#>  du 
#> 0.1 
#> 
#> $ints
#>  [1] 0.0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0 1.1 1.2 1.3 1.4 1.5 1.6 1.7 1.8
#> [20] 1.9 2.0 2.1 2.2 2.3 2.4 2.5 2.6 2.7 2.8 2.9 3.0 3.1 3.2 3.3 3.4 3.5 3.6 3.7
#> [39] 3.8 3.9 4.0 4.1 4.2 4.3 4.4 4.5 4.6 4.7 4.8 4.9 5.0 5.1 5.2 5.3 5.4 5.5 5.6
#> [58] 5.7 5.8 5.9 6.0 6.1 6.2 6.3 6.4 6.5 6.6 6.7 6.8 6.9 7.0 7.1
#> 
#> $smin
#> smin 
#>    0 
#> 
#> $smax
#> smax 
#>  7.1 
#> 
#> $ds
#>  ds 
#> 0.1 
#> 
```
