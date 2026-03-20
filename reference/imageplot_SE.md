# Image Plot of Standard Errors for the 2ts hazard

`imageplot_SE()` plots an image of the SEs of the two time scales hazard
with contour lines.

## Usage

``` r
imageplot_SE(x, y, z, plot_options = list(), ...)
```

## Arguments

- x:

  The coordinates for the x-axis. This is a vector of intervals over the
  `u` axis (default), a matrix with the corner points of the
  parallelograms over the `t` time scale, or a vector of intervals for
  the `t` time scale.

- y:

  The coordinates for the y-axis. This is a vector of intervals over the
  `s` time scale (default), or a matrix with the corner points of the
  parallelograms over the `s` time scale.

- z:

  The values of the surface to plot, organized in a matrix with
  dimensions compatible with those of `x` and `y`. These can be the SEs
  for the hazard, the SEs for the log-hazard or the SEs for the
  log10-hazard.

- plot_options:

  A list of options for the plot:

  - `loghazard` A Boolean. Default is `FALSE`. If `FALSE` the function
    returns a plot of the standard errors of the hazard surface, if
    `TRUE` the function returns a plot of the standard errors of the
    log-hazard surface.

  - `log10hazard` A Boolean. Default is `FALSE`. If `TRUE`, the function
    returns a plot of the standard errors of the log10-hazard surface

  - `original` A Boolean. Default is `TRUE`. Plot the (log-)hazard in
    the (t,s)-plane. If `FALSE`, the (log-)hazard will be plotted in the
    (u,s)-plane.

  - `rectangular_grid` A Boolean. Default is `FALSE`. If `TRUE`, a
    rectangular grid is used for plotting also in the (t,s)-plane as
    opposed to the grid of parallelograms used as default in the
    (t,s)-plane.

  - `col_palette` A function defining the color palette. The default
    palette is `rev(colorspace::sequential_hcl(n = 50, "Red-Purple"))`.

  - `n_shades` The number of color shades to plot, default is 50.

  - `breaks` The vector of breaks for the color legend. If `n_shades` is
    provided, this should be of length `n_shades + 1`. Otherwise,
    `n_shades` will be recalculated accordingly.

  - `show_legend` A Boolean. Default is `TRUE`. If `FALSE` no legend
    will be plotted, useful for multi-panel figures with common legend.
    Works only for plots on rectangular grid!

  - `tmax` The maximum value of `t` that should be plotted.

  - `main` The title of the plot.

  - `xlab` The label of the first time axis (plotted on the x axis).

  - `ylab` The label of the second time axis (plotted on the y axis).

  - `xlim` A vector with two elements defining the limits of the time
    scale on the x axis.

  - `ylim` A vector with two elements defining the limits of the time
    scale on the y axis.

  - `cex_main` The magnification to be used for the main title, default
    is `1.2`.

  - `cex_lab` The magnification to be used for the axis labels, default
    is `1`.

  - `contour_lines` A Boolean. Default is `FALSE`. If `TRUE` white
    contour lines are added to the surfaces.

  - `contour_col` The color for the contour lines. Default is `white`.

  - `contour_cex` The magnification to be used for the contour lines.
    Default is `.8`.

  - `contour_nlev` The number of contour levels. Default is `10`.

- ...:

  Further arguments to image.plot or image

## Value

An image plot of the SEs for the (log-) hazard surface.
