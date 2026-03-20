# Plot slices of the (log-) hazard

`plot_slices()` plots slices of the (log-)hazard with two time scales,
at selected values of one of the two time dimensions.

## Usage

``` r
plot_slices(x, y, direction, plot_options = list())
```

## Arguments

- x:

  A vector of values for the x-axis. This is a vector of values over the
  axis opposite to the one where the sliced are cut.

- y:

  A matrix of (log-)hazard values.

- direction:

  Either `"u"` or `"s"`.

- plot_options:

  A list of options for the plot:

  - `loghazard` A Boolean. Default is `FALSE`. If `FALSE` the function
    returns a plot of cross-sections from the hazard surface, if `TRUE`
    the function returns a plot of cross-sections from the log-hazard
    surface.

  - `log10hazard` A Boolean. Default is `FALSE`. If `TRUE` returns a
    plot of cross-sections from the log10-hazard surface.

  - `col_palette` A function defining the color palette. The default
    palette is
    [`grDevices::gray.colors()`](https://rdrr.io/r/grDevices/gray.colors.html).

  - `n_shades` The number of color shades to plot, default is 50.

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

  - `lwd` The line width. Default is `2`.

  - `lty` The line type.

## Value

A plot of the slices of the hazard cut at selected points.
