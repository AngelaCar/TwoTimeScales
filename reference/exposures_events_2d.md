# Bin data on two time scales

`exposures_events_2d()` computes individual or aggregated matrices of
exposure times and event counts starting from individual records of time
at entry in the process (measured over the first time scale), duration
at entry in the process (measured over the second time scale), duration
at exit from the process (measured over the second time scale), and
event's indicator.

## Usage

``` r
exposures_events_2d(u, s_in = NULL, s_out, ev, bins_list, individual = FALSE)
```

## Arguments

- u:

  A vector of fixed times at entry in the process, measured over the
  first time scale.

- s_in:

  A vector of (possibly left truncated) times at entry. If this is not
  provided by the user, the function will consider a value of 0 for all
  observations.

- s_out:

  A vector of times at event or censoring.

- ev:

  A vector of events' indicators (possible values 0/1).

- bins_list:

  is a list with the following (necessary) elements (usually prepared by
  [`make_bins()`](https://angelacar.github.io/TwoTimeScales/reference/make_bins.md)):

  - `bins_u` a vector of extreme values for the bins over the `u` axis

  - `bins_s` a vector of extreme values for the bins over the `s` axis

- individual:

  A Boolean. Default is `FALSE`: if `FALSE` computes the matrices R and
  Y collectively for all observations; if `TRUE` computes the matrices R
  and Y separately for each individual record.

## Value

A list with the following elements:

- `R` an array of exposure times: if `individual == TRUE`, then `R` is
  an array of dimension \\n_u\\ by \\n_s\\ by \\n\\, otherwise is an
  array of dimension \\n_u\\ by \\n_s\\

- `Y`an array of event counts: if `individual == TRUE`, then `Y` is an
  array of dimension \\n_u\\ by \\n_s\\ by \\n\\, otherwise is an array
  of dimension \\n_u\\ by \\n_s\\

## Details

The fixed-time variable `u` and the second time scale `s` are divided
into \\n_u\\ and \\n_s\\ intervals, respectively. The extremes of these
intervals are provided as input to the function. First, the fixed-time
at entry is located in one of the \\n_u\\ bins that cover the whole
range of `u`. Then, the time-at-risk for each individual is split
according to the \\n_s\\ bins that span the whole range of values for
`s`, and an event indicator is placed in the bin where the exit time is
located. This is done by calling the function `exposure_events_1d()`. If
individual matrices of exposure and events are required, then the
function returns two arrays of dimension \\n_u\\ by \\n_s\\ by \\n\\. If
aggregated results are preferred, the individual contributions are
summed in each bin to provide a matrix of total exposure times and a
matrix of total event counts, both of dimensions \\n_u\\ by \\n_s\\. See
also
[`prepare_data()`](https://angelacar.github.io/TwoTimeScales/reference/prepare_data.md)
to conveniently prepare individual data for the analysis with one, or
two time scales.

## Author

Angela Carollo <carollo@demogr.mpg.de>

## Examples

``` r
# ---- Bin colon cancer data by time at randomization and time since recurrence ----
# First create vectors of bins (using function `make_bins()`)
bins <- make_bins(u = reccolon2ts$timer, s_out = reccolon2ts$timesr,
du = 30, ds = 30)
#> `s_in = NULL`. I will use `s_in = 0` for all observations.
# Now bin data (note: the s_in argument is omitted because data are not left truncated)
bindata2d <- exposures_events_2d(u = reccolon2ts$timer,
s_out = reccolon2ts$timesr, ev = reccolon2ts$status, bins = bins)
#> `s_in = NULL`. I will use `s_in = 0` for all observations.
```
