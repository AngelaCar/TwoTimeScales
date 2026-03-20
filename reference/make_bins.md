# Construct bins over one or more time axes

`make_bins()` constructs the bins over the time axes and saves the
extremes of the bins in a vector.

## Usage

``` r
make_bins(
  t_in = NULL,
  t_out = NULL,
  u = NULL,
  s_in = NULL,
  s_out,
  min_t = NULL,
  max_t = NULL,
  min_u = NULL,
  max_u = NULL,
  min_s = NULL,
  max_s = NULL,
  dt = NULL,
  du = NULL,
  ds
)
```

## Arguments

- t_in:

  (optional) A vector of entry times on the time scale `t`.

- t_out:

  (optional) A vector of exit times on the time scale `t`.

- u:

  (optional) A vector of fixed-times at entry in the process.

- s_in:

  (optional) A vector of entry times on the time scale `s`.

- s_out:

  A vector of exit times on the time scale `s`.

- min_t:

  (optional) A minimum value for the bins over `t`. If `NULL`, the
  minimum of `t_in` will be used.

- max_t:

  (optional) A maximum value for the bins over `t`. If `NULL`, the
  maximum of `t_out` will be used.

- min_u:

  (optional) A minimum value for the bins over `u`. If `NULL`, the
  minimum of `u` will be used.

- max_u:

  (optional) A maximum value for the bins over `u`. If `NULL`, the
  maximum of `u` will be used.

- min_s:

  (optional) A minimum value for the bins over `s`. If `NULL`, the
  minimum of `s_in` will be used.

- max_s:

  (optional) A maximum value for the bins over `s`. If `NULL`, the
  maximum of `s_out` will be used.

- dt:

  (optional) A scalar giving the length of the intervals on the `t` time
  scale.

- du:

  (optional) A scalar giving the length of the intervals on the `u`
  axis.

- ds:

  A scalar giving the length of the intervals on the `s` time scale.

## Value

A list with the following elements:

- `bins_t` if `t_out` is provided, this is a vector of bins extremes for
  the time scale `t`

- `midt` if `t_out` is provided, this is a vector with the midpoints of
  the bins over `t`

- `nt` if `t_out` is provided, this is the number of bins over `t`

- `bins_u` if `u` is provided, this is a vector of bins extremes for `u`
  axis

- `midu` if `u` is provided, this is a vector with the midpoints of the
  bins over `u`

- `nu` if `u` is provided, this is the number of bins over `u`

- `bins_s` is a vector of bins extremes for the time scale `s`

- `mids` is a vector with the midpoints of the bins over `s`

- `ns` is the number of bins over `s`

## Details

It allows construction of bins over the time scales `t` and `s` and/or
over the fixed-time axis `u`. The time scale `s` is always required. See
also
[`prepare_data()`](https://angelacar.github.io/TwoTimeScales/reference/prepare_data.md)
to conveniently prepare individual data for the analysis with one, or
two time scales.

A few words about constructing the grid of bins. There is no 'golden
rule' or optimal strategy for setting the number of bins over each time
axis, or deciding on the bins' width. It very much depends on the data
structure, however, we try to give some directions here. First, in most
cases, more bins is better than less bins. A good number is about 30
bins. However, if data are scarce, the user might want to find a
compromise between having a larger number of bins, and having many bins
empty. Second, the chosen width of the bins (that is `du` and `ds`) does
depend on the time unit over which the time scales are measured. For
example, if the time is recorded in days, as in the example below, and
several years of follow-up are available, the user can split the data in
bins of width 30 (corresponding to about one month), 60 (about two
months), 90 (about three months), etc. If the time scale is measured in
years, then appropriate width could be 0.25 (corresponding to a quarter
of a year), or 0.5 (that is half year). However, in some cases, time
might be measure in completed years, as is often the case for age. In
this scenario, an appropriate bin width is 1.

Finally, it is always a good idea to plot your data first, and explore
the range of values over which the time scale(s) are recorded. This will
give insight about reasonable values for the arguments `min_s`, `min_u`,
`max_s` and `max_u` (that in any case are optional).

## Examples

``` r
# Make bins for colon cancer data by time at randomization and time since recurrence
bins <- make_bins(u = reccolon2ts$timer, s_out = reccolon2ts$timesr,
                 du = 30, ds = 30)
#> `s_in = NULL`. I will use `s_in = 0` for all observations.
# Make bins for colon cancer data only over time since recurrence
bins <- make_bins(s_out = reccolon2ts$timesr, ds = 60)
#> `s_in = NULL`. I will use `s_in = 0` for all observations.
```
