# Bin data on the Lexis diagram

`exposures_events_Lexis()` computes aggregated matrices of exposure
times and event counts over two time scales, on the Lexis diagram.

The time scales are `t` and `s`. This function uses functions from the
package `popEpi` and from the package `Epi`, and code shared by Bendix
Carstensen on the website bendixcarstensen.com. See also
[`prepare_data()`](https://angelacar.github.io/TwoTimeScales/reference/prepare_data.md)
to conveniently prepare individual data for the analysis with one, or
two time scales.

## Usage

``` r
exposures_events_Lexis(t_in = NULL, t_out, s_in = NULL, s_out, ev, bins_list)
```

## Arguments

- t_in:

  (optional) A vector of entry times on the time scale `t`.

- t_out:

  (optional) A vector of exit times on the time scale `t`.

- s_in:

  (optional) A vector of entry times on the time scale `s`.

- s_out:

  A vector of exit times on the time scale `s`.

- ev:

  A vector of event indicators (possible values 0/1).

- bins_list:

  A list with the following (necessary) elements:

  - `bins_t` a vector of extreme values for the bins over the `t` axis.

  - `nt` the number of bins over `t`.

  - `bins_s` a vector of extreme values for the bins over the `t` axis.

  - `ns` the number of bins over `s`.

## Value

A list with the following elements:

- `R` an array of exposure times of dimension \\n_t\\ by \\n_s\\

- `Y` an array of event counts of dimension \\n_t\\ by \\n_s\\

## References

Carstensen B, Plummer M, Laara E, Hills M (2022). Epi: A Package for
Statistical Analysis in Epidemiology. R package version 2.47.1,
https://CRAN.R-project.org/package=Epi.

Miettinen J, Rantanen M, Seppa K (2023). popEpi: Functions for
Epidemiological Analysis using Population Data. R package version
0.4.11, https://cran.r-project.org/package=popEpi.

## Author

Angela Carollo <carollo@demogr.mpg.de>

## Examples

``` r
# ---- Bin colon cancer data by time since randomization and time since recurrence ----
# First create vectors of bins (using function `make_bins()`)
bins <- make_bins(t_out = reccolon2ts$timedc, s_out = reccolon2ts$timesr,
dt = 90, ds = 90)
#> `t_in` not provided. I will use `t_in = t_out - s_in`.
#> `s_in = NULL`. I will use `s_in = 0` for all observations.
# Now bin data (note: the t_in and s_in arguments are omitted because data are not left truncated)
bindata2d <- exposures_events_Lexis(t_out = reccolon2ts$timedc,
s_out = reccolon2ts$timesr, ev = reccolon2ts$status, bins = bins)
#> `t_in` not provided. I will use `t_in = t_out - s_in`.
#> `s_in = NULL`. I will use `s_in = 0` for all observations.
#> NOTE: entry.status has been set to 0 for all.
```
