# Prepare raw data by binning them in 1d or 2d

`prepare_data()` prepares the raw individual time-to-event data for
hazard estimation in 1d or 2d.

Given the raw data, this function first constructs the bins over one or
two time axes and then computes the aggregated (or individual) vectors,
or matrices, of exposure times and events indicators.

## Usage

``` r
prepare_data(
  data = NULL,
  t_in = NULL,
  t_out = NULL,
  u = NULL,
  s_in = NULL,
  s_out,
  events,
  min_t = NULL,
  max_t = NULL,
  min_u = NULL,
  max_u = NULL,
  min_s = NULL,
  max_s = NULL,
  dt = NULL,
  du = NULL,
  ds,
  individual = FALSE,
  covs = NULL
)
```

## Arguments

- data:

  A data frame. If the object provided is a tibble, or a data.table, it
  will be converted to data.frame.

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

- events:

  A vector of event's indicators (possible values 0/1).

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

- individual:

  A Boolean. Default is `FALSE`: if `FALSE` computes the matrices `R`
  and `Y` collectively for all observations; if `TRUE` computes the
  matrices `R` and `Y` separately for each individual record.

- covs:

  A data.frame with the variables to be used as covariates, or a vector
  with the names of the covariates to be included. The function will
  create dummy variables for any factor variable passed as argument in
  `covs`. If a variable of class character is passed as argument, it
  will be converted to factor.

## Value

A list with the following elements:

- `bins` a list:

  - `bins_t` if `t_out` is provided, this is a vector of bins extremes
    for the time scale `t`.

  - `mid_t` if `t_out` is provided, this is a vector with the midpoints
    of the bins over `t`.

  - `nt` if `t_out` is provided, this is the number of bins over `t`.

  - `bins_u` if `u` is provided, this is a vector of bins extremes for
    `u` axis.

  - `midu` if `u` is provided, this is a vector with the midpoints of
    the bins over `u`.

  - `nu` if `u` is provided, this is the number of bins over `u`.

  - `bins_s` is a vector of bins extremes for the time scale `s`.

  - `mids` is a vector with the midpoints of the bins over `s`.

  - `ns` is the number of bins over `s`.

- `bindata`:

  - `r` or `R` an array of exposure times: if binning the data over one
    time scale only this is a vector. If binning the data over two time
    scales and if `individual == TRUE` then `R` is an array of dimension
    nu by ns by n, otherwise it is an array of dimension nu by ns

  - `y` or `Y` an array of event counts: if binning the data over one
    time scale only this is a vector. If binning the data over two time
    scales and if `individual == TRUE` then `Y` is an array of dimension
    nu by ns by n, otherwise it is an array of dimension nu by ns

  - `Z` A matrix of covariates' values to be used in the model, of
    dimension n by p

## Details

A few words about constructing the grid of bins. Bins are containers for
the individual data. There is no 'golden rule' or optimal strategy for
setting the number of bins over each time axis, or deciding on the bins'
width. It very much depends on the data structure, however, we try to
give some directions here. First, in most cases, more bins is better
than less bins. A good number is about 30 bins. It also depends on the
user's expectation about how rapidly the hazard of the event changes
over the time scales. However, if data are scarce, the user might want
to find a compromise between having a larger number of bins, and having
many bins empty. Second, the chosen width of the bins (that is `du` and
`ds`) does depend on the time unit over which the time scales are
measured. For example, if the time is recorded in days, as in the
example below, and several years of follow-up are available, the user
can split the data in bins of width 30 (corresponding to about one
month), 60 (about two months), 90 (about three months), etc. If the time
scale is measured in years, then appropriate width could be 0.25
(corresponding to a quarter of a year), or 0.5 (that is half year).
However, in some cases, time might be measure in completed years, as is
often the case for age. In this scenario, an appropriate bin width is 1.

Finally, it is always a good idea to plot the data first, and explore
the range of values over which the time scale(s) are recorded. This will
give insight about reasonable values for the arguments `min_s`, `min_u`,
`max_s` and `max_u` (that in any case are optional).

Regarding the names of covariates or levels of categorical
covariates/factors: When using "LMMsolver" to fit a model with
covariates that have names (or factor labels) including a symbol such as
"+", "-", "\<" or "\>", these will result in an error. To avoid this,
the responsible names (labels) will be rewritten without mathematical
symbols. For example: "Lev+5FU" (in the colon cancer data) is replaced
by "Lev&5FU".

## Author

Angela Carollo <carollo@demogr.mpg.de>

## Examples

``` r
# Bin data over s = time since recurrence only, with intervals of length 30 days
# aggregated data (no covariates)
# The following example provide the vectors of data directly from the dataset
binned_data <- prepare_data(s_out = reccolon2ts$timesr, events = reccolon2ts$status, ds = 30)
#> `s_in = NULL`. I will use `s_in = 0` for all observations.
# Visualize vector of event counts
print(binned_data$bindata$y)
#>  [1] 14 24 16 11 24 22 20 10 24 26 13 19 10  9 11 11 10 10  9 13 10  7 10  7  4
#> [26]  4  7  5  8  4  1  2  5  5  5  0  0  1  1  1  1  1  1  0  0  0  1  2  1  1
#> [51]  0  1  0  0  0  2  1  0  1  0  1  0  0  0  1  0  0  0  0  0  0  0  0  0  0
#> [76]  0  0  1  0  0  0  0  0  0  0  0  0  0  0  0  0
# Visualize midpoints of the bins
print(binned_data$bins$mids)
#>  [1]   15   45   75  105  135  165  195  225  255  285  315  345  375  405  435
#> [16]  465  495  525  555  585  615  645  675  705  735  765  795  825  855  885
#> [31]  915  945  975 1005 1035 1065 1095 1125 1155 1185 1215 1245 1275 1305 1335
#> [46] 1365 1395 1425 1455 1485 1515 1545 1575 1605 1635 1665 1695 1725 1755 1785
#> [61] 1815 1845 1875 1905 1935 1965 1995 2025 2055 2085 2115 2145 2175 2205 2235
#> [76] 2265 2295 2325 2355 2385 2415 2445 2475 2505 2535 2565 2595 2625 2655 2685
#> [91] 2715
# Visualize number of bins
print(binned_data$bins$ns)
#> [1] 91

# Now, the same thing is done by providing a dataset and the name of all relevant variables
binned_data <- prepare_data(data = reccolon2ts, s_out = "timesr", events = "status", ds = 30)
#> `s_in = NULL`. I will use `s_in = 0` for all observations.
# Visualize vector of event counts
print(binned_data$bindata$y)
#>  [1] 14 24 16 11 24 22 20 10 24 26 13 19 10  9 11 11 10 10  9 13 10  7 10  7  4
#> [26]  4  7  5  8  4  1  2  5  5  5  0  0  1  1  1  1  1  1  0  0  0  1  2  1  1
#> [51]  0  1  0  0  0  2  1  0  1  0  1  0  0  0  1  0  0  0  0  0  0  0  0  0  0
#> [76]  0  0  1  0  0  0  0  0  0  0  0  0  0  0  0  0

# Now using ds = .3 and the same variable measured in years
binned_data <- prepare_data(s_out = reccolon2ts$timesr_y, events = reccolon2ts$status, ds = .3)
#> `s_in = NULL`. I will use `s_in = 0` for all observations.
# Visualize vector of exposure timess
print(binned_data$bindata$r)
#>  [1] 128.9234086 108.6935661  86.6990418  67.7699521  55.8161533  43.9006160
#>  [7]  33.6388090  26.2145791  20.8121150  16.9198494  14.2230664  12.1813142
#> [13]  11.2851472  10.0655715   9.0108145   7.7641342   6.0433949   4.7993155
#> [19]   2.9854209   2.0730322   1.2588638   0.9572895   0.8856947   0.3787817
#> [25]   0.2606434


# Bin data over u = time at recurrence and s = time since recurrence, measured in days
# aggregated data (no covariates)
# Note that if we do not provide du this is taken to be equal to ds
binned_data <- prepare_data(
  u = reccolon2ts$timer, s_out = reccolon2ts$timesr,
  events = reccolon2ts$status, ds = 30
)
#> `s_in = NULL`. I will use `s_in = 0` for all observations.
#> `s_in = NULL`. I will use `s_in = 0` for all observations.

# Visualize matrix of event counts
print(binned_data$bindata$Y)
#>       [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8] [,9] [,10] [,11] [,12] [,13]
#>  [1,]    1    0    0    0    1    0    0    0    0     0     0     0     0
#>  [2,]    1    2    3    0    0    0    0    0    1     0     0     0     0
#>  [3,]    0    5    1    0    2    2    1    0    2     2     0     2     0
#>  [4,]    1    4    0    1    1    3    2    1    3     0     1     0     0
#>  [5,]    1    1    2    1    1    1    1    0    0     0     0     0     0
#>  [6,]    2    1    2    1    1    2    1    0    2     4     2     0     0
#>  [7,]    1    0    1    1    1    1    2    0    2     2     1     0     0
#>  [8,]    1    3    0    0    1    2    3    0    1     1     0     5     2
#>  [9,]    1    1    1    0    1    2    0    0    0     0     2     1     0
#> [10,]    2    0    0    1    2    0    1    1    2     2     2     1     1
#> [11,]    0    2    1    0    0    1    0    2    0     2     0     1     1
#> [12,]    0    0    0    1    1    0    2    0    1     1     1     0     1
#> [13,]    0    1    1    0    4    0    1    0    0     0     0     0     0
#> [14,]    0    0    0    0    0    1    0    0    1     0     1     1     0
#> [15,]    0    1    0    0    0    0    1    0    1     1     0     2     0
#> [16,]    0    0    0    1    1    0    0    1    1     0     0     1     0
#> [17,]    0    0    0    0    1    1    2    0    1     0     1     0     2
#> [18,]    0    0    0    0    3    1    0    2    0     1     0     1     0
#> [19,]    0    0    1    0    0    1    1    0    0     1     0     1     0
#> [20,]    0    0    1    1    0    0    0    0    1     1     1     0     0
#> [21,]    0    1    0    0    0    0    0    0    0     0     0     0     1
#> [22,]    0    0    0    0    1    0    0    0    1     1     0     0     1
#> [23,]    0    0    0    0    0    0    0    0    1     1     0     0     0
#> [24,]    0    1    0    0    1    0    0    0    1     0     0     0     0
#> [25,]    0    0    0    0    0    1    0    1    0     2     0     0     0
#> [26,]    0    0    0    0    0    0    0    0    0     0     0     0     0
#> [27,]    0    0    1    0    0    0    1    0    0     0     0     0     0
#> [28,]    1    0    0    0    0    0    0    0    0     0     0     0     0
#> [29,]    0    0    0    0    0    1    0    0    0     1     0     0     0
#> [30,]    0    1    0    0    0    0    0    1    0     0     1     0     0
#> [31,]    1    0    0    0    0    1    0    0    0     1     0     0     0
#> [32,]    0    0    0    1    0    0    0    0    0     0     0     0     0
#> [33,]    0    0    0    0    0    0    0    0    0     0     0     0     0
#> [34,]    0    0    0    0    0    0    1    0    2     0     0     0     0
#> [35,]    0    0    0    0    0    0    0    0    0     1     0     0     0
#> [36,]    0    0    0    0    0    0    0    0    0     0     0     0     0
#> [37,]    0    0    0    0    0    1    0    0    0     0     0     0     0
#> [38,]    0    0    0    0    0    0    0    0    0     0     0     0     0
#> [39,]    0    0    0    0    0    0    0    0    0     0     0     1     0
#> [40,]    0    0    0    0    0    0    0    0    0     0     0     0     0
#> [41,]    0    0    0    0    0    0    0    0    0     0     0     0     1
#> [42,]    0    0    0    0    0    0    0    0    0     0     0     0     0
#> [43,]    0    0    0    0    0    0    0    1    0     0     0     0     0
#> [44,]    0    0    0    1    0    0    0    0    0     0     0     0     0
#> [45,]    0    0    0    0    0    0    0    0    0     0     0     1     0
#> [46,]    0    0    0    0    0    0    0    0    0     0     0     0     0
#> [47,]    0    0    0    0    0    0    0    0    0     0     0     0     0
#> [48,]    0    0    0    0    0    0    0    0    0     0     0     0     0
#> [49,]    0    0    0    0    0    0    0    0    0     0     0     0     0
#> [50,]    0    0    0    0    0    0    0    0    0     0     0     0     0
#> [51,]    0    0    0    0    0    0    0    0    0     0     0     0     0
#> [52,]    0    0    0    0    1    0    0    0    0     0     0     0     0
#> [53,]    0    0    0    0    0    0    0    0    0     0     0     0     0
#> [54,]    0    0    0    0    0    0    0    0    0     1     0     0     0
#> [55,]    0    0    0    0    0    0    0    0    0     0     0     0     0
#> [56,]    0    0    0    0    0    0    0    0    0     0     0     0     0
#> [57,]    0    0    0    0    0    0    0    0    0     0     0     0     0
#> [58,]    0    0    1    0    0    0    0    0    0     0     0     0     0
#> [59,]    1    0    0    0    0    0    0    0    0     0     0     0     0
#> [60,]    0    0    0    0    0    0    0    0    0     0     0     1     0
#> [61,]    0    0    0    0    0    0    0    0    0     0     0     0     0
#> [62,]    0    0    0    0    0    0    0    0    0     0     0     0     0
#> [63,]    0    0    0    0    0    0    0    0    0     0     0     0     0
#> [64,]    0    0    0    0    0    0    0    0    0     0     0     0     0
#> [65,]    0    0    0    0    0    0    0    0    0     0     0     0     0
#> [66,]    0    0    0    1    0    0    0    0    0     0     0     0     0
#> [67,]    0    0    0    0    0    0    0    0    0     0     0     0     0
#> [68,]    0    0    0    0    0    0    0    0    0     0     0     0     0
#> [69,]    0    0    0    0    0    0    0    0    0     0     0     0     0
#> [70,]    0    0    0    0    0    0    0    0    0     0     0     0     0
#> [71,]    0    0    0    0    0    0    0    0    0     0     0     0     0
#> [72,]    0    0    0    0    0    0    0    0    0     0     0     0     0
#> [73,]    0    0    0    0    0    0    0    0    0     0     0     0     0
#> [74,]    0    0    0    0    0    0    0    0    0     0     0     0     0
#> [75,]    0    0    0    0    0    0    0    0    0     0     0     0     0
#> [76,]    0    0    0    0    0    0    0    0    0     0     0     0     0
#>       [,14] [,15] [,16] [,17] [,18] [,19] [,20] [,21] [,22] [,23] [,24] [,25]
#>  [1,]     2     0     0     0     0     0     0     1     0     0     0     0
#>  [2,]     1     0     0     0     1     0     0     0     0     0     0     0
#>  [3,]     0     0     0     0     0     0     0     0     0     1     0     0
#>  [4,]     1     0     0     0     0     1     1     1     0     0     0     0
#>  [5,]     0     1     0     1     0     1     1     1     0     0     1     0
#>  [6,]     1     1     0     1     1     1     1     0     0     0     1     0
#>  [7,]     0     0     1     1     0     0     0     0     0     0     0     0
#>  [8,]     0     0     0     0     0     0     0     0     0     2     0     0
#>  [9,]     0     1     1     2     0     0     0     1     1     1     0     0
#> [10,]     0     0     1     0     0     0     0     0     0     1     0     0
#> [11,]     0     0     2     0     0     1     1     0     0     0     1     0
#> [12,]     0     0     0     0     0     0     1     1     0     0     0     0
#> [13,]     1     1     0     0     0     1     0     0     0     0     0     0
#> [14,]     0     1     0     1     0     1     0     0     0     0     1     0
#> [15,]     0     0     0     0     0     0     1     0     0     0     0     0
#> [16,]     0     0     0     0     1     0     0     0     0     1     0     1
#> [17,]     1     0     0     1     1     0     0     0     0     1     0     0
#> [18,]     0     1     0     0     0     0     1     0     0     0     0     0
#> [19,]     0     0     1     0     0     0     1     1     0     0     0     1
#> [20,]     0     0     0     0     0     0     0     1     1     0     0     1
#> [21,]     0     0     1     0     3     0     0     0     1     0     0     0
#> [22,]     0     0     1     1     0     0     0     0     0     0     0     0
#> [23,]     0     0     0     0     2     0     1     1     0     0     0     0
#> [24,]     0     1     0     0     0     1     1     0     0     0     0     0
#> [25,]     0     0     1     0     0     0     1     0     0     0     0     0
#> [26,]     1     1     0     0     0     0     0     0     0     0     0     0
#> [27,]     0     0     0     0     0     0     0     0     0     0     0     0
#> [28,]     0     1     0     0     0     0     0     0     0     0     1     0
#> [29,]     0     0     0     0     0     0     1     0     1     0     0     0
#> [30,]     0     0     0     0     0     0     0     0     0     0     0     0
#> [31,]     0     1     1     0     0     0     0     0     0     0     0     0
#> [32,]     0     0     0     0     0     0     0     0     0     0     0     0
#> [33,]     0     0     0     0     0     1     0     0     0     0     0     0
#> [34,]     0     0     0     0     0     0     0     0     1     0     0     0
#> [35,]     0     0     0     0     0     0     1     0     0     0     0     1
#> [36,]     0     0     0     0     0     0     0     0     0     0     0     0
#> [37,]     0     0     0     0     0     0     0     0     1     1     0     0
#> [38,]     0     0     0     0     0     0     0     0     0     1     0     0
#> [39,]     0     0     0     0     0     0     0     0     0     0     0     0
#> [40,]     0     0     0     0     0     0     0     0     0     0     0     0
#> [41,]     0     0     0     0     0     0     0     0     0     0     0     0
#> [42,]     0     0     0     0     0     0     0     0     0     0     0     0
#> [43,]     0     0     0     1     0     0     0     0     0     0     0     0
#> [44,]     0     0     0     0     0     0     0     0     0     0     0     0
#> [45,]     0     0     0     0     0     0     0     0     0     0     0     0
#> [46,]     0     0     0     0     0     0     0     0     0     0     0     0
#> [47,]     0     0     0     0     0     0     0     0     0     0     0     0
#> [48,]     0     1     0     0     0     0     0     0     1     0     0     0
#> [49,]     1     0     0     0     0     0     0     0     0     0     0     0
#> [50,]     0     0     0     1     0     0     0     0     0     0     0     0
#> [51,]     0     0     0     0     0     0     0     0     0     0     0     0
#> [52,]     0     0     0     0     0     0     0     1     0     0     0     0
#> [53,]     0     0     0     0     0     0     0     0     0     0     0     0
#> [54,]     0     0     0     0     1     0     0     0     0     0     0     0
#> [55,]     0     0     0     0     0     1     0     0     0     0     0     0
#> [56,]     0     0     0     0     0     0     0     0     0     0     0     0
#> [57,]     0     0     0     0     0     0     0     0     0     0     0     0
#> [58,]     0     0     0     0     0     0     0     0     0     0     0     0
#> [59,]     0     0     0     0     0     0     0     0     0     0     0     0
#> [60,]     0     0     0     0     0     0     0     0     0     0     0     0
#> [61,]     0     0     0     0     0     0     0     0     0     0     0     0
#> [62,]     0     0     0     0     0     0     0     0     0     0     0     0
#> [63,]     0     0     0     0     0     0     0     1     0     0     0     0
#> [64,]     0     0     0     0     0     0     0     0     0     0     0     0
#> [65,]     0     0     0     0     0     0     0     0     0     0     0     0
#> [66,]     0     0     0     0     0     0     0     0     0     0     1     0
#> [67,]     0     0     0     0     0     0     0     0     0     0     0     0
#> [68,]     0     0     0     0     0     0     0     0     0     0     1     0
#> [69,]     0     0     1     0     0     0     0     0     0     0     0     0
#> [70,]     0     0     0     0     0     0     0     0     0     0     0     0
#> [71,]     0     0     0     0     0     0     0     0     0     0     0     0
#> [72,]     0     0     0     0     0     0     0     0     0     0     0     0
#> [73,]     0     0     0     0     0     0     0     0     0     0     0     0
#> [74,]     0     0     0     0     0     0     0     0     0     0     0     0
#> [75,]     0     0     0     0     0     0     0     0     0     1     0     0
#> [76,]     0     0     0     0     0     0     0     0     0     0     0     0
#>       [,26] [,27] [,28] [,29] [,30] [,31] [,32] [,33] [,34] [,35] [,36] [,37]
#>  [1,]     1     0     0     0     0     0     0     0     0     0     0     0
#>  [2,]     0     0     0     0     0     0     0     0     0     0     0     0
#>  [3,]     0     0     0     0     0     0     0     0     0     0     0     0
#>  [4,]     0     1     0     2     0     0     0     1     0     0     0     0
#>  [5,]     0     1     0     0     0     0     0     0     0     0     0     0
#>  [6,]     0     0     0     1     0     0     0     0     0     0     0     0
#>  [7,]     1     0     0     0     1     0     0     0     1     0     0     0
#>  [8,]     0     0     0     0     1     0     0     0     0     1     0     0
#>  [9,]     0     2     0     0     0     0     0     1     0     1     0     0
#> [10,]     0     0     0     0     0     0     0     0     0     0     0     0
#> [11,]     0     0     0     1     0     0     0     0     0     0     0     0
#> [12,]     0     0     1     1     0     0     0     1     0     0     0     0
#> [13,]     0     0     0     1     0     0     0     0     0     1     0     0
#> [14,]     0     0     0     0     0     0     0     0     0     1     0     0
#> [15,]     0     0     1     1     0     0     0     0     1     0     0     0
#> [16,]     0     0     0     0     0     1     0     1     0     0     0     0
#> [17,]     0     0     0     0     0     0     1     0     0     0     0     0
#> [18,]     0     0     1     0     0     0     0     0     1     0     0     0
#> [19,]     0     0     0     0     0     0     0     0     1     0     0     0
#> [20,]     0     1     0     0     0     0     0     0     0     0     0     0
#> [21,]     0     0     1     0     0     0     0     0     0     0     0     0
#> [22,]     0     0     0     0     0     0     0     0     0     0     0     0
#> [23,]     0     0     0     0     0     0     0     0     0     0     0     0
#> [24,]     1     0     0     0     0     0     0     0     0     0     0     0
#> [25,]     0     0     0     0     0     0     0     0     0     0     0     0
#> [26,]     0     0     0     0     0     0     0     0     0     0     0     0
#> [27,]     0     0     1     0     0     0     0     0     0     0     0     0
#> [28,]     0     0     0     0     0     0     0     0     0     0     0     0
#> [29,]     0     0     0     0     0     0     0     0     0     1     0     0
#> [30,]     0     0     0     1     0     0     0     0     0     0     0     0
#> [31,]     0     0     0     0     0     0     0     1     0     0     0     0
#> [32,]     0     0     0     0     0     0     0     0     0     0     0     0
#> [33,]     0     0     0     0     0     0     0     0     0     0     0     0
#> [34,]     0     0     0     0     0     0     0     0     0     0     0     0
#> [35,]     0     1     0     0     0     0     0     0     0     0     0     0
#> [36,]     0     0     0     0     0     0     0     0     0     0     0     0
#> [37,]     0     0     0     0     0     0     0     0     0     0     0     0
#> [38,]     0     0     0     0     1     0     1     0     0     0     0     0
#> [39,]     0     0     0     0     0     0     0     0     0     0     0     0
#> [40,]     0     0     0     0     0     0     0     0     0     0     0     0
#> [41,]     0     0     0     0     1     0     0     0     0     0     0     0
#> [42,]     0     0     0     0     0     0     0     0     0     0     0     0
#> [43,]     0     0     0     0     0     0     0     0     1     0     0     0
#> [44,]     0     0     0     0     0     0     0     0     0     0     0     0
#> [45,]     0     0     0     0     0     0     0     0     0     0     0     0
#> [46,]     0     0     0     0     0     0     0     0     0     0     0     0
#> [47,]     0     0     0     0     0     0     0     0     0     0     0     0
#> [48,]     0     0     0     0     0     0     0     0     0     0     0     0
#> [49,]     0     1     0     0     0     0     0     0     0     0     0     0
#> [50,]     0     0     0     0     0     0     0     0     0     0     0     0
#> [51,]     0     0     0     0     0     0     0     0     0     0     0     0
#> [52,]     0     0     0     0     0     0     0     0     0     0     0     0
#> [53,]     0     0     0     0     0     0     0     0     0     0     0     0
#> [54,]     0     0     0     0     0     0     0     0     0     0     0     0
#> [55,]     0     0     0     0     0     0     0     0     0     0     0     0
#> [56,]     0     0     0     0     0     0     0     0     0     0     0     0
#> [57,]     0     0     0     0     0     0     0     0     0     0     0     0
#> [58,]     0     0     0     0     0     0     0     0     0     0     0     0
#> [59,]     1     0     0     0     0     0     0     0     0     0     0     0
#> [60,]     0     0     0     0     0     0     0     0     0     0     0     0
#> [61,]     0     0     0     0     0     0     0     0     0     0     0     0
#> [62,]     0     0     0     0     0     0     0     0     0     0     0     0
#> [63,]     0     0     0     0     0     0     0     0     0     0     0     0
#> [64,]     0     0     0     0     0     0     0     0     0     0     0     0
#> [65,]     0     0     0     0     0     0     0     0     0     0     0     0
#> [66,]     0     0     0     0     0     0     0     0     0     0     0     0
#> [67,]     0     0     0     0     0     0     0     0     0     0     0     0
#> [68,]     0     0     0     0     0     0     0     0     0     0     0     0
#> [69,]     0     0     0     0     0     0     0     0     0     0     0     0
#> [70,]     0     0     0     0     0     0     0     0     0     0     0     0
#> [71,]     0     0     0     0     0     0     0     0     0     0     0     0
#> [72,]     0     0     0     0     0     0     0     0     0     0     0     0
#> [73,]     0     0     0     0     0     0     0     0     0     0     0     0
#> [74,]     0     0     0     0     0     0     0     0     0     0     0     0
#> [75,]     0     0     0     0     0     0     0     0     0     0     0     0
#> [76,]     0     0     0     0     0     0     0     0     0     0     0     0
#>       [,38] [,39] [,40] [,41] [,42] [,43] [,44] [,45] [,46] [,47] [,48] [,49]
#>  [1,]     0     0     0     0     0     0     0     0     0     0     0     0
#>  [2,]     0     0     0     0     0     0     0     0     0     0     0     0
#>  [3,]     0     0     0     0     0     0     0     0     0     0     0     0
#>  [4,]     0     0     0     0     0     0     0     0     0     0     0     0
#>  [5,]     0     0     0     0     0     0     0     0     0     0     0     0
#>  [6,]     0     0     0     0     0     0     0     0     0     0     0     0
#>  [7,]     0     0     0     0     0     0     0     0     0     0     0     0
#>  [8,]     0     0     0     0     0     0     0     0     0     0     0     0
#>  [9,]     0     0     0     0     0     0     0     0     0     0     0     0
#> [10,]     0     0     0     0     0     0     0     0     0     0     0     0
#> [11,]     0     0     0     0     0     0     0     0     0     0     0     0
#> [12,]     0     0     0     0     0     0     0     0     0     0     0     0
#> [13,]     0     0     1     0     0     0     0     0     0     1     0     0
#> [14,]     0     0     0     0     0     0     0     0     0     0     0     0
#> [15,]     0     0     0     0     0     0     0     0     0     0     1     0
#> [16,]     0     0     0     0     0     1     0     0     0     0     0     0
#> [17,]     0     0     0     0     0     0     0     0     0     0     0     1
#> [18,]     1     0     0     0     0     0     0     0     0     0     0     0
#> [19,]     0     0     0     0     0     0     0     0     0     0     0     0
#> [20,]     0     0     0     0     0     0     0     0     0     0     0     0
#> [21,]     0     0     0     0     1     0     0     0     0     0     0     0
#> [22,]     0     0     0     1     0     0     0     0     0     0     0     0
#> [23,]     0     0     0     0     0     0     0     0     0     0     0     0
#> [24,]     0     0     0     0     0     0     0     0     0     0     0     0
#> [25,]     0     1     0     0     0     0     0     0     0     0     0     0
#> [26,]     0     0     0     0     0     0     0     0     0     0     0     0
#> [27,]     0     0     0     0     0     0     0     0     0     0     0     0
#> [28,]     0     0     0     0     0     0     0     0     0     0     0     0
#> [29,]     0     0     0     0     0     0     0     0     0     0     0     0
#> [30,]     0     0     0     0     0     0     0     0     0     0     0     0
#> [31,]     0     0     0     0     0     0     0     0     0     0     1     0
#> [32,]     0     0     0     0     0     0     0     0     0     0     0     0
#> [33,]     0     0     0     0     0     0     0     0     0     0     0     0
#> [34,]     0     0     0     0     0     0     0     0     0     0     0     0
#> [35,]     0     0     0     0     0     0     0     0     0     0     0     0
#> [36,]     0     0     0     0     0     0     0     0     0     0     0     0
#> [37,]     0     0     0     0     0     0     0     0     0     0     0     0
#> [38,]     0     0     0     0     0     0     0     0     0     0     0     0
#> [39,]     0     0     0     0     0     0     0     0     0     0     0     0
#> [40,]     0     0     0     0     0     0     0     0     0     0     0     0
#> [41,]     0     0     0     0     0     0     0     0     0     0     0     0
#> [42,]     0     0     0     0     0     0     0     0     0     0     0     0
#> [43,]     0     0     0     0     0     0     0     0     0     0     0     0
#> [44,]     0     0     0     0     0     0     0     0     0     0     0     0
#> [45,]     0     0     0     0     0     0     0     0     0     0     0     0
#> [46,]     0     0     0     0     0     0     0     0     0     0     0     0
#> [47,]     0     0     0     0     0     0     0     0     0     0     0     0
#> [48,]     0     0     0     0     0     0     0     0     0     0     0     0
#> [49,]     0     0     0     0     0     0     0     0     0     0     0     0
#> [50,]     0     0     0     0     0     0     0     0     0     0     0     0
#> [51,]     0     0     0     0     0     0     0     0     0     0     0     0
#> [52,]     0     0     0     0     0     0     0     0     0     0     0     0
#> [53,]     0     0     0     0     0     0     0     0     0     0     0     0
#> [54,]     0     0     0     0     0     0     0     0     0     0     0     0
#> [55,]     0     0     0     0     0     0     0     0     0     0     0     0
#> [56,]     0     0     0     0     0     0     0     0     0     0     0     0
#> [57,]     0     0     0     0     0     0     0     0     0     0     0     0
#> [58,]     0     0     0     0     0     0     0     0     0     0     0     0
#> [59,]     0     0     0     0     0     0     0     0     0     0     0     0
#> [60,]     0     0     0     0     0     0     0     0     0     0     0     0
#> [61,]     0     0     0     0     0     0     0     0     0     0     0     0
#> [62,]     0     0     0     0     0     0     0     0     0     0     0     0
#> [63,]     0     0     0     0     0     0     0     0     0     0     0     0
#> [64,]     0     0     0     0     0     0     0     0     0     0     0     0
#> [65,]     0     0     0     0     0     0     0     0     0     0     0     0
#> [66,]     0     0     0     0     0     0     0     0     0     0     0     0
#> [67,]     0     0     0     0     0     0     0     0     0     0     0     0
#> [68,]     0     0     0     0     0     0     0     0     0     0     0     0
#> [69,]     0     0     0     0     0     0     0     0     0     0     0     0
#> [70,]     0     0     0     0     0     0     0     0     0     0     0     0
#> [71,]     0     0     0     0     0     0     0     0     0     0     0     0
#> [72,]     0     0     0     0     0     0     0     0     0     0     0     0
#> [73,]     0     0     0     0     0     0     0     0     0     0     0     0
#> [74,]     0     0     0     0     0     0     0     0     0     0     0     0
#> [75,]     0     0     0     0     0     0     0     0     0     0     0     0
#> [76,]     0     0     0     0     0     0     0     0     0     0     0     0
#>       [,50] [,51] [,52] [,53] [,54] [,55] [,56] [,57] [,58] [,59] [,60] [,61]
#>  [1,]     0     0     0     0     0     0     0     0     0     1     0     0
#>  [2,]     0     0     0     0     0     0     0     0     0     0     0     0
#>  [3,]     0     0     0     0     0     0     0     0     0     0     0     0
#>  [4,]     0     0     0     0     0     0     0     0     0     0     0     0
#>  [5,]     0     0     0     0     0     0     0     0     0     0     0     0
#>  [6,]     1     0     0     0     0     0     0     0     0     0     0     0
#>  [7,]     0     0     0     0     0     0     0     0     0     0     0     0
#>  [8,]     0     0     0     0     0     0     0     0     0     0     0     0
#>  [9,]     0     0     0     0     0     0     1     0     0     0     0     0
#> [10,]     0     0     0     0     0     0     0     0     0     0     0     0
#> [11,]     0     0     0     0     0     0     0     1     0     0     0     1
#> [12,]     0     0     1     0     0     0     0     0     0     0     0     0
#> [13,]     0     0     0     0     0     0     0     0     0     0     0     0
#> [14,]     0     0     0     0     0     0     0     0     0     0     0     0
#> [15,]     0     0     0     0     0     0     0     0     0     0     0     0
#> [16,]     0     0     0     0     0     0     0     0     0     0     0     0
#> [17,]     0     0     0     0     0     0     1     0     0     0     0     0
#> [18,]     0     0     0     0     0     0     0     0     0     0     0     0
#> [19,]     0     0     0     0     0     0     0     0     0     0     0     0
#> [20,]     0     0     0     0     0     0     0     0     0     0     0     0
#> [21,]     0     0     0     0     0     0     0     0     0     0     0     0
#> [22,]     0     0     0     0     0     0     0     0     0     0     0     0
#> [23,]     0     0     0     0     0     0     0     0     0     0     0     0
#> [24,]     0     0     0     0     0     0     0     0     0     0     0     0
#> [25,]     0     0     0     0     0     0     0     0     0     0     0     0
#> [26,]     0     0     0     0     0     0     0     0     0     0     0     0
#> [27,]     0     0     0     0     0     0     0     0     0     0     0     0
#> [28,]     0     0     0     0     0     0     0     0     0     0     0     0
#> [29,]     0     0     0     0     0     0     0     0     0     0     0     0
#> [30,]     0     0     0     0     0     0     0     0     0     0     0     0
#> [31,]     0     0     0     0     0     0     0     0     0     0     0     0
#> [32,]     0     0     0     0     0     0     0     0     0     0     0     0
#> [33,]     0     0     0     0     0     0     0     0     0     0     0     0
#> [34,]     0     0     0     0     0     0     0     0     0     0     0     0
#> [35,]     0     0     0     0     0     0     0     0     0     0     0     0
#> [36,]     0     0     0     0     0     0     0     0     0     0     0     0
#> [37,]     0     0     0     0     0     0     0     0     0     0     0     0
#> [38,]     0     0     0     0     0     0     0     0     0     0     0     0
#> [39,]     0     0     0     0     0     0     0     0     0     0     0     0
#> [40,]     0     0     0     0     0     0     0     0     0     0     0     0
#> [41,]     0     0     0     0     0     0     0     0     0     0     0     0
#> [42,]     0     0     0     0     0     0     0     0     0     0     0     0
#> [43,]     0     0     0     0     0     0     0     0     0     0     0     0
#> [44,]     0     0     0     0     0     0     0     0     0     0     0     0
#> [45,]     0     0     0     0     0     0     0     0     0     0     0     0
#> [46,]     0     0     0     0     0     0     0     0     0     0     0     0
#> [47,]     0     0     0     0     0     0     0     0     0     0     0     0
#> [48,]     0     0     0     0     0     0     0     0     0     0     0     0
#> [49,]     0     0     0     0     0     0     0     0     0     0     0     0
#> [50,]     0     0     0     0     0     0     0     0     0     0     0     0
#> [51,]     0     0     0     0     0     0     0     0     0     0     0     0
#> [52,]     0     0     0     0     0     0     0     0     0     0     0     0
#> [53,]     0     0     0     0     0     0     0     0     0     0     0     0
#> [54,]     0     0     0     0     0     0     0     0     0     0     0     0
#> [55,]     0     0     0     0     0     0     0     0     0     0     0     0
#> [56,]     0     0     0     0     0     0     0     0     0     0     0     0
#> [57,]     0     0     0     0     0     0     0     0     0     0     0     0
#> [58,]     0     0     0     0     0     0     0     0     0     0     0     0
#> [59,]     0     0     0     0     0     0     0     0     0     0     0     0
#> [60,]     0     0     0     0     0     0     0     0     0     0     0     0
#> [61,]     0     0     0     0     0     0     0     0     0     0     0     0
#> [62,]     0     0     0     0     0     0     0     0     0     0     0     0
#> [63,]     0     0     0     0     0     0     0     0     0     0     0     0
#> [64,]     0     0     0     0     0     0     0     0     0     0     0     0
#> [65,]     0     0     0     0     0     0     0     0     0     0     0     0
#> [66,]     0     0     0     0     0     0     0     0     0     0     0     0
#> [67,]     0     0     0     0     0     0     0     0     0     0     0     0
#> [68,]     0     0     0     0     0     0     0     0     0     0     0     0
#> [69,]     0     0     0     0     0     0     0     0     0     0     0     0
#> [70,]     0     0     0     0     0     0     0     0     0     0     0     0
#> [71,]     0     0     0     0     0     0     0     0     0     0     0     0
#> [72,]     0     0     0     0     0     0     0     0     0     0     0     0
#> [73,]     0     0     0     0     0     0     0     0     0     0     0     0
#> [74,]     0     0     0     0     0     0     0     0     0     0     0     0
#> [75,]     0     0     0     0     0     0     0     0     0     0     0     0
#> [76,]     0     0     0     0     0     0     0     0     0     0     0     0
#>       [,62] [,63] [,64] [,65] [,66] [,67] [,68] [,69] [,70] [,71] [,72] [,73]
#>  [1,]     0     0     0     0     0     0     0     0     0     0     0     0
#>  [2,]     0     0     0     0     0     0     0     0     0     0     0     0
#>  [3,]     0     0     0     0     0     0     0     0     0     0     0     0
#>  [4,]     0     0     0     0     0     0     0     0     0     0     0     0
#>  [5,]     0     0     0     0     0     0     0     0     0     0     0     0
#>  [6,]     0     0     0     0     0     0     0     0     0     0     0     0
#>  [7,]     0     0     0     0     0     0     0     0     0     0     0     0
#>  [8,]     0     0     0     0     0     0     0     0     0     0     0     0
#>  [9,]     0     0     0     0     0     0     0     0     0     0     0     0
#> [10,]     0     0     0     0     0     0     0     0     0     0     0     0
#> [11,]     0     0     0     0     0     0     0     0     0     0     0     0
#> [12,]     0     0     0     0     0     0     0     0     0     0     0     0
#> [13,]     0     0     0     0     0     0     0     0     0     0     0     0
#> [14,]     0     0     0     0     0     0     0     0     0     0     0     0
#> [15,]     0     0     0     0     0     0     0     0     0     0     0     0
#> [16,]     0     0     0     0     0     0     0     0     0     0     0     0
#> [17,]     0     0     0     0     0     0     0     0     0     0     0     0
#> [18,]     0     0     0     0     0     0     0     0     0     0     0     0
#> [19,]     0     0     0     0     0     0     0     0     0     0     0     0
#> [20,]     0     0     0     0     0     0     0     0     0     0     0     0
#> [21,]     0     0     0     0     0     0     0     0     0     0     0     0
#> [22,]     0     0     0     0     0     0     0     0     0     0     0     0
#> [23,]     0     0     0     1     0     0     0     0     0     0     0     0
#> [24,]     0     0     0     0     0     0     0     0     0     0     0     0
#> [25,]     0     0     0     0     0     0     0     0     0     0     0     0
#> [26,]     0     0     0     0     0     0     0     0     0     0     0     0
#> [27,]     0     0     0     0     0     0     0     0     0     0     0     0
#> [28,]     0     0     0     0     0     0     0     0     0     0     0     0
#> [29,]     0     0     0     0     0     0     0     0     0     0     0     0
#> [30,]     0     0     0     0     0     0     0     0     0     0     0     0
#> [31,]     0     0     0     0     0     0     0     0     0     0     0     0
#> [32,]     0     0     0     0     0     0     0     0     0     0     0     0
#> [33,]     0     0     0     0     0     0     0     0     0     0     0     0
#> [34,]     0     0     0     0     0     0     0     0     0     0     0     0
#> [35,]     0     0     0     0     0     0     0     0     0     0     0     0
#> [36,]     0     0     0     0     0     0     0     0     0     0     0     0
#> [37,]     0     0     0     0     0     0     0     0     0     0     0     0
#> [38,]     0     0     0     0     0     0     0     0     0     0     0     0
#> [39,]     0     0     0     0     0     0     0     0     0     0     0     0
#> [40,]     0     0     0     0     0     0     0     0     0     0     0     0
#> [41,]     0     0     0     0     0     0     0     0     0     0     0     0
#> [42,]     0     0     0     0     0     0     0     0     0     0     0     0
#> [43,]     0     0     0     0     0     0     0     0     0     0     0     0
#> [44,]     0     0     0     0     0     0     0     0     0     0     0     0
#> [45,]     0     0     0     0     0     0     0     0     0     0     0     0
#> [46,]     0     0     0     0     0     0     0     0     0     0     0     0
#> [47,]     0     0     0     0     0     0     0     0     0     0     0     0
#> [48,]     0     0     0     0     0     0     0     0     0     0     0     0
#> [49,]     0     0     0     0     0     0     0     0     0     0     0     0
#> [50,]     0     0     0     0     0     0     0     0     0     0     0     0
#> [51,]     0     0     0     0     0     0     0     0     0     0     0     0
#> [52,]     0     0     0     0     0     0     0     0     0     0     0     0
#> [53,]     0     0     0     0     0     0     0     0     0     0     0     0
#> [54,]     0     0     0     0     0     0     0     0     0     0     0     0
#> [55,]     0     0     0     0     0     0     0     0     0     0     0     0
#> [56,]     0     0     0     0     0     0     0     0     0     0     0     0
#> [57,]     0     0     0     0     0     0     0     0     0     0     0     0
#> [58,]     0     0     0     0     0     0     0     0     0     0     0     0
#> [59,]     0     0     0     0     0     0     0     0     0     0     0     0
#> [60,]     0     0     0     0     0     0     0     0     0     0     0     0
#> [61,]     0     0     0     0     0     0     0     0     0     0     0     0
#> [62,]     0     0     0     0     0     0     0     0     0     0     0     0
#> [63,]     0     0     0     0     0     0     0     0     0     0     0     0
#> [64,]     0     0     0     0     0     0     0     0     0     0     0     0
#> [65,]     0     0     0     0     0     0     0     0     0     0     0     0
#> [66,]     0     0     0     0     0     0     0     0     0     0     0     0
#> [67,]     0     0     0     0     0     0     0     0     0     0     0     0
#> [68,]     0     0     0     0     0     0     0     0     0     0     0     0
#> [69,]     0     0     0     0     0     0     0     0     0     0     0     0
#> [70,]     0     0     0     0     0     0     0     0     0     0     0     0
#> [71,]     0     0     0     0     0     0     0     0     0     0     0     0
#> [72,]     0     0     0     0     0     0     0     0     0     0     0     0
#> [73,]     0     0     0     0     0     0     0     0     0     0     0     0
#> [74,]     0     0     0     0     0     0     0     0     0     0     0     0
#> [75,]     0     0     0     0     0     0     0     0     0     0     0     0
#> [76,]     0     0     0     0     0     0     0     0     0     0     0     0
#>       [,74] [,75] [,76] [,77] [,78] [,79] [,80] [,81] [,82] [,83] [,84] [,85]
#>  [1,]     0     0     0     0     0     0     0     0     0     0     0     0
#>  [2,]     0     0     0     0     0     0     0     0     0     0     0     0
#>  [3,]     0     0     0     0     0     0     0     0     0     0     0     0
#>  [4,]     0     0     0     0     0     0     0     0     0     0     0     0
#>  [5,]     0     0     0     0     0     0     0     0     0     0     0     0
#>  [6,]     0     0     0     0     0     0     0     0     0     0     0     0
#>  [7,]     0     0     0     0     0     0     0     0     0     0     0     0
#>  [8,]     0     0     0     0     1     0     0     0     0     0     0     0
#>  [9,]     0     0     0     0     0     0     0     0     0     0     0     0
#> [10,]     0     0     0     0     0     0     0     0     0     0     0     0
#> [11,]     0     0     0     0     0     0     0     0     0     0     0     0
#> [12,]     0     0     0     0     0     0     0     0     0     0     0     0
#> [13,]     0     0     0     0     0     0     0     0     0     0     0     0
#> [14,]     0     0     0     0     0     0     0     0     0     0     0     0
#> [15,]     0     0     0     0     0     0     0     0     0     0     0     0
#> [16,]     0     0     0     0     0     0     0     0     0     0     0     0
#> [17,]     0     0     0     0     0     0     0     0     0     0     0     0
#> [18,]     0     0     0     0     0     0     0     0     0     0     0     0
#> [19,]     0     0     0     0     0     0     0     0     0     0     0     0
#> [20,]     0     0     0     0     0     0     0     0     0     0     0     0
#> [21,]     0     0     0     0     0     0     0     0     0     0     0     0
#> [22,]     0     0     0     0     0     0     0     0     0     0     0     0
#> [23,]     0     0     0     0     0     0     0     0     0     0     0     0
#> [24,]     0     0     0     0     0     0     0     0     0     0     0     0
#> [25,]     0     0     0     0     0     0     0     0     0     0     0     0
#> [26,]     0     0     0     0     0     0     0     0     0     0     0     0
#> [27,]     0     0     0     0     0     0     0     0     0     0     0     0
#> [28,]     0     0     0     0     0     0     0     0     0     0     0     0
#> [29,]     0     0     0     0     0     0     0     0     0     0     0     0
#> [30,]     0     0     0     0     0     0     0     0     0     0     0     0
#> [31,]     0     0     0     0     0     0     0     0     0     0     0     0
#> [32,]     0     0     0     0     0     0     0     0     0     0     0     0
#> [33,]     0     0     0     0     0     0     0     0     0     0     0     0
#> [34,]     0     0     0     0     0     0     0     0     0     0     0     0
#> [35,]     0     0     0     0     0     0     0     0     0     0     0     0
#> [36,]     0     0     0     0     0     0     0     0     0     0     0     0
#> [37,]     0     0     0     0     0     0     0     0     0     0     0     0
#> [38,]     0     0     0     0     0     0     0     0     0     0     0     0
#> [39,]     0     0     0     0     0     0     0     0     0     0     0     0
#> [40,]     0     0     0     0     0     0     0     0     0     0     0     0
#> [41,]     0     0     0     0     0     0     0     0     0     0     0     0
#> [42,]     0     0     0     0     0     0     0     0     0     0     0     0
#> [43,]     0     0     0     0     0     0     0     0     0     0     0     0
#> [44,]     0     0     0     0     0     0     0     0     0     0     0     0
#> [45,]     0     0     0     0     0     0     0     0     0     0     0     0
#> [46,]     0     0     0     0     0     0     0     0     0     0     0     0
#> [47,]     0     0     0     0     0     0     0     0     0     0     0     0
#> [48,]     0     0     0     0     0     0     0     0     0     0     0     0
#> [49,]     0     0     0     0     0     0     0     0     0     0     0     0
#> [50,]     0     0     0     0     0     0     0     0     0     0     0     0
#> [51,]     0     0     0     0     0     0     0     0     0     0     0     0
#> [52,]     0     0     0     0     0     0     0     0     0     0     0     0
#> [53,]     0     0     0     0     0     0     0     0     0     0     0     0
#> [54,]     0     0     0     0     0     0     0     0     0     0     0     0
#> [55,]     0     0     0     0     0     0     0     0     0     0     0     0
#> [56,]     0     0     0     0     0     0     0     0     0     0     0     0
#> [57,]     0     0     0     0     0     0     0     0     0     0     0     0
#> [58,]     0     0     0     0     0     0     0     0     0     0     0     0
#> [59,]     0     0     0     0     0     0     0     0     0     0     0     0
#> [60,]     0     0     0     0     0     0     0     0     0     0     0     0
#> [61,]     0     0     0     0     0     0     0     0     0     0     0     0
#> [62,]     0     0     0     0     0     0     0     0     0     0     0     0
#> [63,]     0     0     0     0     0     0     0     0     0     0     0     0
#> [64,]     0     0     0     0     0     0     0     0     0     0     0     0
#> [65,]     0     0     0     0     0     0     0     0     0     0     0     0
#> [66,]     0     0     0     0     0     0     0     0     0     0     0     0
#> [67,]     0     0     0     0     0     0     0     0     0     0     0     0
#> [68,]     0     0     0     0     0     0     0     0     0     0     0     0
#> [69,]     0     0     0     0     0     0     0     0     0     0     0     0
#> [70,]     0     0     0     0     0     0     0     0     0     0     0     0
#> [71,]     0     0     0     0     0     0     0     0     0     0     0     0
#> [72,]     0     0     0     0     0     0     0     0     0     0     0     0
#> [73,]     0     0     0     0     0     0     0     0     0     0     0     0
#> [74,]     0     0     0     0     0     0     0     0     0     0     0     0
#> [75,]     0     0     0     0     0     0     0     0     0     0     0     0
#> [76,]     0     0     0     0     0     0     0     0     0     0     0     0
#>       [,86] [,87] [,88] [,89] [,90] [,91]
#>  [1,]     0     0     0     0     0     0
#>  [2,]     0     0     0     0     0     0
#>  [3,]     0     0     0     0     0     0
#>  [4,]     0     0     0     0     0     0
#>  [5,]     0     0     0     0     0     0
#>  [6,]     0     0     0     0     0     0
#>  [7,]     0     0     0     0     0     0
#>  [8,]     0     0     0     0     0     0
#>  [9,]     0     0     0     0     0     0
#> [10,]     0     0     0     0     0     0
#> [11,]     0     0     0     0     0     0
#> [12,]     0     0     0     0     0     0
#> [13,]     0     0     0     0     0     0
#> [14,]     0     0     0     0     0     0
#> [15,]     0     0     0     0     0     0
#> [16,]     0     0     0     0     0     0
#> [17,]     0     0     0     0     0     0
#> [18,]     0     0     0     0     0     0
#> [19,]     0     0     0     0     0     0
#> [20,]     0     0     0     0     0     0
#> [21,]     0     0     0     0     0     0
#> [22,]     0     0     0     0     0     0
#> [23,]     0     0     0     0     0     0
#> [24,]     0     0     0     0     0     0
#> [25,]     0     0     0     0     0     0
#> [26,]     0     0     0     0     0     0
#> [27,]     0     0     0     0     0     0
#> [28,]     0     0     0     0     0     0
#> [29,]     0     0     0     0     0     0
#> [30,]     0     0     0     0     0     0
#> [31,]     0     0     0     0     0     0
#> [32,]     0     0     0     0     0     0
#> [33,]     0     0     0     0     0     0
#> [34,]     0     0     0     0     0     0
#> [35,]     0     0     0     0     0     0
#> [36,]     0     0     0     0     0     0
#> [37,]     0     0     0     0     0     0
#> [38,]     0     0     0     0     0     0
#> [39,]     0     0     0     0     0     0
#> [40,]     0     0     0     0     0     0
#> [41,]     0     0     0     0     0     0
#> [42,]     0     0     0     0     0     0
#> [43,]     0     0     0     0     0     0
#> [44,]     0     0     0     0     0     0
#> [45,]     0     0     0     0     0     0
#> [46,]     0     0     0     0     0     0
#> [47,]     0     0     0     0     0     0
#> [48,]     0     0     0     0     0     0
#> [49,]     0     0     0     0     0     0
#> [50,]     0     0     0     0     0     0
#> [51,]     0     0     0     0     0     0
#> [52,]     0     0     0     0     0     0
#> [53,]     0     0     0     0     0     0
#> [54,]     0     0     0     0     0     0
#> [55,]     0     0     0     0     0     0
#> [56,]     0     0     0     0     0     0
#> [57,]     0     0     0     0     0     0
#> [58,]     0     0     0     0     0     0
#> [59,]     0     0     0     0     0     0
#> [60,]     0     0     0     0     0     0
#> [61,]     0     0     0     0     0     0
#> [62,]     0     0     0     0     0     0
#> [63,]     0     0     0     0     0     0
#> [64,]     0     0     0     0     0     0
#> [65,]     0     0     0     0     0     0
#> [66,]     0     0     0     0     0     0
#> [67,]     0     0     0     0     0     0
#> [68,]     0     0     0     0     0     0
#> [69,]     0     0     0     0     0     0
#> [70,]     0     0     0     0     0     0
#> [71,]     0     0     0     0     0     0
#> [72,]     0     0     0     0     0     0
#> [73,]     0     0     0     0     0     0
#> [74,]     0     0     0     0     0     0
#> [75,]     0     0     0     0     0     0
#> [76,]     0     0     0     0     0     0

# Visualize midpoints of bins over u
print(binned_data$bins$midu)
#>  [1]   23   53   83  113  143  173  203  233  263  293  323  353  383  413  443
#> [16]  473  503  533  563  593  623  653  683  713  743  773  803  833  863  893
#> [31]  923  953  983 1013 1043 1073 1103 1133 1163 1193 1223 1253 1283 1313 1343
#> [46] 1373 1403 1433 1463 1493 1523 1553 1583 1613 1643 1673 1703 1733 1763 1793
#> [61] 1823 1853 1883 1913 1943 1973 2003 2033 2063 2093 2123 2153 2183 2213 2243
#> [76] 2273


# Bin data over u = time at recurrence and s = time since recurrence, measured in day
# individual-level data required
# we provide two covariates: nodes (numerical) and rx (factor)
covs <- subset(reccolon2ts, select = c("nodes", "rx"))
binned_data <- prepare_data(
  u = reccolon2ts$timer, s_out = reccolon2ts$timesr,
  events = reccolon2ts$status, ds = 30, individual = TRUE, covs = covs
)
#> `s_in = NULL`. I will use `s_in = 0` for all observations.
#> `s_in = NULL`. I will use `s_in = 0` for all observations.

# Visualize structure of binned data
print(str(binned_data$bindata))
#> List of 3
#>  $ R: num [1:76, 1:91, 1:461] 0 0 0 0 0 0 0 0 0 0 ...
#>  $ Y: num [1:76, 1:91, 1:461] 0 0 0 0 0 0 0 0 0 0 ...
#>  $ Z: num [1:461, 1:3] 5 7 6 22 9 5 1 3 1 6 ...
#>   ..- attr(*, "dimnames")=List of 2
#>   .. ..$ : NULL
#>   .. ..$ : chr [1:3] "nodes" "rx_Lev" "rx_Lev+5FU"
#> NULL

# Alternatevely:
binned_data <- prepare_data(
  data = reccolon2ts,
  u = "timer", s_out = "timesr",
  events = "status", ds = 30, individual = TRUE, covs = c("nodes", "rx")
)
#> `s_in = NULL`. I will use `s_in = 0` for all observations.
#> `s_in = NULL`. I will use `s_in = 0` for all observations.
```
