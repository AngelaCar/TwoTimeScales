# Process data to fit model with LMMsolver

Process data to fit model with LMMsolver

## Usage

``` r
prepare_data_LMMsolver(Y = Y, R = R, Z = NULL, bins = bins)
```

## Arguments

- Y:

  A matrix (or 3d-array) of event counts of dimension nu by ns (or nu by
  ns by n).

- R:

  A matrix (or 3d-array) of exposure times of dimension nu by ns (or nu
  by ns by n).

- Z:

  (optional) A regression matrix of covariates values of dimensions n by
  p.

- bins:

  a list with the specification for the bins. This is created by the
  function `prepare_data`. If a list prepared externally from such
  function if provided, it should contain the following elements: \*
  `bins_u` A vector of bins extremes for the time scale `u`. \* `midu` A
  vector with the midpoints of the bins over `u`. \* `nu` The number of
  bins over `u`. \* `bins_s` A vector of bins extremes for the time
  scale `s`. \* `mids` A vector with the midpoints of the bins over `s`.
  \* `ns` The number of bins over `s`.

## Value

A dataset in long form to fit the model with LMMsolver
