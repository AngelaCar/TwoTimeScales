# Print method for a `data2ts` object

Print method for an object of class `data2ts`

## Usage

``` r
# S3 method for class 'data2ts'
print(x, ...)
```

## Arguments

- x:

  of class `data2ts`, as prepared by
  [`prepare_data`](https://angelacar.github.io/TwoTimeScales/reference/prepare_data.md)

- ...:

  Further arguments to print

## Value

No return value

## Author

Angela Carollo <carollo@demogr.mpg.de>

## Examples

``` r
# Bin the colon cancer data over s (time since recurrence)

dt1ts <- prepare_data(data = reccolon2ts,
                      s_in = "entrys",
                      s_out = "timesr",
                      events = "status",
                      ds = 180)

print(dt1ts)
#> An object of class 'data2ts'
#> 
#> Data:
#> List of 2
#>  $ bins   :List of 3
#>  $ bindata:List of 2
#>  - attr(*, "class")= chr "data2ts"
#> NULL
#> 
#> Range covered by the bins: 
#> $bins_s
#> [1]    0 2880
#> 
#> 
#> Number of bins: 
#> $ns
#> [1] 16
#> 
#> 
#> Overview of the binned data:
#> Total exposure time: 236006
#> Total number of events: 409

# Bin the colon cancer data over u (time at recurrence) and s (time since recurrence)
dt2ts <- prepare_data(data = reccolon2ts,
                      u = "timer",
                      s_in = "entrys",
                      s_out = "timesr",
                      events = "status",
                      ds = 180)
print(dt2ts)
#> An object of class 'data2ts'
#> 
#> Data:
#> List of 2
#>  $ bins   :List of 6
#>  $ bindata:List of 2
#>  - attr(*, "class")= chr "data2ts"
#> NULL
#> 
#> Range covered by the bins: 
#> $bins_u
#> [1]    8 2348
#> 
#> $bins_s
#> [1]    0 2880
#> 
#> 
#> Number of bins: 
#> $nu
#> [1] 13
#> 
#> $ns
#> [1] 16
#> 
#> 
#> Overview of the binned data:
#> Total exposure time: 236006
#> Total number of events: 409
```
