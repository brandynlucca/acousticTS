# Along-matrix summing function

Along-matrix summing function

## Usage

``` r
along_sum(rpos, iterations)
```

## Arguments

- rpos:

  Position vector

- iterations:

  Number of iterations

## Value

A numeric matrix containing the adjacent column sums of `rpos`.

## Examples

``` r
positions <- matrix(1:8, nrow = 2)
along_sum(positions, ncol(positions))
#>      [,1] [,2] [,3]
#> [1,]    4    8   12
#> [2,]    6   10   14
```
