# Calculates the Euclidean norm across each row of a given matrix.

Calculates the Euclidean norm across each row of a given matrix.

## Usage

``` r
vecnorm( x )
```

## Arguments

- x:

  A matrix with numeric, real values.

## Value

Calculates the Euclidean norm of a vector.

## Examples

``` r
values <- matrix(c(1, 2, 3), ncol = 3)
vecnorm(values) # should yield 3.741657
#> [1] 3.741657
```
