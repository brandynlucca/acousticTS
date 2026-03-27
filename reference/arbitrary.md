# Creates arbitrary body shape from user inputs

Creates arbitrary body shape from user inputs

## Usage

``` r
arbitrary(..., length_units = "m")
```

## Arguments

- ...:

  Any coordinate or parameter vectors (e.g., x_body, y_body, z_body,
  radius_body, w_body, zU_body, zL_body, custom\_\*). Inputs are stored
  and padded to a common length; validation happens downstream in
  model-specific code.

- length_units:

  Units for body length. Defaults to meters: "m"

## Value

An
[`Arbitrary`](https://brandynlucca.github.io/acousticTS/reference/Arbitrary-class.md)
shape object containing the padded position matrix and stored shape
metadata.

## See also

[`Arbitrary`](https://brandynlucca.github.io/acousticTS/reference/Arbitrary-class.md)

## Examples

``` r
# Symmetric arbitrary shape using x/y/z + radius
arbitrary(
  x_body = c(0, 0.01), y_body = c(0, 0), z_body = c(0, 0),
  radius_body = c(0, 0.002)
)
#> An object of class "Arbitrary"
#> Slot "position_matrix":
#>         x y z     a    zU     zL
#> [1,] 0.00 0 0 0.000 0.000  0.000
#> [2,] 0.01 0 0 0.002 0.002 -0.002
#> 
#> Slot "shape_parameters":
#> $n_segments
#> [1] 1
#> 
#> $length_units
#> [1] "m"
#> 
#> $diameter_units
#> [1] "m"
#> 
#> $radius
#> [1] 0.000 0.002
#> 
#> $mean_radius
#> [1] 0.001
#> 
#> $max_radius
#> [1] 0.002
#> 
#> 

# Dorsal/ventral style inputs
arbitrary(
  x_body = c(0, 0.015),
  w_body = c(0.005, 0.0075),
  zU_body = c(0.001, 0.002),
  zL_body = c(-0.001, -0.002)
)
#> An object of class "Arbitrary"
#> Slot "position_matrix":
#>          x      w    zU     zL
#> [1,] 0.000 0.0050 0.001 -0.001
#> [2,] 0.015 0.0075 0.002 -0.002
#> 
#> Slot "shape_parameters":
#> $n_segments
#> [1] 1
#> 
#> $length_units
#> [1] "m"
#> 
#> $diameter_units
#> [1] "m"
#> 
#> $radius
#> [1] NA
#> 
#> $mean_radius
#> [1] NA
#> 
#> $max_radius
#> [1] NA
#> 
#> 
```
