# Creates a polynomial deformed cylinder

Creates a polynomial deformed cylinder

## Usage

``` r
polynomial_cylinder(length_body, radius_body, n_segments, polynomial,
length_units)
```

## Arguments

- length_body:

  Length (m).

- radius_body:

  Maximum/uniform radius (m).

- n_segments:

  Number of segments to discretize object shape. Defaults to 1e2
  segments.

- polynomial:

  Polynomial coefficient vector.

- length_units:

  Units (default is meters, "m").

## Value

Creates the position vector for a polynomial deformed cylinder.

## References

Smith, J.N., Ressler, P.H., and Warren, J.D. 2013. A distorted wave Born
approximation target strength model for Bering Sea euphausiids. ICES
Journal of Marine Science, 70(1): 204-214.
https://doi.org/10.1093/icesjms/fss140

## See also

[`PolynomialCylinder`](https://brandynlucca.github.io/acousticTS/reference/PolynomialCylinder-class.md)

## Examples

``` r
# We can use the polynomial coefficients defined in Smith et al. (2013) to
# define the position vector of a sub-Arctic krill.
poly_vec <- c(0.83, 0.36, -2.10, -1.20, 0.63, 0.82, 0.64)
# Create the position vector
# This outputs a list containing "rpos" and "radius"
pos <- polynomial_cylinder(
  length_body = 15e-3, radius_body = 2e-3,
  polynomial = poly_vec
)
str(pos)
#> Formal class 'PolynomialCylinder' [package "acousticTS"] with 2 slots
#>   ..@ position_matrix : num [1:101, 1:3] 0 0.00015 0.0003 0.00045 0.0006 ...
#>   .. ..- attr(*, "dimnames")=List of 2
#>   .. .. ..$ : NULL
#>   .. .. ..$ : chr [1:3] "x" "y" "z"
#>   ..@ shape_parameters:List of 5
#>   .. ..$ length             : num 0.015
#>   .. ..$ radius             : num [1:101] 4.00e-05 8.75e-05 1.35e-04 1.83e-04 2.30e-04 ...
#>   .. ..$ length_radius_ratio: num 7.85
#>   .. ..$ n_segments         : num 100
#>   .. ..$ length_units       : chr "m"
```
