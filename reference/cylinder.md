# Creates a cylinder.

Creates a cylinder.

## Usage

``` r
cylinder(length_body, radius_body, length_radius_ratio,
taper, n_segments, length_units)
```

## Arguments

- length_body:

  Length (m).

- radius_body:

  Maximum/uniform radius (m).

- length_radius_ratio:

  Optional ratio input when radius is not explicitly known.

- taper:

  Optional input that is the degree of taper to round ends of the
  cylinder.

- n_segments:

  Number of segments to discretize object shape. Defaults to 1e2
  segments.

- length_units:

  Units (default is meters, "m").

## Value

Creates the position vector for a tapered or untapered cylinder.

## See also

[`Cylinder`](https://brandynlucca.github.io/acousticTS/reference/Cylinder-class.md)
