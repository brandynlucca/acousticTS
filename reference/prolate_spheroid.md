# Creates a prolate spheroid.

Creates a prolate spheroid.

## Usage

``` r
prolate_spheroid(
  length_body,
  radius_body,
  length_radius_ratio = NULL,
  n_segments = 18,
  length_units = "m",
  theta_units = "radians"
)
```

## Arguments

- length_body:

  Semi-major axis length (m).

- radius_body:

  Semi-minor axis length (m). This can also be stylized as the "maximum
  radius" of the scattering object.

- length_radius_ratio:

  Optional ratio input when radius is not explicitly known.

- n_segments:

  Number of segments to discretize object shape. Defaults to 18
  segments.

- length_units:

  Units for body matrix (defaults to m).

- theta_units:

  Units for body orientation (defaults to radians).

## Value

Creates the position vector for a prolate spheroid object of defined
semi-major and -minor axes.

## See also

[`ProlateSpheroid`](https://brandynlucca.github.io/acousticTS/reference/ProlateSpheroid-class.md)
