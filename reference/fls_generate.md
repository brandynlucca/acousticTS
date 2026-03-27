# Manually generate a FLS object.

Manually generate a FLS object.

## Usage

``` r
fls_generate(
  shape = "arbitrary",
  x_body = NULL,
  y_body = NULL,
  z_body = NULL,
  length_body = NULL,
  radius_body = NULL,
  radius_curvature_ratio = NULL,
  n_segments = 18,
  g_body,
  h_body,
  theta_body = pi/2,
  ID = NULL,
  length_units = "m",
  theta_units = "radians",
  ...
)
```

## Arguments

- shape:

  Optional input argument that dictates shape-type, if desired, for
  generalized shapes.

- x_body:

  Vector containing x-axis body (m) shape data.

- y_body:

  Vector containing y-axis body (m) shape data.

- z_body:

  Vector containing z-axis body (m) shape data.

- length_body:

  Optional input for a generic length value input.

- radius_body:

  Vector containing radii (m).

- radius_curvature_ratio:

  Length-to-curvature ratio (pc/L).

- n_segments:

  Number of body segments.

- g_body:

  Density contrast. This can either be a single value (i.e. homogenous)
  or a vector of values (i.e. inhomogenous).

- h_body:

  Soundspeed contrast. This can either be a single value (i.e.
  homogenous) or a vector of values (i.e. inhomogenous).

- theta_body:

  Orientation of the target relative to the transmit source
  (\\\theta\\). Broadside incidence is considered 90 degrees, or pi/2.
  Default value is pi/2; input should be in radians.

- ID:

  Optional metadata entry.

- length_units:

  Units used for position vector. Defaults to "m".

- theta_units:

  Units used for orientation. Defaults to "radians".

- ...:

  Additional parameters.

## Value

FLS-class object

## See also

[`FLS`](https://brandynlucca.github.io/acousticTS/reference/FLS-class.md)
