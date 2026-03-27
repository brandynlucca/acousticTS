# Manually generate a SBF-class object.

Manually generate a SBF-class object.

## Usage

``` r
sbf_generate(
  x_body,
  w_body,
  zU_body,
  zL_body,
  x_bladder,
  w_bladder,
  zU_bladder,
  zL_bladder,
  sound_speed_body,
  sound_speed_bladder,
  density_body,
  density_bladder,
  theta_body = pi/2,
  theta_bladder = pi/2,
  theta_units = "radians",
  length_units = "m",
  ID = NULL
)
```

## Arguments

- x_body:

  Vector containing along-body axis (m).

- w_body:

  Vector containing across-body axis (m).

- zU_body:

  Vector containing dorsal-body axis (m).

- zL_body:

  Vector containing ventral-body axis (m).

- x_bladder:

  Vector containing along-bladder axis (m).

- w_bladder:

  Vector containing across-bladder axis (m).

- zU_bladder:

  Vector containing dorsal-bladder axis (m).

- zL_bladder:

  Vector containing ventral-bladder axis (m).

- sound_speed_body:

  Flesh sound speed (c;_(body), m s⁻¹).

- sound_speed_bladder:

  Bladder sound speed (c, m \\s^-1\\.

- density_body:

  Flesh density (ρ_(body), kg m³).

- density_bladder:

  Bladder density (\\\rho\\, kg m³).

- theta_body:

  Angle of body relative to wavefront (\\\theta_body\\, radians).

- theta_bladder:

  Angle of body relative to wavefront (\\\theta_bladder\\, radians).

- theta_units:

  Angular units.

- length_units:

  Angular units.

- ID:

  Angular units.

## Value

Generates a SBF-class object.

## See also

[`SBF`](https://brandynlucca.github.io/acousticTS/reference/SBF-class.md)
