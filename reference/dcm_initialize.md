# Initialize object for modeling using the DCM.

Initialize object for modeling using the DCM.

## Usage

``` r
dcm_initialize(
  object,
  frequency,
  radius_cylinder = NULL,
  radius_curvature = NULL,
  radius_curvature_ratio = 3,
  radius_cylinder_fun = "max",
  length = NULL,
  g = NULL,
  h = NULL,
  theta = NULL,
  sound_speed_sw = 1500,
  density_sw = 1026,
  alpha_B = 0.8
)
```

## Arguments

- object:

  FLS-class object.

- frequency:

  Transmit frequency (kHz)

- radius_cylinder:

  Optional input to override current shape radius.

- radius_curvature:

  Numeric input for the radius of curvature

- radius_curvature_ratio:

  Ratio between body length and the radius of curvature. Defaults to
  3.0.

- radius_cylinder_fun:

  Defines which radius value will be used from the radius vector.
  Defaults to "max", but also accepts "mean" and "median".

- length:

  Body length (m).

- g:

  Density contrast.

- h:

  Sound speed contrast.

- theta:

  Body orientation relative to the incident sound wave.

- sound_speed_sw:

  Seawater sound speed (c_(body), m s⁻¹).

- density_sw:

  Seawater density (ρ_(body), kg m³)

- alpha_B:

  Numerical coefficient (α_(B)).
