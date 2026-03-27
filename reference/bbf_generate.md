# Generate a BBF-class object.

Generate a BBF-class object.

## Usage

``` r
bbf_generate(
  body_shape,
  backbone_shape,
  density_body = NULL,
  sound_speed_body = NULL,
  g_body = NULL,
  h_body = NULL,
  density_backbone,
  sound_speed_longitudinal_backbone,
  sound_speed_transversal_backbone,
  theta_body = pi/2,
  theta_backbone = pi/2,
  x_offset_backbone = 0,
  z_offset_backbone = 0,
  theta_units = "radians",
  length_units = "m",
  ID = NULL
)
```

## Arguments

- body_shape:

  Pre-built body shape. Must inherit from `Shape`.

- backbone_shape:

  Pre-built backbone shape. Must be a cylindrical `Shape`.

- density_body:

  Flesh density (ρ_(body), kg m³).

- sound_speed_body:

  Flesh sound speed (c_(body), m s⁻¹).

- g_body:

  Body density contrast.

- h_body:

  Body sound speed contrast.

- density_backbone:

  Backbone density (ρ_(bb), kg m³).

- sound_speed_longitudinal_backbone:

  Longitudinal wave speed in the backbone (m/s).

- sound_speed_transversal_backbone:

  Transversal wave speed in the backbone (m/s).

- theta_body:

  Body orientation relative to the incident wave (radians).

- theta_backbone:

  Backbone orientation relative to the incident wave (radians).

- x_offset_backbone:

  Along-body translation applied to the backbone geometry (m).

- z_offset_backbone:

  Dorsoventral translation applied to the backbone geometry (m).

- theta_units:

  Compatibility argument. Scatterer constructors now assume radians and
  ignore non-SI alternatives.

- length_units:

  Compatibility argument. Scatterer constructors now assume meters and
  ignore non-SI alternatives.

- ID:

  Optional metadata identifier.

## Value

Generates a BBF-class object.

## Details

`bbf_generate()` is intended for swimbladder-less fish workflows where
the flesh and backbone should remain explicit, separately parameterized
components. The body is stored using the same segmented-body
representation used by
[`FLS`](https://brandynlucca.github.io/acousticTS/reference/FLS-class.md),
while the backbone is stored as an explicit cylindrical component
carrying elastic material properties for use by cylinder-based modal
models.

## See also

[`BBF`](https://brandynlucca.github.io/acousticTS/reference/BBF-class.md),
[`FLS`](https://brandynlucca.github.io/acousticTS/reference/FLS-class.md),
[`cylinder`](https://brandynlucca.github.io/acousticTS/reference/cylinder.md)

## Examples

``` r
body_shape <- arbitrary(
  x_body = c(0, 0.04, 0.08),
  zU_body = c(0.001, 0.004, 0.001),
  zL_body = c(-0.001, -0.004, -0.001)
)
backbone_shape <- cylinder(
  length_body = 0.06,
  radius_body = 0.0008,
  n_segments = 40
)
bbf_generate(
  body_shape = body_shape,
  backbone_shape = backbone_shape,
  density_body = 1070,
  sound_speed_body = 1570,
  density_backbone = 1900,
  sound_speed_longitudinal_backbone = 3500,
  sound_speed_transversal_backbone = 1700
)
#> BBF-object
#>  Backboned fish (BBF) 
#>  ID:UID
#> Body dimensions:
#>  Length:0.08 m(n = 2 segments)
#>  Mean radius:0.002 m
#>  Max radius:0.004 m
#> Backbone dimensions:
#>  Length:0.06 m(n = 40 segments)
#>  Mean radius:8e-04 m
#>  Max radius:8e-04 m
#> Body material properties:
#>  Density: 1070 kg m^-3 | Sound speed: 1570 m s^-1
#> Backbone elastic properties:
#>  Density: 1900 kg m^-3 | cL: 3500 m s^-1 | cT: 1700 m s^-1
#> Body orientation:1.571 radians | Backbone orientation:1.571 radians
```
