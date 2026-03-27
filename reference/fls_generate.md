# Generate a FLS-class object.

Generate a FLS-class object.

## Usage

``` r
fls_generate(
  shape = NULL,
  x_body = NULL,
  y_body = NULL,
  z_body = NULL,
  length_body = NULL,
  radius_body = NULL,
  radius_curvature_ratio = NULL,
  n_segments = 18,
  g_body = NULL,
  h_body = NULL,
  density_body = NULL,
  sound_speed_body = NULL,
  theta_body = pi/2,
  ID = NULL,
  length_units = "m",
  theta_units = "radians",
  ...
)
```

## Arguments

- shape:

  Pre-built `Shape` object describing the target geometry. If omitted,
  explicit profile coordinates such as `x_body`, `y_body`, and `z_body`
  are treated as the manual geometry pathway. Legacy character dispatch
  such as `"sphere"` or `"arbitrary"` is retained only for backward
  compatibility and is now deprecated.

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

- density_body:

  Absolute density (kg/m^3) if contrasts are not supplied.

- sound_speed_body:

  Absolute sound speed (m/s) if contrasts are not supplied.

- theta_body:

  Orientation of the target relative to the transmit source
  (\\\theta\\). Broadside incidence is considered 90 degrees, or pi/2.
  Default value is pi/2; input should be in radians.

- ID:

  Optional metadata entry.

- length_units:

  Compatibility argument. Scatterer constructors now assume meters and
  ignore non-SI alternatives.

- theta_units:

  Compatibility argument. Scatterer constructors now assume radians and
  ignore non-SI alternatives.

- ...:

  Additional parameters.

## Value

FLS-class object

## Details

The preferred workflow is to build a geometry first with
[`sphere()`](https://brandynlucca.github.io/acousticTS/reference/sphere.md),
[`cylinder()`](https://brandynlucca.github.io/acousticTS/reference/cylinder.md),
[`prolate_spheroid()`](https://brandynlucca.github.io/acousticTS/reference/prolate_spheroid.md),
or
[`arbitrary()`](https://brandynlucca.github.io/acousticTS/reference/arbitrary.md),
then pass that `Shape` object to `fls_generate()`. Material properties
can be supplied either as contrasts (`g_body`/`h_body`) or as absolute
density/sound-speed values (`density_body`/`sound_speed_body`), but not
both for the same property pair. Downstream models derive contrasts
automatically when only absolute values are supplied.

The only supported public geometry paths are now:

1.  supply a pre-built `Shape` object, or

2.  supply explicit profile coordinates directly to the constructor.

Character-based shape dispatch is retained only as a compatibility
pathway and is now deprecated. Internally, every pathway is resolved to
the same `Shape`-first geometry contract before the `FLS` object is
built.

Scatterer constructors store geometry in meters and orientations in
radians. `length_units` and `theta_units` are retained as compatibility
arguments, but non-SI values are normalized to the package-standard
representation.

## See also

[`FLS`](https://brandynlucca.github.io/acousticTS/reference/FLS-class.md)

## Examples

``` r
shape <- prolate_spheroid(
  length_body = 0.04, radius_body = 0.004, n_segments = 50
)
fls_generate(
  shape = shape, density_body = 1045, sound_speed_body = 1520
)
#> FLS-object
#>  Fluid-like scatterer 
#>  ID:UID
#> Body dimensions:
#>  Length:0.04 m(n = 50 cylinders)
#>  Mean radius:0.0031 m
#>  Max radius:0.004 m
#> Shape parameters:
#>  Defined shape:ProlateSpheroid
#>  L/a ratio:10
#>  Taper order:N/A
#> Material properties:
#>  Density: 1045 kg m^-3 | Sound speed: 1520 m s^-1
#> Body orientation (relative to transducer face/axis):1.571 radians
```
