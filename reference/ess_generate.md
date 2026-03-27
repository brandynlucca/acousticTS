# Generate an ESS-class object

Generate an ESS-class object

## Usage

``` r
ess_generate(
  shape = NULL,
  x_body = NULL,
  y_body = NULL,
  z_body = NULL,
  radius_shell = NULL,
  shell_thickness = NULL,
  g_fluid = NULL,
  density_fluid = NULL,
  h_fluid = NULL,
  sound_speed_fluid = NULL,
  g_shell = NULL,
  density_shell = NULL,
  h_shell = NULL,
  sound_speed_shell = NULL,
  E = NULL,
  G = NULL,
  K = NULL,
  nu = NULL,
  theta_shell = pi/2,
  ID = NULL,
  theta_units = "radians",
  length_units = "m"
)
```

## Arguments

- shape:

  Pre-built `Shape` object describing the outer shell geometry. If
  omitted, explicit shell profile coordinates such as `x_body`,
  `y_body`, and `z_body` are treated as the manual geometry pathway.
  Legacy character dispatch such as `"sphere"` is retained only for
  backward compatibility and is now deprecated.

- x_body:

  Vector containing x-axis body (m) shape data.

- y_body:

  Vector containing y-axis body (m) shape data.

- z_body:

  Vector containing z-axis body (m) shape data.

- radius_shell:

  Radius of shell (m).

- shell_thickness:

  Optional shell thickness (m).

- g_fluid:

  Optional density contrast for fluid-like body.

- density_fluid:

  Optional density for fluid-like body (kg/m³).

- h_fluid:

  Optional sound speed contrast for fluid-like body.

- sound_speed_fluid:

  Optional sound speed for fluid-like body (m/s).

- g_shell:

  Density contrast for the shell.

- density_shell:

  Optional density for the shell (kg/m³).

- h_shell:

  Sound speed contrast for the shell.

- sound_speed_shell:

  Optional sound speed for the shell (m/s).

- E:

  Young's modulus (Pa) of the shell material.

- G:

  Shear modulus (Pa) of the shell material.

- K:

  Bulk modulus (Pa) of the shell material.

- nu:

  Poisson's ratio (Dimensionless) of the shell material.

- theta_shell:

  Object orientation relative to incident sound wave.

- ID:

  Optional metadata entry.

- theta_units:

  Compatibility argument. Scatterer constructors now assume radians and
  ignore non-SI alternatives.

- length_units:

  Compatibility argument. Scatterer constructors now assume meters and
  ignore non-SI alternatives.

## Value

ESS-class object

## Details

The preferred workflow is to build the shell geometry first as a `Shape`
object and pass it through `shape`, or to supply explicit shell profile
coordinates directly. Character-based shape dispatch remains available
only as a compatibility pathway and is now deprecated. Material
properties for both the shell and the inner fluid may be supplied either
as contrasts or as absolute values, but each property pair must use only
one representation.

Scatterer constructors store geometry in meters and orientations in
radians. `length_units` and `theta_units` are retained as compatibility
arguments, but non-SI values are normalized to the package-standard
representation.

## See also

[`ESS`](https://brandynlucca.github.io/acousticTS/reference/ESS-class.md)

## Examples

``` r
shell <- sphere(radius_body = 0.03, n_segments = 80)
ess_generate(
  shape = shell,
  shell_thickness = 0.001,
  density_shell = 1050,
  sound_speed_shell = 2350,
  density_fluid = 1030,
  sound_speed_fluid = 1500,
  E = 3.5e9,
  nu = 0.34
)
#> ESS-object 
#>  Elastic-shelled scatterer 
#>   ID: UID 
#>  Material:  
#>    Shell: 
#>      Density: 1050 kg m^-3
#>      Sound speed: 2350 m s^-1
#>      Young's modulus (E): 3.5e+09 Pa
#>      Poisson's ratio: 0.34
#>      Bulk modulus (K): 3645833333.3333 Pa
#>      Shear modulus (G): 1305970149.2537 Pa  
#>    Internal fluid-like body: 
#>      Density: 1030 kg m^-3
#>      Sound speed: 1500 m s^-1  
#>  Shape: 
#>    Shell: 
#>      Radius: 0.03 m  
#>      Diameter: 0.06 m  
#>      Outer thickness: 0.001 m 
#>    Internal fluid-like body: 
#>      Radius: 0.029 m  
#>      Diameter: 0.058 m  
#>  Propagation direction of the incident sound wave: 1.571 radians
```
