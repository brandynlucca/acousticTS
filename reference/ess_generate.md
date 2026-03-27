# Generate ESS shape

Generate ESS shape

## Usage

``` r
ess_generate(
  shape = "sphere",
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

  Optional input argument that dictates shape-type, if desired, for
  generalized shapes.

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

  Units used for orientation. Defaults to "radians".

- length_units:

  Units used for position vector. Defaults to "m".

## Value

ESS-class object

## See also

[`ESS`](https://brandynlucca.github.io/acousticTS/reference/ESS-class.md)
