# Reforge ESS-class object

Resize an elastic-shelled scatterer by applying an isometric scale
factor or specifying a target maximum shell radius. Shell thickness can
be updated independently, which rescales the fluid body so that the
maximum fluid radius equals `new_max_shell_radius - shell_thickness`
(matching the convention in
[`ess_generate`](https://brandynlucca.github.io/acousticTS/reference/ess_generate.md)).
The underlying shape (sphere, prolate spheroid, cylinder, etc.) is
preserved for both shell and fluid bodies; scaling is applied uniformly
to all axes of each position matrix.

## Usage

``` r
# S4 method for class 'ESS'
reforge(
  object,
  scale = NULL,
  radius_target = NULL,
  shell_thickness = NULL,
  n_segments = NULL
)
```

## Arguments

- object:

  ESS-class object.

- scale:

  Single positive scalar applied isometrically to the shell (and fluid
  body, if present). Mutually exclusive with `radius_target`.

- radius_target:

  Target *maximum* outer shell radius (m). Scale factor derived as
  `radius_target / max(current_shell_radius)`. Mutually exclusive with
  `scale`.

- shell_thickness:

  New shell wall thickness (m). The fluid body is rescaled so its
  maximum radius equals `new_max_shell_radius - shell_thickness`. Can be
  combined with `scale`/`radius_target` or used alone.

- n_segments:

  New number of discrete segments. All columns of both the shell and
  fluid position matrices are re-interpolated along the x-axis.

## Value

Modified ESS-class object.
