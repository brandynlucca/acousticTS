# Reforge GAS-class object

Resize a gas-filled fluid scatterer. `GAS` bodies are single-component
fluid-like targets, so the interface mirrors
[`reforge()`](https://brandynlucca.github.io/acousticTS/reference/reforge.md)
for `FLS`: supply either a direct `body_scale` or a `body_target`
describing the desired `length` and/or `radius`. Because both dimensions
are addressable, elongated canonical bodies (prolate spheroids,
cylinders) can be reshaped *non-isometrically* - for example lengthening
the body while holding its radius fixed.

## Usage

``` r
# S4 method for class 'GAS'
reforge(
  object,
  body_scale = NULL,
  body_target = NULL,
  isometric_body = TRUE,
  n_segments_body = NULL,
  scale = NULL,
  radius_target = NULL,
  n_segments = NULL
)
```

## Arguments

- object:

  GAS-class object.

- body_scale:

  Proportional scaling for the body length and/or radius. A single value
  scales both dimensions isometrically; otherwise supply a named numeric
  vector using `length` and/or `radius`.

- body_target:

  Target dimensions (m) for the body length and/or radius, supplied as a
  named numeric vector using `length` and/or `radius`.

- isometric_body:

  Logical; maintain isometric scaling for the body.

- n_segments_body:

  New number of segments along the body profile.

- scale:

  Legacy isometric scale factor applied to every axis. Mutually
  exclusive with `radius_target`.

- radius_target:

  Legacy target *maximum* body radius (m); the isometric scale factor is
  derived as `radius_target / max(current_radius)`. Mutually exclusive
  with `scale`.

- n_segments:

  Legacy alias for `n_segments_body`.

## Value

Modified GAS-class object.

## Details

When the body is a recognized canonical family (sphere, prolate
spheroid, or cylinder) the geometry is regenerated from the requested
dimensions through the corresponding shape constructor, which keeps the
discretization clean and the stored shape descriptor accurate. Arbitrary
shapes fall back to direct per-axis scaling of the position matrix. In
every case the returned object is a `GAS`, with the internal-gas
material contrasts, orientation, and metadata preserved.

The legacy `scale`, `radius_target`, and `n_segments` arguments remain
available and are isometric: they scale every axis together, so a sphere
stays a sphere. Use `body_scale`/`body_target` for independent length
and radius control.

## See also

[`reforge()`](https://brandynlucca.github.io/acousticTS/reference/reforge.md)
