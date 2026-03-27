# Reforge GAS-class object

Resize a gas-filled scatterer by applying an isometric scale factor or
specifying a target maximum radius. Optionally re-discretize the body
representation to a new segment count. The underlying shape (sphere,
prolate spheroid, cylinder, arbitrary, etc.) is preserved; scaling is
applied uniformly to all axes of the position matrix.

## Usage

``` r
# S4 method for class 'GAS'
reforge(object, scale = NULL, radius_target = NULL, n_segments = NULL)
```

## Arguments

- object:

  GAS-class object.

- scale:

  Single positive scalar applied isometrically to every axis of the
  position matrix. Mutually exclusive with `radius_target`.

- radius_target:

  Target *maximum* body radius (m). The scale factor is derived as
  `radius_target / max(current_radius)`. Mutually exclusive with
  `scale`.

- n_segments:

  New number of discrete segments. All position-matrix columns are
  re-interpolated along the x-axis.

## Value

Modified GAS-class object.
