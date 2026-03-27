# Reforge CAL-class object

Resize a calibration sphere by applying an isometric scale factor or
specifying a target diameter. Optionally re-discretize to a new segment
count. CAL objects are always spheres, so the position matrix follows
the
[`sphere()`](https://brandynlucca.github.io/acousticTS/reference/sphere.md)
convention (n_points x 5: x, y, z, zU, zL).

## Usage

``` r
# S4 method for class 'CAL'
reforge(object, scale = NULL, diameter_target = NULL, n_segments = NULL)
```

## Arguments

- object:

  CAL-class object.

- scale:

  Single positive scale factor applied isometrically. Mutually exclusive
  with `diameter_target`.

- diameter_target:

  Target sphere diameter (m). Derives the scale factor internally.
  Mutually exclusive with `scale`.

- n_segments:

  New number of discrete segments along the major axis.

## Value

Modified CAL-class object.
