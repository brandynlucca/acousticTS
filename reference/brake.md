# Support function for bending scatterer body shape and position matrix

Support function for bending scatterer body shape and position matrix

## Usage

``` r
brake(object, radius_curvature, mode = "ratio")
```

## Arguments

- object:

  Dataframe or scatterer-class object

- radius_curvature:

  Radius of curvature that can be parameterized either as a ratio
  relative to body length or actual measurement

- mode:

  Either "ratio" or "measurement"

## Value

A bent version of `object`, returned as the same broad object type with
updated geometry and curvature metadata.
