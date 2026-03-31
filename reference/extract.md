# Extract nested components, slots, or matrix/vector fields from objects

Convenience accessor for reaching into `Scatterer` and `Shape` objects
without spelling out direct slot access repeatedly. `extract()` can also
walk through nested lists, matrices, and named vectors, which makes it
useful for pulling model outputs, component properties, or
position-matrix fields from a common interface. This is especially
helpful after geometry manipulations such as
[`brake()`](https://brandynlucca.github.io/acousticTS/reference/brake.md)
and
[`reforge()`](https://brandynlucca.github.io/acousticTS/reference/reforge.md),
where inspecting the stored body profile is often the fastest way to
confirm what changed.

## Usage

``` r
extract(object, feature)
```

## Arguments

- object:

  Scatterer-class object.

- feature:

  Feature(s) of interest (e.g. body). This can either be a scalar
  string, or a vector of names. When a vector is supplied, the function
  recursively accesses the Scatterer object using the 'feature' vector
  as a directory. For example, `feature = c("body", "rpos", "x")` would
  extract the 'x' coordinate of the position matrix ('rpos') from the
  'body' scattering parameters.

## Value

The extracted object, slot, list element, matrix row/column, or vector
element identified by `feature`.

## See also

[`brake()`](https://brandynlucca.github.io/acousticTS/reference/brake.md),
[`reforge()`](https://brandynlucca.github.io/acousticTS/reference/reforge.md),
[`translate_shape()`](https://brandynlucca.github.io/acousticTS/reference/translate_shape.md),
[`reanchor_shape()`](https://brandynlucca.github.io/acousticTS/reference/reanchor_shape.md),
[`inflate_shape()`](https://brandynlucca.github.io/acousticTS/reference/inflate_shape.md),
[`smooth_shape()`](https://brandynlucca.github.io/acousticTS/reference/smooth_shape.md),
[`resample_shape()`](https://brandynlucca.github.io/acousticTS/reference/resample_shape.md),
[`flip_shape()`](https://brandynlucca.github.io/acousticTS/reference/flip_shape.md),
[`offset_component()`](https://brandynlucca.github.io/acousticTS/reference/offset_component.md)

## Examples

``` r
obj <- fls_generate(
  shape = sphere(radius_body = 0.01, n_segments = 40),
  density_body = 1045,
  sound_speed_body = 1520
)

extract(obj, "body")
#> $rpos
#>    [,1]         [,2]         [,3]         [,4]   [,5]         [,6]         [,7]
#> x  0.02  0.019500000  0.019000000  0.018500000  0.018  0.017500000  0.017000000
#> y  0.00  0.000000000  0.000000000  0.000000000  0.000  0.000000000  0.000000000
#> z  0.00  0.000000000  0.000000000  0.000000000  0.000  0.000000000  0.000000000
#> zU 0.00  0.003122499  0.004358899  0.005267827  0.006  0.006614378  0.007141428
#> zL 0.00 -0.003122499 -0.004358899 -0.005267827 -0.006 -0.006614378 -0.007141428
#>            [,8]   [,9]        [,10]        [,11]        [,12]        [,13]
#> x   0.016500000  0.016  0.015500000  0.015000000  0.014500000  0.014000000
#> y   0.000000000  0.000  0.000000000  0.000000000  0.000000000  0.000000000
#> z   0.000000000  0.000  0.000000000  0.000000000  0.000000000  0.000000000
#> zU  0.007599342  0.008  0.008351647  0.008660254  0.008930286  0.009165151
#> zL -0.007599342 -0.008 -0.008351647 -0.008660254 -0.008930286 -0.009165151
#>           [,14]        [,15]        [,16]        [,17]       [,18]        [,19]
#> x   0.013500000  0.013000000  0.012500000  0.012000000  0.01150000  0.011000000
#> y   0.000000000  0.000000000  0.000000000  0.000000000  0.00000000  0.000000000
#> z   0.000000000  0.000000000  0.000000000  0.000000000  0.00000000  0.000000000
#> zU  0.009367497  0.009539392  0.009682458  0.009797959  0.00988686  0.009949874
#> zL -0.009367497 -0.009539392 -0.009682458 -0.009797959 -0.00988686 -0.009949874
#>           [,20] [,21]        [,22]        [,23]       [,24]        [,25]
#> x   0.010500000  0.01  0.009500000  0.009000000  0.00850000  0.008000000
#> y   0.000000000  0.00  0.000000000  0.000000000  0.00000000  0.000000000
#> z   0.000000000  0.00  0.000000000  0.000000000  0.00000000  0.000000000
#> zU  0.009987492  0.01  0.009987492  0.009949874  0.00988686  0.009797959
#> zL -0.009987492 -0.01 -0.009987492 -0.009949874 -0.00988686 -0.009797959
#>           [,26]        [,27]        [,28]        [,29]        [,30]
#> x   0.007500000  0.007000000  0.006500000  0.006000000  0.005500000
#> y   0.000000000  0.000000000  0.000000000  0.000000000  0.000000000
#> z   0.000000000  0.000000000  0.000000000  0.000000000  0.000000000
#> zU  0.009682458  0.009539392  0.009367497  0.009165151  0.008930286
#> zL -0.009682458 -0.009539392 -0.009367497 -0.009165151 -0.008930286
#>           [,31]        [,32]  [,33]        [,34]        [,35]        [,36]
#> x   0.005000000  0.004500000  0.004  0.003500000  0.003000000  0.002500000
#> y   0.000000000  0.000000000  0.000  0.000000000  0.000000000  0.000000000
#> z   0.000000000  0.000000000  0.000  0.000000000  0.000000000  0.000000000
#> zU  0.008660254  0.008351647  0.008  0.007599342  0.007141428  0.006614378
#> zL -0.008660254 -0.008351647 -0.008 -0.007599342 -0.007141428 -0.006614378
#>     [,37]        [,38]        [,39]        [,40] [,41]
#> x   0.002  0.001500000  0.001000000  0.000500000     0
#> y   0.000  0.000000000  0.000000000  0.000000000     0
#> z   0.000  0.000000000  0.000000000  0.000000000     0
#> zU  0.006  0.005267827  0.004358899  0.003122499     0
#> zL -0.006 -0.005267827 -0.004358899 -0.003122499     0
#> 
#> $radius
#> [1] 0.01
#> 
#> $theta
#> [1] 1.570796
#> 
#> $g
#> NULL
#> 
#> $h
#> NULL
#> 
#> $density
#> [1] 1045
#> 
#> $sound_speed
#> [1] 1520
#> 
#> $radius_curvature_ratio
#> NULL
#> 
extract(obj, c("body", "density"))
#> [1] 1045
extract(obj, c("shape_parameters", "shape"))
#> [1] "Sphere"

bent_obj <- brake(obj, radius_curvature = 5)
head(extract(bent_obj, c("body", "rpos", "z")))
#> [1] -0.0004995835 -0.0004509107 -0.0004047267 -0.0003610325 -0.0003198294
#> [6] -0.0002811182
extract(bent_obj, c("shape_parameters", "radius_curvature_ratio"))
#> [1] 5
```
