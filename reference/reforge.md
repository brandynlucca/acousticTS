# Resize or reparameterize a scatterer object

Generic function for rescaling or otherwise reparameterizing an existing
scatterer while preserving its class semantics. The class-specific
methods handle the component bookkeeping needed to keep the resulting
object structurally valid after lengths, widths, heights, shell
thicknesses, or related geometry descriptors are changed.

## Usage

``` r
reforge(object, ...)
```

## Arguments

- object:

  A scatterer object.

- ...:

  Additional arguments passed to specific methods.

## Value

A scatterer of the same broad class as `object`, rebuilt with the
requested geometric changes.

## Details

`reforge()` is the package's main post-construction geometry-adjustment
tool. It is useful when a target should be modified in place rather than
rebuilt from scratch. The available method-specific arguments depend on
the scatterer class and typically include direct scale factors and/or
target dimensions for length, width, or height.

## See also

[`extract()`](https://brandynlucca.github.io/acousticTS/reference/extract.md),
[`plot.Scatterer()`](https://brandynlucca.github.io/acousticTS/reference/plot.Scatterer.md),
[`brake()`](https://brandynlucca.github.io/acousticTS/reference/brake.md)

## Examples

``` r
obj <- fls_generate(
  shape = sphere(radius_body = 0.01, n_segments = 40),
  density_body = 1045,
  sound_speed_body = 1520
)

bigger_obj <- reforge(obj, length = 0.03)
extract(bigger_obj, c("shape_parameters", "length"))
#> [1] 0.03
```
