# Generic scattering shape object used throughout this package.

A S4 class that provides slots to contain relevant shape data and
metadata for a variety of arbitrary and canonical shapes and geometries.
See
[Scatterer](https://brandynlucca.github.io/acousticTS/reference/Scatterer-class.md)
for a more detailed description on how this S4 object interacts with
generic
[Scatterer](https://brandynlucca.github.io/acousticTS/reference/Scatterer-class.md)
objects.

## Value

An object from the `Shape` S4 class hierarchy.

## Slots

- `position_matrix`:

  Position matrix that provides the 2D representation of the body shape

- `shape_parameters`:

  A list of additional shape specifications

## Examples

``` r
methods::getClass("Shape")
#> Class "Shape" [package "acousticTS"]
#> 
#> Slots:
#>                                         
#> Name:   position_matrix shape_parameters
#> Class:           matrix             list
#> 
#> Known Subclasses: 
#> Class "Arbitrary", directly
#> Class "Sphere", directly
#> Class "ProlateSpheroid", directly
#> Class "OblateSpheroid", directly
#> Class "Cylinder", directly
#> Class "PolynomialCylinder", by class "Cylinder", distance 2
```
