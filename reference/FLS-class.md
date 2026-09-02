# Fluid-like scatterer (FLS) object/class.

A S4 class that provides slots to contain relevant metadata for
scatterers similar to the surrounding fluid medium (i.e fluid-like)
belonging to FLS-class scatterers. See
[Scatterer](https://brandynlucca.github.io/acousticTS/reference/Scatterer-class.md)
for a more detailed description on how this S4 object is organized.

## Value

An object of the `FLS` S4 class.

## See also

[Scatterer](https://brandynlucca.github.io/acousticTS/reference/Scatterer-class.md)

## Examples

``` r
methods::getClass("FLS")
#> Class "FLS" [package "acousticTS"]
#> 
#> Slots:
#>                                                                           
#> Name:          metadata model_parameters            model             body
#> Class:             list             list             list             list
#>                        
#> Name:  shape_parameters
#> Class:             list
#> 
#> Extends: "Scatterer"
```
