# Generic gas-filled scatterer (GAS) object/class.

A S4 class that provides slots to contain relevant metadata for
gas-bearing scatterers belonging to the GAS-class. This object can
include simple gas-filled bubbles to other scatterers with gas
occlusions, swimbladders, and other internal features, if applicable.
The default behavior for this type of object is to only reference the
gaseous/fluid feature with exceptions that are model-dependent. See
[Scatterer](https://brandynlucca.github.io/acousticTS/reference/Scatterer-class.md)
for a more detailed description on how this S4 object is organized.

## Value

An object of the `GAS` S4 class.

## See also

[`Scatterer`](https://brandynlucca.github.io/acousticTS/reference/Scatterer-class.md)

## Examples

``` r
methods::getClass("GAS")
#> Class "GAS" [package "acousticTS"]
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
#> 
#> Known Subclasses: "SBF"
```
