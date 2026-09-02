# Elastic shelled scatterer (ESS) object/class.

A S4 class that provides slots to contain relevant metadata for elastic
shelled scatterers/objects belonging to the ESS-class. This object can
be created using values for both an outer shell and internal tissues, if
applicable. The default behavior for this type of this object is to only
reference the outer shell with few exceptions that are model-dependent.
`ESS` inherits from the broader elastic-based
[ELA](https://brandynlucca.github.io/acousticTS/reference/ELA-class.md)
class. See
[Scatterer](https://brandynlucca.github.io/acousticTS/reference/Scatterer-class.md)
for a more detailed description on how this S4 object is organized.

## Value

An object of the `ESS` S4 class.

## See also

[Scatterer](https://brandynlucca.github.io/acousticTS/reference/Scatterer-class.md)

## Examples

``` r
methods::getClass("ESS")
#> Class "ESS" [package "acousticTS"]
#> 
#> Slots:
#>                                                                           
#> Name:             shell            fluid            model             body
#> Class:             list             list             list             list
#>                                                          
#> Name:  shape_parameters         metadata model_parameters
#> Class:             list             list             list
#> 
#> Extends: 
#> Class "ELA", directly
#> Class "Scatterer", by class "ELA", distance 2
```
