# Solid and calibration sphere (CAL) object/class.

A S4 class that provides slots to contain relevant metadata for solid
sphere objects belonging to the CAL-class scatterers. This object is
created using parameters specific to the outer shell. The default
behavior of this object is to only reference these outer elastic shell
properties with few exceptions that are model-dependent. `CAL` inherits
from the broader elastic-based
[ELA](https://brandynlucca.github.io/acousticTS/reference/ELA-class.md)
class rather than from
[ESS](https://brandynlucca.github.io/acousticTS/reference/ESS-class.md),
because calibration targets are solid elastic bodies rather than elastic
shells. See
[Scatterer](https://brandynlucca.github.io/acousticTS/reference/Scatterer-class.md)
for a more detailed description on how this S4 object is organized.

## Value

An object of the `CAL` S4 class.

## See also

[Scatterer](https://brandynlucca.github.io/acousticTS/reference/Scatterer-class.md)

## Examples

``` r
methods::getClass("CAL")
#> Class "CAL" [package "acousticTS"]
#> 
#> Slots:
#>                                                                           
#> Name:             model             body shape_parameters         metadata
#> Class:             list             list             list             list
#>                        
#> Name:  model_parameters
#> Class:             list
#> 
#> Extends: 
#> Class "ELA", directly
#> Class "Scatterer", by class "ELA", distance 2
```
