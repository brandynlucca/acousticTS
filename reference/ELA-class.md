# Elastic-based scatterer (ELA) object/class.

A shared parent S4 class for elastic-based scatterers. The `ELA`-class
collects the slots common to solid elastic targets such as calibration
spheres
([CAL](https://brandynlucca.github.io/acousticTS/reference/CAL-class.md))
and elastic-shelled targets
([ESS](https://brandynlucca.github.io/acousticTS/reference/ESS-class.md)),
without implying a particular internal structure such as a shell or
homogeneous elastic solid. See
[Scatterer](https://brandynlucca.github.io/acousticTS/reference/Scatterer-class.md)
for a more detailed description of the shared object organization.

## Value

An object of the `ELA` S4 class.

## See also

[Scatterer](https://brandynlucca.github.io/acousticTS/reference/Scatterer-class.md),
[CAL](https://brandynlucca.github.io/acousticTS/reference/CAL-class.md),
[ESS](https://brandynlucca.github.io/acousticTS/reference/ESS-class.md)

## Examples

``` r
methods::getClass("ELA")
#> Class "ELA" [package "acousticTS"]
#> 
#> Slots:
#>                                                                           
#> Name:             model             body shape_parameters         metadata
#> Class:             list             list             list             list
#>                        
#> Name:  model_parameters
#> Class:             list
#> 
#> Extends: "Scatterer"
#> 
#> Known Subclasses: "ESS", "CAL"
```
