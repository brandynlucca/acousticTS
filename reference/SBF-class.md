# Swimbladdered fish (SBF) object/class.

A S4 class that provides slots to contain relevant animal metadata for
parameterizing models for swimbladdered fish (SBF) that are partitioned
into two sets of discretized cylinders: the body and the swimbladder.
Both shapes comprise independent position matrices, material properties,
orientations, and other relevant shape-related data and metadata. See
[Scatterer](https://brandynlucca.github.io/acousticTS/reference/Scatterer-class.md)
for a more detailed description on how this S4 object is organized.

## Value

An object of the `SBF` S4 class.

## See also

[Scatterer](https://brandynlucca.github.io/acousticTS/reference/Scatterer-class.md)

## Examples

``` r
methods::getClass("SBF")
#> Class "SBF" [package "acousticTS"]
#> 
#> Slots:
#>                                                                           
#> Name:          metadata model_parameters            model             body
#> Class:             list             list             list             list
#>                                         
#> Name:           bladder shape_parameters
#> Class:             list             list
#> 
#> Extends: 
#> Class "GAS", directly
#> Class "Scatterer", by class "GAS", distance 2
```
