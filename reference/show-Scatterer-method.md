# Display a scatterer object

Display a scatterer object

## Usage

``` r
# S4 method for class 'Scatterer'
show(object)
```

## Arguments

- object:

  Scattering object.

## Value

Called for its side effect of printing a formatted summary; invisibly
returns `NULL`.

## Examples

``` r
show(cal_generate(material = "WC", n_segments = 20))
#> CAL-object
#>  Calibration sphere
#>  ID:Calibration sphere
#> Material:WC
#>  Sphere longitudinal sound speed:6853m/s
#>  Sphere transversal sound speed:4171m/s
#>  Sphere density:14900kg/m^3
#> Diameter:0.0381 m
#>  Radius:0.01905 m
#> Propagation direction of the incident sound wave:3.142 radians
```
