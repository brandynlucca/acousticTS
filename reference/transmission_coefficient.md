# Plane wave/plane interface transmission coefficient

Plane wave/plane interface transmission coefficient

## Usage

``` r
transmission_coefficient(interface1, interface2, mode = "DWBA")
```

## Arguments

- interface1:

  Dataframe object containing density (kg/m^3) and sound speed (m/s)
  values for a boundary/interface (1)

- interface2:

  Dataframe object containing density (kg/m^3) and sound speed (m/s)
  values for a boundary/interface (2)

- mode:

  Two options: coefficient calculation for "DWBA" and "KRM"

## Value

Pressure-amplitude transmission coefficient at normal incidence.

## Examples

``` r
seawater <- data.frame(density = 1026, sound_speed = 1480)
animal <- data.frame(density = 1050, sound_speed = 1530)
transmission_coefficient(seawater, animal)
#> [1] 1.028167
```
