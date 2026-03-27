# Calculate the density contrast (\\\rho\\) of a scattering boundary

Calculate the density contrast (\\\rho\\) of a scattering boundary

## Usage

``` r
rho(medium, target)
```

## Arguments

- medium:

  Dataframe object containing density (\\kgm^{-3}\\) and sound speed
  (\\ms^{-1}\\) values for the external medium.

- target:

  Dataframe object containing density (\\kgm^{-3}\\) and sound speed
  (\\ms^{-1}\\) values for the target boundary.

## Value

Density contrast, defined as \\(\rho\_{target} -
\rho\_{medium})/\rho\_{target}\\.
