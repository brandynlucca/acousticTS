# Calculate she compressibility (\\\kappa\\) of a scattering boundary/interface.

Calculates the compressibility contrast (\\\kappa\\) between a
scattering interface and the surrounding medium. Compressibility is
defined as: \$\$ K = \frac{1}{\rho c^2} \$\$ where \\\rho\\ is density
(\\kg~m^{-3}\\) and \\c\\ is sound speed (\\m~s^{-1}\\).

The compressibility contrast is then: \$\$ \kappa = \frac{K_2 -
K_1}{K_1} \$\$

where \\K_1\\ is the compressibility of the medium and \\K_2\\ is that
of the target interface.

## Usage

``` r
compressibility(medium, target)
```

## Arguments

- medium:

  Dataframe object containing density (\\kgm^{-3}\\) and sound speed
  (\\ms^{-1}\\) values for a fluid medium external to a scattering
  interface (e.g., seawater).

- target:

  Dataframe object containing density (\\kgm^{-3}\\) and sound speed
  (\\ms^{-1}\\) values for a target boundary.

## Value

Compressibility contrast (\\\kappa\\), dimensionless.
