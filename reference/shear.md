# Calculate the shear modulus (G)

Calculate the shear modulus (G) from two of the three other elastic
moduli to calculate the Lamé's parameter. When more than two values are
input, the function will default to using the bulk (K) and Young's (E)
moduli. This assumes that the input values represent 3D material
properties.

## Usage

``` r
shear(K = NULL, E = NULL, nu = NULL)
```

## Arguments

- K:

  Bulk modulus (Pa).

- E:

  Young's modulus (Pa).

- nu:

  Poisson's ratio (Dimensionless).

## Value

Returns an estimate for the shear modulus (G).
