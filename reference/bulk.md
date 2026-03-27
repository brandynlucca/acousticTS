# Calculate the bulk modulus (K).

Calculate the bulk modulus (K) from two of the three other elastic
moduli to calculate the Lamé's parameter. When more than two values are
input, the function will default to using Young's (E) and shear (G)
moduli. This assumes that the input values represent 3D material
properties.

## Usage

``` r
bulk(E = NULL, G = NULL, nu = NULL)
```

## Arguments

- E:

  Young's modulus (Pa).

- G:

  Shear modulus (Pa).

- nu:

  Poisson's ratio (Dimensionless).

## Value

Returns an estimate for the bulk modulus (K).
