# Calculate Young's modulus (E).

Calculate the Young's modulus (E) from two of the three other elastic
moduli to calculate the Lamé's parameter. When more than two values are
input, the function will default to using the bulk (K) and shear (G)
moduli. This assumes that the input values represent 3D material
properties.

## Usage

``` r
young(K = NULL, G = NULL, nu = NULL)
```

## Arguments

- K:

  Bulk modulus (Pa).

- G:

  Shear modulus (Pa).

- nu:

  Poisson's ratio (Dimensionless).

## Value

Returns an estimate for the Young's modulus (E).
