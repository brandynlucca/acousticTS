# Calculate Lamé's first parameter (\\\lambda\\)

Calculate Lamé's first parameter (\\\lambda\\) from two of the four
other elastic moduli. When more than two values are input, the function
will default to using the bulk (K) and Young's (E) moduli. This assumes
that the input values represent 3D material properties.

## Usage

``` r
lame(K = NULL, E = NULL, G = NULL, nu = NULL)
```

## Arguments

- K:

  Bulk modulus (Pa).

- E:

  Young's modulus (Pa).

- G:

  Shear modulus (Pa).

- nu:

  Poisson's ratio (Dimensionless).

## Value

Returns Lamé's first parameter (\\\lambda\\).
