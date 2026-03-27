# Calculate the shear modulus (G)

\#' Calculates the shear modulus (G) from two of the three other elastic
moduli: bulk modulus (K), Young's modulus (E), or Poisson's ratio
(\\\nu\\). Assumes 3D material properties.

The relationships used are: \$\$G = \frac{3KE}{9K - E}\$\$ \$\$G =
\frac{3K(1 - 2\nu)}{2(1 + \nu)}\$\$ \$\$G = \frac{E}{2(1 + \nu)}\$\$

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

Shear modulus (G, Pa).
