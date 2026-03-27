# Calculate the bulk modulus (K).

Calculates the bulk modulus (K) from two of the three other elastic
moduli: Young's modulus (E), shear modulus (G), or Poisson's ratio
(\\\nu\\). Assumes 3D material properties.

The relationships used are: \$\$K = \frac{E G}{3(3G - E)}\$\$ \$\$K =
\frac{2G(1 + \nu)}{3(1 - 2\nu)}\$\$ \$\$K = \frac{E}{3(1 - 2\nu)}\$\$

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

Bulk modulus (K, Pa).
