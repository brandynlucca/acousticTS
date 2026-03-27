# Calculate Young's modulus (E).

Calculates Young's modulus (E) from two of the three other elastic
moduli: bulk modulus (K), shear modulus (G), or Poisson's ratio
(\\\nu\\). Assumes 3D material properties.

The relationships used are: \$\$E = \frac{9KG}{3K + G}\$\$ \$\$E =
3K(1 - 2\nu)\$\$ \$\$E = 2G(1 + \nu)\$\$

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

Young's modulus (E, Pa).
