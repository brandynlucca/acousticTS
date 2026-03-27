# Calculate Lamé's first parameter (\\\lambda\\)

Calculates Lamé's first parameter (\\\lambda\\) from two of the four
other elastic moduli: bulk modulus (K), Young's modulus (E), shear
modulus (G), or Poisson's ratio (\\\nu\\). Assumes 3D material
properties.

The relationships used are: \$\$\lambda = K - \frac{2G}{3}\$\$
\$\$\lambda = \frac{E\nu}{(1 + \nu)(1 - 2\nu)}\$\$ \$\$\lambda =
\frac{2G\nu}{1 - 2\nu}\$\$ \$\$\lambda = \frac{3K\nu}{1 + \nu}\$\$
\$\$\lambda = \frac{3K(3K - E)}{9K - E}\$\$ \$\$\lambda = \frac{G(E -
2G)}{3G - E}\$\$

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

Lamé's first parameter (\\\lambda\\, Pa).
