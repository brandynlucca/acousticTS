# Calculate the Poisson's ratio (\\\nu\\)

Calculates Poisson's ratio (\\\nu\\) from two of the three other elastic
moduli: bulk modulus (K), Young's modulus (E), or shear modulus (G).
Assumes 3D material properties.

The relationships used are: \$\$\nu = \frac{E}{2G} - 1\$\$ \$\$\nu =
\frac{3K - 2G}{2(3K + G)}\$\$ \$\$\nu = \frac{3K - E}{6K}\$\$

## Usage

``` r
pois(K = NULL, E = NULL, G = NULL)
```

## Arguments

- K:

  Bulk modulus (K, Pa).

- E:

  Young's modulus (E, Pa).

- G:

  Shear modulus (Pa).

## Value

Poisson's ratio (\\\nu\\), dimensionless.
