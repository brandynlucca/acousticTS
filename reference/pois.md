# Calculate the Poisson's ratio (\\\nu\\)

Calculate the Poisson's ratio (\\\nu\\) from two of the three other
elastic moduli to calculate the Lamé's parameter. When more than two
values are input, the function will default to using the bulk (K) and
Young's (E) moduli. This assumes that the input values represent 3D
material properties.

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

Returns a dimensionless ratio known as Poisson's ratio (\\\nu\\).
