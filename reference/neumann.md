# Compute the Neumann factor \\\nu\_{n}\\

The Neumann factor, denoted \\\nu\_{n}\\, is a simple multiplicative
constant commonly used in spherical or spheroidal function theory. It
accounts for the symmetry of cosine terms or the duplication of
even-order contributions in integrals.

## Usage

``` r
neumann(x)
```

## Arguments

- x:

  An integer iterator.

## Value

A numeric vector of the same length as `x`, containing values of
\\\eta_x\\, each equal to 1 or 2.

## Details

Formally, the Neumann factor is defined as: \$\$ \eta_n = \begin{cases}
1, & n = 0, \\ 2, & n \> 0. \end{cases} \$\$

This factor frequently appears in the normalization of spherical Bessel
functions, Legendre expansions, and spheroidal wave functions. It is not
related to the "von Neumann ordinals" used in set theory.

## Examples

``` r
neumann(0) # should return 1
#> [1] 1
neumann(1) # should return 2
#> [1] 2
neumann(2) # should return 2
#> [1] 2

# Vectorized use:
neumann(0:4)
#> [1] 1 2 2 2 2
```
