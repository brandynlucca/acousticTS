# Legendre Polynomial of the First Kind, \\P\_\nu(x)\\

Computes the Legendre polynomial of the first kind, \\P\_\nu(x)\\, for
real order \\\nu\\ (integer or fractional) and real argument \\x\\.

## Usage

``` r
Pn(n, x)
```

## Arguments

- n:

  Numeric vector. Degree (order) of the Legendre polynomial. Can be
  integer or fractional (e.g., 0, 1, 2.5, 3.7). Can also be negative.

- x:

  Numeric vector. Real argument(s) at which to evaluate the polynomial.
  Valid for all real values including \\\|x\| \> 1\\.

## Value

A numeric matrix of dimension `length(n)` by `length(x)`, where element
`[i, j]` contains \\P\_{n_i}(x_j)\\.

## Details

The Legendre polynomial of the first kind satisfies the differential
equation: \$\$(1 - x^2) \frac{d^2 P\_\nu}{dx^2} - 2x
\frac{dP\_\nu}{dx} + \nu(\nu + 1) P\_\nu = 0\$\$

For integer order \\n\\, the function uses the recurrence relation:
\$\$P_0(x) = 1\$\$ \$\$P_1(x) = x\$\$ \$\$P_n(x) = \frac{(2n-1) x
P\_{n-1}(x) - (n-1) P\_{n-2}(x)}{n}\$\$

For fractional order \\\nu\\, the function uses:

- For \\\|x\| \leq 1\\: The Ferrers function via the hypergeometric
  series \\P\_\nu(x) = {}\_2F_1(-\nu, \nu+1; 1; \frac{1-x}{2})\\

- For \\\|x\| \> 1\\: A numerical contour integral representation

## Note

This function calls underlying \\C++\\ code via `Rcpp` for computational
efficiency and to support different cases for both order and argument
that are not readily available in `R`.

## References

Abramowitz, M. and Stegun, I. A. (1972). *Handbook of Mathematical
Functions with Formulas, Graphs, and Mathematical Tables*. Dover
Publications. Chapter 8: Legendre Functions.

NIST Digital Library of Mathematical Functions.
<https://dlmf.nist.gov/14>

## See also

[`Pndk`](https://brandynlucca.github.io/acousticTS/reference/Pndk.md)
for the k^(th) derivative of the Legendre polynomial of the first kind,
[`Qn`](https://brandynlucca.github.io/acousticTS/reference/Qn.md) for
Legendre functions of the second kind.

## Examples

``` r
# Single values
Pn(1, 1)
#>      [,1]
#> [1,]    1

# Multiple orders, single argument
Pn(c(1, 2, 3), 1)
#>      [,1]
#> [1,]    1
#> [2,]    1
#> [3,]    1

# Single order, multiple arguments
Pn(1, c(1, 2, 3))
#>      [,1] [,2] [,3]
#> [1,]    1    2    3

# Multiple orders and arguments (returns a matrix)
Pn(c(1, 2, 3), c(1, 2, 3))
#>      [,1] [,2] [,3]
#> [1,]    1  2.0    3
#> [2,]    1  5.5   13
#> [3,]    1 17.0   63

# Fractional orders
Pn(c(0.5, 2.2, 3), c(1, 2, 3))
#>      [,1]      [,2]      [,3]
#> [1,]    1  1.329138  1.597387
#> [2,]    1  6.849800 17.720229
#> [3,]    1 17.000000 63.000000

# Negative arguments
Pn(c(0.5, 2, -3), c(1, -2.5, 3))
#>      [,1]         [,2]      [,3]
#> [1,]    1 8.994435e-17  1.597387
#> [2,]    1 8.875000e+00 13.000000
#> [3,]    0 0.000000e+00  0.000000
```
