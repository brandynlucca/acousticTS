# Derivative of the Legendre Polynomial of the First Kind

Computes the \\k\\-th derivative of the Legendre polynomial of the first
kind, \\\frac{d^k}{dx^k} P\_\nu(x)\\, with respect to the argument
\\x\\.

## Usage

``` r
Pndk(n, x, k = 1L)
```

## Arguments

- n:

  Numeric vector. Degree (order) of the Legendre polynomial. Can be
  integer or fractional.

- x:

  Numeric vector. Real argument(s) at which to evaluate the derivative.

- k:

  Integer. Order of the derivative (\\k \geq 0\\). Default is 1 for the
  first derivative.

## Value

A numeric matrix of dimension `length(n)` by `length(x)`, where element
`[i, j]` contains \\\frac{d^k}{dx^k} P\_{n_i}(x_j)\\.

## Details

**For integer order \\n\\:**

The derivative is computed using the relationship with associated
Legendre polynomials: \$\$\frac{d^m P_n(x)}{dx^m} =
\frac{1}{(1-x^2)^{m/2}} P_n^m(x)\$\$

where \\P_n^m(x)\\ is the associated Legendre polynomial (the `Boost`
`C++` uses Condon-Shortley phase convention).

At the endpoints \\x = \pm 1\\, the known closed-form expressions are
used: \$\$P_n^{(k)}(1) = \frac{1}{2^k k!} \prod\_{j=0}^{k-1}
(n-j)(n+1+j)\$\$ \$\$P_n^{(k)}(-1) = (-1)^{n+k} P_n^{(k)}(1)\$\$

If \\k \> n\\ for integer \\n\\, the result is 0 (derivative of a
polynomial of degree \\n\\ taken more than \\n\\ times).

**For fractional order \\\nu\\:**

Derivatives are computed using central finite differences:
\$\$\frac{dP\_\nu}{dx} \approx \frac{P\_\nu(x+h) - P\_\nu(x-h)}{2h}\$\$

Higher-order derivatives use the generalized finite difference stencil.
Note that accuracy may be limited for fractional orders due to the
numerical integration underlying
[`Pn`](https://brandynlucca.github.io/acousticTS/reference/Pn.md).

## Note

For fractional orders, the finite difference approximation may have
reduced accuracy (typically 4-6 significant digits) compared to integer
orders. A warning is issued for higher-order derivatives (\\k \> 1\\) of
fractional orders.

This function calls underlying \\C++\\ code via `Rcpp` for computational
efficiency and to support different cases for both order and argument
that are not readily available in `R`.

## References

Abramowitz, M. and Stegun, I. A. (1972). *Handbook of Mathematical
Functions*. Dover Publications. Section 8.5: Associated Legendre
Functions.

NIST Digital Library of Mathematical Functions.
<https://dlmf.nist.gov/14.10>

## See also

[`Pn`](https://brandynlucca.github.io/acousticTS/reference/Pn.md) for
the Legendre polynomial of the first kind.

## Examples

``` r
# First derivative of P_2(x) at x = 0.5
# P_2(x) = (3x^2 - 1)/2, so P'_2(x) = 3x, P'_2(0.5) = 1.5
Pndk(2, 0.5, 1)
#>      [,1]
#> [1,]  1.5

# Second derivative of P_2(x)
# P''_2(x) = 3
Pndk(2, 0.5, 2)
#>      [,1]
#> [1,]    3

# Multiple orders
Pndk(c(1, 2, 3), 0.5, 1)
#>       [,1]
#> [1,] 1.000
#> [2,] 1.500
#> [3,] 0.375

# First derivative at multiple points
Pndk(2, c(-0.5, 0, 0.5), 1)
#>      [,1] [,2] [,3]
#> [1,] -1.5    0  1.5

# Fractional order (uses finite differences)
Pndk(0.5, 0.5, 1)
#>           [,1]
#> [1,] 0.4503717

# Derivative at endpoint
# P'_n(1) = n(n+1)/2
Pndk(3, 1, 1)  # Should be 3*4/2 = 6
#>      [,1]
#> [1,]    6
```
