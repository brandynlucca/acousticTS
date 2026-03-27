# Derivative of the Legendre Function of the Second Kind

Computes the \\k\\-th derivative of the Legendre function of the second
kind, \\\frac{d^k}{dx^k} Q\_\nu(x)\\, with respect to the argument
\\x\\.

## Usage

``` r
Qndk(n, x, k = 1L)
```

## Arguments

- n:

  Numeric vector. Degree (order) of the Legendre function.

- x:

  Numeric vector. Real argument(s) at which to evaluate the derivative.
  Avoid \\x = \pm 1\\.

- k:

  Integer. Order of the derivative (\\k \geq 0\\). Default is 1.

## Value

A complex matrix of dimension `length(n)` by `length(x)`, where element
`[i, j]` contains \\\frac{d^k}{dx^k} Q\_{n_i}(x_j)\\.

## Details

Derivatives are computed using central finite differences:
\$\$\frac{dQ\_\nu}{dx} \approx \frac{Q\_\nu(x+h) - Q\_\nu(x-h)}{2h}\$\$

Higher-order derivatives use the generalized finite difference stencil
with binomial coefficients.

For \\\|x\| \> 1\\, the result is complex since \\Q\_\nu(x)\\ is complex
in that domain.

**Singularities:** The function \\Q\_\nu(x)\\ has logarithmic
singularities at \\x = \pm 1\\. Derivatives near these points will have
reduced accuracy or may be undefined.

## Note

Derivatives are computed via finite differences with step size \\h =
10^{-6}\\. Accuracy is typically 4-6 significant digits for first
derivatives, less for higher orders.

This function calls underlying \\C++\\ code via `Rcpp` for computational
efficiency and to support different cases for both order and argument
that are not readily available in `R`.

## References

Abramowitz, M. and Stegun, I. A. (1972). *Handbook of Mathematical
Functions*. Dover Publications. Chapter 8: Legendre Functions.

## See also

[`Qn`](https://brandynlucca.github.io/acousticTS/reference/Qn.md) for
Legendre functions of the second kind.

## Examples

``` r
# First derivative of Q_1(x) at x = 0.5
Qndk(1, 0.5, 1)
#>             [,1]
#> [1,] 1.215973+0i

# Compare with numerical derivative
h <- 1e-6
(Qn(1, 0.5 + h) - Qn(1, 0.5 - h)) / (2 * h)
#>             [,1]
#> [1,] 1.215973+0i

# Second derivative
Qndk(2, 0.5, 2)
#>            [,1]
#> [1,] 5.42566+0i

# Multiple orders
Qndk(c(1, 2, 3), 0.5, 1)
#>               [,1]
#> [1,]  1.2159728+0i
#> [2,] -0.8427075+0i
#> [3,] -2.8773435+0i

# Complex result for |x| > 1
Qndk(1, 2.0, 1)
#>                      [,1]
#> [1,] -0.1173605-1.570796i
```
