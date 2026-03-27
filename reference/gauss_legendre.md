# Gauss-Legendre nodes and weights

Compute Gauss-Legendre quadrature nodes and weights on an interval
\\\[a,~b\]\\.

## Usage

``` r
gauss_legendre(n, a = -1, b = 1)
```

## Arguments

- n:

  Number of quadrature nodes (n \>= 1).

- a:

  Left endpoint of the integration interval.

- b:

  Right endpoint of the integration interval (b \> a).

## Value

A list with components:

- nodes:

  Quadrature abscissae \\x_i\\ in \\\[a,~b\]\\.

- weights:

  Quadrature weights \\w_i\\ such that \\\int_a^b f(x)\\dx \approx
  \sum\_{i=1}^n w_i\\f(x_i).\\

## Details

Gauss-Legendre quadrature provides exact integration for polynomials of
degree up to \\2n-1\\ using n nodes and weights chosen as the roots of
the Legendre polynomial \\P_n(x)\\ on the canonical interval
\\\[-1,1\]\\. For a general interval \\\[a,b\]\\ the canonical nodes are
shifted and rescaled onto \\\[a,b\]\\, and the weights are scaled by
\\(b-a)/2\\.

This wrapper performs basic argument validation and calls the C++
routine to obtain nodes and weights with high accuracy for moderate `n`.

## References

Davis, P. J., & Rabinowitz, P. (2007). Methods of Numerical Integration
(2nd ed.).
