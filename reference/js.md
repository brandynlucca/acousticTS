# Spherical Bessel function of the first kind, \\j\_\nu(z)\\, and its respective derivatives

Computes the spherical Bessel function of the first kind (\\j\_\nu(z)\\)
and its k-th derivative.

## Usage

``` r
js(l, n)

jsdk(l, n, k)
```

## Arguments

- l:

  Numeric. The order of the spherical Bessel function. Can be integer or
  fractional.

- n:

  Numeric or matrix. The argument (\\z\\) at which to evaluate the
  function. If a matrix is provided, the function is applied
  column-wise.

- k:

  Non-negative integer. The order of the derivative for `jsdk`.

## Value

A numeric vector or matrix (matching the input structure) containing:

- `js`: \\j\_\nu(z)\\

- `jsdk`: \\j^{(k)'}\_l(z)\\ (k-th derivative)

## Details

The spherical Bessel function of the first kind is related to the
cylindrical Bessel function by: \$\$j\_\nu(z) = \sqrt{\frac{\pi}{2z}}
J\_{\nu+1/2}(z)\$\$

where \\J\_\nu(z)\\ is the cylindrical Bessel function of the first
kind.

The spherical Bessel functions satisfy the differential equation: \$\$
z^2 \frac{d^2 j\_\nu}{dz^2} + 2z \frac{dj\_\nu}{dz} + \[z^2 - l(l+1)\]
j\_\nu = 0 \$\$

**Special cases:**

- \\j\_\nu(0) = 0\\ for all \\\nu\\.

- \\j_0(z) = \frac{\sin(z)}{z}\\

- \\j_1(z) = \frac{\sin(z)}{z^2} - \frac{\cos(z)}{z}\\

**Derivatives:**

- First derivative: \\ j'\_\nu(z) = j\_{\nu-1}(z) - \frac{\nu+1}{z}
  j\_\nu(z) \\

- Second derivative: \\j''\_\nu(z) = \frac{(\nu+1)(\nu+2) - z^2}{z^2}
  j\_\nu(z) - \frac{2}{z} j\_{\nu-1}(z)\\

## References

Abramowitz, M. and Stegun, I.A. (Eds.). (1964). *Handbook of
Mathematical Functions with Formulas, Graphs, and Mathematical Tables*.
National Bureau of Standards, Applied Mathematics Series 55. Chapter 10.

NIST Digital Library of Mathematical Functions. <https://dlmf.nist.gov/>

- Spherical Bessel functions: <https://dlmf.nist.gov/10.47>

- Relation to cylindrical Bessel functions: Eq. 10.47.3 at
  <https://dlmf.nist.gov/10.47>

## See also

[`jc`](https://brandynlucca.github.io/acousticTS/reference/jc.md) for
cylindrical Bessel functions of the first kind,
[`ys`](https://brandynlucca.github.io/acousticTS/reference/ys.md) for
spherical Bessel functions of the second kind,
[`hs`](https://brandynlucca.github.io/acousticTS/reference/hs.md) for
spherical Hankel functions.

## Examples

``` r
# Spherical Bessel function
js(0, 1)
#> [1] 0.841471
js(1, 2.5)
#> [1] 0.416213

# Fractional order
js(0.5, 3)
#> [1] 0.04704

# Vector input
js(0, c(1, 2, 3))
#> [1] 0.8414710 0.4546487 0.0470400

# Matrix input (applied column-wise)
js(1, matrix(1:6, nrow = 2))
#> [1]  0.30116868  0.43539777  0.34567750  0.11611075 -0.09508941 -0.16778992

# First derivative
jsdk(1, 2, 1)
#> [1] 0.01925094

# Second derivative
jsdk(1, 2, 2)
#> [1] -0.2369498
```
