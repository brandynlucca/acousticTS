# Spherical Bessel function of the second kind, \\y\_\nu(z)\\, and its respective derivatives

Computes the spherical Bessel function of the second kind
(\\y\_\nu(z)\\), also known as the spherical Neumann function, and its
k-th derivatives.

## Usage

``` r
ys(l, n)

ysdk(l, n, k)
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

  Non-negative integer. The order of the derivative for `ysdk`.

## Value

A numeric vector or matrix (matching the input structure) containing:

- `ys`: \\y\_\nu(z)\\

- `ysd`: \\y^{(k)'}\_l(z)\\ (k-th derivative)

## Details

The spherical Bessel function of the second kind is related to the
cylindrical Bessel function by:

\$\$y\_\nu(z) = \sqrt{\frac{\pi}{2z}} Y\_{\nu+1/2}(z)\$\$

where \\Y\_\nu(z)\\ is the cylindrical Bessel function of the second
kind.

The spherical Bessel functions satisfy the same differential equation as
\\j\_\nu(z)\\: \$\$ z^2 \frac{d^2 y\_\nu}{dz^2} + 2z
\frac{dy\_\nu}{dz} + \[z^2 - \nu(\nu+1)\] y\_\nu = 0 \$\$

**Special cases:**

- \\y\_\nu(0) = -\infty\\ (singularity at the origin).

- \\y_0(z) = -\frac{\cos(z)}{z}\\

- \\y_1(z) = -\frac{\cos(z)}{z^2} - \frac{\sin(z)}{z}\\

**Derivatives:**

- First derivative: \\ y'\_\nu(z) = \frac{\nu}{z} y\_\nu(z) -
  y\_{\nu+1}(z) \\

- Second derivative: \\y''\_\nu(z) = -\frac{\nu}{z^2} y\_\nu(z) +
  \frac{\nu}{z} y'\_\nu(z) - y'\_{\nu+1}(z)\\

## References

Abramowitz, M. and Stegun, I.A. (Eds.). (1964). *Handbook of
Mathematical Functions with Formulas, Graphs, and Mathematical Tables*.
National Bureau of Standards, Applied Mathematics Series 55. Chapter 10.

NIST Digital Library of Mathematical Functions. <https://dlmf.nist.gov/>

- Spherical Bessel functions: <https://dlmf.nist.gov/10.47>

- Relation to cylindrical Bessel functions: Eq. 10.47.4 at
  <https://dlmf.nist.gov/10.47>

## See also

[`yc`](https://brandynlucca.github.io/acousticTS/reference/yc.md) for
cylindrical Bessel functions of the second kind,
[`js`](https://brandynlucca.github.io/acousticTS/reference/js.md) for
spherical Bessel functions of the first kind,
[`hs`](https://brandynlucca.github.io/acousticTS/reference/hs.md) for
spherical Hankel functions.

## Examples

``` r
# Spherical Bessel function of the second kind
ys(0, 1)
#> [1] -0.5403023
ys(1, 2.5)
#> [1] -0.1112059

# Fractional order
ys(0.5, 3)
#> [1] 0.3299975

# Vector input
ys(0, c(1, 2, 3))
#> [1] -0.5403023  0.2080734  0.3299975

# Singularity at origin
ys(0, 0) # Returns -Inf
#> [1] -Inf

# First derivative
ysdk(1, 2, 1)
#> [1] 0.5586854

# Second derivative
ysdk(1, 2, 2)
#> [1] -0.3833794

# 3rd derivative
ysdk(1, 1, 3)
#> [1] 23.85335

# 4th derivative
ysdk(1, 1, 4)
#> [1] -120.1082
```
