# Cylindrical Bessel function of the second kind, \\Y\_\nu(x)\\, and its respective derivatives

Computes the cylindrical Bessel function of the second kind
(\\Y\_\nu(z)\\), also known as the Neumann function or Weber function,
and its derivatives through `ycdk()`.

## Usage

``` r
yc(l, n)

ycdk(l, n, k)
```

## Arguments

- l:

  Numeric. The order (\\\nu\\) of the Bessel function. Must be purely
  real; complex orders are not supported.

- n:

  Numeric or complex. The argument (\\z\\) at which to evaluate the
  function. Supports purely real or purely imaginary values. General
  complex arguments (\\x + iy\\ with \\x \neq 0\\ and \\y \neq 0\\) are
  not supported.

- k:

  Non-negative integer. The order of the derivative for `ycdk`.

## Value

A complex vector containing:

- `yc`: \\Y\_\nu(z)\\

- `ycdk(..., k = 1)`: \\Y'\_\nu(z)\\ (first derivative)

- `ycdk(..., k = 2)`: \\Y''\_\nu(z)\\ (second derivative)

- `ycdk`: \\Y\_\nu^{(k)}(z)\\ (k-th derivative)

## Details

The cylindrical Bessel function of the second kind satisfies the same
differential equation as \\J\_\nu(z)\\: \$\$z^2 \frac{d^2
Y\_\nu}{dz^2} + z \frac{dY\_\nu}{dz} + (z^2 - \nu^2) Y\_\nu = 0\$\$

but represents the linearly independent second solution.

**Supported argument types:**

- **Purely real arguments** (\\z = x\\, where \\x \in \mathbb{R}\\):
  Fully supported for both positive and negative values.

- **Purely imaginary arguments** (\\z = iy\\, where \\y \in
  \mathbb{R}\\): Computed using the identity \\Y\_\nu(iy) = i
  e^{-i\pi\nu/2} I\_\nu(y) - \frac{2}{\pi} e^{i\pi\nu/2} K\_\nu(y)\\
  where \\I\_\nu\\ and \\K\_\nu\\ are modified Bessel functions.

- **General complex arguments** (\\z = x + iy\\, where \\x \neq 0\\ and
  \\y \neq 0\\): **Not supported**.

**Special cases:**

- \\Y\_\nu(0) = -\infty\\ (singularity at the origin).

- For negative real arguments: \\Y\_\nu(-x) = \cos(\pi\nu) Y\_\nu(x) +
  \sin(\pi\nu) J\_\nu(x)\\.

- For integer order \\n\\: \\Y_n(-x) = (-1)^n Y_n(x)\\.

- k-th derivative (DLMF 10.6.1): \$\$ \frac{d^k}{dz^k} Y\_\nu(z) =
  \frac{1}{2^k} \sum\_{j=0}^{k} (-1)^j \binom{k}{j} Y\_{\nu - k + 2j}(z)
  \$\$

**Derivative:** \$\$Y'\_\nu(z) = Y\_{\nu-1}(z) - \frac{\nu}{z}
Y\_\nu(z)\$\$

## References

Abramowitz, M. and Stegun, I.A. (Eds.). (1964). *Handbook of
Mathematical Functions with Formulas, Graphs, and Mathematical Tables*.
National Bureau of Standards, Applied Mathematics Series 55. Chapter 9.

NIST Digital Library of Mathematical Functions. <https://dlmf.nist.gov/>

- Bessel's equation: <https://dlmf.nist.gov/10.2>

- Negative argument identity (\\Y\_\nu(-z)\\): Eq. 10.4.1 at
  <https://dlmf.nist.gov/10.4>

- Imaginary argument identity (\\Y\_\nu(iz)\\): Eq. 10.27.8 at
  <https://dlmf.nist.gov/10.27>

## See also

[`jc`](https://brandynlucca.github.io/acousticTS/reference/jc.md) for
Bessel functions of the first kind,
[`hc`](https://brandynlucca.github.io/acousticTS/reference/hc.md) for
Hankel functions (third kind),
[`ys`](https://brandynlucca.github.io/acousticTS/reference/ys.md) for
spherical Bessel functions of the second kind.

## Examples

``` r
# Real argument, integer order
yc(0, 1)
#> [1] 0.08825696+0i
yc(1, 2.5)
#> [1] 0.1459181+0i

# Real argument, fractional order
yc(0.5, 3)
#> [1] 0.4560488+0i

# Negative real argument
yc(2, -1.5)
#> [1] -0.9321938+0i

# Purely imaginary argument
yc(1, 1i)
#> [1] 0.5651591-0.383186i
yc(2, 3i)
#> [1] 0.03915877-2.245212i

# Singularity at origin
yc(0, 0) # Returns -Inf
#> [1] -Inf+0i

# First derivative
ycdk(1, 2, 1)
#> [1] 0.5638919+0i
```
