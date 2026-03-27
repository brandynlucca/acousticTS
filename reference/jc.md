# Cylindrical Bessel function of the first kind, \\J\_\nu(z)\\, and its respective derivatives

Computes the cylindrical Bessel function of the first kind
(\\J\_\nu(z)\\) and its k-th derivatives.

## Usage

``` r
jc(l, n)

jcdk(l, n, k)
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

  Non-negative integer. The order of the derivative for `jcdk`.

## Value

A complex vector containing:

- `jc`: \\J\_\nu(z)\\

- `jcdk(..., k = 1)`: \\J'\_\nu(z)\\ (first derivative)

- `jcdk(..., k = 2)`: \\J''\_\nu(z)\\ (second derivative)

- `jcdk`: \\J\_\nu^{(k)}(z)\\ (k-th derivative)

## Details

The cylindrical Bessel function of the first kind satisfies Bessel's
differential equation:

\$\$z^2 \frac{d^2 J\_\nu}{dz^2} + z \frac{dJ\_\nu}{dz} + (z^2 - \nu^2) J
\_\nu = 0\$\$

**Supported argument types:**

- **Purely real arguments** (\\z = x\\, where \\x \in \mathbb{R}\\):
  Fully supported for both positive and negative values.

- **Purely imaginary arguments** (\\z = iy\\, where \\y \in
  \mathbb{R}\\): Computed using the identity \\J\_\nu(iy) =
  e^{i\pi\nu/2} I\_\nu(y)\\ where \\I\_\nu\\ is the modified Bessel
  function of the first kind.

- **General complex arguments** (\\z = x + iy\\, where \\x \neq 0\\ and
  \\y \neq 0\\): **Not supported**.

**Special cases:**

- \\J\_\nu(0) = 1\\ if \\\nu = 0\\, otherwise \\J\_\nu(0) = 0\\.

- For negative real arguments with integer order \\n\\: \\J_n(-x) =
  (-1)^n J_n(x)\\.

- For negative real arguments with non-integer order \\\nu\\:
  \\J\_\nu(-x) = e^{i\pi\nu} J\_\nu(x)\\ (complex result).

**Derivatives:**

- First derivative: \\J'\_\nu(z) = J\_{\nu-1}(z) - \frac{\nu}{z}
  J\_\nu(z)\\

- Second derivative: \\J''\_\nu(z) = \frac{1}{4}\left\[J\_{\nu-2}(z) -
  2J\_\nu(z) + J\_{\nu+2}(z)\right\]\\

- k-th derivative (DLMF 10.6.1): \$\$ \frac{d^k}{dz^k} J\_\nu(z) =
  \frac{1}{2^k} \sum\_{j=0}^{k} (-1)^j \binom{k}{j} J\_{\nu - k + 2j}(z)
  \$\$

## References

Abramowitz, M. and Stegun, I.A. (Eds.). (1964). *Handbook of
Mathematical Functions with Formulas, Graphs, and Mathematical Tables*.
National Bureau of Standards, Applied Mathematics Series 55. Chapter 9.

NIST Digital Library of Mathematical Functions. <https://dlmf.nist.gov/>

- Bessel's equation: <https://dlmf.nist.gov/10.2>

- Negative argument identity (\\J\_\nu(-z)\\): Eq. 10.4.1 at
  <https://dlmf.nist.gov/10.4>

- Imaginary argument identity (\\J\_\nu(iz)\\): Eq. 10.27.6 at
  <https://dlmf.nist.gov/10.27>

## See also

[`yc`](https://brandynlucca.github.io/acousticTS/reference/yc.md) for
Bessel functions of the second kind,
[`hc`](https://brandynlucca.github.io/acousticTS/reference/hc.md) for
Hankel functions (third kind),
[`js`](https://brandynlucca.github.io/acousticTS/reference/js.md) for
spherical Bessel functions of the first kind.

## Examples

``` r
# Real argument, integer order
jc(0, 1)
#> [1] 0.7651977+0i
jc(1, 2.5)
#> [1] 0.4970941+0i

# Real argument, fractional order
jc(0.5, 3)
#> [1] 0.06500818+0i

# Negative real argument (integer order gives real result)
jc(2, -1.5)
#> [1] 0.2320877+0i

# Negative real argument (fractional order gives complex result)
jc(0.5, -2)
#> [1] 3.141318e-17+0.5130161i

# Purely imaginary argument
jc(1, 1i)
#> [1] 0+0.5651591i
jc(2, 3i)
#> [1] -2.245212+0i

# First derivative
jcdk(1, 2, 1)
#> [1] -0.06447162+0i

# Second derivative
jcdk(1, 2, 2)
#> [1] -0.4003078+0i
```
