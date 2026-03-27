# Cylindrical Bessel function of the third kind (Hankel), \\H\_\nu(x)\\, and its respective derivatives

Computes the cylindrical Hankel function of the first kind
(\\H^{(1)}\_\nu(z)\\) and its derivatives through `hcdk()`.

## Usage

``` r
hc(l, n)

hcdk(l, n, k)
```

## Arguments

- l:

  Numeric. The order (\\\nu\\) of the Hankel function. Must be purely
  real; complex orders are not supported.

- n:

  Numeric or complex. The argument (\\z\\) at which to evaluate the
  function. Supports purely real or purely imaginary values. General
  complex arguments are not supported.

- k:

  Non-negative integer. The order of the derivative for `hcdk`.

## Value

A complex vector containing:

- `hc`: \\H^{(1)}\_\nu(z)\\

- `hcdk(..., k = 1)`: \\\frac{d}{dz}H^{(1)}\_\nu(z)\\ (first derivative)

- `hcdk(..., k = 2)`: \\\frac{d^2}{dz^2}H^{(1)}\_\nu(z)\\ (second
  derivative)

- `hcdk`: \\\frac{d^k}{dz^k}H^{(1)}\_\nu(z)\\ (k-th derivative)

## Details

The Hankel function of the first kind is defined as: \$\$H^{(1)}\_\nu(z)
= J\_\nu(z) + i Y\_\nu(z)\$\$

where \\J\_\nu(z)\\ is the Bessel function of the first kind and
\\Y\_\nu(z)\\ is the Bessel function of the second kind.

**Supported argument types:** Since \\H^{(1)}\_\nu(z)\\ is computed from
\\J\_\nu(z)\\ and \\Y\_\nu(z)\\, the same restrictions apply:

- **Purely real arguments** (\\z = x\\): Fully supported.

- **Purely imaginary arguments** (\\z = iy\\): Supported.

- **General complex arguments** (\\z = x + iy\\): **Not supported**.

**Derivatives:**

- First derivative: \\\frac{d}{dz}H^{(1)}\_\nu(z) = \frac{\nu}{z}
  H^{(1)}\_\nu(z) - H^{(1)}\_{\nu+1}(z)\\

- Second derivative: \\\frac{d^2}{dz^2}H^{(1)}\_\nu(z) =
  H^{(1)}\_{\nu-2}(z) - \frac{2\nu-1}{z} H^{(1)}\_{\nu-1}(z) +
  \frac{\nu^2+\nu}{z^2} H^{(1)}\_\nu(z)\\

- k-th derivative (DLMF 10.6.1): \\\frac{d^k}{dz^k}H^{(1)}\_\nu(z) =
  \frac{1}{2^k} \sum\_{j=0}^{k} (-1)^j \binom{k}{j}
  H^{(1)}\_{\nu-k+2j}(z)\\

## References

Abramowitz, M. and Stegun, I.A. (Eds.). (1964). *Handbook of
Mathematical Functions with Formulas, Graphs, and Mathematical Tables*.
National Bureau of Standards, Applied Mathematics Series 55. Chapter 9.

NIST Digital Library of Mathematical Functions. <https://dlmf.nist.gov/>

- Hankel function definition: <https://dlmf.nist.gov/10.2>

- k-th derivative formula: Eq. 10.6.1 at <https://dlmf.nist.gov/10.6>

## See also

[`jc`](https://brandynlucca.github.io/acousticTS/reference/jc.md) for
Bessel functions of the first kind,
[`yc`](https://brandynlucca.github.io/acousticTS/reference/yc.md) for
Bessel functions of the second kind,
[`hs`](https://brandynlucca.github.io/acousticTS/reference/hs.md) for
spherical Hankel functions.

## Examples

``` r
# Hankel function
hc(0, 1)
#> [1] 0.7651977+0.08825696i
hc(1, 2.5)
#> [1] 0.4970941+0.1459181i

# Fractional order
hc(0.5, 3)
#> [1] 0.06500818+0.4560488i

# Purely imaginary argument
hc(1, 1i)
#> [1] 0.383186+1.130318i

# First derivative
hcdk(1, 2, 1)
#> [1] -0.06447162+0.5638919i

# Second derivative
hcdk(1, 2, 2)
#> [1] -0.4003078-0.2016716i

# k-th derivative
hcdk(1, 2, 3) # Third derivative
#> [1] 0.08820851-0.154352i
```
