# Spherical Bessel function of the third kind (Hankel), \\h\_\nu(x)\\, and its respective derivatives

Computes the spherical Hankel function of the first kind
(\\h^{(1)}\_\nu(z)\\) and its first-th derivative.

## Usage

``` r
hs(l, n)

hsdk(l, n, k)
```

## Arguments

- l:

  Numeric. The order of the spherical Hankel function. Can be integer or
  fractional.

- n:

  Numeric. The argument (\\z\\) at which to evaluate the function.

- k:

  Non-negative integer. The order of the derivative for `hsdk`.

## Value

A complex vector containing:

- `hs`: \\h^{(1)}\_\nu(z)\\

- `hsdK`: \\\frac{d}{dz^k}h^{(1)}\_l(z)\\ (k-th derivative)

## Details

The spherical Hankel function of the first kind is defined as:
\$\$h^{(1)}\_\nu(z) = j\_\nu(z) + i y\_\nu(z)\$\$

where \\j\_\nu(z)\\ is the spherical Bessel function of the first kind
and \\y\_\nu(z)\\ is the spherical Bessel function of the second kind.

It is related to the cylindrical Hankel function by: \$\$h^{(1)}\_\nu(z)
= \sqrt{\frac{\pi}{2z}} H^{(1)}\_{\nu+1/2}(z)\$\$

The spherical Hankel functions are used extensively in scattering theory
to represent outgoing spherical waves.

**Derivative:** \$\$ \frac{d}{dz}h^{(1)}\_\nu(z) = \frac{\nu}{z}
h^{(1)}\_\nu(z) - h^{(1)}\_{\nu+1}(z) \$\$

## References

Abramowitz, M. and Stegun, I.A. (Eds.). (1964). *Handbook of
Mathematical Functions with Formulas, Graphs, and Mathematical Tables*.
National Bureau of Standards, Applied Mathematics Series 55. Chapter 10.

NIST Digital Library of Mathematical Functions. <https://dlmf.nist.gov/>

- Spherical Bessel functions: <https://dlmf.nist.gov/10.47>

- Relation to cylindrical Hankel functions: Eq. 10.47.5 at
  <https://dlmf.nist.gov/10.47>

## See also

[`hc`](https://brandynlucca.github.io/acousticTS/reference/hc.md) for
cylindrical Hankel functions,
[`js`](https://brandynlucca.github.io/acousticTS/reference/js.md) for
spherical Bessel functions of the first kind,
[`ys`](https://brandynlucca.github.io/acousticTS/reference/ys.md) for
spherical Bessel functions of the second kind.

## Examples

``` r
# Spherical Hankel function
hs(0, 1)
#> [1] 0.841471-0.5403023i
hs(1, 2.5)
#> [1] 0.416213-0.1112059i

# Fractional order
hs(0.5, 3)
#> [1] 0.04704+0.3299975i

# Vector input
hs(0, c(1, 2, 3))
#> [1] 0.8414710-0.5403023i 0.4546487+0.2080734i 0.0470400+0.3299975i

# First derivative
hsdk(1, 2, 1)
#> [1] 0.01925094+0.5586854i
```
