# Legendre Function of the Second Kind, \\Q\_\nu(x)\\

Computes the Legendre function of the second kind, \\Q\_\nu(x)\\, for
real order \\\nu\\ (integer or fractional) and real argument \\x\\.
Returns complex values when \\\|x\| \> 1\\.

## Usage

``` r
Qn(n, x)
```

## Arguments

- n:

  Numeric vector. Degree (order) of the Legendre function. Can be
  integer or fractional (e.g., 0, 1, 2.5, 3.7).

- x:

  Numeric vector. Real argument(s) at which to evaluate the function.
  Valid for all real values. Note that \\x = \pm 1\\ are singularities.

## Value

A complex matrix of dimension `length(n)` by `length(x)`, where element
`[i, j]` contains \\Q\_{n_i}(x_j)\\. For \\\|x\| \< 1\\, the imaginary
part is zero.

## Details

The Legendre function of the second kind satisfies the same differential
equation as \\P\_\nu(x)\\: \$\$(1 - x^2) \frac{d^2 Q\_\nu}{dx^2} - 2x
\frac{dQ\_\nu}{dx} + \nu(\nu + 1) Q\_\nu = 0\$\$

but represents the linearly independent second solution.

**For \\\|x\| \< 1\\ (real result):**

- For integer order \\n\\: Uses the `C++` `Boost` library implementation

- For fractional order \\\nu\\: Uses the Ferrers identity \$\$Q\_\nu(x)
  = \frac{\pi}{2 \sin(\pi \nu)} \left\[ \cos(\pi \nu) P\_\nu(x) -
  P\_\nu(-x) \right\]\$\$

**For \\x = \pm 1\\:** Returns infinity (singularity).

**For \\\|x\| \> 1\\ (complex result):**

- Real part: Computed via the integral representation
  \$\$\text{Re}\[Q\_\nu(x)\] = \int_0^\infty \frac{dt}{(x + \sqrt{x^2-1}
  \cosh t)^{\nu+1}}\$\$

- Imaginary part: \$\$\text{Im}\[Q\_\nu(x)\] = -\frac{\pi}{2}
  P\_\nu(x)\$\$

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

[`Qndk`](https://brandynlucca.github.io/acousticTS/reference/Qndk.md)
for the k^(th) derivative of the Legendre polynomial of the second kind,
[`Pn`](https://brandynlucca.github.io/acousticTS/reference/Pn.md) for
Legendre functions of the first kind.

## Examples

``` r
# Single values
Qn(1, 0.5)
#>               [,1]
#> [1,] -0.7253469+0i

# Multiple orders, single argument
Qn(c(1, 2, 3), 0.5)
#>               [,1]
#> [1,] -0.7253469+0i
#> [2,] -0.8186633+0i
#> [3,] -0.1986548+0i

# Single order, multiple arguments
Qn(1, c(-0.5, 0, 0.5))
#>               [,1]  [,2]          [,3]
#> [1,] -0.7253469+0i -1+0i -0.7253469+0i

# Multiple orders and arguments (returns a matrix)
Qn(c(1, 2, 3), c(0.25, 0.5, 0.75))
#>               [,1]          [,2]          [,3]
#> [1,] -0.9361468+0i -0.7253469+0i -0.2702837+0i
#> [2,] -0.4787615+0i -0.8186633+0i -0.7905467+0i
#> [3,]  0.4246139+0i -0.1986548+0i -0.8079942+0i

# Fractional orders
Qn(c(0.5, 1.5, 2.7), c(0.5, 0.9))
#>               [,1]          [,2]
#> [1,] -0.2655964+0i  0.7873898+0i
#> [2,] -0.8959028+0i -0.0245780+0i
#> [3,] -0.4202294+0i -0.5761958+0i

# Arguments |x| > 1 return complex values
Qn(c(1, 2, 3.5), c(2, 3.5))
#>                        [,1]                     [,2]
#> [1,] 0.098612289- 3.141593i 0.0286266636-  5.497787i
#> [2,] 0.021183794- 8.639380i 0.0033433176- 28.077984i
#> [3,] 0.002370082-47.971060i 0.0001501177-390.254185i

# Singularity at x = 1
Qn(1, 1)
#>          [,1]
#> [1,] Inf+Infi
```
