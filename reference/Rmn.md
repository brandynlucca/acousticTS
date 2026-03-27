# Prolate Spheroidal Radial Functions

Computes the prolate spheroidal radial functions of the first
(\\R\_{mn}^{(1)}\\), second (\\R\_{mn}^{(2)}\\), third
(\\R\_{mn}^{(3)}\\), and fourth (\\R\_{mn}^{(4)}\\) kinds and their
first derivatives with respect to \\\xi\\ for given order \\m\\, degree
\\n\\, size parameter \\c\\, and radial coordinate \\\xi\\.

This function is an R wrapper for compiled C++ code, which in turn calls
the underlying Fortran library (`prolate_swf`) for high-performance
numerical computation. All heavy computation is performed in compiled
code for speed and accuracy.

## Usage

``` r
Rmn(m, n, c, xi, kind = 1, precision = "double")
```

## Arguments

- m:

  Non-negative integer. The order of the spheroidal function (\\m \geq
  0\\).

- n:

  Non-negative integer. The degree of the spheroidal function (\\n \geq
  m\\).

- c:

  Numeric. The scalar size parameter (also denoted \\\gamma\\ in some
  references).

- xi:

  Numeric. The radial coordinate at which to evaluate the function. Must
  satisfy \\\xi \geq 1\\.

- kind:

  Integer. Specifies which kind of radial function to compute: `1`
  (first kind), `2` (second kind), `3` (third kind/outgoing), or `4`
  (fourth kind/incoming). Default is `1`.

- precision:

  Character. Either `"double"` (default) or `"quad"`. Controls the
  floating-point precision used in the underlying Fortran computation.

  `"double"`

  :   Uses standard double-precision (64-bit) arithmetic. Fastest and
      sufficient for most applications.

  `"quad"`

  :   Uses quadruple-precision (128-bit) arithmetic for higher numerical
      accuracy in challenging parameter regimes (e.g., large \\m\\,
      \\n\\, \\c\\, or near singularities). Computation is significantly
      slower.

## Value

A list containing:

- `value`:

  The function value \\R\_{mn}^{(k)}(c, \xi)\\. Real for `kind = 1` or
  `2`; complex for `kind = 3` or `4`.

- `derivative`:

  The first derivative \\\frac{d}{d\xi}R\_{mn}^{(k)}(c, \xi)\\. Real for
  `kind = 1` or `2`; complex for `kind = 3` or `4`.

## Details

The prolate spheroidal radial functions are solutions to the radial part
of the scalar Helmholtz equation in prolate spheroidal coordinates. They
satisfy the differential equation: \$\$(\xi^2 - 1) \frac{d^2
R\_{mn}}{d\xi^2} + 2\xi \frac{dR\_{mn}}{d\xi} - \left(\lambda\_{mn} -
c^2 \xi^2 + \frac{m^2}{\xi^2 - 1}\right) R\_{mn} = 0\$\$

where \\\lambda\_{mn}\\ is the separation constant (eigenvalue).

**Function kinds:**

- `kind = 1`: Radial function of the first kind \\R\_{mn}^{(1)}(c,
  \xi)\\. Regular at the origin; analogous to spherical Bessel functions
  \\j_n\\.

- `kind = 2`: Radial function of the second kind \\R\_{mn}^{(2)}(c,
  \xi)\\. Singular at the focal points (\\\xi = 1\\); analogous to
  spherical Neumann functions \\y_n\\.

- `kind = 3`: Radial function of the third kind (outgoing Hankel-type)
  \\R\_{mn}^{(3)}(c, \xi) = R\_{mn}^{(1)}(c, \xi) + i R\_{mn}^{(2)}(c,
  \xi)\\. Used for outgoing wave solutions.

- `kind = 4`: Radial function of the fourth kind (incoming Hankel-type)
  \\R\_{mn}^{(4)}(c, \xi) = R\_{mn}^{(1)}(c, \xi) - i R\_{mn}^{(2)}(c,
  \xi)\\. Used for incoming wave solutions.

**Domain restrictions:**

- The radial coordinate must satisfy \\\xi \geq 1\\ (prolate domain).

- The order must satisfy \\m \geq 0\\.

- The degree must satisfy \\n \geq m\\.

**Normalization:** The functions use the Morse-Feshbach normalization by
default, where the radial functions reduce to spherical Bessel/Neumann
functions as \\c \rightarrow 0\\.

**Implementation:** This function is an R wrapper for a compiled C++
interface (`Rmn_cpp`), which itself wraps the Fortran subroutine
`profcn` from the `prolate_swf` library developed by Arnie Lee Van Buren
and Jeffrey Boisvert. The underlying algorithm uses a combination of
forward and backward recursion with the Bouwkamp eigenvalue method for
high accuracy across wide parameter ranges. The C++ layer manages
memory, precision selection, and data conversion between R and Fortran
for robust and efficient computation.

## References

Van Buren, A. L. and Boisvert, J. E. "Prolate Spheroidal Wave
Functions." GitHub repository:
<https://github.com/MathieuandSpheroidalWaveFunctions/Prolate_swf>

Morse, P. M. and Feshbach, H. (1953). *Methods of Theoretical Physics*.
McGraw-Hill, New York. Chapter 21.

Flammer, C. (1957). *Spheroidal Wave Functions*. Stanford University
Press.

NIST Digital Library of Mathematical Functions. Chapter 30: Spheroidal
Wave Functions. <https://dlmf.nist.gov/30>

## See also

[`Smn`](https://brandynlucca.github.io/acousticTS/reference/Smn.md) for
prolate spheroidal angular functions.

## Examples

``` r
# First kind radial function
Rmn(m = 2, n = 3, c = 1, xi = 1.5, kind = 1)
#> $value
#> [1] 0.01634468
#> 
#> $derivative
#> [1] 0.04735489
#> 

# Second kind radial function
Rmn(m = 0, n = 2, c = 5, xi = 2.0, kind = 2)
#> $value
#> [1] -0.1082697
#> 
#> $derivative
#> [1] 0.2477105
#> 

# Third kind (outgoing) radial function
Rmn(m = 1, n = 1, c = 2, xi = 1.2, kind = 3)
#> $value
#> [1] 0.342192-0.6005553i
#> 
#> $derivative
#> [1] 0.5788199+2.304993i
#> 

# Fourth kind (incoming) radial function
Rmn(m = 1, n = 1, c = 2, xi = 1.2, kind = 4)
#> $value
#> [1] 0.342192+0.6005553i
#> 
#> $derivative
#> [1] 0.5788199-2.304993i
#> 

# Double precision (default)
Rmn(m = 2, n = 3, c = 1, xi = 1.5)
#> $value
#> [1] 0.01634468
#> 
#> $derivative
#> [1] 0.04735489
#> 

# Quad precision
Rmn(m = 2, n = 3, c = 1, xi = 1.5, precision = "quad")
#> $value
#> [1] 0.01634468
#> 
#> $derivative
#> [1] 0.04735489
#> 
```
