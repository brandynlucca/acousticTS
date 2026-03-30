# Prolate Spheroidal Angular Function of the First Kind, \\S^{1}\_{mn}(c, \eta)\\

Computes the prolate spheroidal angular function of the first kind,
\\S\_{mn}^{(1)}(c, \eta)\\, and its first derivative with respect to
\\\eta\\ for given order \\m\\, degree \\n\\, size parameter \\c\\, and
angular coordinate \\\eta\\.

This function is an R wrapper for compiled C++ code, which in turn calls
the underlying Fortran library (`prolate_swf`) for high-performance
numerical computation. All heavy computation is performed in compiled
code for speed and accuracy.

## Usage

``` r
Smn(m, n, c, eta, normalize = FALSE, precision = "double")
```

## Arguments

- m:

  Non-negative integer. The order of the spheroidal function (\\m \geq
  0\\).

- n:

  Non-negative integer. The degree of the spheroidal function (\\n \geq
  m\\).

- c:

  Numeric. The scalar size parameter.

- eta:

  Numeric vector. The angular coordinate(s) at which to evaluate the
  function. Must satisfy \\\|\eta\| \leq 1\\.

- normalize:

  Logical. If `TRUE`, the angular functions are normalized to have unity
  norm. If `FALSE` (default), the Meixner-Schäfke normalization is used.

- precision:

  Character. Either `"double"` (default) or `"quad"`. Controls the
  floating-point precision used in the underlying Fortran computation.

  `"double"`

  :   Uses standard double-precision (64-bit) arithmetic. Fastest and
      sufficient for most applications.

  `"quad"`

  :   Uses quadruple-precision (128-bit) arithmetic for higher numerical
      accuracy in challenging parameter regimes (e.g., large \\m\\,
      \\n\\, or near singularities). Computation is significantly
      slower.

## Value

A list containing:

- `value`:

  Numeric vector of function values \\S\_{mn}^{(1)}(c, \eta)\\ at each
  input `eta`.

- `derivative`:

  Numeric vector of first derivatives \\\frac{d}{d\eta}S\_{mn}^{(1)}(c,
  \eta)\\ at each input `eta`.

## Details

The prolate spheroidal angular functions are solutions to the angular
part of the scalar Helmholtz equation in prolate spheroidal coordinates.
They satisfy the differential equation: \$\$(1 - \eta^2) \frac{d^2
S\_{mn}}{d\eta^2} - 2\eta \frac{dS\_{mn}}{d\eta} + \left(\lambda\_{mn} -
c^2 \eta^2 + \frac{m^2}{1 - \eta^2}\right) S\_{mn} = 0\$\$

where \\\lambda\_{mn}\\ is the separation constant (eigenvalue).

**Domain restrictions:**

- The angular coordinate must satisfy \\\|\eta\| \leq 1\\.

- The order must satisfy \\m \geq 0\\.

- The degree must satisfy \\n \geq m\\.

**Normalization:** When `normalize = FALSE` (default), the functions use
the Meixner-Schäfke normalization, which matches the normalization of
the corresponding associated Legendre functions. When
`normalize = TRUE`, the functions are scaled to have unity norm.

**Implementation:** This function is an R wrapper for a compiled C++
interface (`Smn_cpp`), which itself wraps the Fortran subroutine
`profcn` from the `prolate_swf` library developed by Arnie Lee Van Buren
and Jeffrey Boisvert. The underlying algorithm uses a combination of
forward and backward recursion with the Bouwkamp eigenvalue method for
high accuracy across wide parameter ranges. The C++ layer manages
memory, precision selection, and data conversion between R and Fortran
for robust and efficient computation.

## References

Van Buren, A. L. and Boisvert, J. E. "Prolate Spheroidal Wave
Functions." GitHub repository:
<https://github.com/MathieuandSpheroidalWaveFunctions/prolate_swf>

Meixner, J. and Schäfke, F. W. (1954). *Mathieusche Funktionen und
Sphäroidfunktionen*. Springer-Verlag, Berlin.

Flammer, C. (1957). *Spheroidal Wave Functions*. Stanford University
Press.

NIST Digital Library of Mathematical Functions. Chapter 30: Spheroidal
Wave Functions. <https://dlmf.nist.gov/30>

## See also

[`Rmn`](https://brandynlucca.github.io/acousticTS/reference/Rmn.md) for
prolate spheroidal radial functions.

## Examples

``` r
# Single evaluation
Smn(m = 2, n = 3, c = 1, eta = 0.5)
#> $value
#> [1] 5.650368
#> 
#> $derivative
#> [1] 3.454326
#> 

# Multiple eta values
Smn(m = 0, n = 2, c = 5, eta = c(-0.5, 0, 0.5))
#> [[1]]
#> [[1]]$value
#>           [,1]
#> [1,] 0.3284913
#> 
#> [[1]]$derivative
#>           [,1]
#> [1,] -1.926133
#> 
#> 
#> [[2]]
#> [[2]]$value
#>            [,1]
#> [1,] -0.4272929
#> 
#> [[2]]$derivative
#>      [,1]
#> [1,]    0
#> 
#> 
#> [[3]]
#> [[3]]$value
#>           [,1]
#> [1,] 0.3284913
#> 
#> [[3]]$derivative
#>          [,1]
#> [1,] 1.926133
#> 
#> 

# With unity normalization
Smn(m = 1, n = 1, c = 2, eta = 0.3, normalize = TRUE)
#> $value
#> [1] 0.856142
#> 
#> $derivative
#> [1] -0.4724318
#> 

# Double precision (default)
Smn(m = 2, n = 3, c = 1, eta = 0.5, precision = "double")
#> $value
#> [1] 5.650368
#> 
#> $derivative
#> [1] 3.454326
#> 

# Quad precision
Smn(m = 2, n = 3, c = 1, eta = 0.5, precision = "quad")
#> $value
#> [1] 5.650368
#> 
#> $derivative
#> [1] 3.454326
#> 
```
