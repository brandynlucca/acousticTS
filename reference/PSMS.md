# Prolate spheroidal modal series (PSMS) solution

Prolate spheroidal modal series (PSMS) solution

## Value

No value; this help topic documents the PSMS model.

## Details

Calculates the far-field scattering amplitude and related quantities for
a prolate spheroid using the modal series solution, supporting various
boundary conditions (rigid, pressure-release, liquid-filled, and
gas-filled).

## Usage

This model is accessed via:


    target_strength(
      ...,
      model="psms",
      phi_body,
      boundary,
      adaptive,
      simplify_Amn,
      precision,
      n_integration,
      sound_speed_sw,
      density_sw
    )

## Arguments

- `phi_body`:

  Incident roll angle (radians).

- `boundary`:

  Boundary condition at the spheroid surface. One of `"fixed_rigid"`,
  `"pressure_release"`, `"liquid_filled"`, or `"gas_filled"`. See the
  boundary conditions documentation for more details on these different
  boundary conditions.

- `adaptive`:

  Apply adaptive modal support and quadrature selection. Defaults to
  `FALSE`.

- `simplify_Amn`:

  When `FALSE` (default), solve the coupled pressure and normal-velocity
  boundary equations for liquid- or gas-filled spheroids. When `TRUE`,
  retain the diagonal approximation described in Theory. This
  approximation can differ near gas resonances. Ignored for other
  boundary conditions.

- `precision`:

  A literal argument that allows for levels of precision when
  calculating the expansion matrix, \\A\_{mn}\\. There are two choices:
  `"double"` for double-precision (64-bit) and `"quad"` for
  quadruple-precision (128-bit). Use `"quad"` when strong gas contrast
  or high acoustic size makes double precision unstable.

- `n_integration`:

  Number of Gauss-Legendre integration points used to calculate angular
  overlap integrals for full liquid-filled solutions in either
  precision. With `adaptive = FALSE`, an explicit value is used as
  supplied and `NULL` selects automatic integration. With
  `adaptive = TRUE`, quadrature order is selected internally and an
  explicit value is ignored with a warning. Full gas solutions calculate
  overlaps from normalized expansion coefficients and do not use this
  quadrature count. Fixed quadrature must resolve the angular functions.
  Automatic quadrature integrates the retained angular expansion
  products using the Gauss-Legendre polynomial exactness condition.
  Under this condition, orthogonality allows the integral to be
  evaluated algebraically from the complete expansion coefficients.
  Fixed integration counts continue to use their requested quadrature
  rule. See
  [gauss_legendre](https://brandynlucca.github.io/acousticTS/reference/gauss_legendre.md).

- `sound_speed_sw`:

  Seawater sound speed (\\m~s^{-1}\\).

- `density_sw`:

  Seawater density (\\kg~m^{-3}\\).

## Theory

some relevant parameters are used to describe the geometry. A spheroid
surface is given by \\\xi = \xi_0 = \mathrm{constant}\\. Denoting the
major radius \\a\\, the minor radius \\b\\, and the semi-focal-length
\\q\\, then:

\$\$ a = \xi_0 q \\ \xi_0 = \left\[1 -
\left(\frac{b}{a}\right)^2\right\]^{-1/2} \$\$

Densities are expressed by \\\rho\\, sound speeds by \\c\\, and
wavenumbers by \\k\\, each followed by a subscript 0 for the surrounding
medium or 1 for the spheroidal body. For the spheroid, an alternative of
the reduced frequency \\ka\\ is \\h = kq\\.

The far-field scattering amplitude is given by:

\$\$ f\_\infty(\theta, \phi \| \theta', \phi') = \frac{2}{j k_0}
\sum\limits\_{m=0}^\infty \sum\limits\_{n=m}^\infty
\frac{\epsilon_m}{N\_{mn}(h_0)} S\_{mn}(h_0, \cos\theta') \\ A\_{mn}
S\_{mn}(h_0, \cos\theta) \cos m(\phi - \phi') \$\$

where \\\epsilon_m = (-1)^{m/2}\\ for even m, and \\N\_{mn}(h_0)\\ is
the norm. \\S\_{mn}\\ is the angle spheroidal wave function of the first
kind of order \\m\\ and degree \\n\\, and \\A\_{mn}\\ is the expansion
coefficient for the scattered wave, determined by the boundary
conditions. The parameters \\\theta\\ and \\\phi\\ denote the spherical
angle coordinates of an observed point along the body. The \\\theta'\\
and \\\phi'\\ denote similar spherical angle coordinates of the incident
direction. Effectively, \\\phi\\ and \\\theta\\ correspond to the roll
and tilt angles (relative the incident sound wave), respectively. It is
assumed that \\\theta'\\ and \\\phi'\\ are perpendicular to the incident
wave such that:

\$\$ \theta' = \pi - \theta \\ \phi' = \pi + \phi \$\$

For pressure release (or soft) and rigid spheroids:

\$\$ A\_{mn} = \frac{-\Delta R\_{mn}^{(1)}(h_0, \xi_0)}{ \Delta
R\_{mn}^{(3)}(h_0, \xi_0) } \$\$

where \\\Delta = 1\\ for soft and \\\Delta = \partial / \partial \xi\\ f
or rigid spheroid, and \\R\_{mn}^{(i)}\\ is the radial spheroidal wave
function of the i-th kind.

For the case of a fluid-filled spheroid, the following simultaneous
equation must be solved:

\$\$ \sum\limits\_{n'=m}^\infty K\_{nn'}^{(3)} A\_{mn'} +
\sum\limits\_{n'=m}^\infty K\_{nn'}^{(1)} \alpha\_{mn'}^{m}
E\_{n'}^{m(1)} = 0 \$\$

where

\$\$ K\_{nl}^{(i)} = \frac{1}{N\_{mn}(h_0)} \int\limits\_{-1}^{1}
S\_{mn}(h_0, \cos\theta) S\_{ml}(h_1, \eta) d\eta \$\$

and

\$\$ \alpha\_{mn}^{m} = \frac{1}{N\_{mn}(h_1)} \int\limits\_{-1}^{1}
S\_{mn}(h_0, \eta) S\_{mn}(h_1, \eta) d\eta \$\$

In practice, these kernel matrices are solved numerically to determine
\\A\_{mn}\\. The implementation uses compiled dense linear algebra and
keeps the fluid-filled solve in the requested arithmetic so that the
linear-system stage does not become detached from the rest of the chosen
precision pathway.

For full gas solutions, angular overlaps are dot products of expansion
coefficients in an orthonormal associated Legendre basis. This follows
from the spheroidal expansion and Legendre orthogonality identities in
[DLMF 30.8.1](https://dlmf.nist.gov/30.8.E1) and [DLMF
14.17.6](https://dlmf.nist.gov/14.17.E6). The coefficients and overlap
sums retain the requested arithmetic precision.

In the case that \\h_0 \approx h_1\\, \\\alpha\_{mn}^{m} \approx 0\\ for
\\n \neq l\\, and a much simplified expression is derived:

\$\$ A\_{mn} = -\frac{E\_{mn}^{m(1)}}{E\_{mn}^{m(3)}} \$\$

The maximum values of \\m\\ and \\n\\ can be estimated by:

\$\$ m\_{\max} = \[2 k_0 b\] \\ n\_{\max} = m\_{\max} + \[h_0 / 2\] \$\$

Full liquid and gas solutions also use the interior-wave and
incidence-dependent estimates accompanying Khodabandeloo et al. (2025).
These estimates initialize the retained modal support. Adaptive
quad-precision backscatter calculations can enlarge support when
degree-tail or boundary checks require it and stop negligible azimuthal
contributions early.

Full fluid solutions preserve pressure and normal-velocity continuity.
Radial mantissas and exponents remain separate during assembly to avoid
intermediate overflow.

## Implementation

**`C++`**

This model is primarily implemented in `C++` since it leverages the
`Fortran` algorithm developed by Arnie Lee van Buren and Jeffery
Boisvert. Another reason this is primarily written in `C++` is due to
computational and performance reasons. As \\k_0\\, \\b\\, and \\h\\
increase, so do \\m\_{\max}\\ and \\n\_{\max}\\. While this is not as
much of a concern for the pressure-release and fixed rigid cases, this
results in increasingly large and unwieldy kernel matrices required for
solving the fluid-filled expansion matrices (\\A\_{mn}\\). Consequently,
this can result in the model taking an impractical amount of time to
compute \\f\_\infty(\theta, \phi \| \theta', \phi')\\. `C++` helps
reduce this burden compared to a pure `R` implementation. The current
implementation uses compiled matrix algebra throughout, applies
backscatter parity to avoid recomputing one of the two angular
\\S\_{mn}\\ matrices, and solves the fluid-filled kernel system natively
in the requested arithmetic.

**Precision**

Another consideration is the floating point precision inherent to `R`.
`R` uses double-precision, which stores numeric values in a 64-bit
format. At low \\m\_{\max}\\ and \\n\_{\max}\\, there is lower numerical
instability in the prolate spheroidal wave functions used in the model.
However, this is not the case at greater \\m\_{\max}\\ and \\n\_{\max}\\
where the difference between double- (64-bit) and quadruple-precision
(128-bit) can contribute to differences in target strength of more than
1 dB. While `R`-packages like `Rmpfr` provide access to the (`GNU`)
`MPFR C` library, the (`GCC`) `libquadmath C++` library. Quad precision
nevertheless remains much more computationally intensive than double
precision because the spheroidal function evaluations, overlap
integrals, and kernel systems all grow rapidly with \\m\_{\max}\\ and
\\n\_{\max}\\.

## References

Furusawa, M. (1988). Prolate spheroidal models for predicting general
trends of fish target strength. Journal of the Acoustical Society of
Japan, 9: 13-24.

Khodabandeloo, B., Heggelund, Y., Ystad, B., Marx, S. A. B., and
Pedersen, G. (2025). High-precision model and open-source software for
acoustic backscattering by liquid- and gas-filled prolate spheroids
across a wide frequency range and incident angles: Implications for
fisheries acoustics. Journal of Sound and Vibration, 616: 119227.
[doi:10.1016/j.jsv.2025.119227](https://doi.org/10.1016/j.jsv.2025.119227)
.

Sanderson, C., and Curtin, R. (2019). Practical sparse matrices in C++
with hybrid storage and template-based expression optimisation.
Mathematical and Computational Applications, 24: 70.

Spencer, R.D., and Granger, S. (1951). The scattering of sound from a
prolate spheroid. The Journal of the Acoustical Society of America, 23:
701-706.

Van Buren, A. L. and Boisvert, J. E. "Prolate Spheroidal Wave
Functions." GitHub repository:
<https://github.com/MathieuandSpheroidalWaveFunctions/Prolate_swf>

## See also

[`target_strength`](https://brandynlucca.github.io/acousticTS/reference/target_strength.md),
[`FLS`](https://brandynlucca.github.io/acousticTS/reference/FLS-class.md),
[`GAS`](https://brandynlucca.github.io/acousticTS/reference/GAS-class.md),
[`ESS`](https://brandynlucca.github.io/acousticTS/reference/ESS-class.md),
[`ProlateSpheroid`](https://brandynlucca.github.io/acousticTS/reference/ProlateSpheroid-class.md),
[`prolate_spheroid`](https://brandynlucca.github.io/acousticTS/reference/prolate_spheroid.md),
[`Smn`](https://brandynlucca.github.io/acousticTS/reference/Smn.md),
[`Rmn`](https://brandynlucca.github.io/acousticTS/reference/Rmn.md)

## Examples

``` r
subset(available_models(), model == "psms")
#>      model slot  source persistent aliases
#> psms  psms PSMS builtin      FALSE        
```
