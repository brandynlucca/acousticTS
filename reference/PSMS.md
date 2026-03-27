# Prolate spheroidal modal series (PSMS) solution

Prolate spheroidal modal series (PSMS) solution

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

  Boundary condition at the cylinder surface. One of `"fixed_rigid"`,
  `"pressure_release"`, `"liquid_filled"`, or `"gas_filled"`. See the
  [boundary conditions
  documentation](https://brandynlucca.github.io/acousticTS/reference/boundary_conditions.md)
  for more details on these different boundary conditions.

- `adaptive`:

  A boolean argument controlling whether the PSMS backend is allowed to
  stop once the retained modal tail becomes numerically negligible and
  to choose an adaptive quadrature order for full fluid- or gas-filled
  overlap integrals. When `FALSE` (default), the implementation uses the
  published hard truncation rule and a fixed quadrature order unless the
  user overrides it directly. When `TRUE`, the hard \\m\_{\max}\\ and
  \\n\_{\max}\\ rules remain the upper bounds, but the backscatter
  implementation is allowed to stop early once the retained tail is both
  numerically small and gradient-flat.

- `simplify_Amn`:

  A boolean argument that flags whether or not to use the simplified
  calculation for the exansion matrix, \\A\_{mn}\\, for a fluid-filled
  scatterer. See Theory for the mathematical formulation and
  assumptions. **Note:** this argument **only** applies to when
  `boundary = "liquid_filled"` or `boundary = "gas_filled"`. It is
  otherwise left unused for all other boundary conditions.

- `precision`:

  A literal argument that allows for levels of precision when
  calculating the expansion matrix, \\A\_{mn}\\. There are two choices:
  `"double"` for double-precision (64-bit) and `"quad"` for
  quadruple-precision (128-bit). See Details for more double- and
  quadruple-precision usages are implemented.

- `n_integration`:

  An integer argument that informs the model how many integration points
  will be used to calculate \\\alpha\_{mn}^{m}\\, which is numerically
  calculated using Gauss-Legendre quadrature. When left as `NULL`, the
  model uses 96 integration points unless `adaptive = TRUE`, in which
  case a reduced-frequency-based quadrature rule is selected internally
  for the full fluid- or gas-filled solve. See
  [gauss_legendre](https://brandynlucca.github.io/acousticTS/reference/gauss_legendre.md)
  for a full description of how `n_integration` is used. **Note:** this
  argument **only** applies to when `boundary = "liquid_filled"` or
  `boundary = "gas_filled"`. It is otherwise left unused for all other
  boundary conditions.

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

In the case that \\h_0 \approx h_1\\, \\\alpha\_{mn}^{m} \approx 0\\ for
\\n \neq l\\, and a much simplified expression is derived:

\$\$ A\_{mn} = -\frac{E\_{mn}^{m(1)}}{E\_{mn}^{m(3)}} \$\$

The maximum values of \\m\\ and \\n\\ can be estimated by:

\$\$ m\_{\max} = \[2 k_0 b\] \\ n\_{\max} = m\_{\max} + \[h_0 / 2\] \$\$

## Implementation

***`C++`***

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

***Precision***

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
