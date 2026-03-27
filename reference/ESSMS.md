# Elastic-shelled sphere modal series (ESSMS) solution

Calculates the far-field scattering amplitude and related quantities
elastic shelled sphere using the modal series solution, as described by
Goodman and Stern (1962).

## Usage

This model is accessed via:

    target_strength(
      ...,
      model="essms",
      sound_speed_sw,
      density_sw,
      m_limit
    )

## Arguments

- `sound_speed_sw`:

  Seawater sound speed (\\m~s^{-1}\\).

- `density_sw`:

  Seawater density (\\kg~m^{-3}\\).

- `m_limit`:

  Optional model truncation limit used to cap the number of modes in the
  numerical calculation. By default, this is set to the maximum value of
  \\ka + 10\\.

## Theory

The elastic-shelled sphere model solves the acoustic scattering problem
by expanding the incident and scattered fields in terms of spherical
Bessel and Hankel functions and Legendre polynomials. The displacement
field in the shell is represented as the sum of a scalar potential and a
vector potential:

\$\$ \mathbf{u} = \nabla \phi + \nabla \times \mathbf{\Psi} \$\$

where \\\phi\\ and \\\mathbf{\Psi}\\ satisfy the Helmholtz equations
with velocities determined by the elastic properties of the shell:

\$\$ (1/c_L^2) \frac{\partial^2 \phi}{\partial t^2} = \nabla^2 \phi,
\qquad (1/c_T^2) \frac{\partial^2 \mathbf{\Psi}}{\partial t^2} =
\nabla^2 \mathbf{\Psi} \$\$

with \\c_L = \sqrt{(\lambda + 2\mu)/\rho}\\ (longitudinal) and \\c_T =
\sqrt{\mu/\rho}\\ (transverse) wave speeds, where \\\lambda\\ and
\\\mu\\ are the Lamé parameters and \\\rho\\ is the shell density.

The scattered field is expanded as:

\$\$ f\_{bs} = -\frac{i}{k} \sum\_{m=0}^{\infty} (2m+1) b_m
P_m(\cos\theta) \$\$

where \\k\\ is the wavenumber in the surrounding fluid, \\P_m\\ is the
Legendre polynomial of order \\m\\, and \\b_m\\ are the modal
coefficients determined by the shell and fluid properties.

The modal coefficients \\b_m\\ are determined by enforcing boundary
conditions at the shell-fluid interfaces, resulting in a system of six
equations for each mode. These are solved using Cramer's rule as the
ratio of two 6 \\\times\\ 6 determinants.

At mode 0:

\$\$ b_m(0) = \frac{ \left\| \begin{array}{ccccc} a_1 & \alpha\_{12} &
\alpha\_{14} & 0 & 0 \\ a_2 & \alpha\_{22} & \alpha\_{24} & 0 & 0 \\ 0 &
\alpha\_{42} & \alpha\_{44} & \alpha\_{46} & 0 \\ 0 & \alpha\_{52} &
\alpha\_{54} & \alpha\_{56} & 0 \\ 0 & \alpha\_{62} & \alpha\_{64} &
\alpha\_{66} & 0 \end{array} \right\| }{ \left\| \begin{array}{ccccc}
\alpha\_{11} & \alpha\_{12} & \alpha\_{14} & 0 & 0 \\ \alpha\_{21} &
\alpha\_{22} & \alpha\_{24} & 0 & 0 \\ 0 & \alpha\_{42} & \alpha\_{44} &
\alpha\_{46} & 0 \\ 0 & \alpha\_{52} & \alpha\_{54} & \alpha\_{56} & 0
\\ 0 & \alpha\_{62} & \alpha\_{64} & \alpha\_{66} & 0 \end{array}
\right\| } \$\$

At modes greater than 0:

\$\$ b_m(m) = -i^m (2m+1) \frac{ \left\| \begin{array}{cccccc} a_1 &
\alpha\_{12} & \alpha\_{13} & \alpha\_{14} & \alpha\_{15} & 0 \\ a_2 &
\alpha\_{22} & \alpha\_{23} & \alpha\_{24} & \alpha\_{25} & 0 \\ 0 &
\alpha\_{32} & \alpha\_{33} & \alpha\_{34} & \alpha\_{35} & 0 \\ 0 &
\alpha\_{42} & \alpha\_{43} & \alpha\_{44} & \alpha\_{45} & \alpha\_{46}
\\ 0 & \alpha\_{52} & \alpha\_{53} & \alpha\_{54} & \alpha\_{55} &
\alpha\_{56} \\ 0 & \alpha\_{62} & \alpha\_{63} & \alpha\_{64} &
\alpha\_{65} & \alpha\_{66} \end{array} \right\| }{ \left\|
\begin{array}{cccccc} \alpha\_{11} & \alpha\_{12} & \alpha\_{13} &
\alpha\_{14} & \alpha\_{15} & 0 \\ \alpha\_{21} & \alpha\_{22} &
\alpha\_{23} & \alpha\_{24} & \alpha\_{25} & 0 \\ 0 & \alpha\_{32} &
\alpha\_{33} & \alpha\_{34} & \alpha\_{35} & 0 \\ 0 & \alpha\_{42} &
\alpha\_{43} & \alpha\_{44} & \alpha\_{45} & \alpha\_{46} \\ 0 &
\alpha\_{52} & \alpha\_{53} & \alpha\_{54} & \alpha\_{55} & \alpha\_{56}
\\ 0 & \alpha\_{62} & \alpha\_{63} & \alpha\_{64} & \alpha\_{65} &
\alpha\_{66} \end{array} \right\| } \$\$

The elements of these matrices depend on the shell's elastic moduli,
thickness, densities, and the acoustic properties of the interior and
exterior fluids. The exact values for each element can be calcaulted
using equations provided by Goodman and Stern (1962) and Stanton (1990).
Since \\\theta = \pi\\ in the backscattering case, the equation for
\\f\_{bs}\\ becomes:

\$\$ f\_{bs} = -\frac{i}{k} \sum\_{m=0}^{\infty} (2m+1) b_m
P_m(\cos\theta) \\ \phantom{f\_{bs}} = -\frac{i}{k} \sum\_{m=0}^{\infty}
(2m+1) b_m P_m(-1) \\ \phantom{f\_{bs}} = -\frac{i}{k}
\sum\_{m=0}^{\infty} (2m+1) b_m (-1)^m \$\$

## Implementation

***`C++`***

The computation for \\b_m\\ was done in `C++` due to relatively large
computational costs with increasing \\ka\\ (and subsequently larger
limits for \\m\\) and the efficiency in solving the systems of linear
equations containing complex values. Since each matrix is a square, this
implementation specifically utilizes lower-upper (LU) decomposition to
break down the matrices into the products of the resulting lower and
upper triangular matrices. However, this means that the determinant
computations are sensitive to ill-conditioned matrices that can amplify
numerical errors.

There are guards in-place that partially address singularity issues when
the denominator is 0. The algorithm does not handle near-singular
matrices directly, but it will raise a warning when a matrix is
ill-conditioned. This is determined based on the pivot ratio from the
calculated condition number. Partial pivoting and row-scaling are also
incorporated to improve numerical stability and reduce the effect of
high-leverage values in a matrix.

***Modal Truncation***

The maximum number of terms for \\n\\ is chosen as \\k_w a\_{shell} +
10\\ (rounded to the nearest integer), which is sufficient for
convergence in most practical cases.

## References

Anderson, V.C. (1950). Sound scattering from a fluid sphere. The Journal
of The Acoustical Society of America, 22: 426–431.

Gaunaurd, G.C., and Wertman, W. (1991). Transient acoustic scattering by
fluid-loaded elastic shells. International Journal of Solids and
Structures, 27: 699-811.

Stanton, T.K. (1990). Sound scattering by spherical and elongated
shelled bodies. The Journal of the Acoustical Society of America, 88:
1619-1633.

## See also

See the boundary conditions documentation for more details on how
elastic-shelled scattering relates to other boundary conditions,
[`target_strength`](https://brandynlucca.github.io/acousticTS/reference/target_strength.md),
[`ESS`](https://brandynlucca.github.io/acousticTS/reference/ESS-class.md),
[`Sphere`](https://brandynlucca.github.io/acousticTS/reference/Sphere-class.md),
[`sphere`](https://brandynlucca.github.io/acousticTS/reference/sphere.md)
