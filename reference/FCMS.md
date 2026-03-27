# Finite cylinder modal series (FCMS) solution

Calculates the far-field scattering amplitude and related quantities for
a finite cylinder using the modal series solution, supporting various
boundary conditions (rigid, pressure-release, liquid-filled, and
gas-filled).

## Usage

This model is accessed via:

    target_strength(
      ...,
      model="fcms",
      boundary,
      sound_speed_sw,
      density_sw,
      m_limit
    )

## Arguments

- `boundary`:

  Boundary condition at a cylindrical surface. One of `"fixed_rigid"`,
  `"pressure_release"`, `"liquid_filled"`, or `"gas_filled"`. See the
  boundary conditions documentation for more details on these different
  boundary conditions.

- `sound_speed_sw`:

  Seawater sound speed (\\m~s^{-1}\\).

- `density_sw`:

  Seawater density (\\kg~m^{-3}\\).

- `m_limit`:

  Optional model truncation limit used to cap the number of modes in the
  numerical calculation.

## Theory

The modal series solution for a finite cylinder expresses the
backscattering amplitude as:

\$\$ f\_{bs} = -\frac{L}{\pi} \frac{\sin(k L \cos \theta)}{k L \cos
\theta} \sum\_{m=0}^{\infty} i^{m+1} B_m \$\$

where \\L\\ is the cylinder length, \\k\\ is the wavenumber in the
surrounding medium, and \\\theta\\ is the angle between the cylinder
axis and the incident wave direction. The coefficients \\B_m\\ depend on
the boundary condition at the cylinder surface.

**Boundary Conditions and Modal Coefficients**

- **Rigid (fixed) cylinder:** The normal velocity at the surface is
  zero. The modal coefficient is: \$\$ B_m = (-1)^m \epsilon_m
  \frac{J_m'(K a)}{H_m^{(1)'}(K a)} \$\$ where \\J_m\\ and \\H_m^{(1)}\\
  are the cylindrical Bessel and Hankel functions of order \\m\\, the
  prime denotes differentiation with respect to the argument, \\K = k
  \sin \theta\\, and \\a\\ is the cylinder radius. The Neumann factor is
  \\\epsilon_0 = 1\\, \\\epsilon_m = 2\\ for \\m \geq 1\\.

- **Pressure-release cylinder:** The acoustic pressure at the surface is
  zero. The modal coefficient is: \$\$ B_m = (-1)^m \epsilon_m
  \frac{J_m(K a)}{H_m^{(1)}(K a)} \$\$

- **Fluid-filled (or gas-filled) cylinder:** Both pressure and normal
  velocity are nonzero at the surface. The modal coefficient is: \$\$
  B_m = -\epsilon_m / (1 + i C_m) \$\$ where \$\$ C_m = \frac{ \left\[
  J_m'(K' a) Y_m(K a) \right\] / \left\[ J_m(K' a) J_m'(K a) \right\] -
  g h \left\[ Y_m'(K a) / J_m'(K a) \right\] }{ \left\[ J_m'(K' a)
  J_m(K a) \right\] / \left\[ J_m(K' a) J_m'(K a) \right\] - g h } \$\$
  Here, \\Y_m\\ is the cylindrical Bessel function of the second kind,
  \\K' = K / h\\, \\g\\ is the density contrast (target to medium), and
  \\h\\ is the sound speed contrast (target to medium).

**Modal Truncation**

The modal sum is truncated at a maximum order determined by \\m\_{\max}
= \max(\lceil k a \rceil) + 10\\, which is sufficient for convergence in
most practical cases.

## References

Stanton, T.K. (1988). Sound scattering by cylinders of finite length. I.
Fuid cylinders. The Journal of the Acoustical Society of America, 83:
55-63.

Stanton, T.K. (1989). Sound scattering by cylinders of finite length.
III. Deformed cylinders. The Journal of the Acoustical Society of
America, 85: 232-237.

## See also

[`target_strength`](https://brandynlucca.github.io/acousticTS/reference/target_strength.md),
[`FLS`](https://brandynlucca.github.io/acousticTS/reference/FLS-class.md),
[`GAS`](https://brandynlucca.github.io/acousticTS/reference/GAS-class.md),
[`Cylinder`](https://brandynlucca.github.io/acousticTS/reference/Cylinder-class.md),
[`cylinder`](https://brandynlucca.github.io/acousticTS/reference/cylinder.md)
