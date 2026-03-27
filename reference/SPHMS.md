# Spherical modal series (SPHMS) solution

Spherical modal series (SPHMS) backscatter model. Supports rigid,
pressure-release, fluid-filled, shelled (pressure-release/liquid/gas).
Frequencies in Hz; sound speed in m/s; density in kg/m^3. Material
properties may be contrasts or absolute (contrasts derived relative to
seawater). Full equations and boundary-condition details are documented
in the SPHMS vignette.

## Usage

This model is accessed via:

    target_strength(
      ...,
      model="sphms",
      boundary,
      sound_speed_sw,
      density_sw,
      m_limit
    )

## Arguments

- `boundary`:

  Boundary condition at a spherical surface. One of `"fixed_rigid"`,
  `"pressure_release"`, `"liquid_filled"`, `"gas_filled"`,
  `"shelled_pressure_release"`, `"shelled_liquid"`, or `"shelled_gas"`.
  See the boundary conditions documentation for more details on these
  different boundary conditions.

- `sound_speed_sw`:

  Seawater sound speed (\\m~s^{-1}\\).

- `density_sw`:

  Seawater density (\\kg~m^{-3}\\).

- `m_limit`:

  Optional model truncation limit used to cap the number of modes in the
  numerical calculation.

## Theory

The modal series solution for a sphere expands backscatter from a
spherical scatterer in spherical Bessel/Hankel modes. The form function
is:

\\ f\_{bs} = -\frac{i}{k_w} \sum\_{n=0}^{\infty} (-1)^n (2n+1) A_n, \\

## References

Anderson, V.C. (1950). Sound scattering from a fluid sphere. The Journal
of The Acoustical Society of America, 22: 426–431.

## See also

[`target_strength`](https://brandynlucca.github.io/acousticTS/reference/target_strength.md),
[`FLS`](https://brandynlucca.github.io/acousticTS/reference/FLS-class.md),
[`GAS`](https://brandynlucca.github.io/acousticTS/reference/GAS-class.md),
[`ESS`](https://brandynlucca.github.io/acousticTS/reference/ESS-class.md),
[`Sphere`](https://brandynlucca.github.io/acousticTS/reference/Sphere-class.md),
[`sphere`](https://brandynlucca.github.io/acousticTS/reference/sphere.md)
