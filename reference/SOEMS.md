# Solid elastic (calibration) sphere modal series (SOEMS) solution

Calculates the far-field scattering amplitude and related quantities for
a solid elastic (calibration) sphere using a modal series solution.

## Usage

This model is accessed via:

    target_strength(
      ...,
      model="calibration",
      sound_speed_sw,
      density_sw,
      adaptive = TRUE
    )

## Arguments

- `sound_speed_sw`:

  Seawater sound speed (\\m~s^{-1}\\).

- `density_sw`:

  Seawater density (\\kg~m^{-3}\\).

- `adaptive`:

  Logical. If `TRUE`, extend the partial-wave sum beyond the initial
  \\\mathrm{round}(ka)+10\\ modal cap until the tail term falls below
  the internal convergence threshold. If `FALSE`, use the original fixed
  modal cutoff only.

## Theory

The calibration sphere model computes acoustic scattering from a solid
elastic sphere by expanding the incident and scattered fields in terms
of spherical Bessel and Hankel functions and Legendre polynomials. Both
compressional and shear waves within the sphere are included, and the
appropriate boundary conditions are enforced at the sphere-water
interface.

The dimensionless frequency parameter is defined as \\q = ka\\, where
\\k\\ is the wavenumber in water and \\a\\ is the sphere radius. The
longitudinal and transverse wave numbers inside the sphere are

\$\$ q_1 = \frac{qc}{c_1} \\ q_2 = \frac{qc}{c_2}, \$\$

where \\c_1\\ and \\c_2\\ are the longitudinal and transverse sound
speeds in the sphere, respectively.

The far-field backscattering form function is given by:

\$\$ f\_\infty(q) = -\frac{2}{q} \sum\limits\_{\ell=0}^{\infty}
(-1)^\ell (2\ell+1) \sin \eta\_\ell \exp(i \eta\_\ell) \$\$

The phase angle for each mode is then given by:

\$\$ \tan \eta\_\ell = -\frac{B_2 j\_\ell'(q) - B_1 j\_\ell(q)}{B_2
y\_\ell'(q) - B_1 y\_\ell(q)} \$\$

where the phase angle \\\eta\_\ell\\ is determined by the boundary
conditions and the material properties of the sphere and surrounding
fluid. The phase angle \\\eta\_\ell\\ is calculated using a series of
coefficients. The auxilliary quantities used for determining the
boundary conditions and subsequently \\\eta\_\ell\\ are reported in
MacLennan (1981).

## References

Hickling, R. (1962). Analysis of echoes from a solid elastic sphere in
water. The Journal of the Acoustical Society of America, 34: 1582-1592.

MacLennan D. N. (1981). The theory of solid spheres as sonar calibration
targets. Scottish Fisheries Research No. 22, Department of Agriculture
and Fisheries for Scotland.

## See also

[`target_strength`](https://brandynlucca.github.io/acousticTS/reference/target_strength.md),
[`CAL`](https://brandynlucca.github.io/acousticTS/reference/CAL-class.md),
[`Sphere`](https://brandynlucca.github.io/acousticTS/reference/Sphere-class.md),
[`sphere`](https://brandynlucca.github.io/acousticTS/reference/sphere.md)
