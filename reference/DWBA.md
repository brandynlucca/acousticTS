# Distorted wave Born approximation (DWBA) for weak scatterers

Calculates backscatter for fluid-like weak scatterers using the
distorted wave Born approximation (DWBA). Frequencies in Hz; sound speed
in m/s; density in kg/m^3. Material properties may be provided as
contrasts or absolute values (contrasts derived relative to seawater).

## Usage

This model is accessed via:

    target_strength(
      ...,
      model="dwba",
      sound_speed_sw,
      density_sw
    )

## Arguments

- `sound_speed_sw`:

  Seawater sound speed (\\m~s^{-1}\\).

- `density_sw`:

  Seawater density (\\kg~m^{-3}\\).

## Theory

The DWBA approach is derived under the weak scattering assumption,
meaning that the differences in compressibility (\\\kappa\\) and density
(\\\rho\\) between the scatterer and surrounding fluid are sufficiently
small to linearize the acoustic scattering problem. This linearization
allows the scattered field to be expressed as an integral over the
scatterer volume.

The far-field backscattering amplitude is given by:

\$\$ f\_{bs} = \frac{k_1}{4\pi} \int \int \int\limits_v
(\gamma\_\kappa - \gamma\_\rho) e^{2i k_2 \cdot r_v} dv \$\$

where

\$\$ \gamma\_\kappa = \frac{\kappa_2 - \kappa_1}{\kappa_1}, \quad
\gamma\_\rho = \frac{\rho_2 - \rho_1}{\rho_2}, \$\$

and

\$\$ \kappa = (\rho c^2)^{-1}. \$\$

where \\c\\ is the sound speed (m s⁻¹). For elongated bodies, the volume
integral reduces to a line intergral along the body axis:

\$\$ f\_{bs} = \frac{k_1}{4} \int\_{r\_{pos}} (\gamma\_\kappa -
\gamma\_\rho) e^{2i k_2 \cdot r\_{pos}} \frac{\text{J}\_1(2 k_2 a \cos
\beta\_{tilt})}{\cos \beta\_{tilt}} \|dr\_{pos}\|, \$\$

where \\a\\ is the local radius and \\\beta\_{tilt}\\ the local tilt
angle. The cylindrical Bessel function of the first kind of order 1,
\\\text{J}\_1\\, accounts for the circular cross-section of each
segment. The wavenumber \\k_2\\ is evaluated inside the scatterer.

**Assumptions**

The DWBA assumes that the scatterer is weakly inhomogenous where the
material property contrasts for the interior (\\c_2\\, \\\rho_2\\) and
surrounding fluid (\\c_1\\, \\\rho_1\\) where:

\$\$ g = \frac{\rho_2}{\rho_1} \approx 1, \quad h = \frac{c_2}{c_1}
\approx 1. \$\$

In practice, \\c_2\\ and \\\rho_2\\ within 5% of \\c_1\\ and \\\rho_1\\,
respectively, can be considered to be sufficent for the weak scattering
assumption whereby:

\$\$ \|h - 1\| \le 0.05, \quad \|g - 1\| \le 0.05. \$\$

This model also assumes that the scatterer's body has no sharp edges or
irregularities, and its cross-section is symmetric around the
longitudinal axis. This allows the scattering integral to be reduced to
a line integral along the body axis in the first place, simplifying the
computation of phase contributions from different segments. Moreover,
this enables the body to be discretized into along-axis segments that
can approximate arbitrary body shapes. Since the body is axisymmetric
and smooth, this further allows the DWBA to be applied for arbitrary
orientation angles without additional correction terms.

The DWBA provides a first-order approximation that neglects multiple
scattering within the body whereby secondary interactions between
different parts of the body are ignored. This is valid for weakly
scattering objects where the amplitude of scattered waves is small and
is consistent with the Born approximation in wave physics. This model
also disregards any elastic or shelled effects, treating scatterers are
being purely fluid-like. Consequently, the lack of internal elasticity
means there is no support for shear waves or resonances due to solid
boundaries. Mathematically, this means that only the contrasts in
compressibility and density between the scatterer and the surrounding
medium contribute to the scattered field.

## Implementation

The model extracts geometric and acoustic parameters from the input
object, constructs rotation matrices and wavenumber matrices, and
evaluates the integral numerically for each frequency of interest. The
DWBA is computationally efficient and handles elongated, weakly
scattering targets such as zooplankton and small fish.

## References

Morse, P.M., and Ingard, K.U. (1968). Theoretical Acoustics. Princeton
University Press.

Stanton, T.K., Chu, D., and Wiebe, P.H. (1998). Sound scattering by
several zooplankton groups. II. Scattering models. The Journal of the
Acoustical Society of America, 103, 236-253.

## See also

See the boundary conditions documentation for more details on how weak
scattering relates to other boundary conditions,
[`target_strength`](https://brandynlucca.github.io/acousticTS/reference/target_strength.md),
[`FLS`](https://brandynlucca.github.io/acousticTS/reference/FLS-class.md)
