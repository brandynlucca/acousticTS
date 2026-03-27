# High-pass approximation (HPA) scattering model

Computes the far-field scattering amplitude and related quantities
spherical and elongated scatterers using a high-pass approximation
model. The high-pass model provides a simplified analytical expression
for the backscattering cross-section that is valid for all values of
\\ka\\, where \\k\\ is the acoustic wavenumber and \\a\\ is a
characteristic dimension of the scatterer. The model is named after the
analogy to a two-pole high-pass filter in electrical circuit theory,
where the frequency response resembles the backscattering behavior as a
function of \\ka\\. The high-pass model is computationally efficient and
captures the essential physics of acoustic scattering from weakly
scattering bodies, including spheres, prolate spheroids, and straight or
bent cylinders. Two implementations are available: the Johnson (1977)
formulation for fluid spheres, and the Stanton (1989) generalization for
spheres, prolate spheroids, and cylinders.

## Usage

This model is accessed via:

    target_strength(
      ...,
      model = "HPA",
      method,
      deviation_fun,
      null_fun,
      sound_speed_sw,
      density_sw
    )

## Arguments

- `method`:

  Method for computing the high-pass model. One of `"johnson"` (Johnson
  1977 formulation for fluid spheres) or `"stanton"` (Stanton 1989
  formulation for spheres, prolate spheroids, and cylinders).

- `deviation_fun`:

  Function or scalar specifying the expected deviation in the phase of
  the scattered wave due to shape complexity, flexure, or other factors.
  Expressed as a function of \\ka\\, or as a scalar constant. Default is
  1.

- `null_fun`:

  Function or scalar specifying the expected reduction in scattering
  amplitude due to nulls in the scattering pattern. Expressed as a
  function of \\ka\\, or as a scalar constant. Default is 1.

- `sound_speed_sw`:

  Seawater sound speed (\\m~s^{-1}\\).

- `density_sw`:

  Seawater density (\\kg~m^{-3}\\).

## Theory

The high-pass model is based on the observation that the backscattering
cross-section of a scatterer can be approximated by a simple analytical
expression that captures the low-frequency Rayleigh scattering regime
(where \\\sigma\_{bs} \propto (ka)^4\\) and the high-frequency geometric
scattering regime (where \\\sigma\_{bs}\\ approaches a constant). The
transition between these regimes is governed by an auxiliary material
property parameter, \\\alpha\\, which depends on the density contrast
\\g\\ and sound speed contrast \\h\\ between the scatterer and the
surrounding medium. For spheres and prolate spheroids, the auxiliary
parameter is given by

\$\$ \alpha\_{\pi s} = \frac{1 - g h^2}{3 g h^2} + \frac{1 - g}{1 + 2g}
\$\$

and for cylinders and prolate spheroids (in the cylindrical limit), the
auxiliary parameter is

\$\$ \alpha\_{\pi c} = \frac{1 - g h^2}{2 g h^2} + \frac{1 - g}{1 + g}
\$\$

where \\g\\ is the density contrast (target to medium) and \\h\\ is the
sound speed contrast (target to medium). The reflection coefficient at
the scatterer surface is given by

\$\$ \mathcal{R} = \frac{gh - 1}{gh + 1} \$\$

The Johnson (1977) high-pass model for a fluid sphere is expressed as

\$\$ \sigma\_{bs} = \frac{a^2 (ka)^4 \alpha\_{\pi s}^2}{1 + \frac{3}{2}
(ka)^4} \$\$

where \\a\\ is the spherical radius. This formulation is valid for all
\\ka\\ and provides a good approximation for fluid spheres. The Stanton
(1989) high-pass model generalizes this approach to spheres, prolate
spheroids, and cylinders, and incorporates empirical terms to account
for nulls in the scattering pattern and deviations due to shape
complexity. For a sphere, the Stanton formulation is

\$\$ \sigma\_{bs} = \frac{ a^2 (ka)^4 \alpha\_{\pi s}^2 \mathcal{G} }{
1 + \frac{4(ka)^4 \alpha\_{\pi s}^2}{\mathcal{R}^2 \mathcal{F}} } \$\$

where \\\mathcal{G}\\ is a null function that accounts for reductions in
scattering amplitude at certain frequencies, and \\\mathcal{F}\\ is a
deviation function that accounts for phase variability. For a prolate
spheroid, the Stanton formulation is

\$\$ \sigma\_{bs} = \frac{ \frac{1}{9} L^2 (ka)^4 \alpha\_{\pi c}^2
\mathcal{G} }{ 1 + \frac{\frac{16}{9}(ka)^4 \alpha\_{\pi c}^2}
{\mathcal{R}^2 \mathcal{F}} } \$\$

where \\L\\ is the length of the prolate spheroid. For a straight
cylinder at angle \\\theta\\, the Stanton formulation is

\$\$ \sigma\_{bs} = \frac{ \frac{1}{4} L^2 (Ka)^4 \alpha\_{\pi c}^2 s^2
\mathcal{G} }{ 1 + \frac{\pi (Ka)^3 \alpha\_{\pi c}^2}{\mathcal{R}^2
\mathcal{F}} } \$\$

where \\K = k \sin \theta\\, \\a\\ is the cylindrical radius, and \\s =
\sin(kL \cos \theta) / (kL \cos \theta)\\ accounts for the finite length
and orientation of the cylinder. For a bent cylinder with radius of
curvature \\\rho_c\\ relative to the cylinder length \\L\\, the Stanton
formulation is

\$\$ \sigma\_{bs} = \frac{ \frac{1}{4} L^2 (ka)^4 \alpha\_{\pi c}^2
\mathcal{H}^2 \mathcal{G} }{ 1 + \frac{ L^2 (ka)^4 \alpha\_{\pi c}^2
\mathcal{H}^2 }{ \rho_c a \mathcal{R}^2 \mathcal{F} } } \$\$

where \\\mathcal{H} = \frac{1}{2} + \frac{1}{2} (\rho_c / L) \sin(L /
\rho_c)\\ is an effective length factor that accounts for the curvature.
The empirical functions \\\mathcal{F}\\ and \\\mathcal{G}\\ are
specified by the user as functions of \\ka\\, and can be used to fit the
model to measured or numerically simulated scattering data. In all
cases, the wavenumber satisfies \\kr \gg 1\\ and \\L \ll 2\sqrt{r
\lambda}\\, where \\r\\ is the distance from the scatterer to the
receiver and \\\lambda\\ is the acoustic wavelength, ensuring that the
far-field and elongated object assumptions are valid.

## References

Johnson, R.K. (1977). Sound scattering from a fluid sphere revisited.
The Journal of the Acoustical Society of America, 61: 375-377.

Johnson, R.K. (1978). Erratum: Sound scattering from a fluid sphere
revisited. The Journal of the Acoustical Society of America, 63: 626.

Stanton, T.K. (1989). Simple approximate formulas for backscattering of
sound by spherical and elongated objects. The Journal of the Acoustical
Society of America, 86: 1499-1510.
