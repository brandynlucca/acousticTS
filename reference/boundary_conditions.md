# Boundary conditions for scatterers

Details of all boundary types supported by various target strength
models.

## Simple boundaries

**Fixed rigid**

*Function alias: `"fixed_rigid"`*

*Models:
[`FCMS`](https://brandynlucca.github.io/acousticTS/reference/FCMS.md),
[`PSMS`](https://brandynlucca.github.io/acousticTS/reference/PSMS.md),
[`SPHMS`](https://brandynlucca.github.io/acousticTS/reference/SPHMS.md)*

Fixed rigid boundary indicative of a "hard" scatterer. This corresponds
to the mechanical impedance: \$\$ Z_s = \frac{p}{v_n}, \$\$ where
\\Z_s\\ is the mechanical impedance, \\p\\ is the acoustic pressure, and
\\v_n\\ is the normal component of the particle velocity along the
scattering surface. For a rigid boundary, \\v_n\\ vanishes, which
results in a Neumann boundary condition on the surface: \$\$
\frac{\partial p}{\partial n} = \nabla p \cdot \textbf{n} = 0, \$\$
where \\\nabla p\\ is the pressure gradient and \\\textbf{n}\\ is the
unit normal vector to the surface that points outward from the scatterer
into the fluid. The resulting infinite mechanical impedance corresponds
to a scattering surface that does not move in response to the incident
sound field. If the surface is sufficiently smooth and large relative to
the acoustic wavelength, this boundary condition permits specular
reflection in the geometric (high-frequency) limit.

**Pressure release**

*Function alias: `"pressure_release"`*

*Models:
[`FCMS`](https://brandynlucca.github.io/acousticTS/reference/FCMS.md),
[`PSMS`](https://brandynlucca.github.io/acousticTS/reference/PSMS.md),
[`SPHMS`](https://brandynlucca.github.io/acousticTS/reference/SPHMS.md)*

Pressure-release boundary indicative of a "soft" scatterer. This
corresponds to the limiting case of vanishing \\Z_s\\: \$\$ Z_s =
\frac{p}{v_n} \to 0. \$\$ In this limit, the boundary cannot sustain
acoustic pressure fluctuations, and the pressure at the surface
vanishes, resulting in a Dirichlet boundary condition: \$\$ p = 0. \$\$
Physically, this corresponds to a surface that is free to move or
compress in response to the incident sound field, such that pressure
fluctuations are instantaneously relieved. For a sufficiently smooth
surface and in the geometric limit, this boundary condition also permits
specular reflection, though with a phase inversion relative to a
fixed-rigid boundary.

**Fluid-filled**

*Function aliases: `"liquid_filled"`, `"gas_filled"`*

*Models:
[`FCMS`](https://brandynlucca.github.io/acousticTS/reference/FCMS.md),
[`PSMS`](https://brandynlucca.github.io/acousticTS/reference/PSMS.md),
[`SPHMS`](https://brandynlucca.github.io/acousticTS/reference/SPHMS.md)*

A fluid-filled scatterer with boundary conditions indicative of a
compliant cavity. At the seawater-fluid (i.e., external-internal)
interface, the acoustic field satisfies continuity of pressure and
normal velocity: \$\$ p\_{\text{inside}} = p\_{\text{outside}}, \quad
v\_{n,\text{inside}} = v\_{n,\text{outside}}, \$\$ where \\p\\ is the
acoustic pressure and \\v_n\\ is the normal component of particle
velocity. Physically, the surface is free to move in response to the
incident sound field according to the fluid properties within the
scatterer. The resulting scattering depends strongly on the
compressibility and density contrast between the interior fluid and
surrounding medium. For a smooth surface in the geometric limit, this
can produce specular reflection components as well as
resonance-modulated backscatter.

## Inelastic shelled boundaries

**Pressure release interior**

*Function alias: `"shelled_pressure_release"`*

*Models:
[`SPHMS`](https://brandynlucca.github.io/acousticTS/reference/SPHMS.md)*

An inelastic (non-deforming)shell surrounding a pressure-release
interior. The inelastic shell does not support shear or elastic waves.
It is therefore assumed that the mechanical impedance of the shell is
large enough to prevent motion normal to its surface. At the
shell-exterior interface, continuity of normal velocity and pressure is
assumed: \$\$ p\_{\text{shell}} = p\_{\text{outside}}, \quad
v\_{n,\text{shell}} = v\_{n,\text{outside}}. \$\$ At the inner interface
with the pressure-release interior, the pressure vanishes:
\$\$p\_{\text{shell}} = 0, \quad v\_{n,\text{shell}} \neq 0.\$\$ The
shell must accommodate these two conditions simultaneously. The outer
interface is mechanically compliant with the surrounding fluid, and the
inner interface allows free motion consistent with a zero-pressure
cavity. The shell can support motion consistent with its material
properties whereas the interior cannot sustain pressure (i.e., a
compliant cavity).

**Fluid-filled interior**

*Function aliases: `"shelled_liquid", "shelled_gas"`*

*Models:
[`SPHMS`](https://brandynlucca.github.io/acousticTS/reference/SPHMS.md)*

An inelastic shell encasing a fluid-filled interior. The shell is
assumed to have infinite mechanical impedance compared to both the
surrounding and interior fluids, so it does not deform elastically or
support shear waves. At the external interface represented by the
shell-surrounding fluid, the acoustic field satisfies a Neumann boundary
condition for the shell: \$\$ \frac{\partial p\_{\text{shell}}}{\partial
n} = 0, \$\$ where \\p\_{\text{shell}}\\ is the pressure just inside the
shell and \\n\\ is the unit normal pointing outward. At the
shell-interior fluid interface, the acoustic field satisfies continuity
of pressure and normal velocity with the surrounding medium: \$\$
p\_{\text{interior}} = p\_{\text{shell}}, \quad v\_{n,\text{interior}} =
v\_{n,\text{shell}}, \$\$ where \\p\_{\text{interior}}\\ and
\\v\_{n,\text{interior}}\\ are the pressure and normal velocity of the
interior fluid, respectively, and \\v\_{n,\text{shell}} \approx 0\\ due
to the inelasticity of the shell. The inelastic shell thus imposes a
rigid boundary that does not remove in response to the incident sound
field. The interior fluid responds freely to the pressure applied by the
surrounding shell according to its compressibility and material
properties.

## Elastic shelled boundaries

**Elastic shell boundary**

*Models:
[`ESSMS`](https://brandynlucca.github.io/acousticTS/reference/ESSMS.md),
[`VESMS`](https://brandynlucca.github.io/acousticTS/reference/VESMS.md)*

An elastic shell encasing a fluid or gas-filled interior. The shell is
assumed to have finite thickness and elastic properties defined by
Lam&eacute parameters \\\lambda\\ and \\\mu\\ and density \\\rho_s\\.
Unlike inelastic shells, the elastic shell supports both longitudinal
(P) and transverse (S) waves, resulting in coupled motion at the shell
interfaces.

The acoustic field along the exterior boundary (shell-fluid interface)
satisfies: \$\$ p\_{\text{fluid}} = \sigma\_{rr}, \quad
v\_{n,\text{fluid}} = u_r, \$\$ where \\\sigma\_{rr}\\ is the radial
stress in the shell, and \\u_r\\ is the radial displacement at the
exterior shell surface. The acoustic field along the interior boundary
(shell-interior interface) satisfies: \$\$ p\_{\text{interior}} =
\sigma\_{rr}, \quad v\_{n,\text{interior}} = u_r, \$\$ ensuring
continuity of normal displacement and stress across the interface.

The displacement field in the shell is expressed as the sum of scalar
and vector potentials: \$\$ \mathbf{u} = \nabla \phi + \nabla \times
\mathbf{\Psi}, \$\$ where \\\phi\\ and \\\Psi\\ satisfy the Helmholtz
equations with wave numbers \\k_L = \omega / c_L\\ and \\k_T = \omega /
c_T\\, and \\c_L = \sqrt{(\lambda + 2\mu)/\rho_s}\\, \\c_T =
\sqrt{\mu/\rho_s}\\.

The resulting system of equations couples the incident acoustic field
with elastic waves in the shell, producing scattering coefficients
\\b_m\\ that depend on shell thickness, elastic moduli, density, and
interior/exterior fluid properties. This formulation generalizes the
inelastic (rigid) shell case by allowing finite radial motion and
elastic wave propagation in the shell.

## Weak scattering

**Fluid-like boundary**

*Models:
[`DWBA`](https://brandynlucca.github.io/acousticTS/reference/DWBA.md),
[`PCDWBA`](https://brandynlucca.github.io/acousticTS/reference/PCDWBA.md),
[`FCMS`](https://brandynlucca.github.io/acousticTS/reference/FCMS.md),
[`BCMS`](https://brandynlucca.github.io/acousticTS/reference/BCMS.md),
[`KRM`](https://brandynlucca.github.io/acousticTS/reference/KRM.md),
[`PSMS`](https://brandynlucca.github.io/acousticTS/reference/PSMS.md),
[`SDWBA`](https://brandynlucca.github.io/acousticTS/reference/SDWBA.md),
[`SPHMS`](https://brandynlucca.github.io/acousticTS/reference/SPHMS.md)*

Represents a scatterer whose acoustic response is dominated by the
surrounding fluid, i.e., the scatterer does not present a significant
impedance contrast. In this regime:

\$\$ p\_{\text{scatterer}} \approx p\_{\text{fluid}} \$\$

where \\p\_{\text{scatterer}}\\ and \\p\_{\text{fluid}}\\ are the
pressures in the scatterer and surrounding medium, respectively. The
scattered field is weak and primarily arises from small density or
compressibility contrasts.
