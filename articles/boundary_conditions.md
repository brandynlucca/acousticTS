# Scattering boundary conditions

## Introduction

The boundary families summarized here correspond to the classical fluid,
rigid, elastic, shell, and viscous-interface scattering literature
([Anderson 1950](#ref-anderson_sound_1950); [Faran
1951](#ref-faran_sound_1951); [Goodman and Stern
1962](#ref-goodman_reflection_1962); [Feuillade and Nero
1998](#ref-feuillade_nero_1998)).

Boundary conditions are the mathematical statements that turn a general
wave equation into a specific scattering problem. In target-strength
modeling, they specify what the target surface can and cannot do in
response to the incident sound field. A boundary may be perfectly rigid,
perfectly pressure releasing, fluid transmitting, weakly contrasting, or
dynamically elastic. Those are not just different algebraic cases. They
are different physical idealizations, and they produce different
scattering mechanisms even for the same gross geometry.

This page gathers the main boundary families that recur throughout the
theory articles. The aim is to explain both the mathematics and the
physics clearly enough that a reader can understand what is being
assumed before turning to any model-specific derivation.

It is best read alongside the shared [acoustic scattering
primer](https://brandynlucca.github.io/acousticTS/articles/acoustic-scattering-primer/acoustic-scattering-primer.md)
and [notation
guide](https://brandynlucca.github.io/acousticTS/articles/notation-and-symbols/notation-and-symbols.md).
In the package theory pages, medium `1` is always the surrounding
seawater, medium `2` is the first target region encountered from the
outside, and deeper internal regions are numbered inward.

The main boundary families separate most naturally according to the
field quantity constrained most directly at the interface: pressure,
normal velocity, layered acoustic transmission, or elastic traction.

## Shared interface notation

Before specializing to any one interface, it helps to keep the common
two-medium notation explicit. If media `i` and `j` meet at an interface,
then the density and sound-speed contrasts are:

g\_{ij} = \frac{\rho_i}{\rho_j}, \qquad h\_{ij} = \frac{c_i}{c_j},

and the acoustic wavenumbers are:

k_i = \frac{\omega}{c_i}, \qquad k_j = \frac{\omega}{c_j}.

For a one-region target in seawater, the most common contrasts are
therefore g\_{21} and h\_{21}. For a shelled target, the shell-to-water
interface uses g\_{21} and h\_{21}, while an internal fluid relative to
the shell uses g\_{32} and h\_{32}.

## Simple boundaries

### Fixed rigid

A fixed rigid boundary is the limiting case of a scatterer whose surface
cannot move in the normal direction. In an inviscid surrounding fluid,
the acoustic pressure field p(\mathbf{x}, t) and the particle velocity
field \mathbf{v}(\mathbf{x}, t) satisfy the linearized momentum
equation:

\rho \frac{\partial \mathbf{v}}{\partial t} = -\nabla p,

where \rho is the surrounding-fluid density, \mathbf{x} is the position
vector, and \nabla p is the pressure gradient. If \mathbf{n} denotes the
outward unit normal vector to the target surface, then the normal
component of particle velocity is:

v_n = \mathbf{v} \cdot \mathbf{n},

where the dot denotes the Euclidean inner product. The mechanical
impedance at the boundary is then defined by:

Z_s = \frac{p}{v_n},

where Z_s is the surface impedance. For a rigid boundary, the surface
cannot move, so:

v_n = 0.

With harmonic time dependence, that implies the Neumann boundary
condition:

\frac{\partial p}{\partial n} = \nabla p \cdot \mathbf{n} = 0,

where \partial / \partial n denotes the derivative taken along the
outward normal direction. The exterior acoustic field is therefore
solved subject to:

\nabla^2 p + k^2 p = 0,

where k = \omega / c is the surrounding-fluid wavenumber, \omega is
angular frequency, and c is sound speed. Physically, the interface
stores no normal motion. In partial-wave and modal formulations, that
condition fixes the scattering coefficients by eliminating any solution
that would produce radial displacement of the boundary. In the
high-frequency limit, a sufficiently smooth rigid surface supports
specular reflection without the phase inversion associated with a
pressure-release boundary.

### Pressure release

A pressure-release boundary is the opposite limiting case. Here the
surface cannot sustain an acoustic pressure fluctuation, so the
interface behaves like a perfectly soft boundary. In impedance language:

Z_s = \frac{p}{v_n} \to 0.

For finite normal motion, the only consistent condition is:

p = 0.

This is a Dirichlet boundary condition. It removes any admissible
solution of the Helmholtz equation that carries nonzero pressure at the
interface. Physically, this corresponds to a surface that relieves
pressure fluctuations essentially instantaneously, as in an idealized
gas interface or other highly compliant boundary. Relative to the rigid
case, the reflected pressure undergoes a phase reversal. In modal
descriptions, the scattering coefficients are fixed by enforcing zero
pressure rather than zero normal velocity at the surface.

### Fluid-filled

For a fluid-filled scatterer, both the exterior and interior regions
support acoustic pressure fields. Let medium 1 denote the surrounding
fluid and medium 2 the interior fluid. Then:

\nabla^2 p_1 + k_1^2 p_1 = 0, \qquad \nabla^2 p_2 + k_2^2 p_2 = 0,

where p_1 and p_2 are the exterior and interior pressures, and:

k_j = \frac{\omega}{c_j}, \qquad j = 1, 2,

with c_j the sound speed in medium j. At the interface, continuity of
normal traction gives continuity of acoustic pressure:

p_1 = p_2.

Continuity of normal particle velocity gives:

\mathbf{v}\_1 \cdot \mathbf{n} = \mathbf{v}\_2 \cdot \mathbf{n},

where \mathbf{v}\_1 and \mathbf{v}\_2 are the particle-velocity fields
in the two media. Using the linearized momentum relation in each fluid:

\mathbf{v}\_j = -\frac{1}{i \omega \rho_j}\nabla p_j,

the normal-velocity condition can also be written as:

\frac{1}{\rho_1}\frac{\partial p_1}{\partial n} =
\frac{1}{\rho_2}\frac{\partial p_2}{\partial n}.

These are the standard transmission conditions for an acoustic
fluid-fluid boundary. Neither pressure nor normal motion is forced to
vanish. Instead, both media participate in the interface dynamics, and
the scattering is controlled by density and compressibility contrast
across the boundary.

## Inelastic shelled boundaries

The inelastic shelled cases are layered acoustic idealizations rather
than true elastic shell problems. The shell is treated as a
non-deforming, non-shear-supporting layer that constrains the adjacent
acoustic fields kinematically. That makes these boundaries useful as
simplified layered models, but they should be distinguished clearly from
a genuinely elastic shell that supports both compressional and shear
waves.

### Pressure release interior

In this configuration, an inelastic shell surrounds a pressure-release
interior. The shell is assumed not to support elastic wave propagation
and not to deform according to Hooke’s law. It therefore acts as a
kinematic constraint rather than as a dynamic elastic medium. A
convenient way to describe the surrounding-fluid and shell-adjacent
acoustic fields is to introduce an exterior pressure field
p\_{\text{exterior}} and an acoustic field immediately adjacent to the
shell, written here as p\_{\text{shell}}.

At the shell-exterior interface, continuity of pressure and normal
velocity with the surrounding fluid is enforced:

p\_{\text{shell}} = p\_{\text{exterior}}, \qquad v\_{n,\text{shell}} =
v\_{n,\text{exterior}},

where v\_{n,\text{shell}} and v\_{n,\text{exterior}} are the normal
particle velocities just inside and outside the shell interface,
respectively. At the shell-interior interface, the interior is pressure
releasing, so:

p\_{\text{shell}} = 0, \qquad v\_{n,\text{shell}} \neq 0.

This idealization means that the shell must transmit normal motion while
accommodating a zero-pressure cavity on its inner side. The shell itself
does not carry elastic shell-wave physics, but it does modify how the
exterior acoustic field couples to the interior cavity. In scattering
formulations, the resulting layered boundary conditions differ from both
rigid and ordinary fluid-filled cases and can introduce low-frequency
cavity-type resonances.

Mathematically, this case may be viewed as an acoustic transmission
problem with an outer continuity condition and an inner Dirichlet
condition. The shell is therefore acting as a geometric and kinematic
intermediary between the exterior fluid and a pressure-release cavity,
not as an elastic medium in its own right.

### Fluid-filled interior

For an inelastic shell enclosing a fluid-filled interior, the shell is
taken to have effectively infinite mechanical impedance and therefore
negligible normal displacement:

v\_{n,\text{shell}} = 0.

At the shell-exterior interface, that produces a Neumann condition on
the acoustic field:

\frac{\partial p\_{\text{shell}}}{\partial n} = 0.

At the shell-interior fluid interface, the interior fluid must satisfy
continuity of pressure and normal velocity with the shell-adjacent
field:

p\_{\text{interior}} = p\_{\text{shell}}, \qquad v\_{n,\text{interior}}
= v\_{n,\text{shell}},

where p\_{\text{interior}} is the interior-fluid pressure and
v\_{n,\text{interior}} is the interior-fluid normal velocity at the
inner shell surface. Because the shell is inelastic, one typically also
takes:

v\_{n,\text{shell}} \approx 0.

The interior fluid therefore experiences an effectively rigid boundary
while still supporting compressional waves of its own. In that sense,
the system is a layered acoustic cavity rather than an elastic shell.
The resulting scattering is controlled by acoustic layering and cavity
behavior rather than by shell-wave propagation.

To make that structure explicit, the interior acoustic field still
satisfies:

\nabla^2 p\_{\text{interior}} + k\_{\text{interior}}^2
p\_{\text{interior}} = 0,

The associated interior wavenumber is:

k\_{\text{interior}} = \frac{\omega}{c\_{\text{interior}}},

where c\_{\text{interior}} is the sound speed of the interior fluid. The
shell therefore modifies the interior field through boundary constraints
even though it does not itself support elastic motion.

## Elastic shelled boundaries

### Elastic shell boundary

An elastic shell supports both longitudinal and transverse elastic
motion. The displacement field \mathbf{u}(\mathbf{x}, t) within the
shell satisfies the Navier equation:

(\lambda + 2 \mu)\nabla(\nabla \cdot \mathbf{u}) - \mu \nabla \times
(\nabla \times \mathbf{u}) + \rho_s \frac{\partial^2
\mathbf{u}}{\partial t^2} = 0,

where \lambda and \mu are the Lam'e parameters and \rho_s is the shell
density. Using the Helmholtz decomposition:

\mathbf{u} = \nabla \phi + \nabla \times \mathbf{\Psi},

the displacement separates into longitudinal and transverse potentials,
each satisfying its own Helmholtz equation with wavenumbers:

k\_\ell = \frac{\omega}{c\_\ell}, \qquad k\_\tau =
\frac{\omega}{c\_\tau},

and wave speeds:

c\_\ell = \sqrt{\frac{\lambda + 2 \mu}{\rho_s}}, \qquad c\_\tau =
\sqrt{\frac{\mu}{\rho_s}}.

The elastic stress tensor is:

\sigma\_{ij} = \lambda \delta\_{ij}\nabla \cdot \mathbf{u} + 2 \mu
\varepsilon\_{ij},

where \delta\_{ij} is the Kronecker delta and \varepsilon\_{ij} is the
infinitesimal strain tensor:

\varepsilon\_{ij} = \frac{1}{2} \left( \frac{\partial u_i}{\partial
x_j} + \frac{\partial u_j}{\partial x_i} \right).

At a fluid-solid interface, the fluid pressure must balance the normal
traction in the shell, the normal fluid velocity must match the normal
shell velocity, and the tangential traction must vanish because the
fluid is inviscid. If u_n = \mathbf{u} \cdot \mathbf{n} is the shell
displacement in the normal direction, \sigma\_{nn} is the normal
traction, and \sigma\_{nt} is the tangential traction, then the
interface conditions may be written schematically as:

p = -\sigma\_{nn}, \qquad -\frac{1}{i \omega \rho_f}\frac{\partial
p}{\partial n} = -i \omega u_n, \qquad \sigma\_{nt} = 0,

where \rho_f is the fluid density. The exact sign of the
pressure-traction relation depends on the outward-normal convention, but
the physical content is unchanged. The shell carries both normal and
tangential motion, and that coupling makes the scattering problem much
richer than either the rigid or the inelastic-shell cases. The resulting
modal systems support resonance structure associated with shell
thickness, elastic moduli, and both the exterior and interior fluid
properties.

## Weak, fluid-like

The weak fluid-like case is best understood as an approximation regime
rather than as a sharply enforcing boundary type. The main assumption is
that the target differs only slightly from the surrounding medium, so
the acoustic field inside the scatterer remains close to the incident
field. In that limit:

p\_{\text{interior}} \approx p\_{\text{exterior}}.

The corresponding material contrasts are small. If \rho_1 and c_1 denote
the surrounding-fluid density and sound speed, and \rho_2 and c_2 denote
the target values, one typically writes:

g = \frac{\rho_2}{\rho_1} \approx 1, \qquad h = \frac{c_2}{c_1}
\approx 1.

Equivalent density and compressibility contrasts may be written as:

\gamma\_\rho = \frac{\rho_2 - \rho_1}{\rho_2}, \qquad \gamma\_\kappa =
\frac{\kappa_2 - \kappa_1}{\kappa_1},

where \kappa_j = (\rho_j c_j^2)^{-1} is compressibility in medium j. The
scattered field is then treated as a first-order perturbation produced
by these small departures from the surrounding fluid rather than by a
strongly reflecting or traction-supporting surface.

This is why weak-scattering treatments are usually expressed as
perturbative approximations instead of exact boundary-value problems
with a sharply enforcing interface condition. Physically, the target
behaves approximately like the surrounding fluid, and the backscatter is
produced by residual density and compressibility contrasts rather than
by a large impedance jump at the boundary.

## References

Anderson, Victor C. 1950. “Sound Scattering from a Fluid Sphere.” *The
Journal of the Acoustical Society of America* 22 (4): 426–31.
<https://doi.org/10.1121/1.1906621>.

Faran, James J. 1951. “Sound Scattering by Solid Cylinders and Spheres.”
*The Journal of the Acoustical Society of America* 23 (4): 405–18.
<https://doi.org/10.1121/1.1906780>.

Feuillade, C., and R. W. Nero. 1998. “A Viscous-Elastic Swimbladder
Model for Describing Enhanced-Frequency Resonance Scattering from Fish.”
*The Journal of the Acoustical Society of America* 103 (6): 3245–55.
<https://doi.org/10.1121/1.423076>.

Goodman, Ralph R., and Raya Stern. 1962. “Reflection and Transmission of
Sound by Elastic Spherical Shells.” *The Journal of the Acoustical
Society of America* 34 (3): 338–44. <https://doi.org/10.1121/1.1928120>.
