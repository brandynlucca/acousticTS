# Scattering boundary conditions

## Introduction

Boundary conditions turn the governing wave equations into a particular
scattering problem. They state whether an interface can sustain
pressure, whether it can move, and which field quantities pass
continuously between adjacent media. A fixed-rigid surface, a
pressure-release surface, and a penetrable fluid target can therefore
scatter very differently even when they have the same geometry.

The conditions must act on the **total** field at an exterior boundary.
The incident and scattered fields are separated to solve the exterior
problem, but neither field generally satisfies the boundary condition
alone ([Colton and Kress 2013](#ref-Colton_2013)).

See [Notation and
Symbols](https://brandynlucca.github.io/acousticTS/articles/notation-and-symbols/notation-and-symbols.md)
and the [Acoustic Scattering
Primer](https://brandynlucca.github.io/acousticTS/articles/acoustic-scattering-primer/acoustic-scattering-primer.md).

## Shared interface notation

Medium `1` is the exterior fluid. Target regions are numbered inward.
For an unlayered target, \Gamma denotes the interface between media `1`
and `2`. For a shell with outer radius a and inner radius b, the outer
and inner interfaces are \Gamma_a and \Gamma_b.

The exterior pressure is decomposed as:

p_1^{\mathrm{tot}} =p_1^{\mathrm{inc}}+p_1^{\mathrm{scat}}.

At a fluid-fluid interface, \mathbf{n} is one unit normal used on both
sides of the interface. Its orientation does not change the continuity
conditions as long as it is used consistently. The normal derivative is:

\partial_n p =\frac{\partial p}{\partial n}
=\mathbf{n}\mathbin{\cdot}\nabla p.

With the e^{-i\omega t} convention, the linearized momentum equation and
its phasor form are ([Morse and Ingard 1968](#ref-Morse_1968)):

\rho_j\frac{\partial\mathbf{v}\_j}{\partial t} =-\nabla p_j, \qquad
-i\omega\rho_j\mathbf{v}\_j=-\nabla p_j.

The fluid particle velocity is therefore:

\mathbf{v}\_j =\frac{1}{i\omega\rho_j}\nabla p_j, \qquad v\_{n,j}
=\frac{1}{i\omega\rho_j}\partial_n p_j.

Each homogeneous fluid pressure satisfies the Helmholtz equation:

\nabla^2p_j+k_j^2p_j=0, \qquad k_j=\frac{\omega}{c_j}.

The directional density and sound-speed ratios used for adjacent media
are:

g\_{ij}=\frac{\rho_i}{\rho_j}, \qquad h\_{ij}=\frac{c_i}{c_j}.

Thus, g\_{21} and h\_{21} compare the first target region with the
exterior. For a shelled target, g\_{32} and h\_{32} compare the core
with the shell.

## Simple boundaries

### Fixed rigid

A fixed-rigid boundary cannot move in its normal direction. The normal
fluid velocity at the surface must therefore vanish:

v\_{n,1}=\mathbf{v}\_1\mathbin{\cdot}\mathbf{n}=0 \qquad\text{on
}\Gamma.

The same limit can be expressed through the specific surface impedance:

Z_s=\frac{p_1^{\mathrm{tot}}}{v\_{n,1}}, \qquad
\|Z_s\|\longrightarrow\infty.

Substitution of the phasor momentum relation gives the Neumann boundary
condition:

\partial_n p_1^{\mathrm{tot}}=0 \qquad\text{on }\Gamma.

In terms of the incident and scattered fields, the prescribed data for
the scattered field are:

\partial_n p_1^{\mathrm{scat}} =-\partial_n p_1^{\mathrm{inc}}
\qquad\text{on }\Gamma.

The pressure itself need not vanish. Instead, the reflected field
cancels the incident normal velocity at the surface. In partial-wave
formulations, the Neumann condition fixes each scattering coefficient by
excluding any total field that would produce radial surface motion. At
high frequency, a smooth rigid surface supports a locally specular
reflection without the pressure phase reversal of a pressure-release
surface.

### Pressure release

A pressure-release surface cannot sustain an acoustic pressure
fluctuation. It is the zero-impedance limit:

Z_s=\frac{p_1^{\mathrm{tot}}}{v\_{n,1}}, \qquad \|Z_s\|\longrightarrow0.

For finite normal velocity, the total pressure obeys the Dirichlet
condition:

p_1^{\mathrm{tot}}=0 \qquad\text{on }\Gamma.

The corresponding condition on the scattered field is:

p_1^{\mathrm{scat}}=-p_1^{\mathrm{inc}} \qquad\text{on }\Gamma.

Normal motion is not constrained to vanish. The idealization is
appropriate when the material on the other side has negligible acoustic
impedance relative to the exterior fluid, as for an idealized gas
boundary away from effects that require an explicit interior solution. A
locally reflected pressure wave undergoes a phase reversal. In a modal
solution, zero surface pressure replaces the zero-normal-velocity
condition used for a rigid target.

### Fluid-filled

A penetrable fluid target supports acoustic pressure on both sides of
its boundary. The exterior total field and interior field satisfy
([Anderson 1950](#ref-Anderson_1950)):

\nabla^2p_1^{\mathrm{tot}}+k_1^2p_1^{\mathrm{tot}}=0, \qquad
\nabla^2p_2+k_2^2p_2=0.

An inviscid fluid has Cauchy stress
\boldsymbol{\sigma}\_j=-p_j\mathbf{I}. Continuity of normal traction
therefore requires pressure continuity:

p_1^{\mathrm{tot}}=p_2 \qquad\text{on }\Gamma.

The interface also admits neither separation nor interpenetration, so
normal particle velocity is continuous:

\mathbf{v}\_1\mathbin{\cdot}\mathbf{n}
=\mathbf{v}\_2\mathbin{\cdot}\mathbf{n} \qquad\text{on }\Gamma.

Using the momentum equation gives the derivative form of the second
transmission condition:

\frac{1}{\rho_1}\partial_n p_1^{\mathrm{tot}}
=\frac{1}{\rho_2}\partial_n p_2 \qquad\text{on }\Gamma.

Pressure and normal motion are both generally nonzero. Their continuity
couples the exterior and interior solutions, while density and
compressibility contrasts determine the strength and phase of the
reflected and transmitted fields ([Medwin and Clay
1998](#ref-Medwin_1998)). Tangential velocity need not be continuous
across an ideal inviscid fluid-fluid interface because neither fluid
transmits shear traction.

## Acoustic shelled boundaries

An acoustic shell is a finite fluid layer, not an immobile shell and not
an elastic solid. Medium `2` occupies the region between \Gamma_a and
\Gamma_b. It supports compressional pressure waves but no shear stress.
The shell field satisfies:

\nabla^2p_2+k_2^2p_2=0 \qquad\text{between }\Gamma_a\text{ and
}\Gamma_b.

At the outer interface, pressure and normal velocity are continuous:

p_1^{\mathrm{tot}}=p_2, \qquad \frac{1}{\rho_1}\partial_n
p_1^{\mathrm{tot}} =\frac{1}{\rho_2}\partial_n p_2 \qquad\text{on
}\Gamma_a.

For a spherical shell, both independent radial solutions are admissible
because the layer does not include the origin. An axisymmetric shell
field can be expanded as:

p_2(r,\theta) =\sum\_{m=0}^{\infty} \left\[
B_m\\j_m(k_2r)+C_m\\y_m(k_2r) \right\]P_m(\cos\theta), \qquad b\<r\<a.

Here, j_m and y_m are spherical Bessel functions of the first and second
kinds. The singularity of y_m(k_2r) at the origin is irrelevant because
r=0 is outside the shell layer. The condition at \Gamma_b determines the
relative weights of these two branches.

### Pressure release interior

For a pressure-release core, the shell pressure vanishes at the inner
interface:

p_2=0 \qquad\text{on }\Gamma_b.

In a spherical modal solution, a radial combination that satisfies this
inner condition identically is:

R_m^{\mathrm{pr}}(r) =j_m(k_2r)y_m(k_2b)-y_m(k_2r)j_m(k_2b).

Its value at the inner radius is:

R_m^{\mathrm{pr}}(b)=0.

The outer pressure and velocity conditions then couple this admissible
shell field to p_1^{\mathrm{tot}}. The shell can carry normal acoustic
motion even though the pressure is zero at its inner boundary. Its
finite thickness and wavenumber add propagation phase between the two
interfaces. The resulting response can contain cavity-like resonances
that are absent from a single pressure-release surface.

### Fluid-filled interior

If medium `3` is a fluid or gas core, it supports its own Helmholtz
field:

\nabla^2p_3+k_3^2p_3=0 \qquad\text{inside }\Gamma_b.

Pressure and normal velocity are continuous at the inner interface:

p_2=p_3, \qquad \frac{1}{\rho_2}\partial_n p_2
=\frac{1}{\rho_3}\partial_n p_3 \qquad\text{on }\Gamma_b.

For a spherical core containing the origin, regularity excludes the
spherical Bessel function of the second kind. The core expansion is
therefore:

p_3(r,\theta) =\sum\_{m=0}^{\infty}D_m\\j_m(k_3r)P_m(\cos\theta), \qquad
0\leq r\<b.

Together, the two conditions at \Gamma_a and the two at \Gamma_b
determine the exterior scattering field, the two shell branches, and the
regular core field for each mode. The shell is not assigned infinite
impedance, and its normal velocity is not set to zero. Resonance and
interference arise from compressional propagation through the shell and
core rather than from elastic shear or bending waves.

## Elastic shelled boundaries

### Elastic shell boundary

An elastic shell supports longitudinal and transverse motion. Its
displacement \mathbf{u}(\mathbf{x},t) follows from balance of linear
momentum and Hooke’s law. For a homogeneous isotropic shell, the
time-domain Navier equation is ([Achenbach 1973](#ref-Achenbach_1973)):

(\lambda+2\mu)\nabla(\nabla\mathbin{\cdot}\mathbf{u})
-\mu\nabla\mathbin{\times}(\nabla\mathbin{\times}\mathbf{u})
-\rho_s\frac{\partial^2\mathbf{u}}{\partial t^2} =0.

For a displacement phasor with e^{-i\omega t} dependence, this becomes:

(\lambda+2\mu)\nabla(\nabla\mathbin{\cdot}\mathbf{u})
-\mu\nabla\mathbin{\times}(\nabla\mathbin{\times}\mathbf{u})
+\rho_s\omega^2\mathbf{u} =0.

The Helmholtz decomposition separates dilatational and equivoluminal
motion:

\mathbf{u}=\nabla\Phi+\nabla\mathbin{\times}\mathbf{\Psi}, \qquad
\nabla\mathbin{\cdot}\mathbf{\Psi}=0.

The scalar and vector potentials satisfy separate Helmholtz equations:

\nabla^2\Phi+k_L^2\Phi=0, \qquad
\nabla^2\mathbf{\Psi}+k_T^2\mathbf{\Psi}=0.

Their wavenumbers and wave speeds are:

k_L=\frac{\omega}{c_L}, \qquad k_T=\frac{\omega}{c_T}, \qquad
c_L=\sqrt{\frac{\lambda+2\mu}{\rho_s}}, \qquad
c_T=\sqrt{\frac{\mu}{\rho_s}}.

The infinitesimal strain and Cauchy stress tensors are:

\boldsymbol{\varepsilon} =\frac{1}{2}\left\[
\nabla\mathbf{u}+(\nabla\mathbf{u})^{\mathsf T} \right\], \qquad
\boldsymbol{\sigma} =\lambda(\nabla\mathbin{\cdot}\mathbf{u})\mathbf{I}
+2\mu\boldsymbol{\varepsilon}.

At either fluid-solid interface, let \mathbf{n}\_s point from the solid
shell into the adjacent fluid. The fluid stress is
\boldsymbol{\sigma}\_f=-p_f\mathbf{I}. Continuity of traction is then
the vector condition:

\boldsymbol{\sigma}\mathbf{n}\_s =-p_f\mathbf{n}\_s.

Separating that relation into normal and tangential parts gives:

\mathbf{n}\_s\mathbin{\cdot} \boldsymbol{\sigma}\mathbf{n}\_s=-p_f,
\qquad (\mathbf{I}-\mathbf{n}\_s\mathbf{n}\_s^{\mathsf T})
\boldsymbol{\sigma}\mathbf{n}\_s=\mathbf{0}.

The second condition states that an inviscid fluid applies no tangential
traction. It does not require the tangential solid displacement to
vanish. Continuity of normal velocity supplies the remaining kinematic
condition:

\mathbf{v}\_f\mathbin{\cdot}\mathbf{n}\_s
=-i\omega\mathbf{u}\mathbin{\cdot}\mathbf{n}\_s.

Written in terms of fluid pressure, the same condition is:

\frac{1}{i\omega\rho_f}\partial\_{n_s}p_f =-i\omega u_n, \qquad
u_n=\mathbf{u}\mathbin{\cdot}\mathbf{n}\_s.

At the outer shell surface, p_f=p_1^{\mathrm{tot}}. At the inner
surface, p_f=p_3. The normal \mathbf{n}\_s points radially outward at
r=a and radially inward at r=b. Stating the normal this way removes any
ambiguity in the pressure-traction sign.

For a solid elastic target with no core, these conditions are imposed
only at the exterior surface. Regularity at the origin excludes singular
longitudinal and transverse radial solutions. The resulting fluid-solid
coupling is the starting point for classical elastic-sphere and
elastic-cylinder scattering theory ([Faran 1951](#ref-Faran_1951)).

For an axisymmetric shell, normal traction, normal velocity, and one
tangential traction condition are imposed at each radius. These six
scalar conditions couple the exterior acoustic field, both elastic
potentials in the shell, and the interior acoustic field ([Goodman and
Stern 1962](#ref-Goodman_1962)). They admit resonances governed by shell
thickness, density, longitudinal and transverse wave speeds, and fluid
loading. This is the essential distinction from an acoustic shell, which
has no transverse potential and cannot support shear or bending motion
([Stanton 1990](#ref-Stanton_1990)).

## Viscous interfaces

A compressible Newtonian viscous fluid supports deviatoric stress as
well as pressure. If \eta is the shear viscosity, \zeta is the bulk
viscosity, and \mathbf{D} is the rate-of-deformation tensor, its Cauchy
stress is ([Pierce 1989](#ref-Pierce_1989)):

\boldsymbol{\sigma}\_v =-p_v\mathbf{I} +2\eta\mathbf{D}
+\left(\zeta-\frac{2}{3}\eta\right)
(\nabla\mathbin{\cdot}\mathbf{v}\_v)\mathbf{I}.

The rate-of-deformation tensor is:

\mathbf{D} =\frac{1}{2}\left\[
\nabla\mathbf{v}\_v+(\nabla\mathbf{v}\_v)^{\mathsf T} \right\].

At an interface with an inviscid acoustic fluid, let \mathbf{n} point
from the viscous medium into the inviscid fluid. Normal traction and
normal velocity are continuous:

\mathbf{n}\mathbin{\cdot}\boldsymbol{\sigma}\_v\mathbf{n} =-p_f, \qquad
\mathbf{v}\_v\mathbin{\cdot}\mathbf{n}
=\mathbf{v}\_f\mathbin{\cdot}\mathbf{n}.

Because the inviscid fluid cannot apply shear traction, the tangential
traction on the viscous side vanishes:

(\mathbf{I}-\mathbf{n}\mathbf{n}^{\mathsf T})
\boldsymbol{\sigma}\_v\mathbf{n}=\mathbf{0}.

At a bonded viscous-fluid and elastic-solid interface, both velocity and
traction are continuous. With \mathbf{n} used consistently on both
sides, the phasor conditions are:

\mathbf{v}\_v=-i\omega\mathbf{u}\_s, \qquad
\boldsymbol{\sigma}\_v\mathbf{n} =\boldsymbol{\sigma}\_s\mathbf{n}.

These vector conditions couple normal and tangential motion. Viscous
compressional and shear branches are generally attenuating, so the
interface can broaden and damp resonances that would be sharper in a
lossless elastic shell ([Feuillade and Nero 1998](#ref-Feuillade_1998)).

## Circumferential waves and resonance families

Boundary conditions do not introduce circumferential waves as separate
terms. Instead, they determine which wave families the complete field
solution can support and how those waves couple to the exterior fluid.
Geometry, frequency, material properties, and attenuation then control
their phase, propagation distance, and resonance structure ([Überall
1973](#ref-Uberall_1973)).

| Phenomenon | Boundary setting | Physical interpretation |
|----|----|----|
| Franz or creeping waves | Smooth convex rigid or pressure-release boundaries | Exterior diffracted waves that follow the surface into the geometrical shadow while continually radiating into the surrounding fluid. |
| Rayleigh-type waves | Fluid-loaded elastic solids and sufficiently thick elastic shells | Surface-localized elastic motion formed from coupled longitudinal and transverse fields. |
| Lamb-type waves | Thin or moderately thin elastic shells | Dispersive guided shell waves whose motion depends on the conditions at both shell surfaces. |
| Whispering-gallery waves | Penetrable fluid or elastic targets | High-order internal or surface-guided waves retained near the boundary by repeated refraction or reflection. |
| Glory | Spherical or nearly axisymmetric targets | A directional enhancement produced when circumferential contributions converge coherently near a symmetry direction. It is an interference or focusing effect, not a separate wave family. |

These labels describe identifiable contributions within an exact modal
solution. They are not additional boundary conditions. For example, the
longitudinal and transverse potentials of an elastic target already
contain its surface and guided-wave contributions once the traction and
velocity conditions have been enforced. Separating those contributions
generally requires modal phase, resonance-pole, time-domain, or
asymptotic analysis. Surface-elastic-wave interpretations of
elastic-sphere and elastic-cylinder backscatter provide a classical
example ([Marston 1988](#ref-Marston_1988)).

## Weak, fluid-like

Weak fluid-like scattering is an approximation regime, not a replacement
for the fluid-fluid transmission conditions. For spatially varying
density and compressibility, the source-free pressure equation can be
written as:

\nabla\mathbin{\cdot} \left(\frac{1}{\rho}\nabla p\right)
+\omega^2\kappa p=0, \qquad \kappa=\frac{1}{\rho c^2}.

When medium `2` differs only slightly from medium `1`, its total
internal field may be replaced at first order by the incident field:

p_2(\mathbf{x})\approx p_1^{\mathrm{inc}}(\mathbf{x}).

The density and sound-speed ratios must remain close to unity:

g\_{21}=\frac{\rho_2}{\rho_1}\approx1, \qquad
h\_{21}=\frac{c_2}{c_1}\approx1.

The corresponding density and compressibility perturbations are:

\gamma\_\rho =\frac{\rho_2-\rho_1}{\rho_2} =1-\frac{1}{g\_{21}}, \qquad
\gamma\_\kappa =\frac{\kappa_2-\kappa_1}{\kappa_1}
=\frac{1}{g\_{21}h\_{21}^2}-1.

The weak-contrast assumption requires:

\|\gamma\_\rho\|\ll1, \qquad \|\gamma\_\kappa\|\ll1.

The scattered field is then the first-order response to the residual
density and compressibility contrasts ([Chu and Ye
1999](#ref-Chu_1999)). Strong impedance jumps, internal multiple
scattering, and resonant behavior violate this approximation even though
the exact fluid-fluid boundary conditions remain valid.

## References

Achenbach, J. D. 1973. *Wave Propagation in Elastic Solids*.
North-Holland Series in Applied Mathematics and Mechanics, v. 16.
North-Holland Pub. Co. American Elsevier Pub. Co.

Anderson, Victor C. 1950. “Sound Scattering from a Fluid Sphere.” *The
Journal of the Acoustical Society of America* 22 (4): 426–31.
<https://doi.org/10.1121/1.1906621>.

Chu, Dezhang, and Zhen Ye. 1999. “A Phase-Compensated Distorted Wave
Born Approximation Representation of the Bistatic Scattering by Weakly
Scattering Objects: Application to Zooplankton.” *The Journal of the
Acoustical Society of America* 106 (4): 1732–43.
<https://doi.org/10.1121/1.428036>.

Colton, David, and Rainer Kress. 2013. *Inverse Acoustic and
Electromagnetic Scattering Theory*. 3rd ed. Vol. 93. Springer.
<https://doi.org/10.1007/978-1-4614-4942-3>.

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

Marston, Philip L. 1988. “GTD for Backscattering from Elastic Spheres
and Cylinders in Water and the Coupling of Surface Elastic Waves with
the Acoustic Field.” *The Journal of the Acoustical Society of America*
83 (1): 25–37. <https://doi.org/10.1121/1.396428>.

Medwin, Herman, and Clarence S. Clay. 1998. *Fundamentals of Acoustical
Oceanography*. Applications of Modern Acoustics. Academic Press.

Morse, Philip M., and K. Uno Ingard. 1968. *Theoretical Acoustics*.
McGraw-Hill.

Pierce, Allan D. 1989. *Acoustics: An Introduction to Its Physical
Principles and Applications*. Acoustical Society of America.

Stanton, T. K. 1990. “Sound Scattering by Spherical and Elongated
Shelled Bodies.” *The Journal of the Acoustical Society of America* 88
(3): 1619–33. <https://doi.org/10.1121/1.400321>.

Überall, Herbert. 1973. “Surface Waves in Acoustics.” In *Physical
Acoustics*, edited by Warren P. Mason and R. N. Thurston, vol. 10.
Academic Press.
