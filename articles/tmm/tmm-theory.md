# Transition matrix method (TMM) theory

## Introduction

Benchmarked Partially validated Experimental

[Overview](https://brandynlucca.github.io/acousticTS/articles/tmm/index.md)
[Implementation](https://brandynlucca.github.io/acousticTS/articles/tmm/tmm-implementation.md)

The transition matrix method (`TMM`) represents a target by the linear
map from incident-wave coefficients to scattered-wave coefficients
([Waterman 1969](#ref-Waterman_1969)). Once constructed, that map can be
evaluated for different incident and receive directions without
resolving the boundary-value problem.

This coefficient map supports monostatic and bistatic scattering,
rotations, and orientation averages when the retained basis and
numerical solution are valid for those operations.

The useful basis depends on geometry. Spheres and spherical shells use
spherical waves. Prolate spheroids use prolate spheroidal waves. Oblate
spheroids can use spherical waves enforced on r(\theta). Finite
cylinders favor cylindrical modes because the sidewall-endcap junction
converges poorly in a smooth spherical basis ([Varadan et al.
1982](#ref-Varadan_1982); [Waterman 2009](#ref-Waterman_2009)).

Medium indices, field conventions, and reporting quantities follow
[Notation and
symbols](https://brandynlucca.github.io/acousticTS/articles/notation-and-symbols/notation-and-symbols.md).
Boundary types are defined on the [Boundary
conditions](https://brandynlucca.github.io/acousticTS/articles/boundary_conditions.md)
page.

## General single-target T-matrix formulation

### Incident and scattered modal expansions

For time-harmonic pressure with implicit factor e^{-i\omega t}, the
pressure in a homogeneous region satisfies the Helmholtz equation
([Morse and Ingard 1968](#ref-Morse_1968)):

\nabla^2 p + k^2 p = 0.

In a modal T-matrix formulation, the incident and scattered fields are
expanded as:

p^{inc} = \sum\_{\nu} a\_{\nu}\\\psi^{(1)}\_{\nu}, \qquad p^{sca} =
\sum\_{\nu} f\_{\nu}\\\psi^{(3)}\_{\nu}.

Here \psi^{(1)}\_{\nu} is regular and \psi^{(3)}\_{\nu} is outgoing. The
transition matrix is the linear map:

\mathbf{f} = \mathbf{T}\mathbf{a}.

For an axisymmetric target, azimuthal order m decouples and \mathbf T is
block diagonal in m. If a\_{m'n'}(\widehat{\mathbf k}\_i) are the
plane-wave coefficients for incident direction \widehat{\mathbf k}\_i,
the far-field amplitude in receive direction \widehat{\mathbf k}\_s has
the general form:

f(\widehat{\mathbf k}\_s,\widehat{\mathbf k}\_i) = \frac{1}{k_1}
\sum\_{mn}\sum\_{m'n'} \mathcal Y\_{mn}(\widehat{\mathbf k}\_s)
T\_{mn,m'n'} a\_{m'n'}(\widehat{\mathbf k}\_i),

where \mathcal Y\_{mn} includes the angular basis and outgoing-wave
phase. Backscatter is the special case \widehat{\mathbf
k}\_s=-\widehat{\mathbf k}\_i ([Mishchenko et al.
2002](#ref-Mishchenko_2002)).

### Boundary conditions

Consider the four scalar acoustic boundary types: `fixed_rigid`,
`pressure_release`, `liquid_filled`, and `gas_filled`.

The first two use only the exterior basis. For a rigid target, the
normal velocity vanishes at the boundary, which in pressure language
means the normal derivative of pressure vanishes. For a pressure-release
target, the pressure itself vanishes at the surface.

For fluid- or gas-filled targets, an interior field must also be
represented. The boundary conditions are then:

p^{ext} = p^{int}, \qquad \frac{1}{\rho\_{ext}} \frac{\partial
p^{ext}}{\partial n} = \frac{1}{\rho\_{int}} \frac{\partial
p^{int}}{\partial n}.

Density and sound speed therefore enter the modal boundary operator.

## Spherical-coordinate branch

### General axisymmetric surface formulation

Suppose the scatterer surface is axisymmetric and can be written in
spherical coordinates as:

r = r(\theta).

The exterior regular and outgoing basis states are then built from
spherical partial waves:

\begin{align\*} \psi^{(1)}\_{mn}(r,\theta,\phi) &=
j_n(kr)\\P_n^m(\cos\theta)\\e^{im\phi}, \\
\psi^{(3)}\_{mn}(r,\theta,\phi) &=
h_n^{(1)}(kr)\\P_n^m(\cos\theta)\\e^{im\phi}. \end{align\*}

Along the curved meridional profile r(\theta), the outward normal
derivative is:

\frac{\partial}{\partial n} = \frac{1}{\sqrt{1 + \left\[r\_\theta /
r\right\]^2}} \left( \frac{\partial}{\partial r} -
\frac{r\_\theta}{r^2}\frac{\partial}{\partial \theta} \right).

where r\_\theta=dr/d\theta.

A nonspherical boundary therefore couples radial and angular
derivatives. The boundary residual is projected back onto the retained
spherical basis to form each m-block.

### Special geometries in the spherical branch

#### Sphere

For a sphere:

r(\theta) = a.

so r\_\theta = 0 and the normal derivative reduces to the ordinary
radial derivative. The spherical-coordinate T-matrix formulation
therefore collapses to the classical spherical partial-wave problem.

#### Spherical shells

For a concentric shell with outer radius a and inner radius b, spherical
symmetry keeps the T-matrix diagonal in (m,n). The diagonal coefficient
is the exterior modal scattering coefficient obtained from the two
interface systems. A fluid shell couples exterior, shell, and core
acoustic fields. An elastic shell also carries longitudinal and
transverse elastic potentials. The coefficient map is therefore:

f\_{mn}=T_n a\_{mn},

with T_n supplied by the corresponding shell boundary determinant. See
[VESM
theory](https://brandynlucca.github.io/acousticTS/articles/vesm/vesm-theory.md)
for fluid shells and [ESSMS
theory](https://brandynlucca.github.io/acousticTS/articles/essms/essms-theory.md)
for elastic shells.

#### Oblate spheroid

For an oblate spheroid with axial semiaxis c and equatorial semiaxis a,
where c \le a:

r(\theta) = \left( \frac{\cos^2\theta}{c^2} + \frac{\sin^2\theta}{a^2}
\right)^{-1/2}.

Differentiating gives:

r\_\theta = -\sin\theta\cos\theta \left( \frac{1}{a^2} - \frac{1}{c^2}
\right) \left( \frac{\cos^2\theta}{c^2} + \frac{\sin^2\theta}{a^2}
\right)^{-3/2}.

So the oblate branch remains a spherical-basis T-matrix formulation, but
with the actual oblate meridional geometry entering through r(\theta)
and r\_\theta.

#### Finite cylinder

For a right circular finite cylinder with half-length a and radius b,
the meridional surface can be written piecewise in spherical coordinates
as:

r(\theta) = \min\left( \frac{a}{\|\cos\theta\|},
\frac{b}{\|\sin\theta\|} \right).

A ray from the origin first meets either an end cap or the sidewall. The
profile is continuous but not differentiable at their junction, which
slows convergence of a spherical-wave representation ([Waterman
2009](#ref-Waterman_2009)).

### Monostatic reconstruction

Once the retained block coefficients are obtained, the backscatter
amplitude is reconstructed by evaluating the outgoing expansion in the
receive direction opposite to the incident plane wave. The
backscattering cross section then follows from:

\sigma\_{bs} = \|f\_{bs}\|^2

and target strength is ([MacLennan et al. 2002](#ref-MacLennan_2002)):

TS = 10 \log\_{10}\left( \frac{\sigma\_{bs}}{1\\ \mathrm{m}^2}\right).

## Cylindrical interpretation

For a finite cylinder, cylindrical partial waves describe the
cross-section and an axial operator describes finite length. Near
broadside in monostatic scattering, that operator reduces to the
familiar finite-cylinder coherence factor.

The sidewall-endcap corners are the main obstacle. Angular products
require particular care because a near-broadside coherence closure is
not a general bistatic cylinder solution.

## Prolate spheroid branch

### Why spherical coordinates are not the best exact basis

A prolate spheroid is not a constant-r surface. So while spherical-wave
expansions can still be written down, they do not align naturally with
the geometry. This is exactly the regime where the classic
spheroidal-coordinate transition-matrix literature becomes relevant
([Varadan et al. 1982](#ref-Varadan_1982); [Hackman and Todoroff
1984](#ref-Hackman_1984)).

In prolate spheroidal coordinates, the boundary is the coordinate
surface \xi=\xi_1.

### Prolate spheroidal coordinates

Let q be the semifocal length and let (\xi,\eta,\phi) be prolate
spheroidal coordinates. In Cartesian coordinates:

\begin{align\*} x &= q\sqrt{(\xi^2 - 1)(1 - \eta^2)}\cos\phi, \\ y &=
q\sqrt{(\xi^2 - 1)(1 - \eta^2)}\sin\phi, \\ z &= q \xi \eta.
\end{align\*}

The coordinate ranges are:

\xi \ge 1, \qquad -1 \le \eta \le 1, \qquad 0 \le \phi \< 2\pi.

If the body surface is \xi = \xi_1, then the major semi-axis a and minor
semi-axis b satisfy:

a = \xi_1 q, \qquad b = q\sqrt{\xi_1^2 - 1}.

so equivalently:

\xi_1 = \left\[1 - \left(\frac{b}{a}\right)^2\right\]^{-1/2}.

### Spheroidal modal representation

In a homogeneous region, the separated pressure field is written as the
product of a radial spheroidal function in \xi, an angular spheroidal
function in \eta, and an azimuthal factor in \phi.

The resulting field expansion has the same logical structure as the
generic T-matrix expression above, but the basis is geometry-matched.

For rigid and pressure-release prolates, the retained degrees remain
effectively local in the exact spheroidal basis. For liquid- and
gas-filled prolates, the interior and exterior reduced frequencies
differ, so the angular bases no longer match exactly. This introduces
overlap-driven coupling between retained degrees, exactly as in the
exact prolate spheroidal modal-series solution ([Hackman and Todoroff
1984](#ref-Hackman_1984); [Spence and Granger 1951](#ref-Spence_1951);
[Furusawa 1988](#ref-Furusawa_1988)).

### T-matrix interpretation in a geometry-matched basis

The T-matrix concept is independent of a particular coordinate system.
It requires regular and outgoing bases, a boundary operator, and a
converged coefficient map between them.

## Mathematical assumptions and scope

The formulations above assume linear time-harmonic acoustics, a single
target, an axisymmetric geometry, and homogeneous properties within each
region. Accuracy depends on modal truncation, boundary quadrature, and
using a basis suited to the surface. Smooth spherical and spheroidal
surfaces support general angular reconstruction once their blocks
converge. A finite-cylinder near-broadside closure does not by itself
establish a general bistatic cylinder operator.

## References

Furusawa, Masahiko. 1988. “Prolate Spheroidal Models for Predicting
General Trends of Fish Target Strength.” *Journal of the Acoustical
Society of Japan (E)* 9 (1): 13–24. <https://doi.org/10.1250/ast.9.13>.

Hackman, Roger H., and Douglas G. Todoroff. 1984. “An Application of the
Spheroidal-Coordinate-Based Transition Matrix: Acoustic Scattering from
High Aspect Ratio Solids.” *The Journal of the Acoustical Society of
America* 76 (S1): S8–8. <https://doi.org/10.1121/1.2022083>.

MacLennan, David N., Percy G. Fernandes, and John Dalen. 2002. “A
Consistent Approach to Definitions and Symbols in Fisheries Acoustics.”
*ICES Journal of Marine Science* 59 (2): 365–69.
<https://doi.org/10.1006/jmsc.2001.1158>.

Mishchenko, Michael I., Larry D. Travis, and Andrew A. Lacis. 2002.
*Scattering, Absorption, and Emission of Light by Small Particles*.
Cambridge University Press.

Morse, Philip M., and K. Uno Ingard. 1968. *Theoretical Acoustics*.
McGraw-Hill.

Spence, R. D., and Sara Granger. 1951. “The Scattering of Sound from a
Prolate Spheroid.” *The Journal of the Acoustical Society of America* 23
(6): 701–6. <https://doi.org/10.1121/1.1906827>.

Varadan, V. K., V. V. Varadan, Louis R. Dragonette, and Lawrence Flax.
1982. “Computation of Rigid Body Scattering by Prolate Spheroids Using
the *t* -Matrix Approach.” *The Journal of the Acoustical Society of
America* 71 (1): 22–25. <https://doi.org/10.1121/1.387311>.

Waterman, P. C. 1969. “New Formulation of Acoustic Scattering.” *The
Journal of the Acoustical Society of America* 45 (6): 1417–29.
<https://doi.org/10.1121/1.1911619>.

Waterman, P. C. 2009. “T -Matrix Methods in Acoustic Scattering.” *The
Journal of the Acoustical Society of America* 125 (1): 42–51.
<https://doi.org/10.1121/1.3035839>.
