# Phase-compensated distorted wave Born approximation

## Introduction

Validated Experimental

[Overview](https://brandynlucca.github.io/acousticTS/articles/pcdwba/index.md)
[Implementation](https://brandynlucca.github.io/acousticTS/articles/pcdwba/pcdwba-implementation.md)

The phase-compensated distorted wave Born approximation (`PCDWBA`)
extends the ordinary `DWBA` to elongated bodies whose centerlines are
curved rather than straight ([Chu and Ye
1999](#ref-chu_phase-compensated_1999); [Stanton
1989](#ref-stanton_sound_1989)). It preserves the same weak-scattering
local kernel used for slender fluid-like bodies, but it replaces the
straight-body phase bookkeeping by a centerline-dependent phase that
follows the bent geometry explicitly.

The family therefore sits conceptually between two extremes: it is more
physical than treating a curved body as if it were straight, but it
remains a Born-type perturbation model rather than a full curved-body
boundary-value solution.

The notation below follows the shared package conventions: medium `1` is
the surrounding seawater and medium `2` is the weakly scattering body.

## Weak-scattering starting point

### Born contrast term

`PCDWBA` inherits its material-contrast physics from the ordinary
distorted wave Born approximation. For a fluid-like body in seawater,
the density and sound-speed contrasts are:

g\_{21} = \frac{\rho_2}{\rho_1}, \qquad h\_{21} = \frac{c_2}{c_1},

where \rho_1, c_1 are the surrounding-fluid density and sound speed and
\rho_2, c_2 are the corresponding body properties.

The standard fluid-like contrast factor is then:

C\_{21} = \frac{1 - g\_{21} h\_{21}^2}{g\_{21} h\_{21}^2} -
\frac{g\_{21} - 1}{g\_{21}}.

The approximation is most defensible when g\_{21} and h\_{21} remain
close to unity, so that the total field inside the body remains close to
the distorted incident field.

### Volume-integral viewpoint

In the general Born picture, the scattered field is written as a volume
integral over the target:

p^{\mathrm{sca}}(\mathbf{r}) \propto \int_V C\_{21}(\mathbf{r}') \\
G_1(\mathbf{r},\mathbf{r}') \\ p^{\mathrm{ref}}(\mathbf{r}') \\ dV',

where G_1 is the Green’s function of the exterior medium and
p^{\mathrm{ref}} is the chosen reference field. For slender bodies, this
three-dimensional integral is reduced by separating the local
cross-sectional response from the phase accumulated along the body axis.

## Curved centerline geometry

### Centerline parameterization

Let s denote arc length along the body centerline and let
\mathbf{r}\_c(s) denote the corresponding centerline position. The local
body radius is a(s), and the local tangent angle is \beta(s).

For a uniformly bent canonical body, the centerline is wrapped around an
osculating circle of radius \rho_c, with curvature \kappa = 1/\rho_c. A
taper function may also scale the local radius toward the ends, but that
taper only changes the local cross-sectional factor; it does not alter
the phase logic that motivates `PCDWBA`.

### Why curvature matters

In a straight-body `DWBA`, the two-way phase factor depends only on the
axial position along a straight line. Once the body bends, two points
with the same arc-length separation no longer have the same projection
onto the incident or receive directions. That means a curved target
cannot be described correctly by the same straight-axis phase term
unless the curvature is negligible.

## Local cylindrical reduction

At each centerline location, the body is approximated locally by a
circular cylinder with radius a(s). The azimuthal integration over that
local cross-section produces the same Bessel-type factor that appears in
slender-body Born models:

\frac{J_1\\\left(2 k_2 a(s)\\ \chi(s)\right)} {2 k_2 a(s)\\ \chi(s)},

Here J_1 is the cylindrical Bessel function of the first kind, k_2 =
\omega / c_2 = k_1 / h\_{21} is the body wavenumber, and \chi(s) is the
local projection factor set by the incident and receive directions
relative to the local tangent.

For monostatic backscatter, this projection reduces to a local cosine
term involving the body orientation and tangent angle.

## Phase-compensated line integral

### General curved-body form

After the local cross-sectional reduction, the scattered field becomes a
centerline integral. In monostatic form, the backscattering amplitude
can be written schematically as:

f\_{\mathrm{bs}} \propto \int\_{-L/2}^{L/2} C\_{21}(s) \left(\frac{k_2
a(s)}{1}\right)^2 \frac{J_1\\\left(2 k_2 a(s)\\ \chi(s)\right)} {2 k_2
a(s)\\ \chi(s)} \exp\\\left\[ 2 i k_1 \hat{\mathbf{q}}\cdot
\mathbf{r}\_c(s) \right\] ds,

where \hat{\mathbf{q}} denotes the backscatter direction.

The crucial ingredient is the phase factor:

\exp\\\left\[ 2 i k_1 \hat{\mathbf{q}}\cdot \mathbf{r}\_c(s) \right\].

This is what makes the model phase compensated: the phase is evaluated
on the actual curved centerline rather than on a fictitiously straight
axis.

### Discrete form

In segmented form, the same model appears as:

f\_{\mathrm{bs}} \propto \sum_j C\_{21,j} \left(k_2 a_j\right)^2
\frac{J_1\\\left(2 k_2 a_j \chi_j\right)} {2 k_2 a_j \chi_j}
\exp\\\left\[ 2 i k_1 \hat{\mathbf{q}}\cdot \mathbf{r}\_{c,j} \right\]
\Delta s_j,

where the index j labels body segments and \Delta s_j is the local
centerline spacing.

This is the form most directly useful for irregular profiles, because it
treats curvature, taper, and segment spacing in the same
centerline-aligned framework.

## Relation to straight `DWBA`

If the centerline becomes straight, then \mathbf{r}\_c(s) reduces to a
linear function of s and the phase-compensated expression collapses to
the ordinary straight-body `DWBA`.

That limiting behavior is important physically: `DWBA` is the
straight-axis weak-scattering limit, while `PCDWBA` is the curved-axis
weak-scattering extension.

The two models therefore differ in phase bookkeeping, not in the local
cross-sectional contrast physics.

## Backscatter and target strength

Once the complex backscattering amplitude has been assembled from the
curved centerline sum, the standard monostatic outputs are:

\sigma\_{\mathrm{bs}} = \left\|f\_{\mathrm{bs}}\right\|^2, \qquad
\mathrm{TS} = 10 \log\_{10}\left(\sigma\_{\mathrm{bs}}\right).

No new reporting convention is introduced by the curvature correction.
The change is entirely in the underlying phase-sensitive amplitude.

## Assumptions and regime

`PCDWBA` rests on the following assumptions:

1.  weak fluid-like material contrast,
2.  slender-body reduction to a local cylindrical kernel,
3.  single scattering,
4.  curvature enters through centerline phase rather than through a new
    exact cross-sectional boundary solve,
5.  the target is described meaningfully by a centerline and local
    radius profile.

That is why `PCDWBA` should be read as a curved-body extension of
`DWBA`, not as a new exact geometry family. Its main purpose is to keep
the part of the physics that `DWBA` misses first when the body bends:
the geometry-dependent coherent phase.

## References

Chu, Dezhang, and Zhen Ye. 1999. “A Phase-Compensated Distorted Wave
Born Approximation Representation of the Bistatic Scattering by Weakly
Scattering Objects: Application to Zooplankton.” *The Journal of the
Acoustical Society of America* 106 (4): 1732–43.
<https://doi.org/10.1121/1.428036>.

Stanton, T. K. 1989. “Sound Scattering by Cylinders of Finite Length.
III. Deformed Cylinders.” *The Journal of the Acoustical Society of
America* 86 (2): 691–705. <https://doi.org/10.1121/1.398193>.
