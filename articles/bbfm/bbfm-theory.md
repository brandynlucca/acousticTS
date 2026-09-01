# Body-backbone fish model

## Introduction

Unvalidated Experimental

[Overview](https://brandynlucca.github.io/acousticTS/articles/bbfm/index.md)
[Implementation](https://brandynlucca.github.io/acousticTS/articles/bbfm/bbfm-implementation.md)

The body-backbone fish model (`BBFM`) represents a swimbladder-less fish
as two explicit contributors: a weakly scattering flesh body and an
elastic backbone. This separation is motivated by measurements and
models in which skeletal scattering cannot be absorbed reliably into a
homogeneous flesh contrast ([Gorska et al. 2005](#ref-Gorska_2005);
[Stanton et al. 1998](#ref-Stanton_1998_1)).

BBFM is a coherent component model rather than an exact three-medium
boundary solution:

1.  the flesh body is treated as a weakly scattering fluid-like region,
2.  the backbone is treated as an elastic cylindrical structure, and
3.  the two terms are embedded into one body-fixed frame through a phase
    translation before their complex amplitudes are summed.

The complex amplitudes are combined before forming backscattering
cross-section.

Medium indices and target-strength quantities follow [Notation and
symbols](https://brandynlucca.github.io/acousticTS/articles/notation-and-symbols/notation-and-symbols.md).
The component approximations are developed on the [DWBA
theory](https://brandynlucca.github.io/acousticTS/articles/dwba/dwba-theory.md)
and [ECMS
theory](https://brandynlucca.github.io/acousticTS/articles/ecms/ecms-theory.md)
pages.

## Geometry and medium indexing

The family uses the shared package convention: medium `1` is the
surrounding seawater, medium `2` is the flesh body, and medium `3` is
the backbone.

The exterior seawater wavenumber is:

k_1 = \frac{\omega}{c_1},

where c_1 is the seawater sound speed. The flesh density and sound-speed
contrasts are therefore:

g\_{21} = \frac{\rho_2}{\rho_1}, \qquad h\_{21} = \frac{c_2}{c_1}.

For the backbone, the model keeps the absolute elastic properties
explicit:

\rho_3, \qquad c\_{L,3}, \qquad c\_{T,3},

where c\_{L,3} and c\_{T,3} are the longitudinal and transverse wave
speeds of the elastic backbone.

The backbone term is a seawater-referenced elastic-cylinder surrogate.
It is not the exact solution for an elastic region 3 embedded in flesh
region 2.

## Flesh-body contribution

### Weak-fluid assumption

The flesh component follows the distorted-wave Born approximation.
Material contrasts relative to seawater must be weak enough for
first-order scattering, while the phase is retained across the extended
body ([Chu and Ye 1999](#ref-Chu_1999)).

Using the same contrast notation as the `DWBA` theory page, the
compressibility and density perturbations are:

\gamma\_{\kappa,21} = \frac{\kappa_2 - \kappa_1}{\kappa_1}, \qquad
\gamma\_{\rho,21} = \frac{\rho_2 - \rho_1}{\rho_2}.

### Schematic backscattering amplitude

At the volume-integral level, the flesh contribution may be written
schematically as:

f\_{\mathrm{bs}}^{(2)} = \frac{k_1^2}{4\pi} \iiint\_{V_2} \left(
\gamma\_{\kappa,21} - \gamma\_{\rho,21}\cos^2\beta \right) \exp\\\left(2
i \mathbf{k}\_2\cdot \mathbf{r}\right) \\ dV.

Here V_2 is the flesh-body volume, \mathbf{k}\_2 is the distorted
interior propagation vector, and \beta is the local angle between
propagation direction and body tangent.

For an elongated axisymmetric body, analytic transverse integration
reduces this volume expression to the one-dimensional DWBA body
integral.

## Backbone contribution

### Elastic-cylinder surrogate

The backbone is represented as a finite elastic cylinder rather than as
another weak-fluid inclusion. Its local physics therefore follows the
same elastic cylinder modal logic used by `ECMS`.

The elastic interior supports both longitudinal and transverse waves,
with wavenumbers:

k\_{L,3} = \frac{\omega}{c\_{L,3}}, \qquad k\_{T,3} =
\frac{\omega}{c\_{T,3}}.

For each cylindrical modal order m, the elastic boundary conditions
produce an order-dependent phase shift \eta_m. The finite-cylinder
backscattering amplitude may then be written schematically as:

f\_{\mathrm{bs}}^{(3)} = \frac{L_3}{\pi} \frac{\sin(k_1 L_3
\cos\theta_3)} {k_1 L_3 \cos\theta_3} \sum\_{m=0}^{\infty} (-1)^m
\epsilon_m \sin\eta_m e^{-i\eta_m}.

Here L_3 is backbone length, \theta_3 is the backbone incidence angle,
\epsilon_m is the usual Neumann factor, and \eta_m collects the
elastic-cylinder boundary-condition physics.

The backbone is therefore an elastic structure with
longitudinal-to-transverse wave conversion, not another weak contrast
perturbation.

## Spatial placement in the body frame

The flesh and backbone amplitudes cannot be added meaningfully unless
they are referred to the same spatial frame. The family uses the
body-fixed coordinate system for that purpose.

If the representative backbone position is \mathbf{r}\_c, then the
backbone amplitude is translated into the body frame by the monostatic
two-way phase factor:

\exp\\\left(2 i k_1
\hat{\mathbf{q}}\_{\mathrm{bs}}\cdot\mathbf{r}\_c\right).

where \hat{\mathbf{q}}\_{\mathrm{bs}} is the backscatter direction.

In the axisymmetric body-frame convention used here, that projection
becomes:

\hat{\mathbf{q}}\_{\mathrm{bs}}\cdot\mathbf{r}\_c = x_c\cos\theta +
z_c\sin\theta.

with (x_c, z_c) the backbone centroid and \theta the stored body angle.

This translation prevents the two components from being treated as if
they scattered from the same point.

## Coherent composite amplitude

Once the flesh and backbone terms are available in a common frame, the
total backscattering amplitude is:

f\_{\mathrm{bs}}^{(\mathrm{BBFM})} = f\_{\mathrm{bs}}^{(2)} +
f\_{\mathrm{bs}}^{(3)} \exp\\\left(2 i k_1
\hat{\mathbf{q}}\_{\mathrm{bs}}\cdot\mathbf{r}\_c\right).

Adding the amplitudes before squaring retains flesh-backbone
interference.

## Cross-section and interference structure

The linear backscattering cross-section is:

\sigma\_{\mathrm{bs}} =
\left\|f\_{\mathrm{bs}}^{(\mathrm{BBFM})}\right\|^2,

and the target strength is ([MacLennan et al.
2002](#ref-MacLennan_2002)):

\mathrm{TS} = 10 \log\_{10}\left(\frac{\sigma\_{\mathrm{bs}}}{1\\
\mathrm{m}^2}\right).

Expanding the squared magnitude makes the composite physics explicit:

\sigma\_{\mathrm{bs}} = \left\|f\_{\mathrm{bs}}^{(2)}\right\|^2 +
\left\|f\_{\mathrm{bs}}^{(3)}\right\|^2 + 2\\\Re\\\left\\
f\_{\mathrm{bs}}^{(2)} \overline{f\_{\mathrm{bs}}^{(3)}} \exp\\\left( -2
i k_1 \hat{\mathbf{q}}\_{\mathrm{bs}}\cdot\mathbf{r}\_c \right)
\right\\.

The third term is the interference term. It is frequency-dependent and
position-dependent, and it is the reason the composite `TS` does not
reduce to a simple sum of the flesh and backbone `TS` curves.

## What the family does and does not solve

### Included physics

`BBFM` explicitly includes:

1.  a weak-fluid flesh-body contribution,
2.  an elastic backbone contribution,
3.  coherent interference between those two components through a shared
    body frame.

### Excluded physics

`BBFM` does not yet solve:

1.  a true embedded elastic-cylinder-in-flesh transmission problem,
2.  repeated rescattering between flesh and backbone,
3.  shadowing or blockage of one component by the other,
4.  anatomical variability in backbone placement across an ensemble.

These omissions distinguish the component model from a fully coupled
composite-wave solution.

## Why this family is still useful

BBFM is useful when flesh is weakly scattering but the backbone is much
stiffer than the surrounding tissue. Keeping the components separate
exposes their individual amplitudes and their frequency-dependent
interference.

## References

Chu, Dezhang, and Zhen Ye. 1999. “A Phase-Compensated Distorted Wave
Born Approximation Representation of the Bistatic Scattering by Weakly
Scattering Objects: Application to Zooplankton.” *The Journal of the
Acoustical Society of America* 106 (4): 1732–43.
<https://doi.org/10.1121/1.428036>.

Gorska, Natalia, Egil Ona, and Rolf Korneliussen. 2005. “Acoustic
Backscattering by Atlantic Mackerel as Being Representative of Fish That
Lack a Swimbladder. Backscattering by Individual Fish.” *ICES Journal of
Marine Science* 62 (5): 984–95.
<https://doi.org/10.1016/j.icesjms.2005.03.010>.

MacLennan, David N., Percy G. Fernandes, and John Dalen. 2002. “A
Consistent Approach to Definitions and Symbols in Fisheries Acoustics.”
*ICES Journal of Marine Science* 59 (2): 365–69.
<https://doi.org/10.1006/jmsc.2001.1158>.

Stanton, Timothy K., Dezhang Chu, Peter H. Wiebe, Linda V. Martin, and
Robert L. Eastwood. 1998. “Sound Scattering by Several Zooplankton
Groups. I. Experimental Determination of Dominant Scattering
Mechanisms.” *The Journal of the Acoustical Society of America* 103 (1):
225–35. <https://doi.org/10.1121/1.421469>.
