# Body-backbone fish model

## Introduction

Experimental Unvalidated

*Model-family pages:*
[Overview](https://brandynlucca.github.io/acousticTS/articles/bbfm/index.md)
[Implementation](https://brandynlucca.github.io/acousticTS/articles/bbfm/bbfm-implementation.md)
[Theory](https://brandynlucca.github.io/acousticTS/articles/bbfm/bbfm-theory.md)

This family is best read alongside the swimbladder-less fish and
composite-scatterer literature that motivates explicit flesh-body and
backbone terms ([Gorska, Ona, and Korneliussen
2005](#ref-gorska_etal_2005); [Stanton, Chu, and Wiebe
1998](#ref-stanton_sound_1998-1); [Clay and Horne
1994](#ref-clay_horne_1994)).

The body-backbone fish model (`BBFM`) is a composite scattering family
for swimbladder-less targets whose flesh body and backbone should remain
explicit, separately parameterized contributors. The point of the model
is not to solve a fully coupled three-medium boundary-value problem
exactly. The point is to keep the two dominant anatomical components
acoustically visible in one coherent backscatter calculation.

That makes `BBFM` conceptually closer to the hybrid logic of `KRM` than
to a single-region canonical modal solution:

1.  the flesh body is treated as a weakly scattering fluid-like region,
2.  the backbone is treated as an elastic cylindrical structure, and
3.  the two terms are embedded into one body-fixed frame through a phase
    translation before their complex amplitudes are summed.

The result is a transparent composite model rather than a homogenized
single-medium approximation.

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

The important approximation enters immediately here: the current
backbone term is not the exact solution for an elastic region `3`
embedded inside flesh region `2`. Instead, it is a seawater-referenced
elastic-cylinder surrogate that is then positioned inside the same body
frame as the flesh solve.

## Flesh-body contribution

### Weak-fluid assumption

The flesh component is treated as a weakly scattering fluid-like body.
Its physics therefore follows the distorted-wave Born logic: the
material contrasts relative to seawater are small enough that
first-order scattering remains meaningful, while the body is still
extended enough that phase accumulation along the body cannot be
ignored.

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

For the elongated axisymmetric bodies used in practice, this is reduced
to the usual one-dimensional `DWBA` body integral. The important point
for `BBFM` is not the exact quadrature form, but the physical role of
the term: it is the weak-fluid flesh-body amplitude.

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

This term is what distinguishes `BBFM` from a pure body-only
weak-scattering model. The backbone is not just another contrast
perturbation. It is an explicit elastic structure with its own internal
wave-conversion physics.

## Spatial placement in the body frame

The flesh and backbone amplitudes cannot be added meaningfully unless
they are referred to the same spatial frame. The current family uses the
body-fixed coordinate system for that purpose.

If the representative backbone position is \mathbf{r}\_c, then the
backbone amplitude is translated into the body frame by the monostatic
two-way phase factor:

\exp\\\left(2 i k_1
\hat{\mathbf{q}}\_{\mathrm{bs}}\cdot\mathbf{r}\_c\right).

where \hat{\mathbf{q}}\_{\mathrm{bs}} is the backscatter direction.

In the axisymmetric body-frame convention used by the current
implementation, that projection becomes:

\hat{\mathbf{q}}\_{\mathrm{bs}}\cdot\mathbf{r}\_c = x_c\cos\theta +
z_c\sin\theta.

with (x_c, z_c) the backbone centroid and \theta the stored body angle.

This is the step that embeds the backbone inside the same coordinate
frame as the flesh solve. Without it, the model would implicitly assume
that the flesh and backbone scatter from the same effective point.

## Coherent composite amplitude

Once the flesh and backbone terms are available in a common frame, the
total backscattering amplitude is:

f\_{\mathrm{bs}}^{(\mathrm{BBFM})} = f\_{\mathrm{bs}}^{(2)} +
f\_{\mathrm{bs}}^{(3)} \exp\\\left(2 i k_1
\hat{\mathbf{q}}\_{\mathrm{bs}}\cdot\mathbf{r}\_c\right).

This is the core `BBFM` statement. The family is coherent because it
adds the two complex amplitudes before squaring.

That distinction matters. If the model instead added cross-sections
directly, all interference between flesh and backbone would be lost.

## Cross-section and interference structure

The linear backscattering cross-section is:

\sigma\_{\mathrm{bs}} =
\left\|f\_{\mathrm{bs}}^{(\mathrm{BBFM})}\right\|^2,

and the target strength is:

\mathrm{TS} = 10 \log\_{10}\left(\sigma\_{\mathrm{bs}}\right).

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

These omissions matter because they define exactly what `BBFM` is: a
component-resolved coherent model, not a full coupled composite-wave
solver.

## Why this family is still useful

Even with those approximations, `BBFM` fills a real modeling gap. A
swimbladder-less fish often has flesh that behaves broadly like a weakly
scattering body and a backbone that is acoustically much stiffer than
the surrounding tissue.

Folding both of those into one effective region can hide the very
mechanism one is trying to study. `BBFM` therefore earns its place not
by being exact, but by keeping the dominant anatomical contributors
explicit while preserving coherent interference between them.

## References

Clay, Clarence S., and John K. Horne. 1994. “Acoustic Models of Fish:
The Atlantic Cod (*Gadus Morhua*).” *The Journal of the Acoustical
Society of America* 96 (3): 1661–68. <https://doi.org/10.1121/1.410245>.

Gorska, Natalia, Egil Ona, and Rolf Korneliussen. 2005. “Acoustic
Backscattering by Atlantic Mackerel as Being Representative of Fish That
Lack a Swimbladder. Backscattering by Individual Fish.” *ICES Journal of
Marine Science* 62 (5): 984–95.
<https://doi.org/10.1016/j.icesjms.2005.03.010>.

Stanton, Timothy K., Dezhang Chu, and Peter H. Wiebe. 1998. “Sound
Scattering by Several Zooplankton Groups. II. Scattering Models.” *The
Journal of the Acoustical Society of America* 103 (1): 236–53.
<https://doi.org/10.1121/1.421110>.
