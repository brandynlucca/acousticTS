# Choosing a Model

## Start with the physical problem

Model selection is a process of matching assumptions, not ranking
methods from best to worst. Geometry, boundary behavior, material
contrast, acoustic size, and the desired scattering quantity all matter.
A more elaborate formulation is not necessarily more appropriate for a
poorly matched target ([Clay and Horne 1994](#ref-Clay_1994); [Stanton
1996](#ref-Stanton_1996); [Jech et al. 2015](#ref-Jech_2015)).

Before choosing a model, identify:

1.  the target’s gross geometry and any important internal components,
2.  the applicable boundary condition or material description,
3.  the frequency and orientation range,
4.  whether the target is weakly scattering, resonant, or acoustically
    large,
5.  whether the calculation requires monostatic backscatter, bistatic
    scattering, orientation averaging, or retained complex amplitude,
    and
6.  whether the proposed model is validated over the required scope.

The flowchart provides a first pass through those decisions. It narrows
the candidate set but does not replace the assumptions and validation
statements on each model’s pages:

![Revised general model-selection
flowchart](model-selection-flowchart.png)[](https://brandynlucca.github.io/acousticTS/articles/krm/index.md "KRM overview")[](https://brandynlucca.github.io/acousticTS/articles/bbfm/index.md "BBFM overview")[](https://brandynlucca.github.io/acousticTS/articles/hpa/index.md "HPA overview")[](https://brandynlucca.github.io/acousticTS/articles/sphms/index.md "SPHMS overview")[](https://brandynlucca.github.io/acousticTS/articles/calibration/index.md "Calibration overview")[](https://brandynlucca.github.io/acousticTS/articles/fcms/index.md "FCMS overview")[](https://brandynlucca.github.io/acousticTS/articles/essms/index.md "ESSMS overview")[](https://brandynlucca.github.io/acousticTS/articles/trcm/index.md "TRCM overview")[](https://brandynlucca.github.io/acousticTS/articles/bcms/index.md "BCMS overview")[](https://brandynlucca.github.io/acousticTS/articles/ecms/index.md "ECMS overview")[](https://brandynlucca.github.io/acousticTS/articles/vesm/index.md "VESM overview")[](https://brandynlucca.github.io/acousticTS/articles/sdwba/index.md "SDWBA overview")[](https://brandynlucca.github.io/acousticTS/articles/dwba/index.md "DWBA overview")[](https://brandynlucca.github.io/acousticTS/articles/psms/index.md "PSMS overview")[](https://brandynlucca.github.io/acousticTS/articles/pcdwba/index.md "PCDWBA overview")[](https://brandynlucca.github.io/acousticTS/articles/tmm/index.md "TMM overview")[](https://brandynlucca.github.io/acousticTS/articles/tmm/index.md "TMM overview")[](https://brandynlucca.github.io/acousticTS/articles/tmm/index.md "TMM overview")[](https://brandynlucca.github.io/acousticTS/articles/tmm/index.md "TMM overview")

Select a model label to open that model family’s overview.

## Choose the physical representation first

### Canonical shapes

Canonical models are the natural starting point when a sphere, prolate
spheroid, or finite cylinder is a defensible representation. Their
coordinate systems follow the boundary closely and often permit exact or
modal-series solutions ([Morse and Feshbach 1953](#ref-Morse_1953);
[Bowman et al. 1987](#ref-Bowman_1987)). The boundary still has to
match. A rigid sphere, a fluid sphere, an elastic shell, and a solid
elastic calibration sphere share a shape but not the same interface
physics.

| Target representation | Primary candidates | Main distinction |
|:---|:---|:---|
| Homogeneous sphere | [SPHMS](https://brandynlucca.github.io/acousticTS/articles/sphms/index.md), [TMM](https://brandynlucca.github.io/acousticTS/articles/tmm/index.md) | SPHMS directly solves the spherical modal problem. TMM is useful when retained transition-matrix products or angular post-processing are required. |
| Solid elastic calibration sphere | [SOEMS](https://brandynlucca.github.io/acousticTS/articles/calibration/index.md) | Includes longitudinal and shear waves in the solid. |
| Elastic spherical shell | [ESSMS](https://brandynlucca.github.io/acousticTS/articles/essms/index.md) | Resolves shell elasticity and the interior fluid. |
| Prolate spheroid | [PSMS](https://brandynlucca.github.io/acousticTS/articles/psms/index.md), [TMM](https://brandynlucca.github.io/acousticTS/articles/tmm/index.md) | PSMS is the geometry-matched modal series. TMM supports a retained matrix workflow. |
| Finite circular cylinder | [FCMS](https://brandynlucca.github.io/acousticTS/articles/fcms/index.md), [TMM](https://brandynlucca.github.io/acousticTS/articles/tmm/index.md) | FCMS is the direct finite-cylinder modal family. The documented TMM cylindrical scope is narrower. |

[HPA](https://brandynlucca.github.io/acousticTS/articles/hpa/index.md)
can provide a faster asymptotic comparison for several smooth canonical
shapes. It does not reproduce all boundary or resonance physics of the
corresponding modal solution.

### Arbitrary or segmented bodies

For weakly scattering elongated fluid-like bodies,
[DWBA](https://brandynlucca.github.io/acousticTS/articles/dwba/index.md)
is usually the first candidate. It integrates contributions along a body
whose density and sound-speed contrasts remain sufficiently small.
[SDWBA](https://brandynlucca.github.io/acousticTS/articles/sdwba/index.md)
uses the same general physical setting but introduces phase variability
to represent unresolved roughness, posture, or other sources of
incoherence ([Stanton et al. 1998](#ref-Stanton_1998_1); [Demer and
Conti 2003](#ref-Demer_2003_1)). Use DWBA for a specified deterministic
target state. Use SDWBA when the distribution of unresolved phase
variation is part of the intended model.

[PCDWBA](https://brandynlucca.github.io/acousticTS/articles/pcdwba/index.md)
is available for polynomially described cylindrical bodies. Its geometry
and approximation assumptions should be checked against the PCDWBA
theory page before treating it as a general replacement for DWBA.

For a fish-like target in which a gas-filled swimbladder is an explicit
and important component,
[KRM](https://brandynlucca.github.io/acousticTS/articles/krm/index.md)
is often the most natural first model. Its hybrid treatment is intended
for the body-plus-swimbladder problem, not merely for any elongated
outline ([Clay 1992](#ref-Clay_1992); [Foote 1982](#ref-Foote_1982)).

[BBFM](https://brandynlucca.github.io/acousticTS/articles/bbfm/index.md)
represents a body assembled from boundary elements. It is relevant when
explicit boundary geometry is more important than a canonical or
centerline reduction. Review its current validation scope and
computational cost before selecting it only because the geometry is
complex.

## Match the acoustic regime

Geometry alone does not identify the dominant scattering mechanism. Two
models can accept similar shapes while describing different regimes:

- [TRCM](https://brandynlucca.github.io/acousticTS/articles/trcm/index.md)
  describes high-frequency, locally cylindrical scattering through
  coherent reflected and transmitted ray contributions.
- [DWBA](https://brandynlucca.github.io/acousticTS/articles/dwba/index.md)
  describes weak-contrast volume scattering.
- [FCMS](https://brandynlucca.github.io/acousticTS/articles/fcms/index.md)
  resolves finite-cylinder modal behavior.
- [HPA](https://brandynlucca.github.io/acousticTS/articles/hpa/index.md)
  retains broad high-frequency trends with fewer of the details present
  in a geometry-matched modal solution.

Acoustic size should therefore be considered together with contrast and
boundary condition. A cylinder does not become a TRCM problem solely
because it is cylindrical, and an elongated target does not become a
DWBA problem if its contrast violates the weak-scattering premise.

## Match the requested output

Most package models are organized around monostatic backscatter and
reported target strength. That is not the only possible scattering
calculation. The
[TMM](https://brandynlucca.github.io/acousticTS/articles/tmm/index.md)
family can retain transition-matrix information for bistatic scattering,
angular grids, and orientation post-processing. Select a model that
produces the quantity needed by the analysis rather than assuming a
backscatter spectrum can answer every angular-scattering question
([Mishchenko et al. 2002](#ref-Mishchenko_2002); [Waterman
1971](#ref-Waterman_1971)).

Output requirements can also affect whether complex amplitude must be
retained. Target strength is suitable for logarithmic reporting.
Coherent combination requires phase-bearing amplitude, while energetic
averaging is performed in a linear cross-section domain. See [Comparing
models on the same
target](https://brandynlucca.github.io/acousticTS/articles/comparing-models/comparing-models.md)
when the decision depends on how two defensible formulations differ.

## When more than one model is defensible

Model overlap is expected. Useful comparisons include:

- DWBA and SDWBA to assess the effect of phase variability,
- FCMS and TRCM to compare modal and high-frequency ray descriptions,
- SPHMS or PSMS and HPA to assess the loss of detail in an asymptotic
  approximation, and
- a canonical model and TMM when retained angular products are required.

Keep geometry, material properties, medium properties, orientation,
frequency, and output definition fixed during the comparison. A
difference then reflects model structure more cleanly. If those inputs
change at the same time, the exercise becomes a broader sensitivity
study.

## Check scope before committing to a model

A compatible object class only establishes that the software can
dispatch the model. It does not establish physical suitability. Before
using a result:

1.  read the model’s theory page and identify its boundary and
    approximation assumptions,
2.  read the implementation page for supported shapes, options, and
    numerical controls,
3.  inspect the model’s status in [Validation and benchmark
    reproduction](https://brandynlucca.github.io/acousticTS/articles/validation-benchmarks/validation-benchmarks.md),
    and
4.  test sensitivity to the model choice when a neighboring formulation
    is also defensible.

The best initial model is the simplest one that retains the geometry,
boundary physics, scattering configuration, and numerical fidelity
required by the scientific question.

## References

Bowman, J. J., T. B. A. Senior, and P. L. E. Uslenghi. 1987.
*Electromagnetic and Acoustic Scattering by Simple Shapes*. Hemisphere
Publishing Corp.

Clay, Clarence S. 1992. “Composite Ray-Mode Approximations for
Backscattered Sound from Gas-Filled Cylinders and Swimbladders.” *The
Journal of the Acoustical Society of America* 92 (4): 2173–80.
<https://doi.org/10.1121/1.405211>.

Clay, Clarence S., and John K. Horne. 1994. “Acoustic Models of Fish:
The Atlantic Cod (*Gadus Morhua*).” *The Journal of the Acoustical
Society of America* 96 (3): 1661–68. <https://doi.org/10.1121/1.410245>.

Demer, David A., and Stephane G. Conti. 2003. “Reconciling Theoretical
Versus Empirical Target Strengths of Krill: Effects of Phase Variability
on the Distorted-Wave Born Approximation.” *ICES Journal of Marine
Science* 60 (2): 429–34.
<https://doi.org/10.1016/S1054-3139(03)00002-X>.

Foote, Kenneth G. 1982. “Optimizing Copper Spheres for Precision
Calibration of Hydroacoustic Equipment.” *The Journal of the Acoustical
Society of America* 71 (3): 742–47. <https://doi.org/10.1121/1.387497>.

Jech, J. Michael, John K. Horne, Dezhang Chu, et al. 2015. “Comparisons
Among Ten Models of Acoustic Backscattering Used in Aquatic Ecosystem
Research.” *The Journal of the Acoustical Society of America* 138 (6):
3742–64. <https://doi.org/10.1121/1.4937607>.

Mishchenko, Michael I., Larry D. Travis, and Andrew A. Lacis. 2002.
*Scattering, Absorption, and Emission of Light by Small Particles*.
Cambridge University Press.

Morse, Philip M., and Herman Feshbach. 1953. *Methods of Theoretical
Physics*. McGraw-Hill.

Stanton, T. 1996. “Acoustic Scattering Characteristics of Several
Zooplankton Groups.” *ICES Journal of Marine Science* 53 (2): 289–95.
<https://doi.org/10.1006/jmsc.1996.0037>.

Stanton, Timothy K., Dezhang Chu, Peter H. Wiebe, Linda V. Martin, and
Robert L. Eastwood. 1998. “Sound Scattering by Several Zooplankton
Groups. I. Experimental Determination of Dominant Scattering
Mechanisms.” *The Journal of the Acoustical Society of America* 103 (1):
225–35. <https://doi.org/10.1121/1.421469>.

Waterman, Peter C. 1971. “Symmetry, Unitarity, and Geometry in
Electromagnetic Scattering.” *Physical Review D* 3 (4): 825–39.
<https://doi.org/10.1103/PhysRevD.3.825>.
