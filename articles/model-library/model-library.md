# Model Library

## Find a model family

`acousticTS` includes geometry-matched modal solutions, weak-scattering
and high-frequency approximations, composite-target models,
layered-sphere models, and a transition-matrix family. These models are
alternatives only when their geometry, boundary conditions, acoustic
regime, and outputs overlap. The inter-model literature provides useful
comparisons, but does not define one model as the universal default
([Jech et al. 2015](#ref-Jech_2015)).

Use [Choosing a
Model](https://brandynlucca.github.io/acousticTS/articles/model-selection/model-selection.md)
to narrow the library by target and scientific question. Select a family
name below to open its overview.

The badge icons summarize benchmark, validation, and lifecycle status.
Hover over a badge, or move keyboard focus to it, to read the evidence
and scope for that family. See [Validation and Benchmark
Reproduction](https://brandynlucca.github.io/acousticTS/articles/validation-benchmarks/validation-benchmarks.md)
for the full status definitions and evidence registry.

### Browse by family

### Modal-series families

#### [SPHMS](https://brandynlucca.github.io/acousticTS/articles/sphms/index.md)

Benchmarked Validated

Spherical modal-series solution for canonical spherical targets.

#### [FCMS](https://brandynlucca.github.io/acousticTS/articles/fcms/index.md)

Benchmarked Validated

Finite-cylinder modal-series solution for straight cylindrical targets.

#### [PSMS](https://brandynlucca.github.io/acousticTS/articles/psms/index.md)

Benchmarked Validated

Prolate-spheroidal modal-series solution for smooth elongated canonical
bodies.

#### [SOEMS](https://brandynlucca.github.io/acousticTS/articles/calibration/index.md)

Benchmarked Validated

Solid elastic spherical model used mainly for calibration spheres.

#### [ESSMS](https://brandynlucca.github.io/acousticTS/articles/essms/index.md)

Unvalidated

Elastic-shelled spherical family for layered shell targets.

#### [BCMS](https://brandynlucca.github.io/acousticTS/articles/bcms/index.md)

Unvalidated Experimental

Bent-cylinder modal-series family for straight and uniformly bent
cylinders.

#### [ECMS](https://brandynlucca.github.io/acousticTS/articles/ecms/index.md)

Unvalidated Experimental

Elastic-cylinder modal-series family for fully elastic solid cylinders.

### Approximation and ray-based families

#### [DWBA](https://brandynlucca.github.io/acousticTS/articles/dwba/index.md)

Benchmarked Validated

Weak-scattering elongated-body approximation for fluid-like targets.

#### [SDWBA](https://brandynlucca.github.io/acousticTS/articles/sdwba/index.md)

Benchmarked Validated

Stochastic DWBA family for unresolved phase variability.

#### [KRM](https://brandynlucca.github.io/acousticTS/articles/krm/index.md)

Benchmarked Validated

Kirchhoff-ray mode model for segmented fish-like body-plus-inclusion
targets.

#### [HPA](https://brandynlucca.github.io/acousticTS/articles/hpa/index.md)

Benchmarked Validated

High-pass approximation for compact asymptotic screening.

#### [TRCM](https://brandynlucca.github.io/acousticTS/articles/trcm/index.md)

Benchmarked Unvalidated

Two-ray cylindrical family for high-frequency locally cylindrical
targets.

#### [PCDWBA](https://brandynlucca.github.io/acousticTS/articles/pcdwba/index.md)

Validated Experimental

Phase-compensated DWBA for bent weakly scattering targets.

### Composite, layered, and transition-matrix families

#### [BBFM](https://brandynlucca.github.io/acousticTS/articles/bbfm/index.md)

Unvalidated Experimental

Composite flesh-plus-backbone family for swimbladder-less fish-like
targets.

#### [VESM](https://brandynlucca.github.io/acousticTS/articles/vesm/index.md)

Validated Experimental

Viscous-elastic layered-sphere family for gas-core, shell, and
viscous-layer targets.

#### [TMM](https://brandynlucca.github.io/acousticTS/articles/tmm/index.md)

Benchmarked Partially validated Experimental

Single-target transition-matrix family for retained monostatic and
angle-dependent scattering products across supported canonical shapes.

## What each family page contains

Each family has a local set of pages:

- **Overview** identifies the intended target, supported object types,
  and main assumptions.
- **Theory** develops the physical and mathematical formulation.
- **Implementation** documents package inputs, numerical controls,
  examples, figures, and validation details.

Model status can differ across branches of the same family. Read the
family overview and implementation page before treating a compatible
function call as evidence that the model is appropriate for a target.

## References

Jech, J. Michael, John K. Horne, Dezhang Chu, et al. 2015. “Comparisons
Among Ten Models of Acoustic Backscattering Used in Aquatic Ecosystem
Research.” *The Journal of the Acoustical Society of America* 138 (6):
3742–64. <https://doi.org/10.1121/1.4937607>.
