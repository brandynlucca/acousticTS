# SDWBA

## Overview

Benchmarked Validated

[Theory](https://brandynlucca.github.io/acousticTS/articles/sdwba/sdwba-theory.md)
[Implementation](https://brandynlucca.github.io/acousticTS/articles/sdwba/sdwba-implementation.md)

The stochastic distorted-wave Born approximation (`SDWBA`) extends the
weak-fluid `DWBA` by averaging over random phase perturbations along the
body.

### Core idea

Retain the local weak-scattering contributions used by `DWBA`, perturb
their relative phases over repeated realizations, and report the mean
linear backscatter response ([Demer and Conti 2003](#ref-Demer_2003_1);
[Demer and Conti 2005](#ref-Demer_2005)).

### Best for

- Krill-like and zooplankton-like targets with unresolved shape or
  posture variability
- Ensemble predictions in which phase variability is part of the model
- Assessing how deterministic interference structure changes under
  randomization

### Supports

- `FLS` geometry and material inputs accepted by `DWBA`
- Configurable realization count and reference phase-scaling parameters
- Averaged complex amplitude, linear cross-section, and target strength

### Main assumptions

- The same weak-fluid and single-scattering regime as `DWBA`
- Unresolved variability can be represented by the specified phase
  distribution
- Randomization modifies coherence rather than the local scattering
  kernel
- Monte Carlo settings are adequate for the requested summary

### Validation status

- Benchmarked against the canonical spectra stored in benchmark_ts.
- Validated against the CCAMLR, NOAA applet, and echoSMs
  implementations.

### Family pages

- [Implementation](https://brandynlucca.github.io/acousticTS/articles/sdwba/sdwba-implementation.md):
  stochastic controls, reproducibility, output, and comparisons
- [Theory](https://brandynlucca.github.io/acousticTS/articles/sdwba/sdwba-theory.md):
  random phase model and ensemble averaging

## References

Demer, David A., and Stephane G. Conti. 2003. “Reconciling Theoretical
Versus Empirical Target Strengths of Krill: Effects of Phase Variability
on the Distorted-Wave Born Approximation.” *ICES Journal of Marine
Science* 60 (2): 429–34.
<https://doi.org/10.1016/S1054-3139(03)00002-X>.

Demer, David A., and Stéphane G. Conti. 2005. “New Target-Strength Model
Indicates More Krill in the Southern Ocean.” *ICES Journal of Marine
Science* 62 (1): 25–32. <https://doi.org/10.1016/j.icesjms.2004.07.027>.
