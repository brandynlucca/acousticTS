# SDWBA

## Overview

Benchmarked Validated

*Model-family pages:*
[Overview](https://brandynlucca.github.io/acousticTS/articles/sdwba/index.md)
[Implementation](https://brandynlucca.github.io/acousticTS/articles/sdwba/sdwba-implementation.md)
[Theory](https://brandynlucca.github.io/acousticTS/articles/sdwba/sdwba-theory.md)

These pages connect krill-body DWBA models to phase variability,
orientation effects, and practical survey use ([Demer and Conti
2003a](#ref-demer_reconciling_2003); [Demer and Conti
2003b](#ref-demer_validation_2003), [2005](#ref-demer_new_2005); [Conti
and Demer 2006](#ref-conti_improved_2006)).

The stochastic distorted wave Born approximation (`SDWBA`) extends the
deterministic `DWBA` by treating unresolved posture and shape
variability as stochastic phase variability along the body.

### Core idea

Start from the same segmented weak-scattering sum as `DWBA`, then
replace the strictly coherent phase accumulation by a randomized phase
model whose statistics are chosen to mimic unresolved biological
variability.

### Best for

- Krill-like or zooplankton-like targets where deterministic body
  geometry is not known precisely
- Orientation-averaged or ensemble-style weak-scattering predictions
- Situations where fully coherent `DWBA` overpredicts narrow
  interference structure

### Supports

- `FLS` objects with the same geometry support as `DWBA`
- Monostatic target strength based on an averaged linear backscatter
  quantity
- The same local contrast notation as `DWBA`, with seawater as medium
  `1` and the body as medium `2`

### Main assumptions

- Weak-scattering fluid-like body
- Phase variability enters statistically rather than through an explicit
  new boundary-value solve
- Randomization acts on the coherent sum rather than on the local
  scattering kernel

### Validation status

- Benchmarked against published SDWBA weak-scattering comparison cases.
- Validated against the CCAMLR MATLAB and NOAA HTML SDWBA
  implementations.

### Family pages

- [Implementation](https://brandynlucca.github.io/acousticTS/articles/sdwba/sdwba-implementation.md):
  stochastic settings, spectra, and validation workflows
- [Theory](https://brandynlucca.github.io/acousticTS/articles/sdwba/sdwba-theory.md):
  randomized coherent sums and scale-invariant phase statistics

## References

Conti, Stéphane G., and David A. Demer. 2006. “Improved Parameterization
of the SDWBA for Estimating Krill Target Strength.” *ICES Journal of
Marine Science* 63 (5): 928–35.
<https://doi.org/10.1016/j.icesjms.2006.02.007>.

Demer, David A., and Stephane G. Conti. 2003a. “Reconciling Theoretical
Versus Empirical Target Strengths of Krill: Effects of Phase Variability
on the Distorted-Wave Born Approximation.” *ICES Journal of Marine
Science* 60 (2): 429–34.
<https://doi.org/10.1016/S1054-3139(03)00002-X>.

Demer, David A., and Stéphane G. Conti. 2003b. “Validation of the
Stochastic Distorted-Wave Born Approximation Model with Broad Bandwidth
Total Target Strength Measurements of Antarctic Krill.” *ICES Journal of
Marine Science* 60 (3): 625–35.
<https://doi.org/10.1016/S1054-3139(03)00063-8>.

———. 2005. “New Target-Strength Model Indicates More Krill in the
Southern Ocean.” *ICES Journal of Marine Science* 62 (1): 25–32.
<https://doi.org/10.1016/j.icesjms.2004.07.027>.
