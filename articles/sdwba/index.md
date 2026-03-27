# SDWBA

## Overview

Benchmarked Validated

*Model-family pages:*
[Overview](https://brandynlucca.github.io/acousticTS/articles/sdwba/index.md)
[Implementation](https://brandynlucca.github.io/acousticTS/articles/sdwba/sdwba-implementation.md)
[Theory](https://brandynlucca.github.io/acousticTS/articles/sdwba/sdwba-theory.md)

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
