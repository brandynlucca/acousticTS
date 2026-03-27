# DWBA

## Overview

Benchmarked Validated

*Model-family pages:*
[Overview](https://brandynlucca.github.io/acousticTS/articles/dwba/index.md)
[Implementation](https://brandynlucca.github.io/acousticTS/articles/dwba/dwba-implementation.md)
[Theory](https://brandynlucca.github.io/acousticTS/articles/dwba/dwba-theory.md)

These pages follow the weak-scattering elongated-body formulation and
later applied fisheries-acoustics usage of the distorted-wave Born
approximation ([Morse and Ingard 1986](#ref-morse_theoretical_1986);
[Stanton, Chu, and Wiebe 1998](#ref-stanton_sound_1998-1); [Gastauer,
Chu, and Cox 2019](#ref-gastauer_zooscatrspan_2019)).

The distorted wave Born approximation (`DWBA`) is the package’s main
weak-scattering fluid-body model for elongated targets that are no
longer well represented by a single canonical sphere, cylinder, or
spheroid.

### Core idea

Treat the target as a weak perturbation of the surrounding fluid, retain
the two-way phase accumulation along the body, and integrate the local
cross-sectional response over a segmented centerline.

### Best for

- Weakly scattering fluid-like elongated bodies
- Zooplankton-like or swimbladder-less flesh-body calculations where
  density and sound-speed contrasts remain modest
- Arbitrary profiles that do not fit a single canonical exact geometry

### Supports

- `FLS` objects with canonical or arbitrary shapes
- Material contrasts expressed naturally as g\_{21} and h\_{21} relative
  to seawater
- Segmented body-axis integrations in the monostatic backscatter setting

### Main assumptions

- Small density and sound-speed contrasts
- First-order Born linearization
- Single scattering with no strong internal reverberation
- Fluid-like, non-elastic interior response

### Validation status

- Benchmarked against the canonical weakly scattering targets summarized
  by Jech et al. (2015).
- Validated against the published McGehee MATLAB workflow and an
  independent DWBA implementation.

### Family pages

- [Implementation](https://brandynlucca.github.io/acousticTS/articles/dwba/dwba-implementation.md):
  object workflows, spectra, and validation tables
- [Theory](https://brandynlucca.github.io/acousticTS/articles/dwba/dwba-theory.md):
  Born linearization, contrast term, and centerline integral

## References

Gastauer, Sven, Dezhang Chu, and Martin J. Cox. 2019. “ZooScatR—An
\<Span Style="font-Variant:small-Caps;"\>r\</Span\> Package for
Modelling the Scattering Properties of Weak Scattering Targets Using the
Distorted Wave Born Approximation.” *The Journal of the Acoustical
Society of America* 145 (1): EL102–8.
<https://doi.org/10.1121/1.5085655>.

Morse, Philip M., and K. Uno Ingard. 1986. *Theoretical Acoustics*.
Princeton, N.J: Princeton University Press.

Stanton, Timothy K., Dezhang Chu, and Peter H. Wiebe. 1998. “Sound
Scattering by Several Zooplankton Groups. II. Scattering Models.” *The
Journal of the Acoustical Society of America* 103 (1): 236–53.
<https://doi.org/10.1121/1.421110>.
