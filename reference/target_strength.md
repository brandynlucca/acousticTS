# Wrapper function to model acoustic target strength

Wrapper function to model acoustic target strength

## Usage

``` r
target_strength(
  object,
  frequency,
  model,
  verbose = FALSE,
  model_args = NULL,
  ...
)
```

## Arguments

- object:

  Scatterer-class object.

- frequency:

  Frequency (Hz).

- model:

  Model name. Available models currently include `"dwba"`
  ([`DWBA`](https://brandynlucca.github.io/acousticTS/reference/DWBA.md)),
  `"bbfm"`
  ([`BBFM`](https://brandynlucca.github.io/acousticTS/reference/BBFM.md)),
  `"pcdwba"`
  ([`PCDWBA`](https://brandynlucca.github.io/acousticTS/reference/PCDWBA.md)),
  `"sdwba"`
  ([`SDWBA`](https://brandynlucca.github.io/acousticTS/reference/SDWBA.md)),
  `"fcms"`
  ([`FCMS`](https://brandynlucca.github.io/acousticTS/reference/FCMS.md)),
  `"bcms"`
  ([`BCMS`](https://brandynlucca.github.io/acousticTS/reference/BCMS.md)),
  `"ecms"`
  ([`ECMS`](https://brandynlucca.github.io/acousticTS/reference/ECMS.md)),
  `"hpa"`
  ([`HPA`](https://brandynlucca.github.io/acousticTS/reference/HP.md)),
  `"krm"`
  ([`KRM`](https://brandynlucca.github.io/acousticTS/reference/KRM.md)),
  `"psms"`
  ([`PSMS`](https://brandynlucca.github.io/acousticTS/reference/PSMS.md)),
  `"sphms"`
  ([`SPHMS`](https://brandynlucca.github.io/acousticTS/reference/SPHMS.md)),
  `"essms"`
  ([`ESSMS`](https://brandynlucca.github.io/acousticTS/reference/ESSMS.md)),
  `"vesms"`
  ([`VESMS`](https://brandynlucca.github.io/acousticTS/reference/VESMS.md)),
  `"trcm"`
  ([`TRCM`](https://brandynlucca.github.io/acousticTS/reference/TRCM.md)),
  and `"calibration"` / `"soems"`
  ([`SOEMS`](https://brandynlucca.github.io/acousticTS/reference/SOEMS.md)).

- verbose:

  Prints current procedural step occurring from model initialization to
  calculating TS. Defaults to FALSE.

- model_args:

  Optional named list of per-model argument bundles. Each list name
  should match one of the requested model names case-insensitively, and
  each value should be either a named list or a named atomic vector of
  arguments to apply only to that model. When the same argument is
  supplied both through `...` and through `model_args[[model_name]]`,
  the model-specific entry takes precedence.

- ...:

  Additional optional model inputs/parameters.

## Value

The input scatterer object with requested model parameters, model
outputs, and target strength results stored in its model slots.

## Details

This is the main high-level entry point for running target-strength
models in `acousticTS`. The supplied scatterer object is checked, the
requested model or models are initialized, and the resulting outputs are
stored back on the same object.

The available model families span exact modal-series solutions and
approximation-based solutions. Readers should consult the model-specific
help topics for the physical assumptions, valid object types,
boundary-condition options, and model-specific arguments used by each
implementation:

- [`DWBA`](https://brandynlucca.github.io/acousticTS/reference/DWBA.md)
  for the distorted-wave Born approximation applied to weakly scattering
  fluid-like bodies.

- [`BBFM`](https://brandynlucca.github.io/acousticTS/reference/BBFM.md)
  for a composite flesh-plus-backbone model that combines a DWBA flesh
  term with an elastic-cylinder backbone term.

- [`PCDWBA`](https://brandynlucca.github.io/acousticTS/reference/PCDWBA.md)
  for the phase-compensated distorted-wave Born approximation applied to
  elongated fluid-like bodies.

- [`SDWBA`](https://brandynlucca.github.io/acousticTS/reference/SDWBA.md)
  for the stochastic distorted-wave Born approximation.

- [`FCMS`](https://brandynlucca.github.io/acousticTS/reference/FCMS.md)
  for the finite cylinder modal series solution.

- [`BCMS`](https://brandynlucca.github.io/acousticTS/reference/BCMS.md)
  for the bent-cylinder modal series solution.

- [`ECMS`](https://brandynlucca.github.io/acousticTS/reference/ECMS.md)
  for the elastic-cylinder modal series solution.

- [`HPA`](https://brandynlucca.github.io/acousticTS/reference/HP.md) for
  the high-pass approximation family.

- [`KRM`](https://brandynlucca.github.io/acousticTS/reference/KRM.md)
  for the Kirchhoff-Ray Mode model.

- [`PSMS`](https://brandynlucca.github.io/acousticTS/reference/PSMS.md)
  for the prolate spheroidal modal series solution.

- [`SPHMS`](https://brandynlucca.github.io/acousticTS/reference/SPHMS.md)
  for the spherical modal series solution.

- [`ESSMS`](https://brandynlucca.github.io/acousticTS/reference/ESSMS.md)
  for the elastic-shelled spherical modal series solution.

- [`VESMS`](https://brandynlucca.github.io/acousticTS/reference/VESMS.md)
  for the viscous-elastic spherical scattering model applied to
  gas-filled elastic shells with an external viscous layer.

- [`TMM`](https://brandynlucca.github.io/acousticTS/reference/TMM.md)
  for the single-target transition-matrix family.

- [`TRCM`](https://brandynlucca.github.io/acousticTS/reference/TRCM.md)
  for the two-ray cylindrical model.

- [`SOEMS`](https://brandynlucca.github.io/acousticTS/reference/SOEMS.md)
  for the solid elastic calibration-sphere model, accessed through
  `"calibration"` or `"soems"`.

Model-specific inputs are usually passed through `...`. For example,
some models require a `boundary` argument, `HPA` uses a `method`
argument, and several models expose additional numerical controls. When
several models are requested together, shared arguments may be supplied
through `...` and per-model overrides may be supplied through
`model_args`. This is useful when different models should share the same
seawater properties but only one of them needs an extra stochastic or
numerical control:

    target_strength(
      object,
      frequency,
      model = c("dwba", "sdwba"),
      density_sw = 1026,
      sound_speed_sw = 1478,
      model_args = list(
        sdwba = list(phase_sd_init = 0.77)
      )
    )

The legacy curved-entry wrappers `"dwba_curved"` and `"sdwba_curved"`
are deprecated; apply
[`brake`](https://brandynlucca.github.io/acousticTS/reference/brake.md)
to the scatterer first, then run `"dwba"` or `"sdwba"` on the curved
object. Model names are normalized internally, so case-insensitive
inputs such as `"DWBA"` and `"dwba"` resolve to the same family.

## See also

[`DWBA`](https://brandynlucca.github.io/acousticTS/reference/DWBA.md),
[`BBFM`](https://brandynlucca.github.io/acousticTS/reference/BBFM.md),
[`PCDWBA`](https://brandynlucca.github.io/acousticTS/reference/PCDWBA.md),
[`SDWBA`](https://brandynlucca.github.io/acousticTS/reference/SDWBA.md),
[`FCMS`](https://brandynlucca.github.io/acousticTS/reference/FCMS.md),
[`BCMS`](https://brandynlucca.github.io/acousticTS/reference/BCMS.md),
[`ECMS`](https://brandynlucca.github.io/acousticTS/reference/ECMS.md),
[`HPA`](https://brandynlucca.github.io/acousticTS/reference/HP.md),
[`KRM`](https://brandynlucca.github.io/acousticTS/reference/KRM.md),
[`PSMS`](https://brandynlucca.github.io/acousticTS/reference/PSMS.md),
[`SPHMS`](https://brandynlucca.github.io/acousticTS/reference/SPHMS.md),
[`ESSMS`](https://brandynlucca.github.io/acousticTS/reference/ESSMS.md),
[`VESMS`](https://brandynlucca.github.io/acousticTS/reference/VESMS.md),
[`TMM`](https://brandynlucca.github.io/acousticTS/reference/TMM.md),
[`TRCM`](https://brandynlucca.github.io/acousticTS/reference/TRCM.md),
[`SOEMS`](https://brandynlucca.github.io/acousticTS/reference/SOEMS.md)
