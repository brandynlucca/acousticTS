# Diagnostics for stored single-target TMM solutions

Reuses the stored transition-matrix blocks to compute a compact set of
numerical and physics-based diagnostics for one or more stored
frequencies. The summary combines:

- monostatic reconstruction residuals,

- reciprocity residuals under incident/receive-angle exchange,

- an optical-theorem residual based on forward scattering and the
  integrated differential cross section,

- block-level conditioning indicators from the stored transition-matrix
  blocks,

- and, for prolate/oblate targets, an equal-volume sphere-to-spheroid
  continuation path that checks whether the monostatic response deforms
  smoothly away from the exact sphere limit.

These checks are meant to help distinguish "the post-processing is
self-consistent" from "the retained solve was also numerically
comfortable," which is especially helpful for the newer nonspherical TMM
branches. For stored cylinders, the current diagnostics are
intentionally limited to exact monostatic reconstruction because a
validated retained cylinder angular operator is not yet available.

## Usage

``` r
tmm_diagnostics(
  object,
  frequency = NULL,
  theta_body = NULL,
  phi_body = NULL,
  reciprocity_pairs = NULL,
  n_theta = 61,
  n_phi = 121,
  continuation_steps = 6L
)
```

## Arguments

- object:

  Scatterer-object previously evaluated with
  `target_strength(..., model = "TMM", store_t_matrix = TRUE)`.

- frequency:

  Optional stored frequency or vector of stored frequencies (Hz).
  Defaults to all stored frequencies.

- theta_body:

  Incident polar angle (radians) used for the monostatic and
  optical-theorem checks. Defaults to the stored TMM incident angle.

- phi_body:

  Incident azimuth angle (radians) used for the monostatic and
  optical-theorem checks. Defaults to the stored TMM incident angle.

- reciprocity_pairs:

  Optional data frame giving explicit reciprocity test angles. Must
  contain `theta_body`, `phi_body`, `theta_scatter`, and `phi_scatter`
  columns in radians.

- n_theta:

  Number of polar-angle grid points used by the optical-theorem
  integration check.

- n_phi:

  Number of azimuth-angle grid points used by the optical-theorem
  integration check.

- continuation_steps:

  Number of equal-volume aspect-ratio steps used for the
  sphere-to-spheroid continuation check on prolate and oblate targets.
  Set to `0` or `1` to skip the continuation path. Ignored for
  non-spheroidal targets.

## Value

A list with components:

- `summary`:

  Per-frequency diagnostic summary.

- `block_metrics`:

  Per-frequency block-level conditioning and transpose-residual
  summaries.

- `continuation`:

  Equal-volume sphere-to-spheroid continuation path for spheroidal
  targets, or `NULL` for other shapes.

## See also

[`tmm_scattering`](https://brandynlucca.github.io/acousticTS/reference/tmm_scattering.md),
[`tmm_scattering_grid`](https://brandynlucca.github.io/acousticTS/reference/tmm_scattering_grid.md),
[`tmm_products`](https://brandynlucca.github.io/acousticTS/reference/tmm_products.md)
