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

## Examples

``` r
target <- fls_generate(
  shape = sphere(radius_body = 0.005, n_segments = 20),
  g_body = 1,
  h_body = 1
)
stored <- target_strength(
  target,
  frequency = 12e3,
  model = "tmm",
  boundary = "pressure_release",
  store_t_matrix = TRUE
)
tmm_diagnostics(stored, n_theta = 5, n_phi = 9)
#> $summary
#>    shape coordinate_system         boundary frequency monostatic_abs_residual
#> 1 Sphere         spherical pressure_release     12000            8.673617e-19
#>   monostatic_rel_residual reciprocity_rel_residual  sigma_total    sigma_ext
#> 1            1.867274e-16             4.251432e-17 0.0002918541 0.0003078102
#>   optical_theorem_sign optical_theorem_rel_residual min_block_rcond
#> 1                    1                    0.0518376    1.128043e-06
#>   max_block_cond_est max_block_transpose_residual
#> 1           886491.1                 2.066084e-17
#>   continuation_target_aspect_ratio continuation_max_abs_step_TS
#> 1                               NA                           NA
#>   continuation_max_abs_second_diff_TS continuation_any_nonfinite
#> 1                                  NA                         NA
#> 
#> $block_metrics
#> $block_metrics$`12000`
#>   m n_terms        rcond transpose_residual
#> 1 0       6 1.255165e-05       2.066084e-17
#> 2 1       5 1.128043e-06       3.690615e-18
#> 3 2       4 3.340586e-06       4.772712e-18
#> 4 3       3 5.818129e-05       4.471309e-18
#> 5 4       2 3.854585e-03       2.200130e-18
#> 6 5       1 1.000000e+00       0.000000e+00
#> 
#> 
#> $continuation
#> NULL
#> 
```
