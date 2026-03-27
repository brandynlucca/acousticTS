# Collect multiple post-processed products from one stored TMM solve

Provides a higher-level convenience interface for the stored-block TMM
workflow. One call can return the monostatic scattering spectrum, an
orientation average, and a higher-level bistatic summary without
rebuilding the underlying T-matrix solve. For stored cylinders, the
currently supported products are the exact monostatic scattering
spectrum and the corresponding orientation-averaged monostatic outputs;
cylinder bistatic summaries remain unavailable until a validated
retained cylinder angular operator is added.

## Usage

``` r
tmm_products(
  object,
  frequency = NULL,
  theta_body = NULL,
  phi_body = NULL,
  orientation = NULL,
  bistatic_summary = FALSE,
  include_grid = FALSE,
  n_theta = 91,
  n_phi = 181,
  n_psi = 181,
  sectors = NULL,
  drop_dB = 3
)
```

## Arguments

- object:

  Scatterer-object previously evaluated with
  `target_strength(..., model = "TMM", store_t_matrix = TRUE)`.

- frequency:

  Stored frequency (Hz) used when `bistatic_summary = TRUE`.

- theta_body:

  Incident polar angle (radians) for the monostatic and bistatic
  products. Defaults to the stored TMM incident angle.

- phi_body:

  Incident azimuth angle (radians) for the monostatic and bistatic
  products. Defaults to the stored TMM incident angle.

- orientation:

  Optional orientation distribution created by
  [`tmm_orientation_distribution`](https://brandynlucca.github.io/acousticTS/reference/tmm_orientation_distribution.md).

- bistatic_summary:

  Logical; include the output of
  [`tmm_bistatic_summary`](https://brandynlucca.github.io/acousticTS/reference/tmm_bistatic_summary.md).

- include_grid:

  Logical; keep the 2D scattering grid inside the bistatic summary
  output.

- n_theta:

  Number of receive polar-angle samples for the bistatic summary.

- n_phi:

  Number of receive azimuth samples for the bistatic summary.

- n_psi:

  Number of forward-centered angular samples used by the local summary
  slices.

- sectors:

  Optional angular-sector table passed to
  [`tmm_bistatic_summary`](https://brandynlucca.github.io/acousticTS/reference/tmm_bistatic_summary.md).

- drop_dB:

  Positive dB drop used to define the backscatter-lobe width.

## Value

A named list containing the requested post-processed TMM products.

## See also

[`tmm_scattering`](https://brandynlucca.github.io/acousticTS/reference/tmm_scattering.md),
[`tmm_average_orientation`](https://brandynlucca.github.io/acousticTS/reference/tmm_average_orientation.md),
[`tmm_bistatic_summary`](https://brandynlucca.github.io/acousticTS/reference/tmm_bistatic_summary.md)
