# Summarize bistatic products from a stored TMM object

Reuses the stored T-matrix blocks at one stored frequency to compute a
higher-level bistatic summary, including forward- and cross-plane
slices, peak-scattering direction, backscatter-lobe width, and
integrated scattering over coarse angular sectors. In the current
package build, this helper is available for the spherical and spheroidal
stored branches. Stored cylinders intentionally stop at exact monostatic
reuse until a validated retained cylinder angular operator is added.

## Usage

``` r
tmm_bistatic_summary(
  object,
  frequency = NULL,
  theta_body = NULL,
  phi_body = NULL,
  n_theta = 91,
  n_phi = 181,
  n_psi = 181,
  sectors = NULL,
  drop_dB = 3,
  include_grid = FALSE
)
```

## Arguments

- object:

  Scatterer-object previously evaluated with
  `target_strength(..., model = "TMM", store_t_matrix = TRUE)`.

- frequency:

  Stored frequency (Hz) to summarize. Required when the object contains
  more than one stored frequency.

- theta_body:

  Incident polar angle (radians). Defaults to the stored TMM incident
  angle.

- phi_body:

  Incident azimuth angle (radians). Defaults to the stored TMM incident
  angle.

- n_theta:

  Number of receive polar-angle samples used by the summary grid.

- n_phi:

  Number of receive azimuth samples used by the summary grid.

- n_psi:

  Number of forward-centered angular samples used for the local
  great-circle slices.

- sectors:

  Optional data frame with columns `sector`, `psi_min`, and `psi_max`
  (radians). When omitted, three coarse forward/oblique/backward sectors
  are used.

- drop_dB:

  Positive dB drop used to define the backscatter-lobe width on the
  forward-centered slice.

- include_grid:

  Logical; include the underlying scattering grid in the returned list.

## Value

A list containing scalar summary metrics, the named slice data frames,
sector integrals, and optionally the underlying scattering grid.

## See also

[`tmm_scattering_grid`](https://brandynlucca.github.io/acousticTS/reference/tmm_scattering_grid.md),
[`tmm_products`](https://brandynlucca.github.io/acousticTS/reference/tmm_products.md)
