# Orientation-average scattering from a stored TMM object

Reuses the stored TMM blocks to average the differential scattering
cross section over a user-supplied set of incident orientations. By
default the receive direction is taken to be the exact monostatic
direction for each supplied orientation.

## Usage

``` r
tmm_average_orientation(
  object,
  theta_body = NULL,
  weights = NULL,
  phi_body = pi,
  theta_scatter = NULL,
  phi_scatter = NULL,
  distribution = NULL
)
```

## Arguments

- object:

  Scatterer-object previously evaluated with
  `target_strength(..., model = "TMM", store_t_matrix = TRUE)`.

- theta_body:

  Numeric vector of incident polar angles (radians).

- weights:

  Optional numeric vector of averaging weights. When omitted, equal
  weights are used.

- phi_body:

  Incident azimuth angle(s) (radians). Either scalar or the same length
  as `theta_body`. Defaults to `pi`.

- theta_scatter:

  Receive polar angle(s) (radians). Either scalar or the same length as
  `theta_body`. Defaults to the monostatic direction.

- phi_scatter:

  Receive azimuth angle(s) (radians). Either scalar or the same length
  as `theta_body`. Defaults to the monostatic direction.

- distribution:

  Optional orientation distribution created by
  [`tmm_orientation_distribution`](https://brandynlucca.github.io/acousticTS/reference/tmm_orientation_distribution.md).
  When supplied, it overrides the direct `theta_body`, `weights`, and
  `phi_body` inputs.

## Value

A data frame containing the frequency, the orientation-averaged
differential backscattering cross section and the corresponding
orientation-averaged target strength.

## See also

[`tmm_scattering`](https://brandynlucca.github.io/acousticTS/reference/tmm_scattering.md)

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
distribution <- tmm_orientation_distribution(n_theta = 5)
tmm_average_orientation(stored, distribution = distribution)
#>   frequency     sigma_bs        TS
#> 1     12000 2.157667e-05 -46.66015
```
