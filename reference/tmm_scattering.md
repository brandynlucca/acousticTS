# Evaluate scattering from a stored TMM object

Evaluates the far-field scattering amplitude from a previously computed
`TMM` object using the stored transition-matrix blocks. This allows the
same retained modal operator to be reused for arbitrary single-target
incident and receive-angle combinations without rebuilding the boundary
solve.

## Usage

``` r
tmm_scattering(
  object,
  theta_body = NULL,
  phi_body = NULL,
  theta_scatter = NULL,
  phi_scatter = NULL
)
```

## Arguments

- object:

  Scatterer-object previously evaluated with
  `target_strength(..., model = "TMM", store_t_matrix = TRUE)`.

- theta_body:

  Incident polar angle (radians). Defaults to the stored TMM incident
  angle.

- phi_body:

  Incident azimuth angle (radians). Defaults to the stored TMM incident
  angle.

- theta_scatter:

  Receive polar angle (radians). Defaults to the exact monostatic
  direction, `pi - theta_body`.

- phi_scatter:

  Receive azimuth angle (radians). Defaults to the exact monostatic
  direction, `phi_body + pi`.

## Value

A data frame with the frequency, complex scattering amplitude, the
corresponding differential cross section, and its level in dB.

## See also

[`target_strength()`](https://brandynlucca.github.io/acousticTS/reference/target_strength.md),
[`tmm_average_orientation`](https://brandynlucca.github.io/acousticTS/reference/tmm_average_orientation.md)
