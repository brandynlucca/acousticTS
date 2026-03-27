# Evaluate a 2D scattering grid from a stored TMM object

Reuses the stored T-matrix blocks to evaluate the far-field scattering
response over a two-dimensional receive-angle grid at one stored
frequency. This is useful for bistatic scattering maps, heatmaps, and
polar-style visualizations without rebuilding the retained modal solve.
In the current package build, this helper is available for the spherical
and spheroidal stored branches. Stored cylinders intentionally stop at
exact monostatic reuse until a validated retained cylinder angular
operator is added.

## Usage

``` r
tmm_scattering_grid(
  object,
  frequency = NULL,
  theta_body = NULL,
  phi_body = NULL,
  theta_scatter = NULL,
  phi_scatter = NULL,
  n_theta = 91,
  n_phi = 181
)
```

## Arguments

- object:

  Scatterer-object previously evaluated with
  `target_strength(..., model = "TMM", store_t_matrix = TRUE)`.

- frequency:

  Stored frequency (Hz) to evaluate. Required when the object contains
  more than one stored frequency.

- theta_body:

  Incident polar angle (radians). Defaults to the stored TMM incident
  angle.

- phi_body:

  Incident azimuth angle (radians). Defaults to the stored TMM incident
  angle.

- theta_scatter:

  Optional vector of receive polar angles (radians). Defaults to an
  evenly spaced grid on `[0, pi]`.

- phi_scatter:

  Optional vector of receive azimuth angles (radians). Defaults to an
  evenly spaced grid on `[0, 2*pi]`.

- n_theta:

  Number of default polar-angle grid points when `theta_scatter` is not
  supplied.

- n_phi:

  Number of default azimuth grid points when `phi_scatter` is not
  supplied.

## Value

A list containing the stored frequency, the incident angles used to
build the grid, the receive-angle vectors, and matrices for the complex
scattering amplitude, differential scattering cross section, and its
level in dB.

## See also

[`tmm_scattering`](https://brandynlucca.github.io/acousticTS/reference/tmm_scattering.md),
[`tmm_average_orientation`](https://brandynlucca.github.io/acousticTS/reference/tmm_average_orientation.md)
