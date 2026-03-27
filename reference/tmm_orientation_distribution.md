# Build an orientation distribution for stored TMM post-processing

Creates a validated set of incident angles and normalized weights that
can be reused by
[`tmm_average_orientation`](https://brandynlucca.github.io/acousticTS/reference/tmm_average_orientation.md)
or
[`tmm_products`](https://brandynlucca.github.io/acousticTS/reference/tmm_products.md).
The distributions defined here are distributions in `theta_body` itself
rather than isotropic solid-angle distributions.

## Usage

``` r
tmm_orientation_distribution(
  distribution = c("uniform", "normal", "truncated_normal", "quadrature", "pdf"),
  theta_body = NULL,
  weights = NULL,
  pdf = NULL,
  phi_body = pi,
  mean_theta = pi/2,
  sd_theta = pi/12,
  lower = 0,
  upper = pi,
  n_theta = 91
)
```

## Arguments

- distribution:

  Orientation-distribution type. One of `"uniform"`, `"normal"`,
  `"truncated_normal"`, `"quadrature"`, or `"pdf"`.

- theta_body:

  Optional numeric vector of incident polar angles (radians). Required
  for the `"quadrature"` and `"pdf"` pathways.

- weights:

  Optional numeric quadrature weights paired with `theta_body` for
  `distribution = "quadrature"`.

- pdf:

  Optional user-supplied density over `theta_body` for
  `distribution = "pdf"`. This can be either a numeric vector the same
  length as `theta_body` or a function evaluated at `theta_body`.

- phi_body:

  Incident azimuth angle(s) (radians). Either scalar or the same length
  as the resolved `theta_body` grid.

- mean_theta:

  Mean angle (radians) for the normal-family distributions.

- sd_theta:

  Standard deviation (radians) for the normal-family distributions.

- lower:

  Lower bound (radians) for the uniform and truncated-normal
  distributions.

- upper:

  Upper bound (radians) for the uniform and truncated-normal
  distributions.

- n_theta:

  Number of grid points for the analytic distributions.

## Value

A data frame with normalized orientation weights and class
`"TMMOrientationDistribution"`.

## See also

[`tmm_average_orientation`](https://brandynlucca.github.io/acousticTS/reference/tmm_average_orientation.md),
[`tmm_products`](https://brandynlucca.github.io/acousticTS/reference/tmm_products.md)
