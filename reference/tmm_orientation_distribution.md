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

## Examples

``` r
tmm_orientation_distribution(
  distribution = "uniform",
  lower = pi / 3,
  upper = 2 * pi / 3,
  n_theta = 7
)
#>   theta_body phi_body    weights distribution
#> 1   1.047198 3.141593 0.08333333      uniform
#> 2   1.221730 3.141593 0.16666667      uniform
#> 3   1.396263 3.141593 0.16666667      uniform
#> 4   1.570796 3.141593 0.16666667      uniform
#> 5   1.745329 3.141593 0.16666667      uniform
#> 6   1.919862 3.141593 0.16666667      uniform
#> 7   2.094395 3.141593 0.08333333      uniform
```
