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
tmm_scattering_grid(stored, n_theta = 5, n_phi = 9)
#> $frequency
#> [1] 12000
#> 
#> $theta_body
#> [1] 1.570796
#> 
#> $phi_body
#> [1] 3.141593
#> 
#> $theta_scatter
#> [1] 0.0000000 0.7853982 1.5707963 2.3561945 3.1415927
#> 
#> $phi_scatter
#> [1] 0.0000000 0.7853982 1.5707963 2.3561945 3.1415927 3.9269908 4.7123890
#> [8] 5.4977871 6.2831853
#> 
#> $f_scat
#>                           [,1]                      [,2]
#> [1,] -0.004784578+0.001248489i -0.004784578+0.001248489i
#> [2,] -0.004564681+0.001247307i -0.004628729+0.001247653i
#> [3,] -0.004474609+0.001246817i -0.004564681+0.001247307i
#> [4,] -0.004564681+0.001247307i -0.004628729+0.001247653i
#> [5,] -0.004784578+0.001248489i -0.004784578+0.001248489i
#>                           [,3]                      [,4]
#> [1,] -0.004784578+0.001248489i -0.004784578+0.001248489i
#> [2,] -0.004784578+0.001248489i -0.004942167+0.001249324i
#> [3,] -0.004784578+0.001248489i -0.005007954+0.001249670i
#> [4,] -0.004784578+0.001248489i -0.004942167+0.001249324i
#> [5,] -0.004784578+0.001248489i -0.004784578+0.001248489i
#>                           [,5]                      [,6]
#> [1,] -0.004784578+0.001248489i -0.004784578+0.001248489i
#> [2,] -0.005007954+0.001249670i -0.004942167+0.001249324i
#> [3,] -0.005101506+0.001250160i -0.005007954+0.001249670i
#> [4,] -0.005007954+0.001249670i -0.004942167+0.001249324i
#> [5,] -0.004784578+0.001248489i -0.004784578+0.001248489i
#>                           [,7]                      [,8]
#> [1,] -0.004784578+0.001248489i -0.004784578+0.001248489i
#> [2,] -0.004784578+0.001248489i -0.004628729+0.001247653i
#> [3,] -0.004784578+0.001248489i -0.004564681+0.001247307i
#> [4,] -0.004784578+0.001248489i -0.004628729+0.001247653i
#> [5,] -0.004784578+0.001248489i -0.004784578+0.001248489i
#>                           [,9]
#> [1,] -0.004784578+0.001248489i
#> [2,] -0.004564681+0.001247307i
#> [3,] -0.004474609+0.001246817i
#> [4,] -0.004564681+0.001247307i
#> [5,] -0.004784578+0.001248489i
#> 
#> $sigma_scat
#>              [,1]         [,2]         [,3]         [,4]         [,5]
#> [1,] 2.445091e-05 2.445091e-05 2.445091e-05 2.445091e-05 2.445091e-05
#> [2,] 2.239209e-05 2.298177e-05 2.445091e-05 2.598582e-05 2.664128e-05
#> [3,] 2.157667e-05 2.239209e-05 2.445091e-05 2.664128e-05 2.758827e-05
#> [4,] 2.239209e-05 2.298177e-05 2.445091e-05 2.598582e-05 2.664128e-05
#> [5,] 2.445091e-05 2.445091e-05 2.445091e-05 2.445091e-05 2.445091e-05
#>              [,6]         [,7]         [,8]         [,9]
#> [1,] 2.445091e-05 2.445091e-05 2.445091e-05 2.445091e-05
#> [2,] 2.598582e-05 2.445091e-05 2.298177e-05 2.239209e-05
#> [3,] 2.664128e-05 2.445091e-05 2.239209e-05 2.157667e-05
#> [4,] 2.598582e-05 2.445091e-05 2.298177e-05 2.239209e-05
#> [5,] 2.445091e-05 2.445091e-05 2.445091e-05 2.445091e-05
#> 
#> $sigma_scat_dB
#>           [,1]      [,2]      [,3]      [,4]      [,5]      [,6]      [,7]
#> [1,] -46.11705 -46.11705 -46.11705 -46.11705 -46.11705 -46.11705 -46.11705
#> [2,] -46.49905 -46.38616 -46.11705 -45.85264 -45.74445 -45.85264 -46.11705
#> [3,] -46.66015 -46.49905 -46.11705 -45.74445 -45.59276 -45.74445 -46.11705
#> [4,] -46.49905 -46.38616 -46.11705 -45.85264 -45.74445 -45.85264 -46.11705
#> [5,] -46.11705 -46.11705 -46.11705 -46.11705 -46.11705 -46.11705 -46.11705
#>           [,8]      [,9]
#> [1,] -46.11705 -46.11705
#> [2,] -46.38616 -46.49905
#> [3,] -46.49905 -46.66015
#> [4,] -46.38616 -46.49905
#> [5,] -46.11705 -46.11705
#> 
```
