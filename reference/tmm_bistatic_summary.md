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
tmm_bistatic_summary(
  stored,
  n_theta = 5,
  n_phi = 9,
  n_psi = 9
)
#> $metrics
#>   frequency forward_sigma_scat forward_sigma_scat_dB     sigma_bs        TS
#> 1     12000       2.758827e-05             -45.59276 2.157667e-05 -46.66015
#>   peak_theta_scatter peak_phi_scatter peak_psi_scatter peak_sigma_scat
#> 1           1.570796         3.141593                0    2.758827e-05
#>   peak_sigma_scat_dB backscatter_lobe_width
#> 1          -45.59276               3.141593
#> 
#> $slices
#> $slices$forward_scatter
#>                   slice psi_scatter theta_scatter phi_scatter
#> 1 forward_scatter_slice   0.0000000      1.570796    3.141593
#> 2 forward_scatter_slice   0.3926991      1.570796    3.534292
#> 3 forward_scatter_slice   0.7853982      1.570796    3.926991
#> 4 forward_scatter_slice   1.1780972      1.570796    4.319690
#> 5 forward_scatter_slice   1.5707963      1.570796    4.712389
#> 6 forward_scatter_slice   1.9634954      1.570796    5.105088
#> 7 forward_scatter_slice   2.3561945      1.570796    5.497787
#> 8 forward_scatter_slice   2.7488936      1.570796    5.890486
#> 9 forward_scatter_slice   3.1415927      1.570796    6.283185
#>                      f_scat   sigma_scat sigma_scat_dB
#> 1 -0.005101506+0.001250160i 2.758827e-05     -45.59276
#> 2 -0.005077135+0.001250033i 2.733988e-05     -45.63203
#> 3 -0.005007954+0.001249670i 2.664128e-05     -45.74445
#> 4 -0.004905034+0.001249128i 2.561968e-05     -45.91426
#> 5 -0.004784578+0.001248489i 2.445091e-05     -46.11705
#> 6 -0.004665141+0.001247849i 2.332067e-05     -46.32259
#> 7 -0.004564681+0.001247307i 2.239209e-05     -46.49905
#> 8 -0.004497961+0.001246944i 2.178652e-05     -46.61812
#> 9 -0.004474609+0.001246817i 2.157667e-05     -46.66015
#> 
#> $slices$dorsal_ventral
#>                  slice psi_scatter theta_scatter phi_scatter
#> 1 dorsal_ventral_slice   0.0000000     1.5707963    3.141593
#> 2 dorsal_ventral_slice   0.3926991     1.1780972    3.141593
#> 3 dorsal_ventral_slice   0.7853982     0.7853982    3.141593
#> 4 dorsal_ventral_slice   1.1780972     0.3926991    3.141593
#> 5 dorsal_ventral_slice   1.5707963     0.0000000    4.712389
#> 6 dorsal_ventral_slice   1.9634954     0.3926991    6.283185
#> 7 dorsal_ventral_slice   2.3561945     0.7853982    6.283185
#> 8 dorsal_ventral_slice   2.7488936     1.1780972    6.283185
#> 9 dorsal_ventral_slice   3.1415927     1.5707963    6.283185
#>                      f_scat   sigma_scat sigma_scat_dB
#> 1 -0.005101506+0.001250160i 2.758827e-05     -45.59276
#> 2 -0.005077135+0.001250033i 2.733988e-05     -45.63203
#> 3 -0.005007954+0.001249670i 2.664128e-05     -45.74445
#> 4 -0.004905034+0.001249128i 2.561968e-05     -45.91426
#> 5 -0.004784578+0.001248489i 2.445091e-05     -46.11705
#> 6 -0.004665141+0.001247849i 2.332067e-05     -46.32259
#> 7 -0.004564681+0.001247307i 2.239209e-05     -46.49905
#> 8 -0.004497961+0.001246944i 2.178652e-05     -46.61812
#> 9 -0.004474609+0.001246817i 2.157667e-05     -46.66015
#> 
#> 
#> $sector_integrals
#>            sector  psi_min  psi_max integrated_sigma_scat
#> 1  forward_sector 0.000000 1.047198          7.126088e-05
#> 2  oblique_sector 1.047198 2.094395          1.385379e-04
#> 3 backward_sector 2.094395 3.141593          9.800041e-05
#> 
```
