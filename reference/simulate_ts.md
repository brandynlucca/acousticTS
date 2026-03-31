# Simulate target strength (TS) with flexible parameterization and batching

Simulate target strength (TS) with flexible parameterization and
batching

## Usage

``` r
simulate_ts(
  object,
  frequency,
  model,
  n_realizations,
  parameters,
  batch_by = NULL,
  parallel = TRUE,
  n_cores = .default_simulation_cores(),
  verbose = TRUE
)
```

## Arguments

- object:

  Scatterer-class object.

- frequency:

  Frequency (Hz).

- model:

  Model name. If multiple models are specified, the output will be a
  list of data frames, one for each model.

- n_realizations:

  Number of realizations and output TS values.

- parameters:

  List containing the values, distributions, or generating functions of
  parameter values that inform the TS model.

- batch_by:

  Optional. Specifies which parameters in `parameters` to batch over.
  Simulations will be run for all combinations of these parameter
  values. Default is `NULL`.

- parallel:

  Logical; whether to parallelize the simulations. Default is `TRUE`.

- n_cores:

  Optional. Number of CPU cores to use for parallelization. Default is
  the smaller of 2 cores and `parallel::detectCores() - 1`.

- verbose:

  Logical; whether to print progress and status messages to the console.
  Default is `TRUE`.

## Value

A data frame when a single model is requested, or a named list of data
frames when multiple models are requested. Each returned data frame
contains the realized parameter values together with the modeled
acoustic output for each simulated run.

## Details

`simulate_ts()` is a workflow wrapper around repeated
[`target_strength()`](https://brandynlucca.github.io/acousticTS/reference/target_strength.md)
calls. It supports three broad parameter modes inside `parameters`:

- scalars that are recycled across every realization,

- explicit vectors that are either aligned with the full simulation grid
  or with one or more batched dimensions, and

- generating functions that are re-evaluated for each realization, and

- structured values such as named target-dimension vectors used by
  [`reforge()`](https://brandynlucca.github.io/acousticTS/reference/reforge.md)
  (for example `body_target = c(length = 0.03)`).

If `batch_by = "length"` and `parameters[["length"]]` is a vector of
candidate values, then simulations are run for each length value,
repeated `n_realizations` times. When multiple parameters are supplied
through `batch_by`, the function builds the full Cartesian grid of those
parameter values and runs the requested number of realizations inside
each batch cell.

Structured batch values should be wrapped in a list so that each
candidate is preserved as one unit. For example, use
`parameters = list(body_target = list(c(length = 0.02), c(length = 0.03))))`
when batching across multiple explicit
[`reforge()`](https://brandynlucca.github.io/acousticTS/reference/reforge.md)
targets.

Convenience dimension aliases are also supported for compatible
[`reforge()`](https://brandynlucca.github.io/acousticTS/reference/reforge.md)
methods. For example, `length_body = 0.03` is treated the same as
`body_target = c(length = 0.03)` for fluid-like scatterers, while
retaining the original `length_body` column in the returned simulation
output.

Parameter names are interpreted in the same way they would be if
supplied directly to
[`target_strength()`](https://brandynlucca.github.io/acousticTS/reference/target_strength.md)
or to the relevant object constructor /
[`reforge()`](https://brandynlucca.github.io/acousticTS/reference/reforge.md)
path. This means `simulate_ts()` can be used for:

- orientation perturbations,

- material-property perturbations,

- morphology studies that trigger shape rebuilding or reforging, and

- side-by-side comparisons across one or more model families.

## Parallelization

This function uses
[`pbapply::pblapply()`](https://peter.solymos.org/pbapply/reference/pbapply.html)
for parallelized simulation with progress bars. On Windows,
parallelization uses PSOCK clusters, which require all necessary objects
and packages to be exported to worker processes. On Unix-like systems,
forking is used, which is generally simpler.

## Performance Issues

Including too many parameters from `parameters` within `batch_by` may
cause significant performance issues or cause `R` to crash. If intensive
simulations are required, consider breaking them into more manageable
chunks

## Examples

``` r
shape_obj <- cylinder(
  length_body = 0.05,
  radius_body = 0.003,
  n_segments = 40
)

obj <- fls_generate(
  shape = shape_obj,
  density_body = 1045,
  sound_speed_body = 1520
)

res <- simulate_ts(
  object = obj,
  frequency = seq(38e3, 50e3, by = 6e3),
  model = "dwba",
  n_realizations = 2,
  parameters = list(
    theta_body = function() runif(1, min = 0.5 * pi, max = pi),
    density_body = 1045
  ),
  parallel = FALSE,
  verbose = FALSE
)

head(res)
#> $DWBA
#>   model realization theta_body density_body frequency        ka
#> 1  DWBA           1   1.697638         1045     38000 0.4775221
#> 2  DWBA           1   1.697638         1045     44000 0.5529203
#> 3  DWBA           1   1.697638         1045     50000 0.6283185
#> 4  DWBA           2   2.881364         1045     38000 0.4775221
#> 5  DWBA           2   2.881364         1045     44000 0.5529203
#> 6  DWBA           2   2.881364         1045     50000 0.6283185
#>                          f_bs     sigma_bs        TS
#> 1 -0.0001558136-7.493349e-20i 2.427787e-08 -76.14790
#> 2 -0.0002007945-1.118128e-19i 4.031842e-08 -73.94497
#> 3 -0.0002476169-1.566886e-19i 6.131414e-08 -72.12439
#> 4 -0.0001558136-7.493349e-20i 2.427787e-08 -76.14790
#> 5 -0.0002007945-1.118128e-19i 4.031842e-08 -73.94497
#> 6 -0.0002476169-1.566886e-19i 6.131414e-08 -72.12439
#> 
```
