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
  n_cores = parallel::detectCores() - 1,
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
  `parallel::detectCores() - 1`.

- verbose:

  Logical; whether to print progress and status messages to the console.
  Default is `TRUE`.

## Value

A data frame (or list of data frames) with simulation results.

## Details

For example, if `batch_by = "length"` and `parameters["length"]` is a
vector of values, simulations will be run for each value of length
`n_realizations` times. If multiple parameters are specified in
`batch_by`, batching will occur over all combinations of their values.

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

## DEPRECATION WARNING

The `simulate_ts` function will be deprecated in future versions and
will be replaced by the `anneal` function in future versions.
