# Simulation and Parameter Sweeps

## Overview

A single target-strength run answers a question about one specified
target state.
[`simulate_ts()`](https://brandynlucca.github.io/acousticTS/reference/simulate_ts.md)
repeats that workflow over deterministic parameter combinations,
stochastic realizations, or both. This is useful for orientation
distributions, morphology studies, uncertain material properties, and
matched comparisons among models ([Jech et al. 2015](#ref-Jech_2015);
[Demer and Conti 2005](#ref-Demer_2005)).

![A baseline scatterer, parameter rules, the simulation design, and the
collected model output.](simulation-batching-schematic.svg)

A baseline scatterer, parameter rules, the simulation design, and the
collected model output.

[`simulate_ts()`](https://brandynlucca.github.io/acousticTS/reference/simulate_ts.md)
calls
[`target_strength()`](https://brandynlucca.github.io/acousticTS/reference/target_strength.md)
for each job and returns a named list of data frames, one per model.
Each row retains the realized parameters beside the acoustic output.

## A minimal simulation

Begin with a valid `Scatterer`, then specify the frequency, model, and
parameters to vary:

``` r

library(acousticTS)

shape_obj <- cylinder(
  length_body = 0.05,
  radius_body = 0.003,
  n_segments = 80
)

obj <- fls_generate(
  shape = shape_obj,
  density_body = 1045,
  sound_speed_body = 1520
)

sim <- simulate_ts(
  object = obj,
  frequency = seq(38e3, 120e3, by = 6e3),
  model = "dwba",
  n_realizations = 5,
  parameters = list(
    theta_body = function() runif(1, 0.5 * pi, pi),
    density_body = 1045
  ),
  parallel = FALSE
)

head(sim$DWBA)
```

With no `parameters` and no `n_realizations`, the function performs one
run of the unchanged object. Use
[`target_strength()`](https://brandynlucca.github.io/acousticTS/reference/target_strength.md)
directly when that is all the workflow requires.

## How parameters are interpreted

The meaning of a parameter depends on how it is supplied:

![Fixed values, explicit vectors, and generating functions have
different simulation roles.](simulation-parameter-modes.svg)

Fixed values, explicit vectors, and generating functions have different
simulation roles.

- A scalar is fixed across all jobs.
- A deterministic vector with more than one value defines a sweep axis.
- A function is called once for each realization and supplies a new
  value.
- A structured value, such as a named `body_target` vector, is passed to
  the applicable object-rebuilding pathway as one unit.

For example:

``` r

parameters <- list(
  sound_speed_body = 1520,
  theta_body = seq(0.5 * pi, pi, length.out = 5),
  density_body = function() rnorm(1, mean = 1045, sd = 5),
  body_target = c(length = 0.045)
)
```

This defines five deterministic orientation cells. Within each cell,
`density_body` is redrawn for every realization. Sound speed and target
length remain fixed.

A multi-value vector is not matched automatically to `n_realizations`.
It defines a deterministic axis. Use a generating function for
realization-level randomness or `permute = FALSE` to pair several
deterministic axes.

## Design cells and realizations

Every deterministic multi-value parameter becomes an outer design axis.
`n_realizations` controls the number of repeated evaluations inside each
design cell:

![Cartesian parameter cells outside repeated within-cell
realizations.](simulation-batch-cells.svg)

Cartesian parameter cells outside repeated within-cell realizations.

The following design has four orientations, three densities, and three
realizations per cell. It therefore creates 36 jobs for each requested
model:

``` r

sim_grid <- simulate_ts(
  object = obj,
  frequency = seq(38e3, 120e3, by = 6e3),
  model = c("dwba", "sdwba"),
  n_realizations = 3,
  parameters = list(
    theta_body = seq(0.5 * pi, pi, length.out = 4),
    density_body = c(1035, 1045, 1055)
  ),
  batch_by = c("theta_body", "density_body"),
  parallel = FALSE
)
```

Multi-value parameters are detected automatically, so `batch_by` is not
required to make those vectors vary. Naming them in `batch_by` makes the
outer design explicit and is useful when reading or generating a
simulation call.

### Crossing or pairing axes

With `permute = TRUE`, the default, deterministic axes form a Cartesian
grid. Set `permute = FALSE` when their values describe matched states
that should advance together:

``` r

sim_paired <- simulate_ts(
  object = obj,
  frequency = 70e3,
  model = "dwba",
  parameters = list(
    theta_body = seq(0.5 * pi, pi, length.out = 4),
    density_body = c(1035, 1045, 1050, 1055)
  ),
  permute = FALSE,
  parallel = FALSE
)
```

This creates four paired states rather than sixteen crossed states. All
varied axes must have the same length. Length-one inputs remain fixed.

## Parameters that rebuild the target

Some parameters alter the target rather than only the model call.
Compatible arguments are applied through
[`reforge()`](https://brandynlucca.github.io/acousticTS/reference/reforge.md),
and matching component fields can be updated before the model is run.
Each output then represents a distinct target state derived from the
baseline object:

![A baseline scatterer is reforged into distinct geometric states before
each model run.](simulation-object-reforge.svg)

A baseline scatterer is reforged into distinct geometric states before
each model run.

For an `FLS` object, convenience names such as `length_body` use the
same rebuilding pathway as `body_target`:

``` r

sim_length <- simulate_ts(
  object = obj,
  frequency = 120e3,
  model = "dwba",
  n_realizations = 4,
  parameters = list(
    length_body = c(0.04, 0.05, 0.06),
    theta_body = function() runif(1, 0.5 * pi, pi)
  ),
  parallel = FALSE
)
```

Use a list when each sweep level is itself a structured target
specification:

``` r

sim_shape <- simulate_ts(
  object = obj,
  frequency = 120e3,
  model = "dwba",
  parameters = list(
    body_target = list(
      c(length = 0.04, radius = 0.0025),
      c(length = 0.06, radius = 0.0035)
    ),
    isometric_body = FALSE
  ),
  parallel = FALSE
)
```

For an already bent `FLS` object, `length_body` refers to the bent
centerline arc length.
[`simulate_ts()`](https://brandynlucca.github.io/acousticTS/reference/simulate_ts.md)
does not bend a straight baseline automatically. Create the intended
baseline with
[`brake()`](https://brandynlucca.github.io/acousticTS/reference/brake.md)
before the simulation. If bending must vary by realization, construct
those objects outside
[`simulate_ts()`](https://brandynlucca.github.io/acousticTS/reference/simulate_ts.md)
so the order of resizing and bending is explicit.

## Inspect and summarize output

The returned object is organized by model:

``` r

names(sim_grid)
head(sim_grid$DWBA)
head(sim_grid$SDWBA)
```

Group summaries should respect the design columns and `realization`
index. Common summaries include mean or quantile target strength by
frequency, orientation distributions, sensitivity across morphology
cells, and matched differences between models. Perform energetic
averages in a linear cross-section domain rather than averaging target
strength in decibels unless a specific reporting convention requires
otherwise.

## Parallel execution and reproducibility

The parallel pathway distributes the same logical job table over PSOCK
workers:

![Serial and parallel execution use the same simulation jobs with
different scheduling.](simulation-execution-layout.svg)

Serial and parallel execution use the same simulation jobs with
different scheduling.

Use `parallel = FALSE` while developing a design. PSOCK startup and,
when working from a source checkout, preparation of a worker
installation can cost more than the model calls for a small job.
Parallel execution becomes more useful once individual jobs are
sufficiently expensive.

``` r

sim_parallel <- simulate_ts(
  object = obj,
  frequency = seq(38e3, 120e3, by = 2e3),
  model = "dwba",
  n_realizations = 100,
  parameters = list(
    theta_body = function() runif(1, 0.5 * pi, pi)
  ),
  parallel = TRUE,
  n_cores = 2
)
```

Record the parameter definitions, package version, and random-number
settings for reproducible stochastic work. Verify a stochastic design
serially before assuming that worker scheduling will reproduce the same
draws as a sequential run.

## Plan the grid before running it

For crossed axes, the number of jobs per model is the product of the
axis lengths multiplied by `n_realizations`. Frequency points add rows
to each job’s result. Several moderate axes can therefore create a large
output and a large number of model calls.

Before a production run:

1.  calculate the number of design cells and expected output rows,
2.  run a small serial subset and inspect the realized parameters,
3.  confirm that crossed or paired semantics match the scientific
    design,
4.  verify that geometry-changing parameters rebuild the intended
    target, and
5.  split very large studies into recoverable batches or a cached
    project-level pipeline.

[`simulate_ts()`](https://brandynlucca.github.io/acousticTS/reference/simulate_ts.md)
is appropriate for a structured family of related model calls. For a
long-running study with external data dependencies, many intermediate
products, or expensive downstream summaries, a workflow manager such as
`targets` can orchestrate separate
[`simulate_ts()`](https://brandynlucca.github.io/acousticTS/reference/simulate_ts.md)
calls without becoming a dependency of `acousticTS`.

## References

Demer, David A., and Stéphane G. Conti. 2005. “New Target-Strength Model
Indicates More Krill in the Southern Ocean.” *ICES Journal of Marine
Science* 62 (1): 25–32. <https://doi.org/10.1016/j.icesjms.2004.07.027>.

Jech, J. Michael, John K. Horne, Dezhang Chu, et al. 2015. “Comparisons
Among Ten Models of Acoustic Backscattering Used in Aquatic Ecosystem
Research.” *The Journal of the Acoustical Society of America* 138 (6):
3742–64. <https://doi.org/10.1121/1.4937607>.
