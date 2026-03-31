# Simulation and Parameter Sweeps

## Introduction

Parameter sweeps are especially useful when a model has to be
interpreted over the orientation, size, and contrast ranges emphasized
in benchmark and survey studies ([Jech et al.
2015](#ref-jech_etal_2015); [Demer and Conti
2005](#ref-demer_new_2005)).

Single deterministic model runs are useful, but many acoustic questions
are inherently distributional. The package therefore includes simulation
tools for repeated realizations, parameter batching, and ensemble
summaries.

![Batching and realizations
schematic](simulation-batching-schematic.svg)

Batching and realizations schematic

This page is about how to turn a single target description into a
structured family of model runs. The core idea is that the object
remains fixed as the bookkeeping container, while selected inputs are
varied across realizations or across a parameter grid.

## Main entry point

The main current simulation interface is
[`simulate_ts()`](https://brandynlucca.github.io/acousticTS/reference/simulate_ts.md).
It supports repeated draws of user-specified parameters, optional
batching over one or more variables, and optional parallel execution.

In other words,
[`simulate_ts()`](https://brandynlucca.github.io/acousticTS/reference/simulate_ts.md)
is a workflow wrapper around repeated
[`target_strength()`](https://brandynlucca.github.io/acousticTS/reference/target_strength.md)
calls. It does not replace the usual model layer. It automates the
process of rebuilding a set of parameters, optionally reworking the
object, running the model, and collecting the outputs into a set of tidy
result tables.

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

res <- simulate_ts(
  object = obj,
  frequency = seq(38e3, 120e3, by = 6e3),
  model = "DWBA",
  n_realizations = 5,
  parameters = list(
    theta_body = function() runif(1, min = 0.5 * pi, max = pi),
    density_body = 1045
  ),
  parallel = FALSE
)
```

    ## ====================================
    ## Scatterer-class: FLS 
    ## Model(s): DWBA 
    ## Simulated parameters: theta_body, density_body 
    ## Total simulation realizations: 5 
    ## Parallelize TS calculations: FALSE 
    ## ====================================
    ## ====================================
    ## Preparing sequential simulations
    ## ====================================
    ## Simulations complete!
    ## ====================================

The return value is a named list with one data frame per model.

## Parameter specification

The scratch simulation notes are worth keeping in mind here because they
clarify the four main parameter modes supported by
[`simulate_ts()`](https://brandynlucca.github.io/acousticTS/reference/simulate_ts.md):

1.  single fixed values,
2.  explicit vectors,
3.  generating functions that return one draw at a time,
4.  structured
    [`reforge()`](https://brandynlucca.github.io/acousticTS/reference/reforge.md)
    inputs such as named target-dimension vectors.

When `batch_by` is supplied, the selected parameters define a Cartesian
grid and each grid cell is then repeated for `n_realizations`. That
distinction between batch cells and within-cell draws is the core
conceptual point of the interface.

Those four modes answer different questions:

1.  fixed values are for constants that should be repeated in every
    realization,
2.  explicit vectors are for a known set of realization-level values,
3.  functions are for stochastic draws generated on demand,
4.  structured values are for explicit geometry targets that should be
    preserved as one unit.

The distinction matters because it separates deterministic sweeps from
uncertainty propagation.

![Three parameter-specification modes in simulate_ts(): fixed constants,
explicit vectors, and per-realization generating
functions.](simulation-parameter-modes.svg)

Three parameter-specification modes in
[`simulate_ts()`](https://brandynlucca.github.io/acousticTS/reference/simulate_ts.md):
fixed constants, explicit vectors, and per-realization generating
functions.

Fixed values are copied into every run unchanged. Explicit vectors
supply a known sequence of values. Generating functions do not provide a
pre-listed table at all. They create a new draw whenever a realization
is built. Structured values are preserved as list-columns so named
geometry targets remain intact. That difference is what separates a
design grid from on-the-fly uncertainty propagation and from explicit
morphology rebuilding.

``` r
parameters <- list(
  theta_body = seq(0.5 * pi, pi, length.out = 5),
  density_body = function() rnorm(1, mean = 1045, sd = 5),
  sound_speed_body = 1520,
  body_target = c(length = 0.045)
)
```

In this example, `theta_body` is an explicit set of values,
`density_body` is random, `sound_speed_body` is held fixed, and
`body_target` is a structured
[`reforge()`](https://brandynlucca.github.io/acousticTS/reference/reforge.md)
input.

For fluid-like scatterers, convenience aliases such as `length_body` are
also supported. These are normalized internally onto the same
[`reforge()`](https://brandynlucca.github.io/acousticTS/reference/reforge.md)
pathway as `body_target`, while preserving the original simulation
column in the returned output.

``` r
res_length_alias <- simulate_ts(
  object = obj,
  frequency = 120e3,
  model = "DWBA",
  n_realizations = 4,
  parameters = list(
    length_body = function() runif(1, min = 0.04, max = 0.06),
    theta_body = function() runif(1, min = 0.5 * pi, max = pi)
  ),
  parallel = FALSE
)
```

    ## ====================================
    ## Scatterer-class: FLS 
    ## Model(s): DWBA 
    ## Simulated parameters: length_body, theta_body 
    ## Total simulation realizations: 4 
    ## Parallelize TS calculations: FALSE 
    ## ====================================
    ## ====================================
    ## Preparing sequential simulations
    ## ====================================
    ## Simulations complete!
    ## ====================================

``` r
head(res_length_alias$DWBA[, c("length_body", "theta_body", "TS")])
```

    ##   length_body theta_body        TS
    ## 1  0.04932787   2.784270 -67.69139
    ## 2  0.04995555   2.944616 -67.80459
    ## 3  0.04579534   1.845592 -67.42663
    ## 4  0.05465764   1.624582 -69.50127

## Batching versus realizations

The most important conceptual distinction in the interface is the
difference between a batched parameter and a realization-level
parameter.

1.  `batch_by` creates a grid of settings that should each be kept
    distinct in the output.
2.  `n_realizations` controls how many repeated runs are performed
    within each grid cell.

Suppose `theta_body` and `density_body` are both placed in `batch_by`.
Then the simulation does not treat them as uncertain draws around one
target. It treats them as a design grid whose combinations are all run
separately.

``` r
res_grid <- simulate_ts(
  object = obj,
  frequency = seq(38e3, 120e3, by = 6e3),
  model = c("DWBA", "SDWBA"),
  n_realizations = 3,
  parameters = list(
    theta_body = seq(0.5 * pi, pi, length.out = 4),
    density_body = c(1035, 1045, 1055)
  ),
  batch_by = c("theta_body", "density_body"),
  parallel = FALSE
)
```

    ## ====================================
    ## Scatterer-class: FLS 
    ## Model(s): DWBA, SDWBA 
    ## Batching parameter(s): theta_body, density_body 
    ## Simulated parameters: theta_body, density_body 
    ## Total simulation realizations: 36 
    ## Parallelize TS calculations: FALSE 
    ## ====================================
    ## ====================================
    ## Preparing sequential simulations
    ## ====================================
    ## Simulations complete!
    ## ====================================

Here the output contains one row structure per model, per frequency, per
batch cell, and per realization. That is why batched jobs can grow
quickly.

![How batch_by creates Cartesian design cells and n_realizations nests
repeated runs inside each cell.](simulation-batch-cells.svg)

How `batch_by` creates Cartesian design cells and `n_realizations` nests
repeated runs inside each cell.

This is the core simulation structure in the current interface. Batched
parameters define the outer design. Realizations live inside that
design. If a user loses track of that nesting, it becomes very easy to
misread deterministic design cells as random draws or random draws as
deterministic sweep levels.

## Parameters that modify the object

Not every simulated parameter is a simple scalar passed into
[`target_strength()`](https://brandynlucca.github.io/acousticTS/reference/target_strength.md).
Some parameters alter the working object itself. This is especially
important for reshaping workflows.

The current simulation code supports two broad cases:

1.  parameters matching accepted
    [`reforge()`](https://brandynlucca.github.io/acousticTS/reference/reforge.md)
    inputs can trigger per-realization resizing or reparameterization,
2.  parameters matching object-component fields can be inserted into the
    working object before the model call.

This means
[`simulate_ts()`](https://brandynlucca.github.io/acousticTS/reference/simulate_ts.md)
can be used not just for orientation or density perturbations, but also
for controlled morphology studies.

``` r
body_shape <- arbitrary(
  x_body = c(0.00, 0.02, 0.05, 0.08),
  zU_body = c(0.001, 0.004, 0.004, 0.001),
  zL_body = c(-0.001, -0.004, -0.004, -0.001)
)

bladder_shape <- arbitrary(
  x_bladder = c(0.02, 0.04, 0.06),
  zU_bladder = c(0.0012, 0.0018, 0.0012),
  zL_bladder = c(-0.0012, -0.0018, -0.0012)
)

obj_sbf <- sbf_generate(
  body_shape = body_shape,
  bladder_shape = bladder_shape,
  density_body = 1045,
  sound_speed_body = 1520,
  density_bladder = 1.2,
  sound_speed_bladder = 343
)

res_shape <- simulate_ts(
  object = obj_sbf,
  frequency = seq(38e3, 120e3, by = 6e3),
  model = "DWBA",
  n_realizations = 2,
  parameters = list(
    body_scale = c(0.9, 1.1),
    swimbladder_inflation_factor = c(0.8, 1.2)
  ),
  batch_by = c("body_scale", "swimbladder_inflation_factor"),
  parallel = FALSE
)
```

    ## ====================================
    ## Scatterer-class: SBF 
    ## Model(s): DWBA 
    ## Batching parameter(s): body_scale, swimbladder_inflation_factor 
    ## Simulated parameters: body_scale, swimbladder_inflation_factor 
    ## Total simulation realizations: 8 
    ## Parallelize TS calculations: FALSE 
    ## ====================================
    ## ====================================
    ## Preparing sequential simulations
    ## ====================================
    ## Simulations complete!
    ## ====================================

This is one of the cleanest ways to keep a morphology-sensitivity
analysis tied to the same baseline target.

For `FLS` objects, the same idea can be expressed more explicitly
through `body_target` or through the convenience alias `length_body`.

``` r
res_fls_reforge <- simulate_ts(
  object = obj,
  frequency = 120e3,
  model = "DWBA",
  n_realizations = 2,
  parameters = list(
    body_target = list(
      c(length = 0.04, radius = 0.0025),
      c(length = 0.06, radius = 0.0035)
    ),
    isometric_body = FALSE,
    theta_body = function() runif(1, min = 0.5 * pi, max = pi)
  ),
  batch_by = "body_target",
  parallel = FALSE
)
```

    ## ====================================
    ## Scatterer-class: FLS 
    ## Model(s): DWBA 
    ## Batching parameter(s): body_target 
    ## Simulated parameters: body_target, isometric_body, theta_body 
    ## Total simulation realizations: 4 
    ## Parallelize TS calculations: FALSE 
    ## ====================================
    ## ====================================
    ## Preparing sequential simulations
    ## ====================================
    ## Simulations complete!
    ## ====================================

``` r
head(res_fls_reforge$DWBA[, c("body_target", "theta_body", "TS")])
```

    ##    body_target theta_body        TS
    ## 1 0.04, 0.0025   2.490944 -68.15085
    ## 2 0.04, 0.0025   2.235577 -68.15085
    ## 3 0.06, 0.0035   2.730156 -72.25574
    ## 4 0.06, 0.0035   1.845710 -72.25574

When batching across structured geometry targets, wrap the candidates in
a list so each target vector is treated as one batch value.

![Object-modifying simulation workflow, where baseline geometry is
reforged into a family of derived target states before model
execution.](simulation-object-reforge.svg)

Object-modifying simulation workflow, where baseline geometry is
reforged into a family of derived target states before model execution.

Morphology-affecting parameters do not behave like ordinary scalar
metadata. They change the working target itself before the model is run.
The resulting simulations should therefore be interpreted as distinct
geometric states derived from one baseline object, not as repeated
evaluations of a single unchanged shape.

There is one especially important point for curvature-aware workflows:
[`simulate_ts()`](https://brandynlucca.github.io/acousticTS/reference/simulate_ts.md)
does not automatically call
[`brake()`](https://brandynlucca.github.io/acousticTS/reference/brake.md)
for the modern `DWBA` or `SDWBA` paths. If the baseline object is
straight, then `length_body` or `body_target` simply resize that
straight object. If the baseline object is already a bent `FLS`, then
those same parameters resize the existing bent geometry directly.

For bent `FLS` objects, `length_body` and `body_target["length"]` are
interpreted as the target bent centerline arc length, not as the
flattened projected `x` span. That makes the resize semantics consistent
with the object that is actually being modeled.

``` r
obj_bent <- brake(obj, radius_curvature = 5, mode = "ratio")

res_bent <- simulate_ts(
  object = obj_bent,
  frequency = 120e3,
  model = "DWBA",
  n_realizations = 2,
  parameters = list(
    length_body = c(0.04, 0.05),
    theta_body = 0.5 * pi
  ),
  batch_by = "length_body",
  parallel = FALSE
)
```

    ## ====================================
    ## Scatterer-class: FLS 
    ## Model(s): DWBA 
    ## Batching parameter(s): length_body 
    ## Simulated parameters: length_body, theta_body 
    ## Total simulation realizations: 4 
    ## Parallelize TS calculations: FALSE 
    ## ====================================
    ## ====================================
    ## Preparing sequential simulations
    ## ====================================
    ## Simulations complete!
    ## ====================================

``` r
head(res_bent$DWBA[, c("length_body", "theta_body", "TS")])
```

    ##   length_body theta_body        TS
    ## 1        0.04   1.570796 -68.48278
    ## 2        0.04   1.570796 -68.48278
    ## 3        0.05   1.570796 -68.35455
    ## 4        0.05   1.570796 -68.35455

That example is useful because it separates two different ideas cleanly:

1.  the bent state comes from the baseline object passed into
    [`simulate_ts()`](https://brandynlucca.github.io/acousticTS/reference/simulate_ts.md),
2.  the resizing comes from `length_body` or `body_target`.

If the scientific workflow instead needs “resize a straight object and
then bend it differently for each realization,” that should currently be
constructed outside
[`simulate_ts()`](https://brandynlucca.github.io/acousticTS/reference/simulate_ts.md)
or through the deprecated `*_curved` model interfaces that still wire
curvature in internally.

## When simulation is useful

Simulation or batched evaluation is most useful when:

1.  orientation is uncertain,
2.  material properties are better treated as distributions than as
    fixed values,
3.  multiple body lengths or frequencies must be screened
    systematically,
4.  model sensitivity is of greater interest than a single nominal
    prediction.

In practice, this means simulation is less about “making the model
stochastic” and more about mapping the consequences of uncertainty,
variability, or design decisions.

## Interpreting the output

The returned list is model-centered. Each element is a data frame that
carries the realization index, any batched parameters, and the model
outputs. That organization is deliberate: it makes downstream plotting
or aggregation straightforward.

``` r
names(res_grid)
```

    ## [1] "DWBA"  "SDWBA"

``` r
head(res_grid$DWBA)
```

    ##   model realization theta_body_idx density_body_idx theta_body density_body
    ## 1  DWBA           1              1                1   1.570796         1035
    ## 2  DWBA           1              1                1   1.570796         1035
    ## 3  DWBA           1              1                1   1.570796         1035
    ## 4  DWBA           1              1                1   1.570796         1035
    ## 5  DWBA           1              1                1   1.570796         1035
    ## 6  DWBA           1              1                1   1.570796         1035
    ##   frequency        ka                        f_bs     sigma_bs        TS
    ## 1     38000 0.4775221 -0.0001558136-7.493349e-20i 2.427787e-08 -76.14790
    ## 2     44000 0.5529203 -0.0002007945-1.118128e-19i 4.031842e-08 -73.94497
    ## 3     50000 0.6283185 -0.0002476169-1.566886e-19i 6.131414e-08 -72.12439
    ## 4     56000 0.7037168 -0.0002946149-2.087997e-19i 8.679795e-08 -70.61491
    ## 5     62000 0.7791150 -0.0003400723-2.668395e-19i 1.156492e-07 -69.36857
    ## 6     68000 0.8545132 -0.0003822717-3.289791e-19i 1.461317e-07 -68.35256

``` r
head(res_grid$SDWBA)
```

    ##   model realization theta_body_idx density_body_idx theta_body density_body
    ## 1 SDWBA           1              1                1   1.570796         1035
    ## 2 SDWBA           1              1                1   1.570796         1035
    ## 3 SDWBA           1              1                1   1.570796         1035
    ## 4 SDWBA           1              1                1   1.570796         1035
    ## 5 SDWBA           1              1                1   1.570796         1035
    ## 6 SDWBA           1              1                1   1.570796         1035
    ##   frequency                        f_bs     sigma_bs        TS     TS_sd
    ## 1     38000 -9.447098e-05-7.072660e-07i 9.925839e-09 -80.03233 -85.54242
    ## 2     44000 -1.214249e-04-3.378969e-06i 1.590348e-08 -77.98508 -83.24166
    ## 3     50000 -1.522909e-04-6.089763e-06i 2.509366e-08 -76.00436 -81.99001
    ## 4     56000 -1.771822e-04+1.795488e-06i 3.431663e-08 -74.64495 -80.07396
    ## 5     62000 -2.039148e-04-3.418298e-06i 4.576683e-08 -73.39449 -78.60106
    ## 6     68000 -2.284664e-04+7.397460e-06i 5.717119e-08 -72.42823 -77.59347

Common downstream summaries include:

1.  mean and quantiles of `TS` by frequency,
2.  variance decomposition by orientation or material parameter,
3.  pairwise model differences across matched grid cells,
4.  screening plots to identify nonlinear regions before deeper
    analysis.

## Parallel execution on Windows

On Windows, parallel execution uses PSOCK workers rather than forked
processes. That has two practical consequences:

1.  startup overhead is more noticeable for small jobs,
2.  larger simulations benefit more clearly from parallelization than
    tiny exploratory runs.

For this reason, `parallel = FALSE` is often the right default while
debugging a workflow, even if `parallel = TRUE` is the right choice for
production runs.

``` r
res_parallel <- simulate_ts(
  object = obj,
  frequency = seq(38e3, 120e3, by = 2e3),
  model = "DWBA",
  n_realizations = 100,
  parameters = list(theta_body = function() runif(1, 0.5 * pi, pi)),
  parallel = TRUE,
  n_cores = 4
)
```

![Execution layout for serial and parallel simulation runs, emphasizing
the same logical job table with different
scheduling.](simulation-execution-layout.svg)

Execution layout for serial and parallel simulation runs, emphasizing
the same logical job table with different scheduling.

Parallelization changes scheduling, not interpretation. The same batched
cells and realization-level jobs still exist. Parallel execution simply
distributes them across workers. That is why debugging is often easier
in serial mode first, even when the final production run is parallel.

## Practical cautions

The simulation workflow can become expensive quickly, especially when
several parameters are batched simultaneously. The main practical risks
are:

- combinatorial growth in the evaluation grid,
- unnecessary parallel overhead for small jobs,
- confusion between realization-level randomness and batch-level
  parameter sweeps.

Two additional cautions are worth keeping in mind:

1.  if a parameter sequence is meant to represent uncertainty, placing
    it in `batch_by` changes the interpretation from random variation to
    deterministic design cells,
2.  if geometry-affecting parameters are varied, the results should be
    interpreted as new target states rather than as mere perturbations
    of one fixed shape.

For that reason, this page should be read together with [FAQ and
troubleshooting](https://brandynlucca.github.io/acousticTS/articles/faq-troubleshooting/faq-troubleshooting.md).

On Windows, parallel execution uses PSOCK workers rather than forked
processes, so larger simulations benefit from being planned deliberately
instead of relying on a maximal batch grid by default.

## Relationship to `anneal()`

The current documentation and source code indicate that
[`simulate_ts()`](https://brandynlucca.github.io/acousticTS/reference/simulate_ts.md)
is expected to give way to `anneal()` in future releases. The important
point for current users is that the workflow logic on this page remains
useful even if the interface evolves.

The stable ideas are:

1.  keep one object as the baseline target description,
2.  separate batch design from within-cell variation,
3.  make model comparison operate on matched simulation cells,
4.  summarize in linear or logarithmic units according to the scientific
    question.

## Relationship to future interfaces

The current simulation interface is already useful, but the package
documentation indicates that this part of the workflow may continue to
evolve. This article should therefore be read as a workflow guide rather
than as a claim that the current interface is final.

## References

Demer, David A., and Stéphane G. Conti. 2005. “New Target-Strength Model
Indicates More Krill in the Southern Ocean.” *ICES Journal of Marine
Science* 62 (1): 25–32. <https://doi.org/10.1016/j.icesjms.2004.07.027>.

Jech, J. Michael, John K. Horne, Dezhang Chu, David A. Demer, David T.
I. Francis, Natalia Gorska, Benjamin Jones, et al. 2015. “Comparisons
Among Ten Models of Acoustic Backscattering Used in Aquatic Ecosystem
Research.” *The Journal of the Acoustical Society of America* 138 (6):
3742–64. <https://doi.org/10.1121/1.4937607>.
