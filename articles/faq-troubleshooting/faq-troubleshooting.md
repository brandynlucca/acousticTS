# Developer Troubleshooting

## Start from a clean diagnosis

First determine whether R is using the source checkout, an installed
release, or a previously loaded namespace. Many apparent code failures
are version mismatches between those states.

From the repository root, load the working tree and record the session:

``` r

devtools::load_all(".")
packageVersion("acousticTS")
sessionInfo()
```

Reproduce the failure in a fresh R process before changing code. A
result that only fails after other scripts have run usually indicates
leaked options, registry state, attached packages, working-directory
assumptions, or modified objects.

## Installation and compiled code

`acousticTS` contains C++17 and Fortran code. Installing from source
requires a matching compiler toolchain:

- **Windows:** install the Rtools version for the installed R release
  and check that its compiler tools are available to R.
- **macOS:** install the Xcode Command Line Tools and the GNU Fortran
  compiler recommended for the installed R release.
- **Linux:** install `make`, a C++17 compiler, `gfortran`, and the
  development libraries required by the R installation.

Check the toolchain before interpreting a long linker error:

``` r

pkgbuild::check_build_tools(debug = TRUE)
```

Messages mentioning `CXX17`, `gfortran`, `quadmath`, BLAS, undefined
Fortran symbols, or an unavailable linker usually indicate a build
environment problem. Run a source install in a clean process to preserve
the full compiler output:

``` powershell
R CMD INSTALL .
```

Do not edit `R/RcppExports.R` or `src/RcppExports.cpp` by hand. Run
[`Rcpp::compileAttributes()`](https://rdrr.io/pkg/Rcpp/man/compileAttributes.html)
only when an exported C++ interface changes, then review both generated
files. Similarly, use `devtools::document()` after roxygen changes and
review updates to `NAMESPACE` and `man/`.

## A registered model is missing

Inspect the registry before debugging dispatch:

``` r

available_models()
```

For a session model, confirm that
[`register_model()`](https://brandynlucca.github.io/acousticTS/reference/register_model.md)
ran after `acousticTS` was loaded and that its canonical name or alias
appears in the table. Function names do not have to match the model
name. The registry entry itself must point to resolvable initializer and
solver functions.

Persistent registrations require package-qualified references such as
`"myPackage::initialize_tsl"`. Raw function objects cannot be restored
in a new session. Name conflicts are also rejected, including aliases
that collide with built-in or user registrations.

Remove a single user entry with:

``` r

unregister_model("tsl")
```

To test without any session registrations, clear them while retaining
the on-disk configuration:

``` r

reset_model_registry(remove_persisted = FALSE)
```

Use `remove_persisted = TRUE` only when the saved user configuration
should also be deleted. See [Creating
Models](https://brandynlucca.github.io/acousticTS/articles/creating-models-from-scratch/creating-models-from-scratch.md)
for the registry contract.

## Initializer and solver failures

Separate dispatch from calculation. Call the initializer directly with a
small canonical object, inspect its slots, then call the solver:

``` r

prepared <- tsl_initialize(
  object = target,
  frequency = c(38e3, 70e3)
)

extract(prepared, "model_parameters")$TSL
extract(prepared, "model")$TSL

solved <- TSL(prepared)
extract(solved, "model")$TSL
```

Check these contracts:

- both functions return the updated `Scatterer`,
- the initializer and solver use the registry’s same slot name,
- model-specific arguments appear in the initializer formals,
- output rows align with the documented frequency or angle grid,
- `TS`, `sigma_bs`, and `f_bs` agree with their definitions, and
- unsupported classes, boundaries, and parameter values fail explicitly.

If direct calls work but
[`target_strength()`](https://brandynlucca.github.io/acousticTS/reference/target_strength.md)
fails, inspect the registry entry, the requested alias, and
`model_args`. If deterministic runs work but
[`simulate_ts()`](https://brandynlucca.github.io/acousticTS/reference/simulate_ts.md)
fails, test one realization without parallel execution. This usually
exposes an argument-shape or worker-export problem more clearly.

## Tests pass alone but fail in the suite

Run the smallest relevant test first:

``` r

testthat::test_file("tests/testthat/test-model_registry.R")
```

Then run the package tests in a fresh process with `devtools::test()`.
Tests that modify the model registry, options, environment variables,
random-number state, or working directory must restore the original
state with [`on.exit()`](https://rdrr.io/r/base/on.exit.html). Avoid
assertions that depend on test order or on objects created by another
file.

For numerical failures, compare the first differing intermediate
quantity rather than only the final `TS`. Record the geometry, boundary,
material values, frequency, angle, numerical controls, and comparison
domain. The procedure in [Validation
Methods](https://brandynlucca.github.io/acousticTS/articles/validation-benchmarks/validation-benchmarks.md)
separates numerical verification from regression testing.

## Pkgdown and vignette failures

The site builds from the source package in one R process in continuous
integration. Reproduce that mode locally after loading the checkout:

``` r

devtools::load_all(".")

pkgdown::build_site_github_pages(
  new_process = FALSE,
  install = FALSE
)
```

For one page, use its pkgdown article identifier:

``` r

pkgdown::build_article(
  "creating-models-from-scratch/creating-models-from-scratch",
  lazy = FALSE,
  new_process = FALSE
)
```

An error that occurs only with `new_process = TRUE` often means the
child process loaded an installed copy rather than uninstalled source
changes. Install the working tree first, or use the source-loading mode
above to match this repository’s workflow.

When a page fails, also check:

- chunk working-directory assumptions,
- paths to local images, SVG files, CSS hooks, and article links,
- use of unexported helpers from an older installed namespace,
- duplicate chunk labels, and
- code that performs expensive scientific calculations during rendering.

Implementation figures are precomputed. Ordinary vignette and pkgdown
builds should consume committed outputs rather than execute the builders
under `tools/implementation-figures/`. To investigate generated drift,
run the family through `tools/implementation-figures/run_all.R` and
inspect the provenance file under `.tmp/implementation-figures/`.

## Reporting a reproducible problem

A useful issue contains:

1.  a minimal object and function call,
2.  the complete error and warnings,
3.  [`sessionInfo()`](https://rdrr.io/r/utils/sessionInfo.html) and the
    package version or commit,
4.  the operating system and compiler versions for installation
    failures,
5.  the model, boundary, geometry, frequency, and numerical options, and
6.  whether the failure reproduces in a clean R process.

For unexpected physical output, first confirm units, geometry, material
properties, boundary condition, and model domain. Then report the
smallest case that separates the unexpected result from the expected
reference. Questions about model suitability belong with [Choosing a
Model](https://brandynlucca.github.io/acousticTS/articles/model-selection/model-selection.md),
while phase and component-addition questions belong with [Combining
Scattering
Components](https://brandynlucca.github.io/acousticTS/articles/combining-components/combining-components.md).
