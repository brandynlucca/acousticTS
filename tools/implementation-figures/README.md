# Implementation Figure Builders

This directory is now limited to lightweight, package-facing figure rebuilds
and small committed inputs used directly by the public `*-implementation.Rmd`
vignettes.

## Package-side scope

What belongs here:

- vignette-facing figure rebuild scripts
- tiny committed tables consumed directly by vignettes
- lightweight build helpers needed during package development

What does not belong here anymore:

- BEM/FEM solver runners
- external-validation registries
- mesh-generation tooling for external solvers
- heavy compare tables, residual grids, and exploratory validation labs
- raw regeneration inputs that only exist to rebuild external validation
  compare suites

Those heavier workflows now live in:

- [acousticTSValidation](C:/Users/Brandyn/Desktop/acousticTSValidation)

## Current layout

- `run_all.R`
  - lightweight package-side figure rebuild entry point
- `manifest.csv`
  - index of the package-side implementation builders that still remain here
- family directories such as `bbfm/`, `fcms/`, `psms/`, `tmm/`
  - package-facing implementation figure builders
- `helpers/`
  - lightweight shared helpers used by the remaining package-side builders
- `data/`
  - small committed numeric inputs still needed by package-facing builders

## What moved out

The following validation infrastructure has been externalized:

- canonical `bem/` runners and their data roots
- canonical `fem/` runners and their data roots
- heavy BEM/FEM Python helper workflows
- family-specific external-solver validation runners mirrored into the
  validation repo
- large legacy `tmm/` validation artifacts beyond the small pure-R package-side
  outputs
- raw compare/regeneration scripts and source timing/reference tables for
  `bcms`, `ecms`, `pcdwba`, `vesm`, and calibration diagnostics

What intentionally remains here:

- package-facing implementation figure builders
- committed compare outputs that are read directly by those builders or by the
  `*-implementation.Rmd` vignettes
- the small pure-R `tmm` continuation/exact-validation inputs still needed by
  the vignette-facing TMM builders

See:

- [MIGRATION.md](C:/Users/Brandyn/Desktop/acousticTS/tools/implementation-figures/MIGRATION.md)

## Package workflow

Run the lightweight package-side figure refresh from the repository root with:

```r
source("tools/implementation-figures/run_all.R")
```

To rebuild only selected families, set `ACOUSTICTS_IMPL_FAMILIES` before
sourcing the runner:

```r
Sys.setenv(ACOUSTICTS_IMPL_FAMILIES = "dwba,sdwba")
source("tools/implementation-figures/run_all.R")
```

Heavy validation profiles are no longer supported from this repo-level runner.
Use the external validation project for BEM/FEM/external-solver workflows.
