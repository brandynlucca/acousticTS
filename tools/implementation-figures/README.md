# Implementation figure builders

This directory contains the reproducible, package-facing builders for figures
and compact tables used by the public implementation articles.

The builders are development tooling. They are excluded from the package
bundle, and neither package installation nor a routine pkgdown build executes
them. Generated outputs are committed so the documentation can render from a
fresh checkout without recomputing scientific results.

## Layout

- `manifest.csv` declares each family, builder, package-code inputs, and
  expected outputs.
- `run_all.R` selects families, launches one sequential R subprocess per
  family, validates declared outputs, and records provenance.
- `run_family.R` is the subprocess worker used by `run_all.R`.
- Family directories contain the supported figure builders.
- `helpers/` contains shared plotting and repository helpers.
- `data/` contains compact, committed numeric inputs consumed by the builders
  or implementation articles.

Heavy external-solver validation, exploratory analyses, machine-specific
inputs, and scratch-table generators do not belong in this directory.

## Running builders

From the repository root, rebuild every package-facing family with:

```r
source("tools/implementation-figures/run_all.R")
```

To rebuild selected families:

```r
Sys.setenv(ACOUSTICTS_IMPL_FAMILIES = "dwba,sdwba")
source("tools/implementation-figures/run_all.R")
```

An empty family selection or `ACOUSTICTS_IMPL_FAMILIES=all` selects every
family. Only the `light` package-facing profile is supported.

The default provenance report is written to:

```text
.tmp/implementation-figures/provenance.csv
```

It records output sizes and MD5 hashes, elapsed family time, and the package
commit.

## Continuous integration

The pull-request workflow uses the manifest's `inputs` column to select
families affected by model-specific changes. An unmapped change under `R/` or
`src/` conservatively rebuilds every family. Missing outputs or generated drift
fail the relevant pull-request job. Manual runs may keep drift informational by
leaving `fail_on_drift` disabled.
