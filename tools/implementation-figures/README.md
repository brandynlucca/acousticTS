# Implementation Figure Builders

This directory is the permanent home for scripts that regenerate the static
figure and comparison assets used by the `*-implementation.Rmd` vignettes.

## Layout

- `run_all.R`
  - entry point for the lightweight implementation rebuild pass
- `manifest.csv`
  - flat index of the current implementation builders and their outputs
- `bbfm/`
  - `BBFM`-specific figure builders
- `calibration/`
  - calibration-sphere comparison and diagnostic builders
- `helpers/`
  - shared BEMPP and mesh helpers used by the TMM validation figures
- `tmm/`
  - retained-operator and BEM comparison figure builders for `TMM`
- top-level `*.R`
  - family-specific comparison/data builders used by implementation pages

## Associated data

Implementation-page CSV inputs that should not live in `scratch/` are stored in:

- `tools/implementation-figures/data/`

Implementation pages should read from that directory rather than from
`scratch/`.

## Intended workflow

Run the default lightweight figure/data refresh from the repository root with:

```r
source("tools/implementation-figures/run_all.R")
```

To rebuild only selected families, set `ACOUSTICTS_IMPL_FAMILIES` to a
comma-separated list before sourcing the runner:

```r
Sys.setenv(ACOUSTICTS_IMPL_FAMILIES = "dwba,sdwba")
source("tools/implementation-figures/run_all.R")
```

To include the heavier BEM-backed `TMM` validation builders as well, set:

```r
Sys.setenv(ACOUSTICTS_IMPL_PROFILE = "all")
source("tools/implementation-figures/run_all.R")
```

The default `light` profile is the one intended for routine figure refreshes
and CI checks. The `all` profile additionally runs the heavier BEM-backed TMM
validation figure builders.

## CI recommendation

Yes, this is a good fit for GitHub Actions, but I would not run the heaviest
validation scripts on every pull request by default.

The clean split is:

- lightweight PR check:
  - resolve the affected implementation families from the PR diff
  - run `tools/implementation-figures/run_all.R` only for those families
  - fail if tracked figure/data assets change
- optional/manual validation job:
  - set `ACOUSTICTS_IMPL_PROFILE=all`
  - optionally set `ACOUSTICTS_IMPL_FAMILIES`
  - run the same runner
  - use `workflow_dispatch`, a label, or a path filter

That keeps routine PR feedback fast while still making the expensive figure
rebuild path reproducible and visible.

## Rule going forward

Any new static asset added to an implementation vignette should come with:

1. a rebuild script in this directory,
2. any committed CSV inputs under `tools/implementation-figures/data/`, and
3. an entry in `manifest.csv`.
