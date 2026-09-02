# Contributing to acousticTS

Thank you for considering a contribution to `acousticTS`.

## Before opening an issue

Search the existing issues and consult the [package
documentation](https://brandynlucca.github.io/acousticTS/). Use the
appropriate issue form and include a minimal reproducible example. For
numerical problems, identify the model, boundary condition, geometry,
parameters, expected result, and any independent reference calculation.

## Development setup

Install a current R release and the platform compiler toolchain required
for packages containing C++, C, and Fortran source code. On Windows,
install the matching version of Rtools.

Install the package dependencies with:

``` r

pak::pak()
```

## Pull requests

Create a focused branch and keep each pull request limited to one
coherent change. Before submitting it:

1.  Add or update tests for behavior changes.
2.  Add a `NEWS.md` entry for user-visible changes.
3.  Run `devtools::document()` when roxygen comments change.
4.  Run `devtools::test()`.
5.  Run `devtools::check()`.

Do not commit compiled objects, temporary plots, local check
directories, or generated files that are already excluded by
`.gitignore` or `.Rbuildignore`.

## R documentation and style

Keep R source lines within 80 characters where practical. Use `TRUE` and
`FALSE`, never the abbreviations `T` and `F`.

Every function must have a meaningful roxygen title and documentation.
Internal functions must use `@noRd`. Exported documentation and examples
must not refer to internal functions or internal implementation helpers.

Examples must be deterministic, quick, and safe during `R CMD check`.
Examples that intentionally require unsupported precision, substantial
computation, or interactive plotting should be enclosed in `\dontrun{}`
where appropriate.

## Compiled code

Changes under `src/` must compile without new diagnostics on supported
platforms. Preserve upstream copyright and licensing notices. Avoid
adding platform-specific compiler flags to package `Makevars` files.

Numerical changes should include regression tests and clearly state the
reference result, tolerance, and applicable parameter range.

## Reporting security concerns

Do not include credentials, tokens, private datasets, or other sensitive
information in an issue or pull request. Contact the package maintainer
privately when public disclosure would create a security risk.
