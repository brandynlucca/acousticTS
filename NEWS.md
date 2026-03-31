# v.2026.03.0

## What's Changed

### ⚠️ Breaking Changes
* The `DCM` model has been fully deprecated due to its 1) misleading name and 2) limited scope. It has been superseded by `TRCM`, which uses an adjusted implementation that allows for both straight and bent cylinder geometries. 
* The `DWBA_curved` and `SDWBA_curved` models have been deprecated. It is now expected that users explicitly apply `brake` to their shapes (i.e., adds curvature) prior to using the standard `DWBA` and `SDWBA` models. 
* The `MSS_anderson` model has been deprecated and replaced with `SPHMS`. 
* The `high_pass_stanton` model has been deprecated and replaced with `HPA`.
* The `MSS_goodman_stern` model has been deprecated and replaced with `ESSMS`.
* Additional geometry and class guardrails have been added across the package. Some legacy workflows may now error if the supplied shape, scatterer class, or boundary assumptions are incompatible with the selected model.

### New features and enhancements :sparkles:

#### Software :desktop_computer:

##### New & updated models :pencil:
* `SPHMS`: Modal-series solution for canonical spherical targets across multiple supported boundary conditions.
* `FCMS`: Modal-series solution for finite cylinders across multiple supported boundary conditions.
* `PSMS`: Modal-series solution for prolate spheroids across multiple supported boundary conditions.
* `SOEMS` / `calibration`: The solid elastic sphere workflow has been broadened and documented as a general solid-elastic-sphere family rather than only a calibration-sphere workflow. `calibration` remains available as a compatibility alias.
* `ESSMS`: Modal-series solution for elastic-shelled spheres.
* `BCMS`: Modal-series solution for finite, bent cylinders.
* `ECMS`: Modal-series solution for solid elastic finite cylinders.
* `KRM`: Now incorporates the breathing-mode modal-series calculation for `ka < 0.15` and has been updated to better align with alternative body sound-speed assumptions.
* `HPA`: Expanded usable geometries and now includes the Johnson (1977) formulation.
* `TRCM`: The two-ray cylinder model now supports both straight and bent cylinders.
* `PCDWBA`: Phase-compensated distorted-wave Born approximation for bent weakly scattering targets.
* `BBFM`: Composite flesh-plus-backbone model for swimbladder-less fish-like targets, treating the backbone as a solid elastic scatterer.
* `VESMS` (`VESM` family): Viscous-elastic layered-sphere model for gas-core, shell, and viscous-layer targets.
* `TMM`: Transition-matrix framework for monostatic and angle-dependent scattering of spheres, prolate spheroids, oblate spheroids, and guarded finite-cylinder branches.

##### Numerical methods  :calling:
* Cylindrical and spherical Bessel functions of the first, second, and third kinds were overhauled and expanded, including broader support for orders, purely imaginary arguments, and public $k^{\text{th}}$ derivative helpers such as `jcdk()`, `hcdk()`, `jsdk()`, and `ysdk()`.
* Full public support for prolate spheroidal angular and radial wave functions via `Smn()` and `Rmn()`, including derivatives.
* Public Legendre support now includes both the first and second kinds, including `Pn()`, `Pndk()`, `Qn()`, and `Qndk()`.
* New numerical helpers such as `gauss_legendre()` and related modal/scattering utilities were added to support higher-level model implementations.

##### Workflow, geometry, and extensibility :gear:
* Shape- and scatterer-generation workflows were substantially streamlined and reorganized.
* `simulate_ts()` provides parameter-sweep and repeated-realization workflows for TS simulation.
* `canonicalize_shape()` converts arbitrary shapes into canonical surrogates such as spheres, prolate spheroids, oblate spheroids, and cylinders.
* `reforge()` supports resizing and reparameterizing supported scatterer objects.
* `extract()` provides a consistent way to pull nested components and parameters from supported objects.
* New model-registry helpers make it possible to register user-defined models for direct use with the broader `acousticTS` workflow:
  * `available_models()`
  * `register_model()`
  * `unregister_model()`
  * `reset_model_registry()`
* The `TMM` workflow now includes additional post-processing and diagnostics helpers, including:
  * `tmm_scattering()`
  * `tmm_scattering_grid()`
  * `tmm_products()`
  * `tmm_bistatic_summary()`
  * `tmm_average_orientation()`
  * `tmm_orientation_distribution()`
  * `tmm_diagnostics()`
 
##### Other enhancements :pig2: :lipstick:
* New elastic-property helpers were added, including `lame()`, `pois()`, `bulk()`, `young()`, and `shear()`.
* Core acoustic utilities such as `wavenumber()`, `compressibility()`, `transmission_coefficient()`, `degrees()`, and `radians()` were expanded and cleaned up.
* Example-data support was improved, including bundled fish and benchmark objects used throughout the documentation and validation workflows.

#### Validation and benchmarking :white_check_mark:
* A formal validation and benchmarking registry was added to track benchmarked, validated, partially validated, unvalidated, and experimental model families.
* A new `benchmark_ts` dataset was added to support canonical comparison ladders across modal-series families.
* Reproducible implementation-figure and comparison pipelines were added for documented model families.
* Cross-software and cross-reference validation coverage was expanded across models such as `SPHMS`, `FCMS`, `PSMS`, `SOEMS`, `KRM`, `HPA`, `PCDWBA`, `VESMS`, and `TMM`.
* The pkgdown model library now exposes family-level validation badges and supporting summaries directly in the site.

#### Performance :rocket:
* Major compiled backends were added in `C++` and `Fortran`, substantially improving performance across special functions and core models. This affects the Bessel functions, Legendre functions, spheroidal wave functions, `DWBA` / `SDWBA`, `PSMS`, `ESSMS`, `TMM`, and `VESMS`.
* Adaptive methods were added for `SOEMS` / `calibration`, giving users more control over the balance between speed and precision.
* Optional modal-truncation controls such as `m_limit` were added or expanded across modal-series solutions, allowing users to trade precision for speed more explicitly at higher `ka`.

#### Documentation :open_book:
* The package documentation was fully reorganized into a pkgdown site with model-family overview, theory, and implementation pages.
* A large new vignette set was added, including package primers, model selection, numerical foundations, notation, material properties, boundary conditions, example data, shape/scatterer construction, running models, comparing models, plotting results, simulation sweeps, validation benchmarks, and model creation workflows.
* Model-family landing pages and the model library were added to make the package structure and validation status easier to navigate.
* Figure-generation and vignette-rendering workflows were standardized and expanded to support the implementation and validation documentation.

#### Bug-fixes :bug:
* `DWBA` / `SDWBA`: End-on incidence for some shapes, especially cylinders, could previously return `Inf` values for $\sigma_\text{bs}`. This behavior was fixed.
* Numerous fixes were applied across shape generation, scatterer handling, plotting, model dispatch, documentation rendering, and pkgdown build behavior.
* Several numerical edge cases and validation pathways were tightened across the special-function and modal-series infrastructure.

#### Testing and QA :test_tube:
* The automated test suite was expanded substantially across models, special functions, geometry helpers, simulation workflows, validation utilities, and pkgdown helpers.
* Additional validation-specific and implementation-specific regression tests were added to support the new model families and compiled backends.

#### Infrastructure :building_construction:
* GitHub Actions workflows were expanded to cover `R CMD check`, pkgdown builds, documentation generation, implementation figures, pre-CRAN checks, `rhub`, and test coverage.
* Package metadata and release-support files were improved, including `NEWS.md` and `CITATION.cff`.
* Internal package organization was significantly restructured, including model dispatch, utilities, plotting helpers, validation helpers, and compiled source layout.

**Full Changelog**: https://github.com/brandynlucca/acousticTS/compare/v.1.0...v.2026.03.0

# v1.0.1

* Added a `NEWS.md` file to track changes to the package.
