# Package index

## Core Workflows

### Running and simulating models

Primary package entry points for deterministic target-strength modeling
and repeated simulation workflows.

- [`target_strength()`](https://brandynlucca.github.io/acousticTS/reference/target_strength.md)
  : Wrapper function to model acoustic target strength
- [`simulate_ts()`](https://brandynlucca.github.io/acousticTS/reference/simulate_ts.md)
  : Simulate target strength (TS) with flexible parameterization and
  batching

### Model registry and extensions

Helpers for listing built-in models and registering custom model
implementations for the current session or across sessions.

- [`available_models()`](https://brandynlucca.github.io/acousticTS/reference/available_models.md)
  : List available target-strength models
- [`register_model()`](https://brandynlucca.github.io/acousticTS/reference/register_model.md)
  : Register a user-defined target-strength model
- [`unregister_model()`](https://brandynlucca.github.io/acousticTS/reference/unregister_model.md)
  : Remove a user-defined target-strength model registration
- [`reset_model_registry()`](https://brandynlucca.github.io/acousticTS/reference/reset_model_registry.md)
  : Clear user-defined model registrations

### Inspection, plotting, and reshaping

Helpers for pulling components out of package objects, visualizing
stored results, and modifying geometry after construction.

- [`extract()`](https://brandynlucca.github.io/acousticTS/reference/extract.md)
  : Extract nested components, slots, or matrix/vector fields from
  objects
- [`plot(`*`<Scatterer>`*`)`](https://brandynlucca.github.io/acousticTS/reference/plot.Scatterer.md)
  : Plot scatterer geometry, stored model results, or stored TMM
  scattering views
- [`reforge()`](https://brandynlucca.github.io/acousticTS/reference/reforge.md)
  : Resize or reparameterize a scatterer object
- [`brake()`](https://brandynlucca.github.io/acousticTS/reference/brake.md)
  : Support function for bending scatterer body shape and position
  matrix

### Bundled benchmark and example data

Built-in objects used throughout the package examples, validation
workflows, and benchmark reproductions.

- [`benchmark_ts`](https://brandynlucca.github.io/acousticTS/reference/benchmark_ts.md)
  : Benchmark model outputs from Jech et al. (2015)
- [`cod`](https://brandynlucca.github.io/acousticTS/reference/cod.md) :
  Sample cod shape with fully inflated swimbladder
- [`krill`](https://brandynlucca.github.io/acousticTS/reference/krill.md)
  : Sample krill (Euphausia superba) shape taken from McGehee et al.
  (1998)
- [`sardine`](https://brandynlucca.github.io/acousticTS/reference/sardine.md)
  : Sample sardine shape with fully inflated swimbladder

## T-Matrix Post-Processing

### Stored-state scattering products

Post-processing helpers that reuse stored T-matrix state for angular
scattering, orientation averaging, diagnostics, and higher-level summary
products.

- [`tmm_scattering()`](https://brandynlucca.github.io/acousticTS/reference/tmm_scattering.md)
  : Evaluate scattering from a stored TMM object
- [`tmm_scattering_grid()`](https://brandynlucca.github.io/acousticTS/reference/tmm_scattering_grid.md)
  : Evaluate a 2D scattering grid from a stored TMM object
- [`tmm_average_orientation()`](https://brandynlucca.github.io/acousticTS/reference/tmm_average_orientation.md)
  : Orientation-average scattering from a stored TMM object
- [`tmm_orientation_distribution()`](https://brandynlucca.github.io/acousticTS/reference/tmm_orientation_distribution.md)
  : Build an orientation distribution for stored TMM post-processing
- [`tmm_bistatic_summary()`](https://brandynlucca.github.io/acousticTS/reference/tmm_bistatic_summary.md)
  : Summarize bistatic products from a stored TMM object
- [`tmm_diagnostics()`](https://brandynlucca.github.io/acousticTS/reference/tmm_diagnostics.md)
  : Diagnostics for stored single-target TMM solutions
- [`tmm_products()`](https://brandynlucca.github.io/acousticTS/reference/tmm_products.md)
  : Collect multiple post-processed products from one stored TMM solve

## Scatterer Types

### Scatterer classes

General scattering classes that dictate expected parameters and model
compatibility.

- [`Scatterer`](https://brandynlucca.github.io/acousticTS/reference/Scatterer-class.md)
  [`Scatterer-class`](https://brandynlucca.github.io/acousticTS/reference/Scatterer-class.md)
  : Scatterer-class object for target strength estimation
- [`CSC`](https://brandynlucca.github.io/acousticTS/reference/CSC-class.md)
  [`CSC-class`](https://brandynlucca.github.io/acousticTS/reference/CSC-class.md)
  : Composite scatterer (CSC) object/class.
- [`ELA`](https://brandynlucca.github.io/acousticTS/reference/ELA-class.md)
  [`ELA-class`](https://brandynlucca.github.io/acousticTS/reference/ELA-class.md)
  : Elastic-based scatterer (ELA) object/class.
- [`BBF`](https://brandynlucca.github.io/acousticTS/reference/BBF-class.md)
  [`BBF-class`](https://brandynlucca.github.io/acousticTS/reference/BBF-class.md)
  : Backboned fish (BBF) object/class.
- [`CAL`](https://brandynlucca.github.io/acousticTS/reference/CAL-class.md)
  [`CAL-class`](https://brandynlucca.github.io/acousticTS/reference/CAL-class.md)
  : Solid and calibration sphere (CAL) object/class.
- [`ESS`](https://brandynlucca.github.io/acousticTS/reference/ESS-class.md)
  [`ESS-class`](https://brandynlucca.github.io/acousticTS/reference/ESS-class.md)
  : Elastic shelled scatterer (ESS) object/class.
- [`GAS`](https://brandynlucca.github.io/acousticTS/reference/GAS-class.md)
  [`GAS-class`](https://brandynlucca.github.io/acousticTS/reference/GAS-class.md)
  : Generic gas-filled scatterer (GAS) object/class.
- [`FLS`](https://brandynlucca.github.io/acousticTS/reference/FLS-class.md)
  [`FLS-class`](https://brandynlucca.github.io/acousticTS/reference/FLS-class.md)
  : Fluid-like scatterer (FLS) object/class.
- [`SBF`](https://brandynlucca.github.io/acousticTS/reference/SBF-class.md)
  [`SBF-class`](https://brandynlucca.github.io/acousticTS/reference/SBF-class.md)
  : Swimbladdered fish (SBF) object/class.

### Scatterer generation

Constructors for the supported scatterer classes.

- [`cal_generate()`](https://brandynlucca.github.io/acousticTS/reference/cal_generate.md)
  : Generate a CAL-class object.
- [`ess_generate()`](https://brandynlucca.github.io/acousticTS/reference/ess_generate.md)
  : Generate an ESS-class object
- [`fls_generate()`](https://brandynlucca.github.io/acousticTS/reference/fls_generate.md)
  : Generate a FLS-class object.
- [`gas_generate()`](https://brandynlucca.github.io/acousticTS/reference/gas_generate.md)
  : Generate a GAS-class object
- [`sbf_generate()`](https://brandynlucca.github.io/acousticTS/reference/sbf_generate.md)
  : Generate a SBF-class object.
- [`bbf_generate()`](https://brandynlucca.github.io/acousticTS/reference/bbf_generate.md)
  : Generate a BBF-class object.

## Scatterer Geometry

### Shape classes

Canonical and arbitrary shape classes used throughout the package.

- [`Shape`](https://brandynlucca.github.io/acousticTS/reference/Shape-class.md)
  [`Shape-class`](https://brandynlucca.github.io/acousticTS/reference/Shape-class.md)
  : Generic scattering shape object used throughout this package.
- [`Arbitrary`](https://brandynlucca.github.io/acousticTS/reference/Arbitrary-class.md)
  [`Arbitrary-class`](https://brandynlucca.github.io/acousticTS/reference/Arbitrary-class.md)
  : Arbitrary body shape
- [`Cylinder`](https://brandynlucca.github.io/acousticTS/reference/Cylinder-class.md)
  [`Cylinder-class`](https://brandynlucca.github.io/acousticTS/reference/Cylinder-class.md)
  : Cylindrical body shape
- [`PolynomialCylinder`](https://brandynlucca.github.io/acousticTS/reference/PolynomialCylinder-class.md)
  [`PolynomialCylinder-class`](https://brandynlucca.github.io/acousticTS/reference/PolynomialCylinder-class.md)
  : Cylindrical body shape deformed using a polynomial
- [`OblateSpheroid`](https://brandynlucca.github.io/acousticTS/reference/OblateSpheroid-class.md)
  [`OblateSpheroid-class`](https://brandynlucca.github.io/acousticTS/reference/OblateSpheroid-class.md)
  : Oblate spheroidal body shape
- [`ProlateSpheroid`](https://brandynlucca.github.io/acousticTS/reference/ProlateSpheroid-class.md)
  [`ProlateSpheroid-class`](https://brandynlucca.github.io/acousticTS/reference/ProlateSpheroid-class.md)
  : Prolate spheroidal body shape
- [`Sphere`](https://brandynlucca.github.io/acousticTS/reference/Sphere-class.md)
  [`Sphere-class`](https://brandynlucca.github.io/acousticTS/reference/Sphere-class.md)
  : Spherical body shape

### Shape generation

Shape constructors, canonicalization helpers, and the convenience
wrapper for canonical geometry.

- [`canonicalize_shape()`](https://brandynlucca.github.io/acousticTS/reference/canonicalize_shape.md)
  : Canonicalize one shape into a canonical surrogate
- [`create_shape()`](https://brandynlucca.github.io/acousticTS/reference/create_shape.md)
  : A wrapper function that automatically creates generalized and/or
  canonical shapes for TS modeling.
- [`arbitrary()`](https://brandynlucca.github.io/acousticTS/reference/arbitrary.md)
  : Creates arbitrary body shape from user inputs
- [`cylinder()`](https://brandynlucca.github.io/acousticTS/reference/cylinder.md)
  : Creates a cylinder.
- [`polynomial_cylinder()`](https://brandynlucca.github.io/acousticTS/reference/polynomial_cylinder.md)
  : Creates a polynomial deformed cylinder
- [`oblate_spheroid()`](https://brandynlucca.github.io/acousticTS/reference/oblate_spheroid.md)
  : Creates an oblate spheroid.
- [`prolate_spheroid()`](https://brandynlucca.github.io/acousticTS/reference/prolate_spheroid.md)
  : Creates a prolate spheroid.
- [`sphere()`](https://brandynlucca.github.io/acousticTS/reference/sphere.md)
  : Creates a sphere.

## Acoustic and Elastic Utilities

Utility functions for basic acoustic quantities, contrasts, and elastic
material-property conversions.

### Acoustic quantities and coefficients

- [`wavenumber()`](https://brandynlucca.github.io/acousticTS/reference/wavenumber.md)
  : Calculate the acoustic wavenumber (\\k\\) for a given frequency and
  sound speed.
- [`linear()`](https://brandynlucca.github.io/acousticTS/reference/linear.md)
  [`db()`](https://brandynlucca.github.io/acousticTS/reference/linear.md)
  : Convert between logarithmic (dB) and linear domains for backscatter
  values.
- [`transmission_coefficient()`](https://brandynlucca.github.io/acousticTS/reference/transmission_coefficient.md)
  : Plane wave/plane interface transmission coefficient
- [`compressibility()`](https://brandynlucca.github.io/acousticTS/reference/compressibility.md)
  : Calculate she compressibility (\\\kappa\\) of a scattering
  boundary/interface.
- [`rho()`](https://brandynlucca.github.io/acousticTS/reference/rho.md)
  : Calculate the density contrast (\\\rho\\) of a scattering boundary

### Elastic moduli for homogeneous isotropic materials

- [`pois()`](https://brandynlucca.github.io/acousticTS/reference/pois.md)
  : Calculate the Poisson's ratio (\\\nu\\)
- [`bulk()`](https://brandynlucca.github.io/acousticTS/reference/bulk.md)
  : Calculate the bulk modulus (K).
- [`young()`](https://brandynlucca.github.io/acousticTS/reference/young.md)
  : Calculate Young's modulus (E).
- [`shear()`](https://brandynlucca.github.io/acousticTS/reference/shear.md)
  : Calculate the shear modulus (G)
- [`lame()`](https://brandynlucca.github.io/acousticTS/reference/lame.md)
  : Calculate Lamé's first parameter (\\\lambda\\)

## Mathematical Functions

Special functions and numerical helpers used throughout the package.

### Summation, quadrature, and vector helpers

- [`along_sum()`](https://brandynlucca.github.io/acousticTS/reference/along_sum.md)
  : Along-matrix summing function
- [`gauss_legendre()`](https://brandynlucca.github.io/acousticTS/reference/gauss_legendre.md)
  : Gauss-Legendre nodes and weights
- [`neumann()`](https://brandynlucca.github.io/acousticTS/reference/neumann.md)
  : Compute the Neumann factor \\\nu\_{n}\\
- [`vecnorm()`](https://brandynlucca.github.io/acousticTS/reference/vecnorm.md)
  : Calculates the Euclidean norm across each row of a given matrix.

### Bessel functions (cylindrical)

- [`jc()`](https://brandynlucca.github.io/acousticTS/reference/jc.md)
  [`jcdk()`](https://brandynlucca.github.io/acousticTS/reference/jc.md)
  : Cylindrical Bessel function of the first kind, \\J\_\nu(z)\\, and
  its respective derivatives
- [`yc()`](https://brandynlucca.github.io/acousticTS/reference/yc.md)
  [`ycdk()`](https://brandynlucca.github.io/acousticTS/reference/yc.md)
  : Cylindrical Bessel function of the second kind, \\Y\_\nu(x)\\, and
  its respective derivatives
- [`hc()`](https://brandynlucca.github.io/acousticTS/reference/hc.md)
  [`hcdk()`](https://brandynlucca.github.io/acousticTS/reference/hc.md)
  : Cylindrical Bessel function of the third kind (Hankel),
  \\H\_\nu(x)\\, and its respective derivatives

### Bessel functions (spherical)

- [`js()`](https://brandynlucca.github.io/acousticTS/reference/js.md)
  [`jsdk()`](https://brandynlucca.github.io/acousticTS/reference/js.md)
  : Spherical Bessel function of the first kind, \\j\_\nu(z)\\, and its
  respective derivatives
- [`ys()`](https://brandynlucca.github.io/acousticTS/reference/ys.md)
  [`ysdk()`](https://brandynlucca.github.io/acousticTS/reference/ys.md)
  : Spherical Bessel function of the second kind, \\y\_\nu(z)\\, and its
  respective derivatives
- [`hs()`](https://brandynlucca.github.io/acousticTS/reference/hs.md)
  [`hsdk()`](https://brandynlucca.github.io/acousticTS/reference/hs.md)
  : Spherical Bessel function of the third kind (Hankel), \\h\_\nu(x)\\,
  and its respective derivatives

### Legendre functions

- [`Pn()`](https://brandynlucca.github.io/acousticTS/reference/Pn.md) :
  Legendre Polynomial of the First Kind, \\P\_\nu(x)\\
- [`Pndk()`](https://brandynlucca.github.io/acousticTS/reference/Pndk.md)
  : Derivative of the Legendre Polynomial of the First Kind
- [`Qn()`](https://brandynlucca.github.io/acousticTS/reference/Qn.md) :
  Legendre Function of the Second Kind, \\Q\_\nu(x)\\
- [`Qndk()`](https://brandynlucca.github.io/acousticTS/reference/Qndk.md)
  : Derivative of the Legendre function of the second kind

### Spheroidal wave functions

- [`Smn()`](https://brandynlucca.github.io/acousticTS/reference/Smn.md)
  : Prolate Spheroidal Angular Function of the First Kind,
  \\S^{1}\_{mn}(c, \eta)\\
- [`Rmn()`](https://brandynlucca.github.io/acousticTS/reference/Rmn.md)
  : Prolate Spheroidal Radial Functions

### Angle conversion

- [`degrees()`](https://brandynlucca.github.io/acousticTS/reference/degrees.md)
  : Convert angular measurements from radians to degrees
- [`radians()`](https://brandynlucca.github.io/acousticTS/reference/radians.md)
  : Convert angular measurements from degrees to radians.
