# Package index

## Scatterer types

### Scatterer classes

General scattering classes that dictate expected parameters and
compatibility with different classes of scattering models.

- [`Scatterer-class`](https://brandynlucca.github.io/acousticTS/reference/Scatterer-class.md)
  [`Scatterer`](https://brandynlucca.github.io/acousticTS/reference/Scatterer-class.md)
  : Scatterer-class object for target strength estimation
- [`CAL-class`](https://brandynlucca.github.io/acousticTS/reference/CAL-class.md)
  [`CAL`](https://brandynlucca.github.io/acousticTS/reference/CAL-class.md)
  : Solid and calibration sphere (CAL) object/class.
- [`ESS-class`](https://brandynlucca.github.io/acousticTS/reference/ESS-class.md)
  [`ESS`](https://brandynlucca.github.io/acousticTS/reference/ESS-class.md)
  : Elastic shelled scatterer (ESS) object/class.
- [`FLS-class`](https://brandynlucca.github.io/acousticTS/reference/FLS-class.md)
  [`FLS`](https://brandynlucca.github.io/acousticTS/reference/FLS-class.md)
  : Fluid-like scatterer (FLS) object/class.
- [`GAS-class`](https://brandynlucca.github.io/acousticTS/reference/GAS-class.md)
  [`GAS`](https://brandynlucca.github.io/acousticTS/reference/GAS-class.md)
  : Generic gas-filled scatterer (GAS) object/class.
- [`SBF-class`](https://brandynlucca.github.io/acousticTS/reference/SBF-class.md)
  [`SBF`](https://brandynlucca.github.io/acousticTS/reference/SBF-class.md)
  : Swimbladdered fish (SBF) object/class.

### Scatterer generation

Functions that generate each of the supported scatterering classes.

- [`cal_generate()`](https://brandynlucca.github.io/acousticTS/reference/cal_generate.md)
  : Generate a CAL-class object.
- [`ess_generate()`](https://brandynlucca.github.io/acousticTS/reference/ess_generate.md)
  : Generate ESS shape
- [`fls_generate()`](https://brandynlucca.github.io/acousticTS/reference/fls_generate.md)
  : Manually generate a FLS object.
- [`gas_generate()`](https://brandynlucca.github.io/acousticTS/reference/gas_generate.md)
  : Create GAS object
- [`sbf_generate()`](https://brandynlucca.github.io/acousticTS/reference/sbf_generate.md)
  : Manually generate a SBF-class object.

## Scatterer shapes

### Shape classes

Generic shape classes used by various functions that provide the
appropriate information for plotting, model parameterization, and other
features.

- [`Shape-class`](https://brandynlucca.github.io/acousticTS/reference/Shape-class.md)
  [`Shape`](https://brandynlucca.github.io/acousticTS/reference/Shape-class.md)
  : Generic scattering shape object used throughout this package.
- [`Arbitrary-class`](https://brandynlucca.github.io/acousticTS/reference/Arbitrary-class.md)
  [`Arbitrary`](https://brandynlucca.github.io/acousticTS/reference/Arbitrary-class.md)
  : Arbitrary body shape
- [`Cylinder-class`](https://brandynlucca.github.io/acousticTS/reference/Cylinder-class.md)
  [`Cylinder`](https://brandynlucca.github.io/acousticTS/reference/Cylinder-class.md)
  : Cylindrical body shape
- [`PolynomialCylinder-class`](https://brandynlucca.github.io/acousticTS/reference/PolynomialCylinder-class.md)
  [`PolynomialCylinder`](https://brandynlucca.github.io/acousticTS/reference/PolynomialCylinder-class.md)
  : Cylindrical body shape deformed using a polynomial
- [`ProlateSpheroid-class`](https://brandynlucca.github.io/acousticTS/reference/ProlateSpheroid-class.md)
  [`ProlateSpheroid`](https://brandynlucca.github.io/acousticTS/reference/ProlateSpheroid-class.md)
  : Prolate spheroidal body shape
- [`Sphere-class`](https://brandynlucca.github.io/acousticTS/reference/Sphere-class.md)
  [`Sphere`](https://brandynlucca.github.io/acousticTS/reference/Sphere-class.md)
  : Spherical body shape

### Shape generation

Generator functions for producing arbitrary or canonical shapes.

- [`create_shape()`](https://brandynlucca.github.io/acousticTS/reference/create_shape.md)
  : A wrapper function that automatically creates generalized and/or
  canonical shapes for TS modeling.
- [`arbitrary()`](https://brandynlucca.github.io/acousticTS/reference/arbitrary.md)
  : Creates arbitrary body shape from user inputs
- [`cylinder()`](https://brandynlucca.github.io/acousticTS/reference/cylinder.md)
  : Creates a cylinder.
- [`polynomial_cylinder()`](https://brandynlucca.github.io/acousticTS/reference/polynomial_cylinder.md)
  : Creates a polynomial deformed cylinder
- [`prolate_spheroid()`](https://brandynlucca.github.io/acousticTS/reference/prolate_spheroid.md)
  : Creates a prolate spheroid.
- [`sphere()`](https://brandynlucca.github.io/acousticTS/reference/sphere.md)
  : Creates a sphere.

## Scattering properties

Functions used for computing different scattering properties
(e.g. contrasts, impedances, coefficients) used in various scattering
models.

### Elastic moduli for homogoneous isotropic materials

- [`bulk()`](https://brandynlucca.github.io/acousticTS/reference/bulk.md)
  : Calculate the bulk modulus (K).
- [`lame()`](https://brandynlucca.github.io/acousticTS/reference/lame.md)
  : Calculate Lamé's first parameter (\\\lambda\\)
- [`pois()`](https://brandynlucca.github.io/acousticTS/reference/pois.md)
  : Calculate the Poisson's ratio (\\\nu\\)
- [`shear()`](https://brandynlucca.github.io/acousticTS/reference/shear.md)
  : Calculate the shear modulus (G)
- [`young()`](https://brandynlucca.github.io/acousticTS/reference/young.md)
  : Calculate Young's modulus (E).

## Mathematical functions

Mathematical functions used throughout the `acousticTS` package.

### Bessel functions

- [`jc()`](https://brandynlucca.github.io/acousticTS/reference/jc.md)
  [`jcd()`](https://brandynlucca.github.io/acousticTS/reference/jc.md)
  [`jcdd()`](https://brandynlucca.github.io/acousticTS/reference/jc.md)
  : Cylindrical Bessel function of the first kind and its respective
  derivatives
- [`js()`](https://brandynlucca.github.io/acousticTS/reference/js.md)
  [`jsd()`](https://brandynlucca.github.io/acousticTS/reference/js.md)
  [`jsdd()`](https://brandynlucca.github.io/acousticTS/reference/js.md)
  : Spherical Bessel function of the first kind and its respective
  derivatives
- [`yc()`](https://brandynlucca.github.io/acousticTS/reference/yc.md)
  [`ycd()`](https://brandynlucca.github.io/acousticTS/reference/yc.md) :
  Cylindrical Bessel (Neumann) function of the second kind and its
  derivatives
- [`ys()`](https://brandynlucca.github.io/acousticTS/reference/ys.md)
  [`ysd()`](https://brandynlucca.github.io/acousticTS/reference/ys.md)
  [`ysdd()`](https://brandynlucca.github.io/acousticTS/reference/ys.md)
  : Spherical Bessel function of the second kind and its respective
  derivative
- [`hc()`](https://brandynlucca.github.io/acousticTS/reference/hc.md)
  [`hcd()`](https://brandynlucca.github.io/acousticTS/reference/hc.md)
  [`hcdd()`](https://brandynlucca.github.io/acousticTS/reference/hc.md)
  : Cylindrical Bessel (Hankel) function of the third kind and its
  derivatives
- [`hcdk()`](https://brandynlucca.github.io/acousticTS/reference/hcdk.md)
  : k-th derivative of the cylindrical Hankel function of the first kind
- [`hs()`](https://brandynlucca.github.io/acousticTS/reference/hs.md)
  [`hsd()`](https://brandynlucca.github.io/acousticTS/reference/hs.md) :
  Spherical Bessel function of the third kind and its respective
  derivative

### Legendre functions

- [`Pn()`](https://brandynlucca.github.io/acousticTS/reference/Pn.md) :
  Legendre Polynomial function (Pn) of the first kind.

### Vector operations

- [`vecnorm()`](https://brandynlucca.github.io/acousticTS/reference/vecnorm.md)
  : Calculates the Euclidean norm across each row of a given matrix.

## Conversion functions

### Acoustics

- [`linear()`](https://brandynlucca.github.io/acousticTS/reference/linear.md)
  [`db()`](https://brandynlucca.github.io/acousticTS/reference/linear.md)
  : Convert backscatter values from log- to linear-domain.

### Trigonometry

- [`degrees()`](https://brandynlucca.github.io/acousticTS/reference/degrees.md)
  : Convert angular measurements from radians to degrees
- [`radians()`](https://brandynlucca.github.io/acousticTS/reference/radians.md)
  : Convert angular measurements from degrees to radians.
