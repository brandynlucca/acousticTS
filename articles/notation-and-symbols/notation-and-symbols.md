# Notation and symbols

## Introduction

The notation on this page follows standard scattering texts for
Helmholtz, elastic-wave, and basis-expansion formulations ([Morse and
Ingard 1968](#ref-Morse_1968); [Flammer 1957](#ref-Flammer_1957);
[Waterman 2009](#ref-Waterman_2009)).

This page is the compact notation guide for the package theory articles.
The model-family pages still define symbols locally when needed, but the
symbols listed here are intended to carry the same meaning everywhere
unless a page explicitly says otherwise.

## Medium indexing

The package theory pages use one exterior-to-interior convention:

| Symbol | Meaning                                            |
|--------|----------------------------------------------------|
| `1`    | surrounding seawater or ambient exterior fluid     |
| `2`    | first target region encountered from the exterior  |
| `3`    | second target region encountered from the exterior |
| `4`    | third target region encountered from the exterior  |

Examples: an unshelled sphere uses media `1` and `2`, an elastic shell
with internal fluid uses media `1`, `2`, and `3`, and a viscous-elastic
layered sphere can use media `1`, `2`, `3`, and `4`.

## Shared field variables

| Symbol             | Meaning                                                          | Typical units |
|--------------------|------------------------------------------------------------------|---------------|
| p_1^{\mathrm{inc}} | incident pressure in seawater                                    | Pa            |
| p_1^{\mathrm{sca}} | scattered pressure in seawater                                   | Pa            |
| p_1^{\mathrm{tot}} | total exterior pressure, p_1^{\mathrm{inc}} + p_1^{\mathrm{sca}} | Pa            |
| p_j                | pressure in medium j                                             | Pa            |
| \mathbf{u}         | elastic displacement vector                                      | m             |
| \mathbf{v}         | fluid particle velocity                                          | m s^{-1}      |
| \mathbf{n}         | outward unit normal                                              | dimensionless |
| \sigma\_{ij}       | stress-tensor components                                         | Pa            |

## Frequencies, speeds, and wavenumbers

| Symbol   | Meaning                                     | Typical units |
|----------|---------------------------------------------|---------------|
| f        | acoustic frequency                          | Hz            |
| \omega   | angular frequency, 2\pi f                   | rad s^{-1}    |
| c_j      | sound speed in medium j                     | m s^{-1}      |
| k_j      | acoustic wavenumber in medium j, \omega/c_j | m^{-1}        |
| k\_{L,j} | longitudinal elastic wavenumber in medium j | m^{-1}        |
| k\_{T,j} | transverse or shear wavenumber in medium j  | m^{-1}        |
| a        | characteristic radius                       | m             |
| L        | characteristic length                       | m             |
| ka       | acoustic size based on a chosen radius      | dimensionless |

## Material properties and contrasts

| Symbol    | Meaning                                    | Typical units |
|-----------|--------------------------------------------|---------------|
| \rho_j    | density in medium j                        | kg m^{-3}     |
| \kappa_j  | compressibility in medium j                | Pa^{-1}       |
| \lambda_j | Lamé’s first parameter in elastic medium j | Pa            |
| \mu_j     | shear modulus in elastic medium j          | Pa            |
| g\_{ij}   | density contrast, \rho_i / \rho_j          | dimensionless |
| h\_{ij}   | sound-speed contrast, c_i / c_j            | dimensionless |

Important examples: body relative to seawater is expressed with g\_{21}
and h\_{21}, a shell relative to seawater is likewise expressed with
g\_{21} and h\_{21} when the shell is the first interior region, an
internal fluid relative to the shell is expressed with g\_{32} and
h\_{32}, and an internal fluid relative directly to seawater is
expressed with g\_{31} and h\_{31}.

## Scattering quantities

| Symbol                                  | Meaning                                               | Typical units |
|-----------------------------------------|-------------------------------------------------------|---------------|
| f(\theta_s,\phi_s \mid \theta_i,\phi_i) | far-field scattering amplitude                        | m             |
| f\_{\mathrm{bs}}                        | backscattering amplitude                              | m             |
| \sigma\_{\mathrm{bs}}                   | backscattering cross-section, \|f\_{\mathrm{bs}}\|^2  | m^2           |
| \mathrm{TS}                             | target strength, 10 \log\_{10}(\sigma\_{\mathrm{bs}}) | dB            |

## Angles and directions

| Symbol           | Meaning                                                                       | Typical units |
|------------------|-------------------------------------------------------------------------------|---------------|
| \theta_i, \phi_i | incident polar and azimuthal angles                                           | rad           |
| \theta_s, \phi_s | scattered or receive polar and azimuthal angles                               | rad           |
| \theta           | body or target orientation angle when a page uses a 2D axisymmetric reduction | rad           |
| \beta            | local body-tilt angle along a segmented or curved centerline                  | rad           |

## Coordinate-specific notation

Some theory pages use geometry-matched coordinates that introduce
additional symbols:

| Geometry         | Symbols                      | Meaning                                                     |
|------------------|------------------------------|-------------------------------------------------------------|
| sphere           | (r,\theta,\phi)              | ordinary spherical coordinates                              |
| cylinder         | (r,\phi,z)                   | ordinary cylindrical coordinates                            |
| prolate spheroid | (\xi,\eta,\phi)              | prolate spheroidal coordinates                              |
| oblate spheroid  | (\xi,\eta,\phi) or r(\theta) | oblate geometry in a spheroidal or spherical representation |

When a page uses scale factors such as h\_\xi, h\_\eta, or h\_\phi,
those are coordinate metric coefficients, not sound-speed contrasts. The
contrast notation always keeps the medium subscripts explicitly:
h\_{21}, h\_{32}, and so on.

## Companion pages

- [Acoustic scattering
  primer](https://brandynlucca.github.io/acousticTS/articles/acoustic-scattering-primer/acoustic-scattering-primer.md)
- [Scattering boundary
  conditions](https://brandynlucca.github.io/acousticTS/articles/boundary_conditions.md)
- [Material
  properties](https://brandynlucca.github.io/acousticTS/articles/material-properties/material-properties.md)

## References

Flammer, Carson. 1957. *Spheroidal Wave Functions*.
<https://ui.adsabs.harvard.edu/abs/1957spwf.book.....F>.

Morse, Philip M., and K. Uno Ingard. 1968. *Theoretical Acoustics*. New
York, NY: McGraw-Hill.

Waterman, P. C. 2009. “T -Matrix Methods in Acoustic Scattering.” *The
Journal of the Acoustical Society of America* 125 (1): 42–51.
<https://doi.org/10.1121/1.3035839>.
