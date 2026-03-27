# Single-target transition matrix method (TMM)

Computes monostatic backscatter from a single axisymmetric target using
a transition-matrix formulation. The current implementation targets
smooth bodies of revolution and finite cylinders already represented in
the package as a `Sphere`, `OblateSpheroid`, `ProlateSpheroid`, or
`Cylinder`, and supports rigid, pressure-release, and homogeneous
penetrable fluid/gas interiors.

## Details

This implementation is intentionally scoped to **single targets** and
the monostatic backscatter quantity used by
[`target_strength()`](https://brandynlucca.github.io/acousticTS/reference/target_strength.md).

For spheres and oblate spheroids, the current implementation uses a
spherical-wave basis with an explicit projected boundary solve over the
target surface. For prolate spheroids, it instead uses a
spheroidal-coordinate transition-matrix-equivalent backend, which is the
more natural coordinate system for that geometry and is consistent with
the scalar spheroidal T-matrix literature for single-target scattering.
For finite cylinders, the default monostatic branch uses a
cylindrical-coordinate modal/T-matrix-equivalent backend so that the
backscatter benchmark remains aligned with the exact finite-cylinder
family. When `store_t_matrix = TRUE`, cylinders retain lightweight
cylindrical-family state that supports exact monostatic reuse and
orientation-averaged monostatic products, while general-angle cylinder
bistatic post-processing remains outside the current validated scope.
Because of that narrower validation status, cylinder calls emit a
warning by default; see `options(acousticTS.warn_tmm_cylinder = FALSE)`
to silence it in controlled test or benchmarking workflows.

The present solver is therefore a practical single-target acoustic
T-matrix method motivated by the classic transition-matrix literature,
but it is **not yet** a full implementation of the far-field two-part
T-matrix workflow of Ganesh and Hawkins (2008, 2022), which assumes an
external far-field solver as the first stage.

## Usage

This model is accessed via:

    target_strength(
      ...,
      model = "TMM",
      boundary,
      sound_speed_sw,
      density_sw,
      n_max,
      store_t_matrix
    )

## Arguments

- `boundary`:

  Boundary condition at the target surface. One of `"fixed_rigid"`,
  `"pressure_release"`, `"liquid_filled"`, or `"gas_filled"`.

- `sound_speed_sw`:

  Surrounding-medium sound speed (\\m~s^{-1}\\).

- `density_sw`:

  Surrounding-medium density (\\kg~m^{-3}\\).

- `n_max`:

  Optional truncation limit. For spheres and oblate spheroids, this is
  the maximum spherical-wave degree used in the truncated T-matrix
  solve. For the default monostatic cylinder branch, it is the
  cylindrical modal cutoff used in the geometry-matched backend. When
  left as `NULL`, a geometry-aware rule is used frequency-by-frequency.
  This argument is currently ignored for prolate spheroids, which use
  the spheroidal-coordinate branch.

- `store_t_matrix`:

  Logical flag controlling whether the frequency-specific retained state
  is stored under `object@model_parameters$TMM$parameters$t_matrix`. The
  default is `FALSE` to avoid large object sizes. Explicit block
  retention is available for the spherical and spheroidal branches. For
  cylinders, the stored state keeps the geometry-matched cylindrical
  monostatic family available for exact monostatic reuse and
  orientation-averaged monostatic products; full general-angle cylinder
  bistatic post-processing is not yet provided.

## Theory

For a single target, the incident and scattered fields are expanded in
regular and outgoing modal bases, respectively:

\$\$ p^{inc} = \sum\_{\nu} a\_{\nu} \\ \psi\_{\nu}^{(1)}, \qquad p^{sca}
= \sum\_{\nu} f\_{\nu} \\ \psi\_{\nu}^{(3)}, \$\$

where the transition matrix \\\mathbf{T}\\ maps incident coefficients to
scattered coefficients:

\$\$ \mathbf{f} = \mathbf{T}\mathbf{a}. \$\$

For the axisymmetric single-target case used here, the azimuthal orders
decouple and each block is recovered by enforcing the boundary
conditions on the target surface for the retained modal basis functions.
The backscatter amplitude is then obtained by evaluating the outgoing
expansion in the monostatic receive direction opposite to the incident
plane wave.

## References

Waterman, P. C. (1969). New formulation of acoustic scattering. *The
Journal of the Acoustical Society of America*, **45**, 1417-1429.

Varadan, V. K., Varadan, V. V., Bringi, V. N., and Waterman, P. C.
(1982). Computation of rigid body scattering by prolate spheroids using
the T-matrix approach. *The Journal of the Acoustical Society of
America*, **71**, 22-25.

Hackman, R. H. (1984). An application of the spheroidal-coordinate-based
transition matrix: The acoustic scattering from high aspect ratio
solids. *The Journal of the Acoustical Society of America*, **76**,
1058-1070. Waterman, P. C. (2009). T-matrix methods in acoustic
scattering. *The Journal of the Acoustical Society of America*, **125**,
42-51.

Ganesh, M., and Hawkins, S. C. (2008). A far-field based T-matrix method
for three dimensional acoustic scattering. *Wave Motion*, **45**,
1441-1460.

Ganesh, M., and Hawkins, S. C. (2022). A numerically stable T-matrix
method for acoustic scattering by nonspherical particles with large
aspect ratios and size parameters. *The Journal of the Acoustical
Society of America*, **151**, 1978-1988.

## See also

[`target_strength`](https://brandynlucca.github.io/acousticTS/reference/target_strength.md),
[`FLS`](https://brandynlucca.github.io/acousticTS/reference/FLS-class.md),
[`GAS`](https://brandynlucca.github.io/acousticTS/reference/GAS-class.md),
[`Sphere`](https://brandynlucca.github.io/acousticTS/reference/Sphere-class.md),
[`OblateSpheroid`](https://brandynlucca.github.io/acousticTS/reference/OblateSpheroid-class.md),
[`ProlateSpheroid`](https://brandynlucca.github.io/acousticTS/reference/ProlateSpheroid-class.md),
[`Cylinder`](https://brandynlucca.github.io/acousticTS/reference/Cylinder-class.md),
[`sphere`](https://brandynlucca.github.io/acousticTS/reference/sphere.md),
[`oblate_spheroid`](https://brandynlucca.github.io/acousticTS/reference/oblate_spheroid.md),
[`prolate_spheroid`](https://brandynlucca.github.io/acousticTS/reference/prolate_spheroid.md),
[`cylinder`](https://brandynlucca.github.io/acousticTS/reference/cylinder.md)
