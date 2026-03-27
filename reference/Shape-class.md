# Generic scattering shape object used throughout this package.

A S4 class that provides slots to contain relevant shape data and
metadata for a variety of arbitrary and canonical shapes and geometries.
See
[Scatterer](https://brandynlucca.github.io/acousticTS/reference/Scatterer-class.md)
for a more detailed description on how this S4 object interacts with
generic
[Scatterer](https://brandynlucca.github.io/acousticTS/reference/Scatterer-class.md)
objects.

## Slots

- `position_matrix`:

  Position matrix that provides the 2D representation of the body shape

- `shape_parameters`:

  A list of additional shape specifications
