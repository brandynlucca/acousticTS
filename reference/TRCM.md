# Two-ray cylinder model (TRCM) for elongated scatterers

Computes the far-field scattering amplitude and related quantities for
elongated fluid-like scatterers using the two-ray approximation model,
as described by Stanton et al. (1993, 1998). The two-ray model is a
high-frequency ray-based approximation that accounts for interference
between reflections from the front and back interfaces of the scatterer,
along with a directivity pattern that depends on the scatterer's
orientation and curvature. The model is computationally efficient and
captures the essential physics of acoustic scattering from elongated
bodies such as zooplankton, particularly at high frequencies where the
acoustic wavelength is much smaller than the organism size. The model
supports both straight and bent cylinders, with the bent cylinder
formulation incorporating the radius of curvature to account for body
shape effects on the directivity pattern. For more details, see the
[expanded documentation on the two-ray cylinder
model](https://brandynlucca.github.io/acousticTS/articles/trcm/trcm-theory.html).

## Usage

This model is accessed via:

    target_strength(
      ...,
      model = "TRCM",
      radius_curvature,
      radius_curvature_ratio,
      radius_cylinder_fun,
      sound_speed_sw,
      density_sw
    )

## Arguments

- `radius_curvature`:

  Radius of curvature for bent cylinders (\\m\\). If `NULL`, the model
  assumes a straight cylinder unless `radius_curvature_ratio` is
  specified.

- `radius_curvature_ratio`:

  Ratio of radius of curvature to body length (\\\rho_c / L\\). Used to
  compute `radius_curvature` if not explicitly provided. Default is
  `NULL`.

- `radius_cylinder_fun`:

  Method for selecting the representative radius when the cylinder has
  variable radius. One of `"center"` (default), `"mean"`, `"median"`, or
  `"max"`.

- `sound_speed_sw`:

  Seawater sound speed (\\m~s^{-1}\\).

- `density_sw`:

  Seawater density (\\kg~m^{-3}\\).

## References

Stanton, T.K., Chu, D., Wiebe, P.H., and Clay, C.S. (1993a). Average
echoes from randomly oriented random-length finite cylinders:
Zooplankton models. The Journal of the Acoustical Society of America,
94: 3463-3472.

Stanton, T.K., Chu, D., Wiebe, P.H., Martin, L.V., and Eastwood, R.L.
(1998). Sound scattering by several zooplankton groups. I. Experimental
determination of dominant scattering mechanisms. The Journal of the
Acoustical Society of America, 103: 225-235.

Stanton, T.K., Chu, D., and Wiebe, P.H. (1998). Sound scattering by
several zooplankton groups. II. Scattering models. The Journal of the
Acoustical Society of America, 103: 236-253.

## See also

[`target_strength`](https://brandynlucca.github.io/acousticTS/reference/target_strength.md),
[`FLS`](https://brandynlucca.github.io/acousticTS/reference/FLS-class.md),
[`Cylinder`](https://brandynlucca.github.io/acousticTS/reference/Cylinder-class.md),
[`cylinder`](https://brandynlucca.github.io/acousticTS/reference/cylinder.md)
