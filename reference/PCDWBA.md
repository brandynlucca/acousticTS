# Phase-compensated distorted wave Born approximation (PCDWBA)

Computes backscatter from weakly scattering elongated fluid-like targets
using the phase-compensated distorted wave Born approximation (PCDWBA)
described by Chu and Ye (1999). This model is particularly suited to
uniformly bent or gently curved cylindrical bodies with tapered ends,
but it can also operate on arbitrary fluid-like centerline profiles
stored in `FLS` objects.

## Usage

This model is accessed via:

    target_strength(
      ...,
      model = "PCDWBA",
      sound_speed_sw,
      density_sw,
      radius_curvature,
      radius_curvature_ratio
    )

## Arguments

- `sound_speed_sw`:

  Seawater sound speed (\\m~s^{-1}\\).

- `density_sw`:

  Seawater density (\\kg~m^{-3}\\).

- `radius_curvature`:

  Optional radius of curvature (\\m\\) used to rebuild a uniformly bent
  centerline for canonical cylinders.

- `radius_curvature_ratio`:

  Optional radius of curvature divided by body length (\\\rho_c / L\\)
  used to rebuild a uniformly bent centerline for canonical cylinders.

## References

Chu, D., and Ye, Z. (1999). A phase-compensated distorted wave Born
approximation representation of the bistatic scattering by weakly
scattering objects: Application to zooplankton. *The Journal of the
Acoustical Society of America*, 106, 1732-1743.

Stanton, T.K. (1989). Sound scattering by cylinders of finite length.
III. Deformed cylinders. *The Journal of the Acoustical Society of
America*, 86, 691-705.
