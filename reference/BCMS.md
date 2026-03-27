# Bent cylinder modal series (BCMS) solution

Computes backscatter from straight and uniformly bent finite cylinders
by combining the exact finite-cylinder modal-series kernel with the
equivalent coherent-length correction used for uniformly bent cylinders
near normal incidence.

## Usage

This model is accessed via:

    target_strength(
      ...,
      model = "BCMS",
      boundary,
      sound_speed_sw,
      density_sw,
      m_limit
    )

## Arguments

- `boundary`:

  Boundary condition at the cylinder surface. One of `"fixed_rigid"`,
  `"pressure_release"`, `"liquid_filled"`, or `"gas_filled"`.

- `sound_speed_sw`:

  Seawater sound speed (\\m~s^{-1}\\).

- `density_sw`:

  Seawater density (\\kg~m^{-3}\\).

- `m_limit`:

  Optional model truncation limit used to cap the number of retained
  cylindrical modes.

## References

Stanton, T.K. (1988). Sound scattering by cylinders of finite length. I.
Fluid cylinders. *The Journal of the Acoustical Society of America*, 83:
55-63.

Stanton, T.K. (1989). Sound scattering by cylinders of finite length.
III. Deformed cylinders. *The Journal of the Acoustical Society of
America*, 85: 232-237.

Stanton, T.K., Chu, D., Wiebe, P.H., and Clay, C.S. (1993). Average
echoes from randomly oriented random-length finite cylinders:
zooplankton models. *The Journal of the Acoustical Society of America*,
94: 3463-3472.
