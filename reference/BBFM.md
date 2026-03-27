# Body-backbone composite model (BBFM)

Computes a coherent composite backscatter prediction for
swimbladder-less fish represented as a weakly scattering flesh body plus
an explicit elastic backbone. The flesh contribution is evaluated with
the distorted-wave Born approximation (DWBA), while the backbone
contribution is evaluated with the elastic cylinder modal-series
solution (ECMS). The backbone amplitude is then translated into the
stored body coordinate frame using a phase factor based on the backbone
centroid before the two complex amplitudes are summed.

## Usage

This model is accessed via:

    target_strength(
      ...,
      model = "BBFM",
      sound_speed_sw,
      density_sw,
      m_limit
    )

## Arguments

- `sound_speed_sw`:

  Seawater sound speed (\\m~s^{-1}\\).

- `density_sw`:

  Seawater density (\\kg~m^{-3}\\).

- `m_limit`:

  Optional modal truncation limit passed to the backbone ECMS solve.

## Theory

The model follows the same structural decomposition used for
swimbladder-less mackerel by Gorska, Ona, and Korneliussen (2005): the
flesh is treated as a weakly scattering fluid-like body, and the
backbone is treated as an elastic cylindrical structure. In this
implementation, the total complex backscattering amplitude is \$\$
f\_{bs} = f\_{\mathrm{flesh}} + f\_{\mathrm{backbone}} \exp\left\\ 2 i
k\_{sw} v_c \right\\, \$\$ where \\f\_{\mathrm{flesh}}\\ is obtained
from
[`DWBA`](https://brandynlucca.github.io/acousticTS/reference/DWBA.md),
\\f\_{\mathrm{backbone}}\\ is obtained from
[`ECMS`](https://brandynlucca.github.io/acousticTS/reference/ECMS.md),
and \\v_c = x_c \cos \theta + z_c \sin \theta\\ is the projection of the
stored backbone centroid onto the backscatter direction.

## References

Gorska, N., Ona, E., and Korneliussen, R. (2005). Acoustic
backscattering by Atlantic mackerel as being representative of fish that
lack a swimbladder. Backscattering by individual fish. *ICES Journal of
Marine Science*, 62: 984-995.

Stanton, T.K. (1988). Sound scattering by cylinders of finite length.
II. Elastic cylinders. *The Journal of the Acoustical Society of
America*, 83: 64-67.

Stanton, T.K. (1989). Sound scattering by cylinders of finite length.
III. Deformed cylinders. *The Journal of the Acoustical Society of
America*, 86: 691-705.

Stanton, T.K., Chu, D., and Wiebe, P.H. (1998). Sound scattering by
several zooplankton groups. II. Scattering models. *The Journal of the
Acoustical Society of America*, 103: 236-253.
