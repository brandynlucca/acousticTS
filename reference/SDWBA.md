# Stochastic distorted wave Born approximation (SDWBA) for weak scatterers

Calculates the far-field scattering amplitude and related quantities for
fluid-like, weak scatterers using the stochastic distorted wave Born
approximation (SDWBA), as described by Demer and Conti (2003). The SDWBA
extends the deterministic DWBA by incorporating stochastic phase
variability to account for unresolved structural complexity and dynamic
variability in biological scatterers.

## Value

No value; this help topic documents the SDWBA model.

## Usage

This model is accessed via:


    target_strength(
      ...,
      model = "sdwba",
      n_iterations,
      n_segments_init,
      phase_sd_init,
      length_init,
      frequency_init,
      sound_speed_sw,
      density_sw
    )

## Arguments

- `n_iterations`:

  Number of stochastic realizations for averaging target strength
  predictions.

- `n_segments_init`:

  Reference number of body segments.

- `phase_sd_init`:

  Reference phase deviation (radians).

- `length_init`:

  Reference body length (m).

- `frequency_init`:

  Reference frequency (Hz).

- `sound_speed_sw`:

  Seawater sound speed (\\m~s^{-1}\\).

- `density_sw`:

  Seawater density (\\kg~m^{-3}\\).

## Theory

The SDWBA is derived under the weak scattering assumption, where the
differences in compressibility (\\\kappa\\) and density (\\\rho\\)
between the scatterer and the surrounding fluid are small enough to
linearize the acoustic scattering problem (see
[`DWBA`](https://brandynlucca.github.io/acousticTS/reference/DWBA.md)).
in this regime, multiple scattering within the body is neglected, and
the total scattered field is approximated as the coherent sum of
first-order contributions from individual body segments.

The key extension introduced by the SDWBA is the inclusion of stochastic
phase variability to represent unresolved morphological complexity,
internal inhomogeneity, and dynamic effects such as body flexure and
orientation variability. The linear scattering coefficient is written
as:

\$\$ f\_{bs}(\theta) = \sum\limits\_{j=1}^N f\_{bs}^{(j)}(\theta) \exp(i
\varphi_j), \$\$

where \\N\\ is the number of body segments, \\f\_{bs}^{(j)}\\ is the
contribution from segment \\j\\, and \\\varphi_j\\ is a random phase
perturbation drawn independently for each segment.

The phase perturbations are assumed to follow a zero-mean Gaussian
distribution with variance related to the effective signal-to-noise
ratio (SNR) of the scattering process. The minimum expected phase
variance due to noise is given by:

\$\$ \mathbb{V}(\varphi_j) = \frac{1}{2 \mathrm{SNR}}, \$\$

though in practice larger variances are used to account for additional
physical sources of phase decorrelation not explicitly modeled.

The expected backscattering cross-section is obtained by ensemble
averaging over multiple stochastic realizations: \$\$ \langle
\sigma\_{bs}(\theta) \rangle = \mathbb{E}\\\left\[ \left\|
f\_{bs}(\theta) \right\|^2 \right\] \approx \frac{1}{M} \sum\_{m=1}^{M}
\left\| f\_{bs}^{(m)}(\theta) \right\|^2, \$\$

The reported `TS` is \\10\log\_{10}\langle\sigma\_{bs}\rangle\\ in dB re
1 square metre.

Above the minimum segment count, phase variability scales approximately
as:

\$\$ \mathrm{sd}\_{\varphi}(f)\\ f = \mathrm{sd}\_{\varphi_0}\\ f_0,
\$\$

The segment count scales with frequency and length:

\$\$ N(f, L) = \max\left(N_0, \left\lceil N_0 \frac{f L}{f_0
L_0}\right\rceil\right). \$\$

The phase standard deviation at arbitrary frequency and length is then:

\$\$ \mathrm{sd}\_{\varphi}(f, L) = \mathrm{sd}\_{\varphi_0} \frac{N_0
L}{N(f, L) L_0}. \$\$

Here \\L\\ is the stored shape length, inferred from axial span for
coordinate-only objects. Match its convention to `length_init`.

## Implementation

Independent random phases are applied to segment DWBA amplitudes. The
resulting cross sections are averaged before conversion to TS.

Resampling preserves nodes at unchanged resolution; otherwise, it
linearly interpolates positions and radii on a uniform axial grid,
preserving endpoints and direction. This requires monotonic axial
coordinates and changes the independent phase intervals.

`TS_sd` is the sample standard deviation of realization TS in dB. It is
`NA` for fewer than two realizations or any non-finite TS.

## References

Conti, D.A., and Conti, S.G. (2006). Improved parameterization of the
SDWBA for estimating krill target strength. ICES Journal of Marine
Science, 63: 928-935.

Demer, D.A., and Conti, S.G. (2003). Reconciling theoretical versus
empirical target strengths of krill: effects of phase variability on the
distorted-wave Born approximation. ICES Journal of Marine Science, 60:
429-434.

Stanton, T.K., Chu, D., and Wiebe, P.H. (1998). Sound scattering by
several zooplankton groups. II. Scattering models. The Journal of the
Acoustical Society of America, 103, 236-253.

## See also

See the boundary conditions documentation for more details on weak
scattering assumptions,
[`target_strength`](https://brandynlucca.github.io/acousticTS/reference/target_strength.md),
[`FLS`](https://brandynlucca.github.io/acousticTS/reference/FLS-class.md),
[`DWBA`](https://brandynlucca.github.io/acousticTS/reference/DWBA.md)

## Examples

``` r
subset(available_models(), model == "sdwba")
#>       model  slot  source persistent aliases
#> sdwba sdwba SDWBA builtin      FALSE        
```
