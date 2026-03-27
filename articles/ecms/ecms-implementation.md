# Elastic cylinder modal series solution

## acousticTS implementation

Experimental Unvalidated

*Model-family pages:*
[Overview](https://brandynlucca.github.io/acousticTS/articles/ecms/index.md)
[Implementation](https://brandynlucca.github.io/acousticTS/articles/ecms/ecms-implementation.md)
[Theory](https://brandynlucca.github.io/acousticTS/articles/ecms/ecms-theory.md)

These pages sit between the classical elastic-cylinder literature and
later finite-length cylinder approximations used in fisheries acoustics
([Faran 1951](#ref-faran_sound_1951); [Stanton
1988](#ref-stanton_sound_1988)).

The elastic-cylinder modal-series solution is available through
`target_strength(..., model = "ecms")`. The preferred geometry carrier
is an elastic-cylinder `ESS` object, while the elastic material
parameters are supplied through:

- `density_body`
- `sound_speed_longitudinal_body`
- `sound_speed_transversal_body`

The current implementation check uses an independent direct
transcription of the Faran-Stanton algebra rather than a
package-to-package comparison. That is useful for confirming that the
package code path reproduces the intended elastic-cylinder algebra on a
shared grid, but it is not a substitute for an external benchmark ladder
or a separate public software implementation.

`ECMS` is still marked unvalidated because the current check is an
independent algebra transcription rather than an external benchmark
ladder or separate public software implementation.

### Reference case

The reference cylinder uses:

- length `40 mm`
- radius `5 mm`
- body density `2800 kg m^-3`
- longitudinal speed `6398 m s^-1`
- transverse speed `3122 m s^-1`
- surrounding water density `1026.8 kg m^-3`
- surrounding water sound speed `1477.3 m s^-1`
- broadside incidence
- `12-200 kHz` in `2 kHz` steps

In acousticTS, the call is:

``` r
library(acousticTS)

elastic_cylinder <- fls_generate(
  shape = cylinder(
    length_body = 0.04,
    radius_body = 0.005,
    n_segments = 201
  ),
  density_body = 2800,
  sound_speed_body = 1500,
  theta_body = pi / 2
)

elastic_cylinder <- target_strength(
  elastic_cylinder,
  frequency = seq(12e3, 200e3, by = 2e3),
  model = "ecms",
  density_sw = 1026.8,
  sound_speed_sw = 1477.3,
  sound_speed_longitudinal_body = 6398,
  sound_speed_transversal_body = 3122
)

head(extract(elastic_cylinder, "model")$ECMS)
```

    ##   frequency        ka                       f_bs     sigma_bs        TS
    ## 1     12000 0.2551893 -0.001166938+1.400062e-05i 1.361941e-06 -58.65842
    ## 2     14000 0.2977208 -0.001556457+2.480975e-05i 2.423172e-06 -56.15616
    ## 3     16000 0.3402524 -0.001985893+4.061151e-05i 3.945421e-06 -54.03907
    ## 4     18000 0.3827839 -0.002447214+6.270084e-05i 5.992788e-06 -52.22371
    ## 5     20000 0.4253155 -0.002931593+9.261615e-05i 8.602815e-06 -50.65359
    ## 6     22000 0.4678470 -0.003429548+1.322009e-04i 1.177928e-05 -49.28881

### Implementation check

| Diagnostic                       | Value |
|:---------------------------------|------:|
| Max abs. delta TS (dB)           |  0.00 |
| Mean abs. delta TS (dB)          |  0.00 |
| Frequency at max delta (kHz)     | 12.00 |
| acousticTS elapsed (s)           |  0.08 |
| Direct transcription elapsed (s) |  0.03 |

This is not presented as a benchmark. It is an implementation identity
check: the package output and the independent algebra transcription
coincide on the shared grid, so the current `ECMS` code path is
reproducing the stated elastic-cylinder series rather than drifting
numerically from it.

#### Spectrum overlay

![Pre-rendered ECMS comparison showing the direct reference spectrum,
the acousticTS spectrum, and the residual across
frequency.](ecms-reference-comparison.png)

### Closing note

The point of this page is therefore narrower than the benchmarked
modal-series families. `ECMS` is not being claimed here as externally
validated. What is being documented is that the package reproduces the
elastic-cylinder algebra it claims to implement, across a full frequency
band rather than only at a few checkpoint frequencies.

## References

Faran, James J. 1951. “Sound Scattering by Solid Cylinders and Spheres.”
*The Journal of the Acoustical Society of America* 23 (4): 405–18.
<https://doi.org/10.1121/1.1906780>.

Stanton, T. K. 1988. “Sound Scattering by Cylinders of Finite Length. I.
Fluid Cylinders.” *The Journal of the Acoustical Society of America* 83
(1): 55–63. <https://doi.org/10.1121/1.396184>.
