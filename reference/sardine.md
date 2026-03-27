# Sample sardine shape with fully inflated swimbladder

A pre-generated SBF scatterer containing all information required for
target strength modeling. The packaged geometry follows the sardine
entry distributed through the NOAA Fisheries KRM reference collection
and archived in the `echoSMs` resources. The object metadata identifies
the target as *Sardinops sagax caerulea* following Conti and Demer
(2003).

## Usage

``` r
data(sardine)
```

## Format

A named list with the following components:

- metadata:

  Relevant and identifying metadata (`list`).

- model_parameters:

  Specified model parameters (`list`).

- model:

  Model outputs and results (`list`).

- body:

  A list with:

  - `rpos`: Position matrix (x, yw, zU, zL; m).

  - `sound_speed`: Flesh sound speed (\\c\_{body}\\, m/s).

  - `density`: Flesh density (\\\rho\_{body}\\, kg/m\\^3\\).

  - `theta`: Orientation relative to transducer (\\\theta\_{body}\\,
    radians).

- bladder:

  A list with:

  - `rpos`: Position matrix (x, yw, zU, zL; m).

  - `sound_speed`: Bladder sound speed (\\c\_{bladder}\\, m/s).

  - `density`: Bladder density (\\\rho\_{bladder}\\, kg/m\\^3\\).

  - `theta`: Orientation relative to transducer (\\\theta\_{bladder}\\,
    radians).

- shape_parameters:

  A named list with:

  - `body`: A list describing the body:

    - `length`: Body length (m).

    - `ncyl`: Number of discrete cylinders.

    - `theta_units`: Units for orientation angle.

    - `length_units`: Units for length.

  - `bladder`: A list describing the swimbladder:

    - `length`: Bladder length (m).

    - `ncyl`: Number of discrete cylinders.

    - `theta_units`: Units for orientation angle.

    - `length_units`: Units for length.

## Source

NOAA Fisheries KRM model reference collection
(<https://www.fisheries.noaa.gov/data-tools/krm-model>) and the archived
`echoSMs` resource set (<https://github.com/ices-tools-dev/echoSMs>).

## References

Conti, S.G., and Demer, D.A. (2003). Wide-bandwidth acoustical
characterization of anchovy and sardine from reverberation measurements
in an echoic tank. *ICES Journal of Marine Science*, 60, 617-624.
[doi:10.1016/S1054-3139(03)00056-0](https://doi.org/10.1016/S1054-3139%2803%2900056-0)
