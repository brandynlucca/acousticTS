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

## Value

A pre-generated `SBF` scatterer object.

## References

Conti, S.G., and Demer, D.A. (2003). Wide-bandwidth acoustical
characterization of anchovy and sardine from reverberation measurements
in an echoic tank. *ICES Journal of Marine Science*, 60, 617-624.
[doi:10.1016/S1054-3139(03)00056-0](https://doi.org/10.1016/S1054-3139%2803%2900056-0)

## Examples

``` r
data("sardine", package = "acousticTS")
sardine
#> SBF-object
#>  Swimbladdered fish (SBF) 
#>  ID:Sardinops sagax caerulea (Conti and Demer, 2003)
#> Body dimensions:
#>  Length:0.21 m(n = 379 cylinders)
#>  Mean radius:0.0155 m
#>  Max radius:0.0214 m
#> Bladder dimensions:
#>  Length:0.085 m(n = 154 cylinders)
#>  Mean radius:0.0048 m
#>  Max radius:0.0078 m
#> Body material properties:
#>  Density: 1070 kg m^-3 | Sound speed: 1570 m s^-1
#> Bladder fluid material properties:
#>  Density: 1.24 kg m^-3 | Sound speed: 345 m s^-1
#> Body orientation (relative to transducer face/axis):1.571 radians
```
