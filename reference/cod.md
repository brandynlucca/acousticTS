# Sample cod shape with fully inflated swimbladder

A pre-generated SBF scatterer containing all information required for
target strength modeling. The packaged object corresponds to the
historical Atlantic cod example used in the KRM literature and matches
the Cod D case archived in the `echoSMs` KRM shape collection.

## Usage

``` r
data(cod)
```

## Format

A pre-generated SBF scatterer containing all information required for
target strength modeling.

- metadata:

  Relevant and identifying metadata (`list`).

- model_parameters:

  Container for specified model parameters (`list`).

- model:

  Model outputs and results (`list`).

- body:

  List with:

  - `rpos`: Position matrix (x, yw, zU, zL; m).

  - `sound_speed`: Flesh sound speed (\\c\_{body}\\, m/s).

  - `density`: Flesh density (\\\rho\_{body}\\, kg/m\\^3\\).

  - `theta`: Orientation relative to transducer (\\\theta\_{body}\\,
    radians).

- bladder:

  List with:

  - `rpos`: Position matrix (x, yw, zU, zL; m).

  - `sound_speed`: Flesh sound speed (\\c\_{bladder}\\, m/s).

  - `density`: Flesh density (\\\rho\_{bladder}\\, kg/m\\^3\\).

  - `theta`: Orientation relative to transducer (\\\theta\_{bladder}\\,
    radians).

- shape_parameters:

  Named list with:

  - `body`: List with length (m), ncyl (int), theta_units (str),
    length_units (str).

  - `bladder`: List with length (m), ncyl (int), theta_units (str),
    length_units (str).

## Source

Historical KRM cod shape collection archived in `echoSMs`
(<https://github.com/ices-tools-dev/echoSMs>) and associated with the
Clay and Horne Atlantic cod example.

## References

Clay, C.S., and Horne, J.K. (1994). Acoustic models of fish: The
Atlantic cod (*Gadus morhua*). *The Journal of the Acoustical Society of
America*, 96, 1661-1668.
[doi:10.1121/1.410245](https://doi.org/10.1121/1.410245)
