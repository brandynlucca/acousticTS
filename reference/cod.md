# Sample code shape with fully inflated swimbladder.

Sample code shape with fully inflated swimbladder.

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
