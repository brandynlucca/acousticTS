# Sample sardine shape with fully inflated swimbladder

A pre-generated SBF scatterer containing all information required for
target strength modeling.

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
