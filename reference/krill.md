# Sample krill (Euphausia superba) shape taken from McGehee et al. (1998)

A dataset containing a sample krill (Euphausia superba) body shape
proposed by McGehee et al. (1998).

## Usage

``` r
data(krill)
```

## Format

A pre-generated FLS scatterer containing all information required for
target strength modeling.

- metadata:

  Relevant and identifying metadata (`list`).

- model_parameters:

  Container for specified model parameters (`list`).

- model:

  Model outputs and results (`list`).

- body:

  A list with elements:

  - `rpos`: Position matrix (x, y, z; m).

  - `radius`: Radius of each discrete cylinder (m).

  - `g`: Body density contrast relative to medium.

  - `h`: Sound speed contrast relative to medium.

  - `theta`: Orientation angle (\\\theta\_{body}\\, radians).

- shape_parameters:

  A list with:

  - `length`: Body length (m).

  - `ncyl`: Number of cylinders.

  - `theta_units`: Unit of orientation.

  - `length_units`: Unit of length.
