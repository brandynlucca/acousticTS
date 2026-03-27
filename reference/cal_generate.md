# Generate a CAL-class object.

Generate a CAL-class object.

## Usage

``` r
cal_generate(
  material = "WC",
  diameter = 0.0381,
  sound_speed_longitudinal = NULL,
  sound_speed_transversal = NULL,
  density_sphere = NULL,
  theta_sphere = pi,
  ID = NULL,
  diameter_units = "m",
  theta_units = "radians",
  n_segments = 100
)
```

## Arguments

- material:

  Material-type for the soldi sphere. See 'Details' built-in material
  options.

- diameter:

  Spherical diameter (m).

- sound_speed_longitudinal:

  Longitudinal sound speed (m/s).

- sound_speed_transversal:

  Transversal sound speed (m/s).

- density_sphere:

  Density (kg/m^3).

- theta_sphere:

  Backscattering direction (Default: pi radians).

- ID:

  Optional metadata ID input.

- diameter_units:

  Units for diameter. Defaults to "m".

- theta_units:

  Units for direction. Defaults to "radians".

- n_segments:

  Number of segments to discretize object shape.

## Value

Generates a CAL-class object.

## Details

There are several options for the **material** argument:

|              |              |          |         |                   |                    |
|--------------|--------------|----------|---------|-------------------|--------------------|
| **Material** | **Argument** | **c1**   | **c2**  | **\\\rho1\\**     | *Tungsten carbide* |
| "WC"         | 6853         | 4171     | 14900   | *Stainless steel* | "steel"            |
| 5980         | 3297         | 7970     | *Brass* | "brass"           | 4372               |
| 2100         | 8360         | *Copper* | "Cu"    | 4760              | 2288.5             |
| 8947         | *Aluminum*   | "Al"     | 6260    | 3080              | 2700               |

## See also

[`CAL`](https://brandynlucca.github.io/acousticTS/reference/CAL-class.md)
