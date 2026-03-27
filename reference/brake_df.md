# Support function for bending scatterer position matrix dataframe

Support function for bending scatterer position matrix dataframe

## Usage

``` r
brake_df(body_df, radius_curvature, mode = "ratio")
```

## Arguments

- body_df:

  Dataframe object containing body shape information

- radius_curvature:

  Radius of curvature that can be parameterized either as a ratio
  relative to body length or actual measurement

- mode:

  Either "ratio" or "measurement"
