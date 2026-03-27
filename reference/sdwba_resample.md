# Resample shape for SDWBA model with piecewise constant radius

This function resamples the shape of a fluid-like scatterer (FLS) object
for use in stochastic distorted wave Born approximation (SDWBA)
calculations. The resampling preserves the overall shape of the
scatterer while creating a new representation with the specified number
of segments. The radius assignment uses a stepwise algorithm to maintain
piecewise constant radius values across segments.

## Usage

``` r
sdwba_resample(object, n_segments)
```

## Arguments

- object:

  FLS-class object to resample

- n_segments:

  Number of segments in the resampled shape

## Value

FLS object with resampled shape
