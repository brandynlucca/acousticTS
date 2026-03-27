# Wrapper function incorporating phase deviation into contour integration

Wrapper function incorporating phase deviation into contour integration

## Usage

``` r
phase_integrate(x, y, n_iterations, integral, phase_sd)
```

## Arguments

- x:

  Indexing argument for multi-row objects

- y:

  Indexing argument for multi-column objects

- n_iterations:

  Number of phase deviations to average and summarize

- integral:

  Integral function used for numerical integration via adaptive
  quadrature

- phase_sd:

  Phase standard deviation
