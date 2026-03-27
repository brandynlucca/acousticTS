# Convert backscatter values from log- to linear-domain.

The `linear(...)` function converts a given value into the linear
domain, while the `db(...)` function converts inputs into the log
domain.

## Usage

``` r
linear(value, coefficient = 10)

db(value, coefficient = 10)
```

## Arguments

- value:

  Logarithmic (e.g. TS) or linear (\\\sigma_bs\\) value

- coefficient:

  Optional. Numeric coefficient preceding the logarithm. Default is 10.

## Value

Transforms the backscattering response into either the log (`db`) or
linear (`linear`) domains.
