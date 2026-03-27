# Convert between logarithmic (dB) and linear domains for backscatter values.

The `linear` function converts a value from the logarithmic (dB) domain
to the linear domain, while the `db` function converts a value from the
linear domain to the logarithmic (dB) domain. These are commonly used
for target strength (TS) and backscattering coefficient
(\\\sigma\_{bs}\\) conversions.

The conversions are defined as: \$\$\text{linear}(x) = c^{x / c}\$\$
\$\$\text{db}(x) = c \log_c(x)\$\$ where \\c\\ is the coefficient
(default 10).

## Usage

``` r
linear(value, coefficient = 10)

db(value, coefficient = 10)
```

## Arguments

- value:

  Numeric value to convert. For `linear`, this is a logarithmic value
  (e.g., dB TS); for `db`, this is a linear value (e.g.,
  \\\sigma\_{bs}\\).

- coefficient:

  Optional. Numeric coefficient (base) for the logarithm. Default is 10.

## Value

For `linear`, returns the value converted to the linear domain. For
`db`, returns the value converted to the logarithmic (dB) domain.
