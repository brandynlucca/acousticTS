# Register a user-defined target-strength model

Register a user-defined target-strength model

## Usage

``` r
register_model(
  name,
  initialize,
  solver,
  slot = NULL,
  aliases = character(),
  persist = FALSE,
  overwrite = FALSE
)
```

## Arguments

- name:

  Canonical model name passed to
  [`target_strength()`](https://brandynlucca.github.io/acousticTS/reference/target_strength.md).

- initialize:

  Initializer function or function reference.

- solver:

  Solver function or function reference.

- slot:

  Result-slot name stored on the scatterer. Defaults to the usual
  upper-case model label derived from `name`.

- aliases:

  Optional alternative model names.

- persist:

  Logical scalar. When `TRUE`, the model registration is written to the
  user's `R_user_dir()` config path and reloaded in later sessions.
  Persistent registrations require package-qualified function references
  such as `"mypkg::tsl_initialize"` and `"mypkg::TSL"`.

- overwrite:

  Logical scalar. When `TRUE`, an existing user registration with the
  same canonical name can be replaced. Built-in models cannot be
  overwritten.

## Value

Invisibly returns the normalized registry entry.
