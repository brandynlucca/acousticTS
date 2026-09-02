# Remove a user-defined target-strength model registration

Remove a user-defined target-strength model registration

## Usage

``` r
unregister_model(name, remove_persisted = TRUE)
```

## Arguments

- name:

  Canonical model name or alias.

- remove_persisted:

  Logical scalar. When `TRUE`, any persisted registry entry is also
  removed from the user's config file.

## Value

Invisibly returns the removed canonical model name.

## Examples

``` r
if (interactive()) {
  initialize_demo <- function(object, frequency, ...) object
  solve_demo <- function(object) object
  register_model(
    "demo",
    initialize = initialize_demo,
    solver = solve_demo
  )
  unregister_model("demo")
}
```
