# Reforge CAL-class object

Resize a calibration sphere by applying an isometric scale factor or
specifying a target diameter. Optionally re-discretize to a new segment
count. CAL objects are always spheres, so the position matrix follows
the
[`sphere()`](https://brandynlucca.github.io/acousticTS/reference/sphere.md)
convention (n_points x 5: x, y, z, zU, zL).

## Usage

``` r
# S4 method for class 'CAL'
reforge(object, scale = NULL, diameter_target = NULL, n_segments = NULL)
```

## Arguments

- object:

  CAL-class object.

- scale:

  Single positive scale factor applied isometrically. Mutually exclusive
  with `diameter_target`.

- diameter_target:

  Target sphere diameter (m). Derives the scale factor internally.
  Mutually exclusive with `scale`.

- n_segments:

  New number of discrete segments along the major axis.

## Value

Modified CAL-class object.

## Examples

``` r
methods::selectMethod("reforge", "CAL")
#> new("MethodDefinition", .Data = function (object, ...) 
#> {
#>     .local <- function (object, scale = NULL, diameter_target = NULL, 
#>         n_segments = NULL) 
#>     {
#>         if (is.null(scale) && is.null(diameter_target) && is.null(n_segments)) {
#>             stop("Must specify at least one of: scale, diameter_target, or n_segments.", 
#>                 call. = FALSE)
#>         }
#>         if (!is.null(scale) && !is.null(diameter_target)) {
#>             stop("Specify only one of scale or diameter_target, not both.", 
#>                 call. = FALSE)
#>         }
#>         if (!is.null(scale) && (!is.numeric(scale) || length(scale) != 
#>             1 || scale <= 0)) {
#>             stop("'scale' must be a single positive number.", 
#>                 call. = FALSE)
#>         }
#>         if (!is.null(diameter_target) && (!is.numeric(diameter_target) || 
#>             length(diameter_target) != 1 || diameter_target <= 
#>             0)) {
#>             stop("'diameter_target' must be a single positive number.", 
#>                 call. = FALSE)
#>         }
#>         if (!is.null(n_segments) && (!is.numeric(n_segments) || 
#>             length(n_segments) != 1 || n_segments < 1)) {
#>             stop("'n_segments' must be a single positive integer.", 
#>                 call. = FALSE)
#>         }
#>         body <- acousticTS::extract(object, "body")
#>         shape <- acousticTS::extract(object, "shape_parameters")
#>         rpos <- body$rpos
#>         current_radius <- shape$radius_body %||% shape$radius
#>         if (!is.null(diameter_target)) {
#>             scale <- (diameter_target/2)/current_radius
#>         }
#>         if (!is.null(n_segments)) {
#>             rpos <- .resample_rpos(rpos, as.integer(n_segments) + 
#>                 1L)
#>             methods::slot(object, "shape_parameters")$n_segments <- as.integer(n_segments)
#>         }
#>         if (!is.null(scale)) {
#>             rpos <- rpos * scale
#>             new_radius <- current_radius * scale
#>             methods::slot(object, "body")$radius <- new_radius
#>             methods::slot(object, "body")$diameter <- new_radius * 
#>                 2
#>             methods::slot(object, "shape_parameters")$radius <- new_radius
#>             methods::slot(object, "shape_parameters")$diameter <- new_radius * 
#>                 2
#>         }
#>         methods::slot(object, "body")$rpos <- rpos
#>         return(object)
#>     }
#>     .local(object, ...)
#> }, target = new("signature", .Data = "CAL", names = "object", 
#>     package = "acousticTS"), defined = new("signature", .Data = "CAL", 
#>     names = "object", package = "acousticTS"), generic = "reforge")
#> <bytecode: 0x55e96d1420d8>
#> <environment: namespace:acousticTS>
#> attr(,"target")
#> An object of class “signature”
#> object 
#>  "CAL" 
#> attr(,"defined")
#> An object of class “signature”
#> object 
#>  "CAL" 
#> attr(,"generic")
#> [1] "reforge"
#> attr(,"generic")attr(,"package")
#> [1] "acousticTS"
#> attr(,"class")
#> [1] "MethodDefinition"
#> attr(,"class")attr(,"package")
#> [1] "methods"
```
