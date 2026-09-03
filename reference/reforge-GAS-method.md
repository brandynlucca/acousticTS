# Reforge GAS-class object

Resize a gas-filled fluid scatterer. `GAS` bodies are single-component
fluid-like targets, so the interface mirrors
[`reforge()`](https://brandynlucca.github.io/acousticTS/reference/reforge.md)
for `FLS`: supply either a direct `body_scale` or a `body_target`
describing the desired `length` and/or `radius`. Because both dimensions
are addressable, elongated canonical bodies (prolate spheroids,
cylinders) can be reshaped *non-isometrically* - for example lengthening
the body while holding its radius fixed.

## Usage

``` r
# S4 method for class 'GAS'
reforge(
  object,
  body_scale = NULL,
  body_target = NULL,
  isometric_body = TRUE,
  n_segments_body = NULL,
  scale = NULL,
  radius_target = NULL,
  n_segments = NULL
)
```

## Arguments

- object:

  GAS-class object.

- body_scale:

  Proportional scaling for the body length and/or radius. A single value
  scales both dimensions isometrically; otherwise supply a named numeric
  vector using `length` and/or `radius`.

- body_target:

  Target dimensions (m) for the body length and/or radius, supplied as a
  named numeric vector using `length` and/or `radius`.

- isometric_body:

  Logical; maintain isometric scaling for the body.

- n_segments_body:

  New number of segments along the body profile.

- scale:

  Legacy isometric scale factor applied to every axis. Mutually
  exclusive with `radius_target`.

- radius_target:

  Legacy target *maximum* body radius (m); the isometric scale factor is
  derived as `radius_target / max(current_radius)`. Mutually exclusive
  with `scale`.

- n_segments:

  Legacy alias for `n_segments_body`.

## Value

Modified GAS-class object.

## Details

When the body is a recognized canonical family (sphere, prolate
spheroid, or cylinder) the geometry is regenerated from the requested
dimensions through the corresponding shape constructor, which keeps the
discretization clean and the stored shape descriptor accurate. Arbitrary
shapes fall back to direct per-axis scaling of the position matrix. In
every case the returned object is a `GAS`, with the internal-gas
material contrasts, orientation, and metadata preserved.

The legacy `scale`, `radius_target`, and `n_segments` arguments remain
available and are isometric: they scale every axis together, so a sphere
stays a sphere. Use `body_scale`/`body_target` for independent length
and radius control.

## See also

[`reforge()`](https://brandynlucca.github.io/acousticTS/reference/reforge.md)

## Examples

``` r
methods::selectMethod("reforge", "GAS")
#> new("MethodDefinition", .Data = function (object, ...) 
#> {
#>     .local <- function (object, body_scale = NULL, body_target = NULL, 
#>         isometric_body = TRUE, n_segments_body = NULL, scale = NULL, 
#>         radius_target = NULL, n_segments = NULL) 
#>     {
#>         if (!is.null(scale) && !is.null(radius_target)) {
#>             stop("Specify only one of scale or radius_target, not both.", 
#>                 call. = FALSE)
#>         }
#>         if (!is.null(scale) && (!is.numeric(scale) || length(scale) != 
#>             1 || scale <= 0)) {
#>             stop("'scale' must be a single positive number.", 
#>                 call. = FALSE)
#>         }
#>         if (!is.null(radius_target) && (!is.numeric(radius_target) || 
#>             length(radius_target) != 1 || radius_target <= 0)) {
#>             stop("'radius_target' must be a single positive number.", 
#>                 call. = FALSE)
#>         }
#>         if (!is.null(n_segments) && !is.null(n_segments_body)) {
#>             stop("Specify only one of n_segments or n_segments_body, not both.", 
#>                 call. = FALSE)
#>         }
#>         n_seg <- n_segments_body %||% n_segments
#>         if (!is.null(n_seg) && (!is.numeric(n_seg) || length(n_seg) != 
#>             1 || n_seg < 1)) {
#>             stop("'n_segments' must be a single positive integer.", 
#>                 call. = FALSE)
#>         }
#>         if (!is.null(scale) || !is.null(radius_target)) {
#>             if (!is.null(body_scale) || !is.null(body_target)) {
#>                 stop(paste0("Use either the legacy scale/radius_target arguments or the new ", 
#>                   "body_scale/body_target arguments, not both."), 
#>                   call. = FALSE)
#>             }
#>         }
#>         if (!is.null(body_scale) && !is.null(body_target)) {
#>             stop("Specify only one of body_scale or body_target, not both.", 
#>                 call. = FALSE)
#>         }
#>         if (is.null(body_scale) && is.null(body_target) && is.null(scale) && 
#>             is.null(radius_target) && is.null(n_seg)) {
#>             stop(paste0("Must specify at least one of: body_scale, body_target, scale, ", 
#>                 "radius_target, n_segments, or n_segments_body."), 
#>                 call. = FALSE)
#>         }
#>         body <- acousticTS::extract(object, "body")
#>         shape <- acousticTS::extract(object, "shape_parameters")
#>         rpos <- body$rpos
#>         shape_type <- shape$shape
#>         current_length <- shape$length %||% (max(rpos[, 1]) - 
#>             min(rpos[, 1]))
#>         current_max_r <- max(shape$radius, na.rm = TRUE)
#>         current_n_seg <- shape$n_segments %||% (nrow(rpos) - 
#>             1L)
#>         if (!is.null(radius_target)) {
#>             body_scale <- radius_target/current_max_r
#>         }
#>         else if (!is.null(scale)) {
#>             body_scale <- scale
#>         }
#>         body_dims <- c(length = current_length, radius = current_max_r)
#>         body_target <- .validate_dimensions_target(body_target, 
#>             "body_target", c("length", "radius"))
#>         body_scale_lst <- .reforge_scale_vector(body_scale, body_target, 
#>             body_dims)
#>         body_scales <- .validate_dimension_scaling(dims = body_scale_lst$scale, 
#>             dims_name = paste0("body", body_scale_lst$suffix), 
#>             valid_dims = c("length", "radius"), isometry = isometric_body, 
#>             iso_name = "isometric_body")
#>         new_length <- if (is.null(body_scales)) {
#>             current_length
#>         }
#>         else {
#>             current_length * unname(body_scales["length"])
#>         }
#>         new_radius <- if (is.null(body_scales)) {
#>             current_max_r
#>         }
#>         else {
#>             current_max_r * unname(body_scales["radius"])
#>         }
#>         new_n_seg <- if (!is.null(n_seg)) 
#>             as.integer(n_seg)
#>         else as.integer(current_n_seg)
#>         new_shape <- .reforge_gas_canonical_shape(shape_type, 
#>             new_length, new_radius, new_n_seg)
#>         if (!is.null(new_shape)) {
#>             rebuilt <- gas_generate(shape = new_shape, g_fluid = body$g, 
#>                 h_fluid = body$h, theta_body = body$theta)
#>             methods::slot(object, "body")$rpos <- acousticTS::extract(rebuilt, 
#>                 "body")$rpos
#>             methods::slot(object, "body")$radius <- acousticTS::extract(rebuilt, 
#>                 "body")$radius
#>             methods::slot(object, "shape_parameters") <- acousticTS::extract(rebuilt, 
#>                 "shape_parameters")
#>         }
#>         else {
#>             if (!is.null(n_seg)) {
#>                 rpos <- .resample_rpos(rpos, new_n_seg + 1L)
#>             }
#>             length_ratio <- new_length/current_length
#>             radius_ratio <- new_radius/current_max_r
#>             x_anchor <- min(rpos[, 1], na.rm = TRUE)
#>             rpos[, 1] <- x_anchor + (rpos[, 1] - x_anchor) * 
#>                 length_ratio
#>             radius_cols <- intersect(colnames(rpos), c("y", "zU", 
#>                 "zL"))
#>             for (col in radius_cols) {
#>                 rpos[, col] <- rpos[, col] * radius_ratio
#>             }
#>             methods::slot(object, "body")$rpos <- rpos
#>             methods::slot(object, "body")$radius <- shape$radius * 
#>                 radius_ratio
#>             methods::slot(object, "shape_parameters")$radius <- shape$radius * 
#>                 radius_ratio
#>             methods::slot(object, "shape_parameters")$length <- new_length
#>             methods::slot(object, "shape_parameters")$n_segments <- new_n_seg
#>         }
#>         return(object)
#>     }
#>     .local(object, ...)
#> }, target = new("signature", .Data = "GAS", names = "object", 
#>     package = "acousticTS"), defined = new("signature", .Data = "GAS", 
#>     names = "object", package = "acousticTS"), generic = "reforge")
#> <bytecode: 0x5571fd844910>
#> <environment: namespace:acousticTS>
#> attr(,"target")
#> An object of class “signature”
#> object 
#>  "GAS" 
#> attr(,"defined")
#> An object of class “signature”
#> object 
#>  "GAS" 
#> attr(,"generic")
#> [1] "reforge"
#> attr(,"generic")attr(,"package")
#> [1] "acousticTS"
#> attr(,"class")
#> [1] "MethodDefinition"
#> attr(,"class")attr(,"package")
#> [1] "methods"
```
