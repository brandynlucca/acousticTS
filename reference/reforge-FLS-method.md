# Reforge FLS-class object

Resize a fluid-like scatterer's single body by scaling its `length`
and/or `radius`, supplied either as a direct `body_scale` factor or a
`body_target` dimension (m). The cross-section is circular, so `radius`
drives both lateral and vertical extent together.

## Usage

``` r
# S4 method for class 'FLS'
reforge(
  object,
  body_scale = NULL,
  body_target = NULL,
  isometric_body = TRUE,
  n_segments_body = NULL,
  length = NULL,
  radius = NULL,
  length_radius_ratio_constant = TRUE,
  n_segments = NULL
)
```

## Arguments

- object:

  FLS-class object.

- body_scale:

  Proportional scaling to the body length and radius. When a single
  value is supplied, both dimensions are scaled together. Otherwise,
  this input must be a named numeric vector using `length` and/or
  `radius`.

- body_target:

  Target dimensions (m) for the body length and/or radius. This input
  must be a named numeric vector.

- isometric_body:

  Logical; maintain isometric scaling for body.

- n_segments_body:

  New number of segments along the body profile.

- length:

  Legacy alias for a new body length resize.

- radius:

  Legacy alias for a new maximum body radius.

- length_radius_ratio_constant:

  Legacy toggle controlling whether a length-only resize also rescales
  radius.

- n_segments:

  Legacy alias for `n_segments_body`.

## Value

An updated object of class `FLS`.

## Details

For bent bodies (see
[`brake()`](https://brandynlucca.github.io/acousticTS/reference/brake.md))
a `length` resize follows the true centerline arc length and rescales
the curved path, while a `radius` resize changes only the tube thickness
and leaves the centerline in place. The legacy `length`, `radius`,
`length_radius_ratio_constant`, and `n_segments` arguments are retained
as thin wrappers over the `body_scale`/`body_target` pathway.

## Examples

``` r
methods::selectMethod("reforge", "FLS")
#> new("MethodDefinition", .Data = function (object, ...) 
#> {
#>     .local <- function (object, body_scale = NULL, body_target = NULL, 
#>         isometric_body = TRUE, n_segments_body = NULL, length = NULL, 
#>         radius = NULL, length_radius_ratio_constant = TRUE, n_segments = NULL) 
#>     {
#>         if (!is.null(body_scale) && !is.null(body_target)) {
#>             stop("Specify only one of body_scale or body_target, not both.", 
#>                 call. = FALSE)
#>         }
#>         if (!is.null(n_segments) && !is.null(n_segments_body)) {
#>             stop("Specify only one of n_segments or n_segments_body, not both.", 
#>                 call. = FALSE)
#>         }
#>         if ((!is.null(length) || !is.null(radius)) && (!is.null(body_scale) || 
#>             !is.null(body_target))) {
#>             stop(paste0("Use either the legacy length/radius arguments or the new ", 
#>                 "body_scale/body_target arguments, not both."), 
#>                 call. = FALSE)
#>         }
#>         if (is.null(body_scale) && is.null(body_target) && is.null(length) && 
#>             is.null(radius) && is.null(n_segments) && is.null(n_segments_body)) {
#>             stop(paste0("Must specify at least one of: body_scale, body_target, length, ", 
#>                 "radius, n_segments, or n_segments_body."), call. = FALSE)
#>         }
#>         shape <- acousticTS::extract(object, "shape_parameters")
#>         body <- acousticTS::extract(object, "body")
#>         rpos <- body$rpos
#>         if (!is.null(n_segments)) {
#>             n_segments_body <- n_segments
#>         }
#>         if (!is.null(length) || !is.null(radius)) {
#>             if (!is.null(length) && !is.null(radius)) {
#>                 body_target <- c(length = length, radius = radius)
#>                 isometric_body <- FALSE
#>             }
#>             else if (!is.null(length)) {
#>                 body_target <- c(length = length)
#>                 isometric_body <- isTRUE(length_radius_ratio_constant)
#>             }
#>             else if (!is.null(radius)) {
#>                 body_target <- c(radius = radius)
#>                 isometric_body <- FALSE
#>             }
#>         }
#>         use_arc_length <- !is.null(body$arc_length) || (!is.null(body$radius_curvature_ratio) && 
#>             is.finite(body$radius_curvature_ratio))
#>         body_dims <- c(length = if (use_arc_length) {
#>             .shape_arc_length(position_matrix = rpos, row_major = TRUE)
#>         } else {
#>             shape$length %||% .shape_length(position_matrix = rpos, 
#>                 row_major = TRUE)
#>         }, radius = shape$radius %||% max(.shape_radius_profile(rpos, 
#>             row_major = TRUE), na.rm = TRUE))
#>         body_target <- .validate_dimensions_target(body_target, 
#>             "body_target", c("length", "radius"))
#>         body_scale_lst <- .reforge_scale_vector(body_scale, body_target, 
#>             body_dims)
#>         body_scales <- .validate_dimension_scaling(dims = body_scale_lst$scale, 
#>             dims_name = paste0("body", body_scale_lst$suffix), 
#>             valid_dims = c("length", "radius"), isometry = isometric_body, 
#>             iso_name = "isometric_body")
#>         if (!is.null(n_segments_body)) {
#>             rpos <- .reforge_resample_rows(rpos, as.integer(n_segments_body) + 
#>                 1L)
#>             methods::slot(object, "shape_parameters")$n_segments <- as.integer(n_segments_body)
#>         }
#>         if (!is.null(body_scales)) {
#>             axis_scales <- c(length = unname(body_scales["length"]), 
#>                 width = unname(body_scales["radius"]), height = unname(body_scales["radius"]))
#>             rpos <- .reforge_apply_axis_scaling(rpos, axis_scales, 
#>                 scale_centerline = use_arc_length)
#>         }
#>         radii <- .shape_radius_profile(position_matrix = rpos, 
#>             row_major = TRUE, error_context = "FLS body")
#>         methods::slot(object, "body")$rpos <- rpos
#>         methods::slot(object, "body")$radius <- radii
#>         if (use_arc_length) {
#>             methods::slot(object, "body")$arc_length <- .shape_arc_length(position_matrix = rpos, 
#>                 row_major = TRUE)
#>             methods::slot(object, "shape_parameters")$length <- .shape_arc_length(position_matrix = rpos, 
#>                 row_major = TRUE)
#>         }
#>         else {
#>             methods::slot(object, "shape_parameters")$length <- .shape_length(position_matrix = rpos, 
#>                 row_major = TRUE)
#>         }
#>         methods::slot(object, "shape_parameters")$radius <- max(radii, 
#>             na.rm = TRUE)
#>         return(object)
#>         if (!is.null(length)) {
#>             new_scale <- length/shape$length
#>             if (length_radius_ratio_constant) {
#>                 rpos <- rpos * new_scale
#>                 if (is.null(radius)) {
#>                   radii <- radii * new_scale
#>                 }
#>                 else {
#>                   r_scale <- radius/shape$radius
#>                   radii <- radii * r_scale
#>                   correction <- r_scale/new_scale
#>                   if (nrow(rpos) >= 4) {
#>                     rpos[seq(4L, nrow(rpos)), ] <- rpos[seq(4L, 
#>                       nrow(rpos)), ] * correction
#>                   }
#>                 }
#>                 methods::slot(object, "shape_parameters")$radius <- max(radii)
#>             }
#>             else {
#>                 rpos[1L, ] <- rpos[1L, ] * new_scale
#>             }
#>             methods::slot(object, "shape_parameters")$length <- abs(diff(range(rpos[1L, 
#>                 ])))
#>         }
#>         if (!is.null(radius) && is.null(length)) {
#>             r_scale <- radius/shape$radius
#>             radii <- radii * r_scale
#>             if (nrow(rpos) >= 4) {
#>                 rpos[seq(4L, nrow(rpos)), ] <- rpos[seq(4L, nrow(rpos)), 
#>                   ] * r_scale
#>             }
#>             methods::slot(object, "shape_parameters")$radius <- max(radii)
#>         }
#>         methods::slot(object, "body")$rpos <- rpos
#>         methods::slot(object, "body")$radius <- radii
#>         return(object)
#>     }
#>     .local(object, ...)
#> }, target = new("signature", .Data = "FLS", names = "object", 
#>     package = "acousticTS"), defined = new("signature", .Data = "FLS", 
#>     names = "object", package = "acousticTS"), generic = "reforge")
#> <bytecode: 0x558783417fd0>
#> <environment: namespace:acousticTS>
#> attr(,"target")
#> An object of class “signature”
#> object 
#>  "FLS" 
#> attr(,"defined")
#> An object of class “signature”
#> object 
#>  "FLS" 
#> attr(,"generic")
#> [1] "reforge"
#> attr(,"generic")attr(,"package")
#> [1] "acousticTS"
#> attr(,"class")
#> [1] "MethodDefinition"
#> attr(,"class")attr(,"package")
#> [1] "methods"
```
