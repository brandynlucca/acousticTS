# Reforge BBF-class object.

Resize a backboned-fish scatterer by rescaling the flesh body and
elastic backbone independently or together. The body follows the same
profile-based length/width/height scaling used by `SBF`, while the
backbone retains the same cylinder-style component bookkeeping and
preserves its relative axial start within the body when body length
changes.

## Usage

``` r
# S4 method for class 'BBF'
reforge(
  object,
  body_scale = NULL,
  body_target = NULL,
  backbone_scale = NULL,
  backbone_target = NULL,
  isometric_body = TRUE,
  isometric_backbone = TRUE,
  maintain_ratio = TRUE,
  n_segments_body = NULL,
  n_segments_backbone = NULL,
  containment = c("warn", "error", "ignore")
)
```

## Arguments

- object:

  BBF-class object.

- body_scale:

  Proportional scaling to the body length, width, and height dimensions.
  When a single value is supplied, all dimensions are scaled uniformly.
  Otherwise, this input must be a named numeric vector.

- body_target:

  Target dimensions (m) for the body length, width, and height
  dimensions. This input must be a named numeric vector.

- backbone_scale:

  Proportional scaling to the backbone length, width, and height
  dimensions. When a single value is supplied, all dimensions are scaled
  uniformly. Otherwise, this input must be a named numeric vector.

- backbone_target:

  Target dimensions (m) for the backbone length, width, and height
  dimensions. This input must be a named numeric vector.

- isometric_body:

  Logical; maintain isometric scaling for body.

- isometric_backbone:

  Logical; maintain isometric scaling for backbone.

- maintain_ratio:

  Maintain size ratio between body and backbone.

- n_segments_body:

  Number of points along the body profile.

- n_segments_backbone:

  Number of points along the backbone profile.

- containment:

  Containment policy for internal geometry checks. Use `"warn"` to keep
  the current warning behavior, `"error"` to fail fast for invalid
  internal geometries, or `"ignore"` to skip containment checks.

## Value

Modified BBF-class object.

## Examples

``` r
methods::selectMethod("reforge", "BBF")
#> new("MethodDefinition", .Data = function (object, ...) 
#> {
#>     .local <- function (object, body_scale = NULL, body_target = NULL, 
#>         backbone_scale = NULL, backbone_target = NULL, isometric_body = TRUE, 
#>         isometric_backbone = TRUE, maintain_ratio = TRUE, n_segments_body = NULL, 
#>         n_segments_backbone = NULL, containment = c("warn", "error", 
#>             "ignore")) 
#>     {
#>         containment <- match.arg(containment)
#>         if (is.null(body_scale) && is.null(backbone_scale) && 
#>             is.null(body_target) && is.null(backbone_target) && 
#>             is.null(n_segments_body) && is.null(n_segments_backbone)) {
#>             stop("Must specify at least one scaling, target, or segment count parameter.")
#>         }
#>         if (!is.null(body_scale) && !is.null(body_target)) {
#>             stop("Specify only one of body_scale or body_target, not both.")
#>         }
#>         if (!is.null(backbone_scale) && !is.null(backbone_target)) {
#>             stop("Specify only one of backbone_scale or backbone_target, not both.")
#>         }
#>         if ((!is.null(body_scale) || !is.null(body_target)) && 
#>             (!is.null(backbone_scale) || !is.null(backbone_target)) && 
#>             maintain_ratio) {
#>             maintain_ratio <- FALSE
#>             message("Hidden State Change Warning:\n ", "Multiple axes specified for the body and backbone: ", 
#>                 "'maintain_ratio' will be ignored for those axes.")
#>         }
#>         body <- acousticTS::extract(object, "body")
#>         backbone <- acousticTS::extract(object, "backbone")
#>         shape <- acousticTS::extract(object, "shape_parameters")
#>         rpos_b <- body$rpos
#>         rpos_bb <- backbone$rpos
#>         body_dims <- .reforge_component_dimensions(rpos_b, length = shape$body$length)
#>         backbone_dims <- .reforge_component_dimensions(rpos_bb, 
#>             length = shape$backbone$length %||% .shape_length(position_matrix = rpos_bb, 
#>                 row_major = TRUE))
#>         backbone_relative_start <- .reforge_relative_start(rpos_bb, 
#>             rpos_b)
#>         backbone_relative_vertical <- .reforge_relative_vertical_offset(rpos_bb, 
#>             rpos_b)
#>         body_target <- .validate_dimensions_target(body_target, 
#>             "body_target", c("length", "width", "height"))
#>         body_scale_lst <- .reforge_scale_vector(body_scale, body_target, 
#>             body_dims)
#>         backbone_target <- .validate_dimensions_target(backbone_target, 
#>             "backbone_target", c("length", "width", "height"))
#>         backbone_scale_lst <- .reforge_scale_vector(backbone_scale, 
#>             backbone_target, backbone_dims)
#>         body_scales <- .validate_dimension_scaling(dims = body_scale_lst$scale, 
#>             dims_name = paste0("body", body_scale_lst$suffix), 
#>             valid_dims = c("length", "width", "height"), isometry = isometric_body, 
#>             iso_name = "isometric_body")
#>         backbone_scales <- .validate_dimension_scaling(dims = backbone_scale_lst$scale, 
#>             dims_name = paste0("backbone", backbone_scale_lst$suffix), 
#>             valid_dims = c("length", "width", "height"), isometry = isometric_backbone, 
#>             iso_name = "isometric_backbone")
#>         if (maintain_ratio) {
#>             if (!is.null(body_scales) && is.null(backbone_scales)) {
#>                 backbone_scales <- body_scales
#>             }
#>             else if (is.null(body_scales) && !is.null(backbone_scales)) {
#>                 body_scales <- backbone_scales
#>             }
#>         }
#>         if (!is.null(backbone_scales)) {
#>             width_scale <- unname(backbone_scales["width"])
#>             height_scale <- unname(backbone_scales["height"])
#>             width_changed <- !isTRUE(all.equal(width_scale, 1))
#>             height_changed <- !isTRUE(all.equal(height_scale, 
#>                 1))
#>             if (width_changed && height_changed && !isTRUE(all.equal(width_scale, 
#>                 height_scale))) {
#>                 stop("Backbone reforge must preserve a circular cross-section: ", 
#>                   "width and height scaling must match.", call. = FALSE)
#>             }
#>             if (width_changed && !height_changed) {
#>                 backbone_scales["height"] <- width_scale
#>             }
#>             else if (!width_changed && height_changed) {
#>                 backbone_scales["width"] <- height_scale
#>             }
#>         }
#>         rpos_b <- .reforge_resample_rows(rpos_b, n_segments_body)
#>         rpos_bb <- .reforge_resample_rows(rpos_bb, n_segments_backbone)
#>         rpos_b <- .reforge_apply_axis_scaling(rpos_b, body_scales)
#>         rpos_bb <- .reforge_apply_axis_scaling(rpos_bb, backbone_scales)
#>         if (!is.null(body_scales) && body_scales["length"] != 
#>             1) {
#>             rpos_bb <- .reforge_shift_to_relative_start(rpos_bb, 
#>                 backbone_relative_start, rpos_b)
#>         }
#>         rpos_bb <- .reforge_shift_to_relative_vertical_offset(rpos_bb, 
#>             backbone_relative_vertical, rpos_b)
#>         .reforge_check_internal_containment(rpos_b, rpos_bb, 
#>             component_label = "Backbone", action = containment)
#>         body_radius <- .shape_radius_profile(position_matrix = rpos_b, 
#>             row_major = TRUE, error_context = "BBF body")
#>         backbone_radius <- .shape_radius_profile(position_matrix = rpos_bb, 
#>             row_major = TRUE, error_context = "BBF backbone")
#>         methods::slot(object, "body")$rpos <- rpos_b
#>         methods::slot(object, "body")$radius <- body_radius
#>         methods::slot(object, "backbone")$rpos <- rpos_bb
#>         methods::slot(object, "backbone")$radius <- backbone_radius
#>         methods::slot(object, "components")$backbone <- methods::slot(object, 
#>             "backbone")
#>         methods::slot(object, "shape_parameters")$body$length <- .shape_length(position_matrix = rpos_b, 
#>             row_major = TRUE)
#>         methods::slot(object, "shape_parameters")$body$radius <- if (all(is.na(body_radius))) {
#>             NA_real_
#>         }
#>         else {
#>             max(body_radius, na.rm = TRUE)
#>         }
#>         methods::slot(object, "shape_parameters")$body$n_segments <- ncol(rpos_b)
#>         methods::slot(object, "shape_parameters")$backbone$length <- .shape_length(position_matrix = rpos_bb, 
#>             row_major = TRUE)
#>         methods::slot(object, "shape_parameters")$backbone$radius <- max(backbone_radius, 
#>             na.rm = TRUE)
#>         methods::slot(object, "shape_parameters")$backbone$n_segments <- ncol(rpos_bb)
#>         return(object)
#>     }
#>     .local(object, ...)
#> }, target = new("signature", .Data = "BBF", names = "object", 
#>     package = "acousticTS"), defined = new("signature", .Data = "BBF", 
#>     names = "object", package = "acousticTS"), generic = "reforge")
#> <bytecode: 0x5587833d41f0>
#> <environment: namespace:acousticTS>
#> attr(,"target")
#> An object of class “signature”
#> object 
#>  "BBF" 
#> attr(,"defined")
#> An object of class “signature”
#> object 
#>  "BBF" 
#> attr(,"generic")
#> [1] "reforge"
#> attr(,"generic")attr(,"package")
#> [1] "acousticTS"
#> attr(,"class")
#> [1] "MethodDefinition"
#> attr(,"class")attr(,"package")
#> [1] "methods"
```
