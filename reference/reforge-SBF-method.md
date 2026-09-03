# Reforge SBF-class object

Resize a swimbladdered fish scatterer by rescaling the flesh body and
the gas-filled swimbladder together or independently. Body and
swimbladder geometry each support `length`, `width`, and `height`
scaling supplied as a direct `*_scale` factor or a `*_target` dimension
(m).

## Usage

``` r
# S4 method for class 'SBF'
reforge(
  object,
  body_scale = NULL,
  body_target = NULL,
  swimbladder_scale = NULL,
  swimbladder_target = NULL,
  isometric_body = TRUE,
  isometric_swimbladder = TRUE,
  maintain_ratio = TRUE,
  swimbladder_inflation_factor = 1,
  n_segments_body = NULL,
  n_segments_swimbladder = NULL,
  containment = c("warn", "error", "ignore")
)
```

## Arguments

- object:

  SBF-class object.

- body_scale:

  Proportional scaling to the body length, width, and height dimensions.
  When a single value is supplied, all dimensions are scaled using the
  same scaling factor. Otherwise, this input must be a named numeric
  vector.

- body_target:

  Target dimensions (m) for the body length, width, and height
  dimensions. This input must be a named numeric vector.

- swimbladder_scale:

  Proportional scaling to the swimbladder length, width, and height
  dimensions. When a single value is supplied, all dimensions are scaled
  using the same scaling factor. Otherwise, this input must be a named
  numeric vector.

- swimbladder_target:

  Target dimensions (m) for the swimbladder length, width, and height
  dimensions. This input must be a named numeric vector.

- isometric_body:

  Logical; maintain isometric scaling for body.

- isometric_swimbladder:

  Logical; maintain isometric scaling for bladder.

- maintain_ratio:

  Maintain size ratio between body and swimbladder.

- swimbladder_inflation_factor:

  Proportional swimbladder volume where the swimbladder x-axis origin
  and terminus are both held constant.

- n_segments_body:

  Number of segments along the body.

- n_segments_swimbladder:

  Number of segments along the bladder.

- containment:

  Containment policy for internal geometry checks. Use `"warn"` to keep
  the current warning behavior, `"error"` to fail fast for invalid
  internal geometries, or `"ignore"` to skip containment checks.

## Value

An updated object of class `SBF`.

## Details

Vertical scaling is applied about each component's own centerline so
that curved bodies keep their proportions when resized rather than
warping - a resized body remains geometrically similar whether it is
enlarged or shrunk. After scaling, the swimbladder is re-nested inside
the body: its relative along-body start and its relative vertical offset
within the body envelope are restored, so the bladder tracks the body
instead of drifting out of it. The `swimbladder_inflation_factor` scales
only the bladder width and height (its length endpoints are held fixed)
about the bladder's own center. The `containment` policy governs what
happens when the swimbladder exceeds the body envelope.

## Examples

``` r
methods::selectMethod("reforge", "SBF")
#> new("MethodDefinition", .Data = function (object, ...) 
#> {
#>     .local <- function (object, body_scale = NULL, body_target = NULL, 
#>         swimbladder_scale = NULL, swimbladder_target = NULL, 
#>         isometric_body = TRUE, isometric_swimbladder = TRUE, 
#>         maintain_ratio = TRUE, swimbladder_inflation_factor = 1, 
#>         n_segments_body = NULL, n_segments_swimbladder = NULL, 
#>         containment = c("warn", "error", "ignore")) 
#>     {
#>         containment <- match.arg(containment)
#>         if (is.null(body_scale) && is.null(swimbladder_scale) && 
#>             is.null(body_target) && is.null(swimbladder_target) && 
#>             is.null(n_segments_body) && is.null(n_segments_swimbladder) && 
#>             swimbladder_inflation_factor == 1) {
#>             stop("Must specify at least one scaling, target, inflation factor, ", 
#>                 "or segment count parameter.")
#>         }
#>         if (!is.null(body_scale) && !is.null(body_target)) {
#>             stop("Specify only one of body_scale or body_target, not both.")
#>         }
#>         if (!is.null(swimbladder_scale) && !is.null(swimbladder_target)) {
#>             stop("Specify only one of swimbladder_scale or swimbladder_target, not both.")
#>         }
#>         if ((!is.null(body_scale) || !is.null(body_target)) && 
#>             (!is.null(swimbladder_scale) || !is.null(swimbladder_target)) && 
#>             maintain_ratio) {
#>             maintain_ratio <- FALSE
#>             message("Hidden State Change Warning:\n ", "Multiple axes specified for the body and swimbladder: ", 
#>                 "'maintain_ratio' will be ignored for those axes.")
#>         }
#>         body <- acousticTS::extract(object, "body")
#>         bladder <- acousticTS::extract(object, "bladder")
#>         shape <- acousticTS::extract(object, "shape_parameters")
#>         rpos_b <- body$rpos
#>         rpos_sb <- bladder$rpos
#>         body_dims <- .reforge_component_dimensions(rpos_b, length = shape$body$length)
#>         bladder_dims <- .reforge_component_dimensions(rpos_sb, 
#>             length = shape$bladder$length %||% .shape_length(position_matrix = rpos_sb, 
#>                 row_major = TRUE))
#>         bladder_relative_start <- .reforge_relative_start(rpos_sb, 
#>             rpos_b)
#>         bladder_relative_vertical <- .reforge_relative_vertical_offset(rpos_sb, 
#>             rpos_b)
#>         body_target <- .validate_dimensions_target(body_target, 
#>             "body_target", c("length", "width", "height"))
#>         body_scale_lst <- .reforge_scale_vector(body_scale, body_target, 
#>             body_dims)
#>         swimbladder_target <- .validate_dimensions_target(swimbladder_target, 
#>             "swimbladder_target", c("length", "width", "height"))
#>         swimbladder_scale_lst <- .reforge_scale_vector(swimbladder_scale, 
#>             swimbladder_target, bladder_dims)
#>         body_scales <- .validate_dimension_scaling(dims = body_scale_lst$scale, 
#>             dims_name = paste0("body", body_scale_lst$suffix), 
#>             valid_dims = c("length", "width", "height"), isometry = isometric_body, 
#>             iso_name = "isometric_body")
#>         bladder_scales <- .validate_dimension_scaling(dims = swimbladder_scale_lst$scale, 
#>             dims_name = paste0("swimbladder", swimbladder_scale_lst$suffix), 
#>             valid_dims = c("length", "width", "height"), isometry = isometric_swimbladder, 
#>             iso_name = "isometric_swimbladder")
#>         if (maintain_ratio) {
#>             if (!is.null(body_scales) && is.null(bladder_scales)) {
#>                 bladder_scales <- body_scales
#>             }
#>             else if (is.null(body_scales) && !is.null(bladder_scales)) {
#>                 body_scales <- bladder_scales
#>             }
#>         }
#>         rpos_b <- .reforge_resample_rows(rpos_b, n_segments_body)
#>         rpos_sb <- .reforge_resample_rows(rpos_sb, n_segments_swimbladder)
#>         rpos_b <- .reforge_apply_axis_scaling(rpos_b, body_scales)
#>         rpos_sb <- .reforge_apply_axis_scaling(rpos_sb, bladder_scales)
#>         if (!is.null(body_scales) && body_scales["length"] != 
#>             1) {
#>             rpos_sb <- .reforge_shift_to_relative_start(rpos_sb, 
#>                 bladder_relative_start, rpos_b)
#>         }
#>         rpos_sb <- .reforge_shift_to_relative_vertical_offset(rpos_sb, 
#>             bladder_relative_vertical, rpos_b)
#>         if (swimbladder_inflation_factor != 1) {
#>             x_bladder_origin <- bladder$rpos[1, 1]/max(body$rpos[1, 
#>                 ])
#>             xsb_start <- x_bladder_origin * max(rpos_b[1, ])
#>             xsb_offset <- rpos_sb[1, 1] - xsb_start
#>             rpos_sb[1, ] <- rpos_sb[1, ] - xsb_offset
#>             rpos_sb[2, ] <- rpos_sb[2, ] * swimbladder_inflation_factor
#>             rpos_sb[3, ] <- rpos_sb[3, ] * swimbladder_inflation_factor
#>             if (nrow(rpos_sb) >= 4) {
#>                 rpos_sb[4, ] <- rpos_sb[4, ] * swimbladder_inflation_factor
#>             }
#>             rpos_sb <- .reforge_shift_to_relative_vertical_offset(rpos_sb, 
#>                 bladder_relative_vertical, rpos_b)
#>         }
#>         .reforge_check_internal_containment(rpos_b, rpos_sb, 
#>             action = containment)
#>         methods::slot(object, "body")$rpos <- rpos_b
#>         methods::slot(object, "bladder")$rpos <- rpos_sb
#>         methods::slot(object, "shape_parameters")$body$length <- .shape_length(position_matrix = rpos_b, 
#>             row_major = TRUE)
#>         methods::slot(object, "shape_parameters")$bladder$length <- .shape_length(position_matrix = rpos_sb, 
#>             row_major = TRUE)
#>         methods::slot(object, "shape_parameters")$body$n_segments <- ncol(rpos_b)
#>         methods::slot(object, "shape_parameters")$bladder$n_segments <- ncol(rpos_sb)
#>         return(object)
#>     }
#>     .local(object, ...)
#> }, target = new("signature", .Data = "SBF", names = "object", 
#>     package = "acousticTS"), defined = new("signature", .Data = "SBF", 
#>     names = "object", package = "acousticTS"), generic = "reforge")
#> <bytecode: 0x5571fd86c2f0>
#> <environment: namespace:acousticTS>
#> attr(,"target")
#> An object of class “signature”
#> object 
#>  "SBF" 
#> attr(,"defined")
#> An object of class “signature”
#> object 
#>  "SBF" 
#> attr(,"generic")
#> [1] "reforge"
#> attr(,"generic")attr(,"package")
#> [1] "acousticTS"
#> attr(,"class")
#> [1] "MethodDefinition"
#> attr(,"class")attr(,"package")
#> [1] "methods"
```
