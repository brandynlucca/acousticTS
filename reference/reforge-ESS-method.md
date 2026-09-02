# Reforge ESS-class object

Resize an elastic-shelled scatterer by applying an isometric scale
factor or specifying a target maximum shell radius. Shell thickness can
be updated independently, which rescales the fluid body so that the
maximum fluid radius equals `new_max_shell_radius - shell_thickness`
(matching the convention in
[`ess_generate`](https://brandynlucca.github.io/acousticTS/reference/ess_generate.md)).
The underlying shape (sphere, prolate spheroid, cylinder, etc.) is
preserved for both shell and fluid bodies; scaling is applied uniformly
to all axes of each position matrix.

## Usage

``` r
# S4 method for class 'ESS'
reforge(
  object,
  scale = NULL,
  radius_target = NULL,
  shell_thickness = NULL,
  n_segments = NULL
)
```

## Arguments

- object:

  ESS-class object.

- scale:

  Single positive scalar applied isometrically to the shell (and fluid
  body, if present). Mutually exclusive with `radius_target`.

- radius_target:

  Target *maximum* outer shell radius (m). Scale factor derived as
  `radius_target / max(current_shell_radius)`. Mutually exclusive with
  `scale`.

- shell_thickness:

  New shell wall thickness (m). The fluid body is rescaled so its
  maximum radius equals `new_max_shell_radius - shell_thickness`. Can be
  combined with `scale`/`radius_target` or used alone.

- n_segments:

  New number of discrete segments. All columns of both the shell and
  fluid position matrices are re-interpolated along the x-axis.

## Value

Modified ESS-class object.

## Examples

``` r
methods::selectMethod("reforge", "ESS")
#> new("MethodDefinition", .Data = function (object, ...) 
#> {
#>     .local <- function (object, scale = NULL, radius_target = NULL, 
#>         shell_thickness = NULL, n_segments = NULL) 
#>     {
#>         if (is.null(scale) && is.null(radius_target) && is.null(shell_thickness) && 
#>             is.null(n_segments)) {
#>             stop(paste0("Must specify at least one of: scale, radius_target, ", 
#>                 "shell_thickness, or n_segments."), call. = FALSE)
#>         }
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
#>         if (!is.null(shell_thickness) && (!is.numeric(shell_thickness) || 
#>             length(shell_thickness) != 1 || shell_thickness <= 
#>             0)) {
#>             stop("'shell_thickness' must be a single positive number.", 
#>                 call. = FALSE)
#>         }
#>         if (!is.null(n_segments) && (!is.numeric(n_segments) || 
#>             length(n_segments) != 1 || n_segments < 1)) {
#>             stop("'n_segments' must be a single positive integer.", 
#>                 call. = FALSE)
#>         }
#>         shell <- acousticTS::extract(object, "shell")
#>         fluid <- acousticTS::extract(object, "fluid")
#>         shape <- acousticTS::extract(object, "shape_parameters")
#>         rpos_shell <- shell$rpos
#>         rpos_fluid <- fluid$rpos
#>         curr_shell_r <- shape$shell$radius
#>         curr_shell_max_r <- max(curr_shell_r, na.rm = TRUE)
#>         curr_fluid_r <- shape$fluid$radius
#>         curr_fluid_max_r <- if (!is.null(curr_fluid_r) && !all(is.na(curr_fluid_r))) {
#>             max(curr_fluid_r, na.rm = TRUE)
#>         }
#>         else {
#>             NA_real_
#>         }
#>         if (!is.null(radius_target)) 
#>             scale <- radius_target/curr_shell_max_r
#>         if (!is.null(n_segments)) {
#>             n_new <- as.integer(n_segments) + 1L
#>             rpos_shell <- .resample_rpos(rpos_shell, n_new)
#>             if (!is.null(rpos_fluid)) {
#>                 rpos_fluid <- .resample_rpos(rpos_fluid, n_new)
#>             }
#>             methods::slot(object, "shape_parameters")$n_segments <- as.integer(n_segments)
#>         }
#>         if (!is.null(scale)) {
#>             rpos_shell <- rpos_shell * scale
#>             new_shell_r <- curr_shell_r * scale
#>             new_shell_max_r <- curr_shell_max_r * scale
#>             methods::slot(object, "shell")$radius <- new_shell_r
#>             methods::slot(object, "shape_parameters")$shell$radius <- new_shell_r
#>             methods::slot(object, "shape_parameters")$shell$diameter <- new_shell_r * 
#>                 2
#>             if (!is.null(rpos_fluid) && !is.na(curr_fluid_max_r)) {
#>                 if (!is.null(shell_thickness)) {
#>                   new_fluid_max_r <- new_shell_max_r - shell_thickness
#>                   if (new_fluid_max_r <= 0) {
#>                     stop("shell_thickness exceeds new shell radius.", 
#>                       call. = FALSE)
#>                   }
#>                   fluid_scale <- new_fluid_max_r/curr_fluid_max_r
#>                   rpos_fluid <- rpos_fluid * fluid_scale
#>                   new_fluid_r <- curr_fluid_r * fluid_scale
#>                   methods::slot(object, "shell")$shell_thickness <- shell_thickness
#>                 }
#>                 else {
#>                   rpos_fluid <- rpos_fluid * scale
#>                   new_fluid_r <- curr_fluid_r * scale
#>                 }
#>                 methods::slot(object, "fluid")$radius <- new_fluid_r
#>                 methods::slot(object, "shape_parameters")$fluid$radius <- new_fluid_r
#>                 methods::slot(object, "shape_parameters")$fluid$diameter <- new_fluid_r * 
#>                   2
#>             }
#>         }
#>         else if (!is.null(shell_thickness)) {
#>             new_fluid_max_r <- curr_shell_max_r - shell_thickness
#>             if (new_fluid_max_r <= 0) {
#>                 stop("shell_thickness exceeds shell radius.", 
#>                   call. = FALSE)
#>             }
#>             if (!is.null(rpos_fluid) && !is.na(curr_fluid_max_r)) {
#>                 fluid_scale <- new_fluid_max_r/curr_fluid_max_r
#>                 rpos_fluid <- rpos_fluid * fluid_scale
#>                 new_fluid_r <- curr_fluid_r * fluid_scale
#>                 methods::slot(object, "fluid")$radius <- new_fluid_r
#>                 methods::slot(object, "shape_parameters")$fluid$radius <- new_fluid_r
#>                 methods::slot(object, "shape_parameters")$fluid$diameter <- new_fluid_r * 
#>                   2
#>             }
#>             methods::slot(object, "shell")$shell_thickness <- shell_thickness
#>         }
#>         methods::slot(object, "shell")$rpos <- rpos_shell
#>         if (!is.null(rpos_fluid)) {
#>             methods::slot(object, "fluid")$rpos <- rpos_fluid
#>         }
#>         return(object)
#>     }
#>     .local(object, ...)
#> }, target = new("signature", .Data = "ESS", names = "object", 
#>     package = "acousticTS"), defined = new("signature", .Data = "ESS", 
#>     names = "object", package = "acousticTS"), generic = "reforge")
#> <bytecode: 0x564373f0b178>
#> <environment: namespace:acousticTS>
#> attr(,"target")
#> An object of class “signature”
#> object 
#>  "ESS" 
#> attr(,"defined")
#> An object of class “signature”
#> object 
#>  "ESS" 
#> attr(,"generic")
#> [1] "reforge"
#> attr(,"generic")attr(,"package")
#> [1] "acousticTS"
#> attr(,"class")
#> [1] "MethodDefinition"
#> attr(,"class")attr(,"package")
#> [1] "methods"
```
