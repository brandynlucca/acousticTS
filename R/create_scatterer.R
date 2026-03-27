################################################################################
# CREATE SCATTERER FUNCTIONS
################################################################################
################################################################################
# Create SBF-class object
################################################################################
#' Generate a SBF-class object.
#'
#' @param x_body Vector containing along-body axis (m).
#' @param w_body Vector containing across-body axis (m).
#' @param zU_body Vector containing dorsal-body axis (m).
#' @param zL_body Vector containing ventral-body axis (m).
#' @param x_bladder Vector containing along-bladder axis (m).
#' @param w_bladder Vector containing across-bladder axis (m).
#' @param zU_bladder Vector containing dorsal-bladder axis (m).
#' @param zL_bladder Vector containing ventral-bladder axis (m).
#' @param density_body Flesh density
#' (\ifelse{html}{\out{&rho;<sub>body</sub>}}{\eqn{\rho_{body}}},
#' kg \ifelse{html}{\out{m<sup>3</sup>}}{\eqn{m^3}}).
#' @param sound_speed_body Flesh sound speed
#' (\ifelse{html}{\out{c;<sub>body</sub>}}{\eqn{c_{body}}},
#' m \ifelse{html}{\out{s<sup>-1</sup>}}{\eqn{s^{-1}}}).
#' @param density_bladder Bladder density (\eqn{\rho},
#' kg \ifelse{html}{\out{m<sup>3</sup>}}{\eqn{m^3}}).
#' @param sound_speed_bladder Bladder sound speed (c, m \eqn{s^-1}).
#' @param g_body Body density contrast.
#' @param h_body Body sound speed contrast.
#' @param g_bladder Bladder density contrast.
#' @param h_bladder Bladder sound speed contrast.
#' @param theta_body Angle of body relative to wavefront
#' (\eqn{\theta_body}, radians).
#' @param theta_bladder Angle of bladder relative to wavefront
#' (\eqn{\theta_bladder}, radians).
#' @param theta_units Compatibility argument. Scatterer constructors now assume
#'   radians and ignore non-SI alternatives.
#' @param length_units Compatibility argument. Scatterer constructors now assume
#'   meters and ignore non-SI alternatives.
#' @param ID Optional metadata identifier.
#' @param body_shape Optional pre-built Shape for the body; overrides
#'   x_body/w_body/zU_body/zL_body if supplied.
#' @param bladder_shape Optional pre-built Shape for the bladder; overrides
#'   x_bladder/w_bladder/zU_bladder/zL_bladder if supplied.
#' @details
#' The recommended interface is to supply pre-built `Shape` objects through
#' `body_shape` and `bladder_shape`, then specify material properties through
#' either contrasts (`g_*`, `h_*`) or absolute density/sound-speed values.
#' Legacy coordinate-vector inputs are retained for backward compatibility and
#' are converted internally to `Shape` objects before the scatterer is built.
#' Character-based shape dispatch is deprecated; build the component geometry
#' first with a `Shape` constructor and then pass the resulting `Shape` object.
#'
#' Body/bladder width vectors (`w_body`, `w_bladder`) default to zeros when
#' missing. Resulting `rpos` matrices always carry the row names
#' x_body, w_body, zU_body, zL_body for the body and x_bladder, w_bladder,
#' zU_bladder, zL_bladder for the bladder.
#' @examples
#' # Manual body/bladder coordinates
#' sbf_generate(
#'   x_body = seq(0, 0.1, length.out = 5),
#'   w_body = rep(0, 5),
#'   zU_body = seq(0.001, 0.002, length.out = 5),
#'   zL_body = -seq(0.001, 0.002, length.out = 5),
#'   x_bladder = seq(0.02, 0.09, length.out = 4),
#'   w_bladder = rep(0, 4),
#'   zU_bladder = rep(0.0015, 4),
#'   zL_bladder = rep(-0.0015, 4),
#'   sound_speed_body = 1500,
#'   sound_speed_bladder = 340,
#'   density_body = 1040,
#'   density_bladder = 1.2
#' )
#'
#' # Using pre-built shapes
#' body_shape <- arbitrary(
#'   x_body = c(0, 0.1), zU_body = c(0.001, 0.002),
#'   zL_body = c(-0.001, -0.002)
#' )
#' bladder_shape <- arbitrary(
#'   x_bladder = c(0.02, 0.09), w_bladder = c(0, 0),
#'   zU_bladder = c(0.0015, 0.0015),
#'   zL_bladder = c(-0.0015, -0.0015)
#' )
#' sbf_generate(
#'   body_shape = body_shape, bladder_shape = bladder_shape,
#'   sound_speed_body = 1500, sound_speed_bladder = 340,
#'   density_body = 1040, density_bladder = 1.2
#' )
#'
#' @return
#' Generates a SBF-class object.
#'
#' @seealso \code{\link{SBF}}
#'
#' @importFrom methods new
#' @keywords scatterer_type_generation
#' @export
sbf_generate <- function(x_body = NULL,
                         w_body = NULL,
                         zU_body = NULL,
                         zL_body = NULL,
                         x_bladder = NULL,
                         w_bladder = NULL,
                         zU_bladder = NULL,
                         zL_bladder = NULL,
                         sound_speed_body = NULL,
                         sound_speed_bladder = NULL,
                         g_body = NULL,
                         h_body = NULL,
                         g_bladder = NULL,
                         h_bladder = NULL,
                         density_body = NULL,
                         density_bladder = NULL,
                         theta_body = pi / 2,
                         theta_bladder = pi / 2,
                         theta_units = "radians",
                         length_units = "m",
                         ID = NULL,
                         body_shape = NULL,
                         bladder_shape = NULL) {
  # Normalize constructor units ===============================================
  units <- .normalize_scatterer_units(
    theta_units = theta_units,
    length_units = length_units,
    context = "SBF"
  )
  # Validate body and bladder material inputs =================================
  .validate_component_material_inputs(
    g = g_body,
    density = density_body,
    h = h_body,
    sound_speed = sound_speed_body,
    component_label = "body"
  )
  .validate_component_material_inputs(
    g = g_bladder,
    density = density_bladder,
    h = h_bladder,
    sound_speed = sound_speed_bladder,
    component_label = "bladder"
  )
  # Resolve body and bladder shape inputs =====================================
  body_shape_obj <- .resolve_profile_component_shape(
    shape_obj = body_shape,
    suffix = "body",
    x = x_body,
    w = w_body,
    zU = zU_body,
    zL = zL_body
  )
  bladder_shape_obj <- .resolve_profile_component_shape(
    shape_obj = bladder_shape,
    suffix = "bladder",
    x = x_bladder,
    w = w_bladder,
    zU = zU_bladder,
    zL = zL_bladder
  )
  # Create body slot ===========================================================
  body <- .build_fluid_profile_component(
    shape_obj = body_shape_obj,
    suffix = "body",
    theta = theta_body,
    density = density_body,
    sound_speed = sound_speed_body,
    g = g_body,
    h = h_body
  )
  # Create bladder slot ========================================================
  bladder <- .build_fluid_profile_component(
    shape_obj = bladder_shape_obj,
    suffix = "bladder",
    theta = theta_bladder,
    density = density_bladder,
    sound_speed = sound_speed_bladder,
    g = g_bladder,
    h = h_bladder
  )
  # Define shape parameters ====================================================
  shape_parameters <- list(
    body = .shape_common_parameters(
      shape_input = body_shape_obj,
      error_context = "SBF body"
    ),
    bladder = .shape_common_parameters(
      shape_input = bladder_shape_obj,
      error_context = "SBF bladder"
    ),
    length_units = units$length_units,
    theta_units = units$theta_units
  )
  # Create metadata field ======================================================
  metadata <- .scatterer_metadata(ID = ID)
  # Create SBF-class object ====================================================
  return(methods::new("SBF",
    metadata = metadata,
    model_parameters = list(),
    model = list(),
    body = body,
    bladder = bladder,
    shape_parameters = shape_parameters
  ))
}
################################################################################
# Create BBF-class object
################################################################################
#' Generate a BBF-class object.
#'
#' @param body_shape Pre-built body shape. Must inherit from \code{Shape}.
#' @param backbone_shape Pre-built backbone shape. Must be a cylindrical
#'   \code{Shape}.
#' @param density_body Flesh density
#'   (\ifelse{html}{\out{&rho;<sub>body</sub>}}{\eqn{\rho_{body}}},
#'   kg \ifelse{html}{\out{m<sup>3</sup>}}{\eqn{m^3}}).
#' @param sound_speed_body Flesh sound speed
#'   (\ifelse{html}{\out{c<sub>body</sub>}}{\eqn{c_{body}}},
#'   m \ifelse{html}{\out{s<sup>-1</sup>}}{\eqn{s^{-1}}}).
#' @param g_body Body density contrast.
#' @param h_body Body sound speed contrast.
#' @param density_backbone Backbone density
#'   (\ifelse{html}{\out{&rho;<sub>bb</sub>}}{\eqn{\rho_{bb}}},
#'   kg \ifelse{html}{\out{m<sup>3</sup>}}{\eqn{m^3}}).
#' @param sound_speed_longitudinal_backbone Longitudinal wave speed in the
#'   backbone (m/s).
#' @param sound_speed_transversal_backbone Transversal wave speed in the
#'   backbone (m/s).
#' @param theta_body Body orientation relative to the incident wave (radians).
#' @param theta_backbone Backbone orientation relative to the incident wave
#'   (radians).
#' @param x_offset_backbone Along-body translation applied to the backbone
#'   geometry (m).
#' @param z_offset_backbone Dorsoventral translation applied to the backbone
#'   geometry (m).
#' @param theta_units Compatibility argument. Scatterer constructors now assume
#'   radians and ignore non-SI alternatives.
#' @param length_units Compatibility argument. Scatterer constructors now assume
#'   meters and ignore non-SI alternatives.
#' @param ID Optional metadata identifier.
#'
#' @details
#' \code{bbf_generate()} is intended for swimbladder-less fish workflows where
#' the flesh and backbone should remain explicit, separately parameterized
#' components. The body is stored using the same segmented-body representation
#' used by \code{\link{FLS}}, while the backbone is stored as an explicit
#' cylindrical component carrying elastic material properties for use by
#' cylinder-based modal models.
#'
#' @examples
#' body_shape <- arbitrary(
#'   x_body = c(0, 0.04, 0.08),
#'   zU_body = c(0.001, 0.004, 0.001),
#'   zL_body = c(-0.001, -0.004, -0.001)
#' )
#' backbone_shape <- cylinder(
#'   length_body = 0.06,
#'   radius_body = 0.0008,
#'   n_segments = 40
#' )
#' bbf_generate(
#'   body_shape = body_shape,
#'   backbone_shape = backbone_shape,
#'   density_body = 1070,
#'   sound_speed_body = 1570,
#'   density_backbone = 1900,
#'   sound_speed_longitudinal_backbone = 3500,
#'   sound_speed_transversal_backbone = 1700
#' )
#'
#' @return
#' Generates a BBF-class object.
#'
#' @seealso \code{\link{BBF}}, \code{\link{FLS}}, \code{\link{cylinder}}
#'
#' @importFrom methods new
#' @keywords scatterer_type_generation
#' @export
bbf_generate <- function(body_shape,
                         backbone_shape,
                         density_body = NULL,
                         sound_speed_body = NULL,
                         g_body = NULL,
                         h_body = NULL,
                         density_backbone,
                         sound_speed_longitudinal_backbone,
                         sound_speed_transversal_backbone,
                         theta_body = pi / 2,
                         theta_backbone = pi / 2,
                         x_offset_backbone = 0,
                         z_offset_backbone = 0,
                         theta_units = "radians",
                         length_units = "m",
                         ID = NULL) {
  # Normalize constructor units ================================================
  units <- .normalize_scatterer_units(
    theta_units = theta_units,
    length_units = length_units,
    context = "BBF"
  )
  # Validate body and backbone shape inputs ====================================
  if (!methods::is(body_shape, "Shape")) {
    stop("'body_shape' must be a pre-built Shape object.", call. = FALSE)
  }
  if (!methods::is(backbone_shape, "Cylinder")) {
    stop(
      "'backbone_shape' must be a cylindrical Shape object generated by ",
      "`cylinder()` or a compatible subclass.",
      call. = FALSE
    )
  }
  # Validate body material properties ==========================================
  .validate_component_material_inputs(
    g = g_body,
    density = density_body,
    h = h_body,
    sound_speed = sound_speed_body,
    component_label = "body"
  )
  # Create body slot ===========================================================
  body_pos <- acousticTS::extract(body_shape, "position_matrix")
  body_shape_params <- acousticTS::extract(body_shape, "shape_parameters")
  body_radius_profile <- .shape_radius_profile(
    position_matrix = body_pos,
    shape_parameters = body_shape_params,
    error_context = "BBF body"
  )
  body <- .build_row_major_fluid_component(
    shape_obj = body_shape,
    theta = theta_body,
    density = density_body,
    sound_speed = sound_speed_body,
    g = g_body,
    h = h_body,
    radius = body_radius_profile
  )
  # Create backbone slot =======================================================
  backbone_pos <- .translate_shape_position_matrix(
    acousticTS::extract(backbone_shape, "position_matrix"),
    x_offset = x_offset_backbone,
    z_offset = z_offset_backbone
  )
  backbone_shape_params <- acousticTS::extract(
    backbone_shape,
    "shape_parameters"
  )
  backbone <- .build_row_major_elastic_component(
    shape_obj = backbone_shape,
    theta = theta_backbone,
    density = density_backbone,
    sound_speed_longitudinal = sound_speed_longitudinal_backbone,
    sound_speed_transversal = sound_speed_transversal_backbone,
    radius = backbone_shape_params$radius,
    position_matrix_override = backbone_pos
  )
  # Define shape parameters ====================================================
  shape_parameters <- list(
    body = .append_shape_specific_parameters(
      list(
        length = .shape_length(position_matrix = body_pos),
        radius = if (all(is.na(body_radius_profile))) {
          NA_real_
        } else {
          max(body_radius_profile, na.rm = TRUE)
        },
        n_segments = .shape_segment_count(body_pos)
      ),
      body_shape
    ),
    backbone = .append_shape_specific_parameters(
      list(
        length = .shape_length(position_matrix = backbone_pos),
        radius = max(backbone_shape_params$radius, na.rm = TRUE),
        n_segments = .shape_segment_count(backbone_pos)
      ),
      backbone_shape
    ),
    length_units = units$length_units,
    theta_units = units$theta_units
  )
  # Create metadata field ======================================================
  metadata <- .scatterer_metadata(ID = ID)
  # Create BBF-class object ====================================================
  return(methods::new("BBF",
    metadata = metadata,
    model_parameters = list(),
    model = list(),
    body = body,
    backbone = backbone,
    shape_parameters = shape_parameters,
    components = list(backbone = backbone)
  ))
}
################################################################################
# Create CAL-class object
################################################################################
#' Generate a CAL-class object.
#' @param material Material-type for a solid sphere. See `Details` for available
#' options. Default is tungsten carbide (WC).
#' @param diameter Spherical diameter (m).
#' @param n_segments Number of segments to discretize object shape.
#' @param sound_speed_longitudinal Longitudinal sound speed (m/s).
#' @param sound_speed_transversal Transversal sound speed (m/s).
#' @param density_sphere Density (kg/m^3).
#' @param theta_sphere Backscattering direction (Default: pi radians).
#' @param ID Optional metadata ID input.
#' @param diameter_units Compatibility argument. `cal_generate()` now assumes
#'   meters and ignores non-SI alternatives.
#' @param theta_units Compatibility argument. `cal_generate()` now assumes
#'   radians and ignores non-SI alternatives.
#' @param material Material-type for the soldi sphere. See 'Details' built-in
#' material options.
#' @examples
#' cal_generate(material = "WC", diameter = 38.1e-3, n_segments = 120)
#'
#' @details
#' There are several options for the \strong{material} argument:
#' \tabular{rlllll}{
#'  \strong{Material} \tab \strong{Argument} \tab \strong{c1} \tab \strong{c2}
#'  \tab \strong{\eqn{\rho1}}\cr
#'  \emph{Tungsten carbide} \tab "WC" \tab 6853 \tab 4171 \tab 14900\cr
#'  \emph{Stainless steel} \tab "steel" \tab 5980 \tab 3297 \tab 7970\cr
#'  \emph{Brass} \tab "brass" \tab 4372 \tab 2100 \tab 8360\cr
#'  \emph{Copper} \tab "Cu" \tab 4760 \tab 2288.5 \tab 8947\cr
#'  \emph{Aluminum} \tab "Al" \tab 6260 \tab 3080 \tab 2700\cr
#' }
#' @return
#' Generates a CAL-class object.
#'
#' @seealso \code{\link{CAL}}
#'
#' @importFrom methods new
#' @keywords scatterer_type_generation
#' @export
cal_generate <- function(material = "WC",
                         diameter = 38.1e-3,
                         sound_speed_longitudinal = NULL,
                         sound_speed_transversal = NULL,
                         density_sphere = NULL,
                         theta_sphere = pi,
                         ID = NULL,
                         diameter_units = "m",
                         theta_units = "radians",
                         n_segments = 1e2) {
  units <- .normalize_scatterer_units(
    theta_units = theta_units,
    diameter_units = diameter_units,
    context = "CAL"
  )
  metadata <- .scatterer_metadata(
    ID = ID,
    default_id = "Calibration sphere",
    extra = list(Material = material)
  )
  # Create sphere object to define definitions =================================
  sphere_shape <- sphere(
    radius_body = diameter / 2,
    n_segments = n_segments,
    diameter_units = units$diameter_units
  )
  # Define calibration sphere body shape =======================================
  body <- list(
    rpos = sphere_shape@position_matrix,
    diameter = diameter,
    radius = diameter / 2,
    theta = theta_sphere
  )
  # Define material properties =================================================
  material_properties <- switch(material,
    Cu = list(
      sound_speed_longitudinal = 4760,
      sound_speed_transversal = 2288.5,
      density = 8947
    ),
    WC = list(
      sound_speed_longitudinal = 6853,
      sound_speed_transversal = 4171,
      density = 14900
    ),
    Al = list(
      sound_speed_longitudinal = 6260,
      sound_speed_transversal = 3080,
      density = 2700
    ),
    steel = list(
      sound_speed_longitudinal = 5610,
      sound_speed_transversal = 3120,
      density = 7800
    ),
    brass = list(
      sound_speed_longitudinal = 4372,
      sound_speed_transversal = 2100,
      density = 8360
    )
  )
  if (!is.null(sound_speed_longitudinal)) {
    material_properties$sound_speed_longitudinal <- sound_speed_longitudinal
  }
  if (!is.null(sound_speed_transversal)) {
    material_properties$sound_speed_transversal <- sound_speed_transversal
  }
  if (!is.null(density_sphere)) {
    material_properties$density <- density_sphere
  }
  # Append material properties to the shape body ===============================
  body <- append(
    body,
    material_properties
  )
  # Define shape parameters ====================================================
  shape_parameters <- list(
    diameter = diameter,
    radius_body = diameter / 2,
    n_segments = n_segments,
    diameter_units = units$diameter_units,
    theta_units = units$theta_units
  )
  # Generate calibration sphere object =========================================
  return(methods::new("CAL",
    metadata = metadata,
    model_parameters = list(),
    model = list(),
    body = body,
    shape_parameters = shape_parameters
  ))
}
################################################################################
# Create FLS-class object
################################################################################
#' Generate a FLS-class object.
#' @param shape Pre-built `Shape` object describing the target geometry. If
#'   omitted, explicit profile coordinates such as `x_body`, `y_body`, and
#'   `z_body` are treated as the manual geometry pathway. Legacy character
#'   dispatch such as `"sphere"` or `"arbitrary"` is retained only for backward
#'   compatibility and is now deprecated.
#' @param x_body Vector containing x-axis body (m) shape data.
#' @param y_body Vector containing y-axis body (m) shape data.
#' @param z_body Vector containing z-axis body (m) shape data.
#' @param length_body Optional input for a generic length value input.
#' @param radius_body Vector containing radii (m).
#' @param n_segments Number of body segments.
#' @param radius_curvature_ratio Length-to-curvature ratio (pc/L).
#' @param g_body Density contrast. This can either be a single value (i.e.
#' homogenous) or a vector of values (i.e. inhomogenous).
#' @param density_body Absolute density (kg/m^3) if contrasts are not supplied.
#' @param h_body Soundspeed contrast. This can either be a single value (i.e.
#' homogenous) or a vector of values (i.e. inhomogenous).
#' @param sound_speed_body Absolute sound speed (m/s) if contrasts are not
#' supplied.
#' @param theta_body Orientation of the target relative to the transmit source
#' (\eqn{\theta}). Broadside incidence is considered 90 degrees, or pi/2.
#' Default value is pi/2; input should be in radians.
#' @param theta_units Compatibility argument. Scatterer constructors now assume
#'   radians and ignore non-SI alternatives.
#' @param length_units Compatibility argument. Scatterer constructors now assume
#'   meters and ignore non-SI alternatives.
#' @param ID Optional metadata entry.
#' @param ... Additional parameters.
#' @examples
#' shape <- prolate_spheroid(
#'   length_body = 0.04, radius_body = 0.004, n_segments = 50
#' )
#' fls_generate(
#'   shape = shape, density_body = 1045, sound_speed_body = 1520
#' )
#' @return
#' FLS-class object
#'
#' @details
#' The preferred workflow is to build a geometry first with `sphere()`,
#' `cylinder()`, `prolate_spheroid()`, or `arbitrary()`, then pass that `Shape`
#' object to `fls_generate()`. Material properties can be supplied either as
#' contrasts (`g_body`/`h_body`) or as absolute density/sound-speed values
#' (`density_body`/`sound_speed_body`), but not both for the same property
#' pair. Downstream models derive contrasts automatically when only absolute
#' values are supplied.
#'
#' The only supported public geometry paths are now:
#' 1. supply a pre-built `Shape` object, or
#' 2. supply explicit profile coordinates directly to the constructor.
#'
#' Character-based shape dispatch is retained only as a compatibility pathway
#' and is now deprecated. Internally, every pathway is resolved to the same
#' `Shape`-first geometry contract before the `FLS` object is built.
#'
#' Scatterer constructors store geometry in meters and orientations in radians.
#' `length_units` and `theta_units` are retained as compatibility arguments, but
#' non-SI values are normalized to the package-standard representation.
#'
#' @seealso \code{\link{FLS}}
#'
#' @importFrom methods new
#' @keywords scatterer_type_generation
#' @export
fls_generate <- function(shape = NULL,
                         x_body = NULL,
                         y_body = NULL,
                         z_body = NULL,
                         length_body = NULL,
                         radius_body = NULL,
                         radius_curvature_ratio = NULL,
                         n_segments = 18,
                         g_body = NULL,
                         h_body = NULL,
                         density_body = NULL,
                         sound_speed_body = NULL,
                         theta_body = pi / 2,
                         ID = NULL,
                         length_units = "m",
                         theta_units = "radians", ...) {
  arg_pull <- lapply(as.list(match.call()), eval, parent.frame())
  # Normalize constructor units ================================================
  units <- .normalize_scatterer_units(
    theta_units = theta_units,
    length_units = length_units,
    context = "FLS"
  )
  # Validate shape input pathway ===============================================
  if (!methods::is(shape, "Shape") &&
    (is.null(shape) || identical(shape, "arbitrary")) &&
    !.has_explicit_profile_coordinates(arg_pull)) {
    stop(
      "Supply 'shape' as a pre-built Shape object, or provide explicit ",
      "profile coordinates such as x_body/y_body/z_body.",
      call. = FALSE
    )
  }
  # Resolve body shape input ===================================================
  shape_input <- .resolve_scatterer_shape_input(shape, arg_pull)
  pos_mat <- acousticTS::extract(shape_input, "position_matrix")
  shape_params <- acousticTS::extract(shape_input, "shape_parameters")
  shape_parameters <- .shape_common_parameters(
    shape_input = shape_input,
    requested_shape = shape,
    error_context = "FLS body",
    extra_units = list(
      length_units = units$length_units,
      theta_units = units$theta_units
    )
  )
  # Validate material properties ==============================================
  .validate_component_material_inputs(
    g = g_body,
    density = density_body,
    h = h_body,
    sound_speed = sound_speed_body,
    component_label = "body"
  )
  .validate_fls_segmented_material_lengths(
    shape_obj = shape_input,
    g = g_body,
    density = density_body,
    h = h_body,
    sound_speed = sound_speed_body
  )
  # Create body slot ===========================================================
  body <- .build_row_major_fluid_component(
    shape_obj = shape_input,
    theta = theta_body,
    density = density_body,
    sound_speed = sound_speed_body,
    g = g_body,
    h = h_body,
    radius = shape_params$radius,
    extra = list(radius_curvature_ratio = radius_curvature_ratio)
  )
  # Create metadata field ======================================================
  metadata <- .scatterer_metadata(ID = ID)
  # Create FLS-class object ====================================================
  return(methods::new("FLS",
    metadata = metadata,
    model_parameters = list(),
    model = list(),
    body = body,
    shape_parameters = shape_parameters
  ))
}
################################################################################
# Create GAS-class object
################################################################################
#' Generate a GAS-class object
#'
#' @param shape Pre-built `Shape` object describing the gas-filled geometry. If
#'   omitted, explicit profile coordinates such as `x_body`, `y_body`, and
#'   `z_body` are treated as the manual geometry pathway. Legacy character
#'   dispatch such as `"sphere"` is retained only for backward compatibility
#'   and is now deprecated.
#' @param radius_body Radius (m). For non-canonical shapes, this would be the
#' maximum or mean radius at the scatterer midsection.
#' @param h_fluid Sound speed contrast of fluid relative to surrounding
#' medium (h).
#' @param g_fluid Density contrast of fluid relative to surrounding density (g).
#' @param sound_speed_fluid Optional fluid sound speed (m/s).
#' @param density_fluid Optional fluid density (m/s).
#' @param theta_body Orientation of the target relative to the incident wave
#'   (radians).
#' @param ID Optional metadata identifier.
#' @param radius_units Compatibility argument. `gas_generate()` now assumes
#'   meters and ignores non-SI alternatives.
#' @param theta_units Compatibility argument. Scatterer constructors now assume
#'   radians and ignore non-SI alternatives.
#' @param n_segments Number of body segments.
#' @usage gas_generate(shape = NULL, radius_body = NULL, h_fluid = 0.22,
#'   g_fluid = 0.0012, sound_speed_fluid = NULL, density_fluid = NULL,
#'   theta_body = pi/2, ID = NULL, radius_units = "m",
#'   theta_units = "radians", n_segments = 100, ...)
#' @param ... Additional manual profile arguments or legacy canonical shape
#'   arguments used by the compatibility geometry pathway, such as `x_body`,
#'   `y_body`, `z_body`, `radius_body`, `length_body`, or `taper`.
#' @examples
#' shape_gas <- sphere(radius_body = 0.01, n_segments = 60)
#' gas_generate(shape = shape_gas, g_fluid = 0.0012, h_fluid = 0.22)
#' @details
#' The preferred workflow is to supply a pre-built `Shape` object or explicit
#' profile coordinates and then describe the internal gas by either contrasts
#' (`g_fluid`, `h_fluid`) or absolute density/sound-speed values
#' (`density_fluid`, `sound_speed_fluid`). Character-based shape dispatch
#' remains available only as a compatibility pathway and is now deprecated.
#'
#' Scatterer constructors store geometry in meters and orientations in radians.
#' `radius_units` and `theta_units` are retained as compatibility arguments, but
#' non-SI values are normalized to the package-standard representation.
#'
#' @return
#' GAS-class object
#'
#' @seealso \code{\link{GAS}}
#'
#' @importFrom methods new
#' @keywords scatterer_type_generation
#' @export
gas_generate <- function(shape = NULL,
                         radius_body = NULL,
                         h_fluid = 0.2200,
                         g_fluid = 0.0012,
                         sound_speed_fluid = NULL,
                         density_fluid = NULL,
                         theta_body = pi / 2,
                         ID = NULL,
                         radius_units = "m",
                         theta_units = "radians",
                         n_segments = 100, ...) {
  # Normalize constructor units ================================================
  units <- .normalize_scatterer_units(
    theta_units = theta_units,
    radius_units = radius_units,
    context = "GAS"
  )
  # Validate fluid material inputs =============================================
  call_args <- names(as.list(match.call(expand.dots = FALSE)))
  h_supplied <- "h_fluid" %in% call_args
  g_supplied <- "g_fluid" %in% call_args
  c_supplied <- "sound_speed_fluid" %in%
    call_args && !is.null(sound_speed_fluid)
  rho_supplied <- "density_fluid" %in% call_args && !is.null(density_fluid)
  .validate_component_material_inputs(
    g = if (g_supplied) g_fluid else NULL,
    density = if (rho_supplied) density_fluid else NULL,
    h = if (h_supplied) h_fluid else NULL,
    sound_speed = if (c_supplied) sound_speed_fluid else NULL,
    component_label = "fluid",
    require_one = FALSE
  )
  # Resolve body shape input ===================================================
  arg_pull <- lapply(as.list(match.call()), eval, parent.frame())
  shape_input <- .resolve_scatterer_shape_input(shape, arg_pull)
  # Create metadata field ======================================================
  metadata <- .scatterer_metadata(ID = ID)
  # Create body shape field ====================================================
  body <- .build_column_major_fluid_component(
    shape_obj = shape_input,
    theta = theta_body,
    density = if (rho_supplied) density_fluid else NULL,
    sound_speed = if (c_supplied) sound_speed_fluid else NULL,
    g = if (rho_supplied) NULL else g_fluid,
    h = if (c_supplied) NULL else h_fluid
  )
  # Define shape parameters ====================================================
  shape_parameters <- .shape_common_parameters(
    shape_input = shape_input,
    requested_shape = shape,
    error_context = "GAS body",
    extra_units = list(
      radius_units = units$radius_units,
      theta_units = units$theta_units
    )
  )
  # Create GAS-class object ====================================================
  methods::new("GAS",
    metadata = metadata,
    model_parameters = list(),
    model = list(),
    body = body,
    shape_parameters = shape_parameters
  )
}
################################################################################
# Create GAS-class object
################################################################################
#' Generate an ESS-class object
#' @inheritParams fls_generate
#' @param shape Pre-built `Shape` object describing the outer shell geometry. If
#'   omitted, explicit shell profile coordinates such as `x_body`, `y_body`,
#'   and `z_body` are treated as the manual geometry pathway. Legacy character
#'   dispatch such as `"sphere"` is retained only for backward compatibility and
#'   is now deprecated.
#' @param radius_shell Radius of shell (m).
#' @param shell_thickness Optional shell thickness (m).
#' @param g_fluid Optional density contrast for fluid-like body.
#' @param density_fluid Optional density for fluid-like body (kg/m³).
#' @param h_fluid Optional sound speed contrast for fluid-like body.
#' @param sound_speed_fluid Optional sound speed for fluid-like body (m/s).
#' @param g_shell Density contrast for the shell.
#' @param density_shell Optional density for the shell (kg/m³).
#' @param h_shell Sound speed contrast for the shell.
#' @param sound_speed_shell Optional sound speed for the shell (m/s).
#' @param E Young's modulus (Pa) of the shell material.
#' @param nu Poisson's ratio (Dimensionless) of the shell material.
#' @param G Shear modulus (Pa) of the shell material.
#' @param K Bulk modulus (Pa) of the shell material.
#' @param theta_shell Object orientation relative to incident sound wave.
#' @details
#' The preferred workflow is to build the shell geometry first as a `Shape`
#' object and pass it through `shape`, or to supply explicit shell profile
#' coordinates directly. Character-based shape dispatch remains available only
#' as a compatibility pathway and is now deprecated. Material properties for
#' both the shell and the inner fluid may be supplied either as contrasts or as
#' absolute values, but each property pair must use only one representation.
#'
#' Scatterer constructors store geometry in meters and orientations in radians.
#' `length_units` and `theta_units` are retained as compatibility arguments, but
#' non-SI values are normalized to the package-standard representation.
#'
#' @return ESS-class object
#' @examples
#' shell <- sphere(radius_body = 0.03, n_segments = 80)
#' ess_generate(
#'   shape = shell,
#'   shell_thickness = 0.001,
#'   density_shell = 1050,
#'   sound_speed_shell = 2350,
#'   density_fluid = 1030,
#'   sound_speed_fluid = 1500,
#'   E = 3.5e9,
#'   nu = 0.34
#' )
#'
#' @seealso \code{\link{ESS}}
#'
#' @importFrom methods new
#' @keywords scatterer_type_generation
#' @export
ess_generate <- function(shape = NULL,
                         x_body = NULL,
                         y_body = NULL,
                         z_body = NULL,
                         radius_shell = NULL,
                         shell_thickness = NULL,
                         g_fluid = NULL,
                         density_fluid = NULL,
                         h_fluid = NULL,
                         sound_speed_fluid = NULL,
                         g_shell = NULL,
                         density_shell = NULL,
                         h_shell = NULL,
                         sound_speed_shell = NULL,
                         E = NULL,
                         G = NULL,
                         K = NULL,
                         nu = NULL,
                         theta_shell = pi / 2,
                         ID = NULL,
                         theta_units = "radians",
                         length_units = "m") {
  # Normalize constructor units ================================================
  units <- .normalize_scatterer_units(
    theta_units = theta_units,
    length_units = length_units,
    context = "ESS"
  )
  # Create metadata field ======================================================
  metadata <- .scatterer_metadata(ID = ID)
  # Resolve shell and fluid shape inputs =======================================
  arg_pull <- lapply(as.list(match.call()), eval, parent.frame())
  rpos <- Filter(
    Negate(is.null),
    .resolve_ess_shape_components(
      shape = shape,
      arguments = arg_pull,
      radius_shell = radius_shell,
      shell_thickness = shell_thickness,
      x_body = x_body,
      y_body = y_body,
      z_body = z_body
    )
  )
  # Assign elastic properties ==================================================
  elastic_params <- .complete_elastic_moduli(E = E, G = G, K = K, nu = nu)
  # Initialize material property containers ====================================
  material_properties <- list()
  # Validate shell and fluid material inputs ===================================
  .validate_component_material_inputs(
    g = g_shell,
    density = density_shell,
    h = h_shell,
    sound_speed = sound_speed_shell,
    component_label = "shell",
    require_one = FALSE
  )
  .validate_component_material_inputs(
    g = g_fluid,
    density = density_fluid,
    h = h_fluid,
    sound_speed = sound_speed_fluid,
    component_label = "fluid",
    require_one = FALSE
  )
  # Define shell material properties ===========================================
  material_properties[["shell"]] <- c(
    Filter(
      Negate(is.null),
      list(
        g = g_shell, density = density_shell,
        h = h_shell, sound_speed = sound_speed_shell
      )
    ),
    elastic_params
  )
  # Define fluid material properties ===========================================
  material_properties[["fluid"]] <- Filter(
    Negate(is.null),
    list(
      g = g_fluid, density = density_fluid,
      h = h_fluid, sound_speed = sound_speed_fluid
    )
  )
  # Finalize shell slot ========================================================
  shell <- c(
    list(
      rpos = rpos[["shell"]]@position_matrix,
      radius = rpos[["shell"]]@shape_parameters$radius,
      shell_thickness = ifelse(!is.null(shell_thickness),
        shell_thickness,
        NA
      ),
      theta = theta_shell
    ),
    material_properties[["shell"]]
  )
  # Finalize fluid slot ========================================================
  fluid <- c(
    list(
      rpos = if ("fluid" %in% names(rpos)) {
        rpos[["fluid"]]@position_matrix
      } else {
        NULL
      },
      radius = if ("fluid" %in% names(rpos)) {
        rpos[["fluid"]]@shape_parameters$radius
      } else {
        NULL
      },
      theta = theta_shell
    ),
    material_properties[["fluid"]]
  )
  # Define shape parameters ====================================================
  normalize_ess_shape_params <- function(params) {
    .normalize_shape_diameter(params)
  }

  shell_shape_params <- c(
    rpos[["shell"]]@shape_parameters,
    list(length_units = units$length_units)
  )
  shell_shape_params <- normalize_ess_shape_params(shell_shape_params)
  fluid_shape_params <- if ("fluid" %in% names(rpos)) {
    c(
      rpos[["fluid"]]@shape_parameters,
      list(length_units = units$length_units)
    )
  } else {
    list(
      diameter = NA,
      radius = NA,
      length_units = units$length_units
    )
  }
  fluid_shape_params <- normalize_ess_shape_params(fluid_shape_params)
  shape_parameters <- list(
    shell = shell_shape_params,
    fluid = fluid_shape_params,
    shape = class(rpos[["shell"]]),
    theta_units = units$theta_units
  )
  # Create ESS-class object ====================================================
  return(methods::new("ESS",
    metadata = metadata,
    model_parameters = list(),
    model = list(),
    shell = shell,
    fluid = fluid,
    shape_parameters = shape_parameters
  ))
}
################################################################################
# Helper function
################################################################################
#' Resolve shape-specific class assignment
#' @param shape_fun Shape-generating function.
#' @param arguments Input functions that will parameterize the shape-generating
#' function via 'shape_fun'.
#' @return Shape-class object.
#'
#' @keywords internal
#' @noRd
.resolve_shape <- function(shape_fun, arguments) {
  # Pass if already Shape-class ================================================
  if (methods::is(shape_fun, "Shape")) {
    return(shape_fun)
  }
  # Resolve function object if character supplied [LEGACY] =====================
  fun_obj <- if (is.character(shape_fun)) {
    get(shape_fun, mode = "function")
  } else {
    shape_fun
  }
  shape_name <- if (is.character(shape_fun)) {
    shape_fun
  } else {
    deparse(substitute(shape_fun))
  }
  # Validate kwargs ============================================================
  .validate_shape_requirements(shape_name, arguments)
  if (shape_fun != "arbitrary") {
    true_args <- .filter_shape_args(fun_obj, arguments)
  } else {
    true_args <- arguments
  }

  do.call(fun_obj, true_args)
}

#' Validate Shape-class-specific kwargs
#' @param shape_name Name of proposed Shape-class.
#' @param arguments Input functions that will parameterize the shape-generating
#' function via 'shape_fun'.
#' @return Shape-class object.
#'
#' @keywords internal
#' @noRd
.validate_shape_requirements <- function(shape_name, arguments) {
  # Format inputs to get Shape-class generator =================================
  name <- tolower(shape_name)
  missing_arg <- function(key) {
    if (!(key %in% names(arguments))) {
      return(TRUE)
    }
    value <- arguments[[key]]
    is.null(value) || (is.atomic(value) && all(is.na(value)))
  }
  stop_missing <- function(msg) stop(msg, call. = FALSE)
  # Toggle through the valid Shape-classes =====================================
  switch(name,
    sphere = {
      if (missing_arg("radius_body")) {
        stop_missing("Sphere requires 'radius_body'.")
      }
    },
    cylinder = {
      if (missing_arg("length_body")) {
        stop_missing("Cylinder requires 'length_body'.")
      }
      if (missing_arg("radius_body") && missing_arg("length_radius_ratio")) {
        stop_missing(
          "Cylinder requires either 'radius_body' or 'length_radius_ratio'."
        )
      }
    },
    prolate_spheroid = {
      length_missing <- missing_arg("length_body") &&
        missing_arg("semimajor_length")
      radius_missing <- missing_arg("radius_body") &&
        missing_arg("semiminor_length")
      if (length_missing) {
        stop_missing(
          "Prolate spheroid requires 'length_body' or 'semimajor_length'."
        )
      }
      if (radius_missing && missing_arg("length_radius_ratio")) {
        stop_missing(
          "Prolate spheroid requires 'radius_body', 'semiminor_length', or ",
          "'length_radius_ratio'."
        )
      }
    },
    oblate_spheroid = {
      length_missing <- missing_arg("length_body") &&
        missing_arg("semiminor_length")
      radius_missing <- missing_arg("radius_body") &&
        missing_arg("semimajor_length")
      if (length_missing) {
        stop_missing(
          "Oblate spheroid requires 'length_body' or 'semiminor_length'."
        )
      }
      if (radius_missing && missing_arg("length_radius_ratio")) {
        stop_missing(
          "Oblate spheroid requires 'radius_body', 'semimajor_length', or ",
          "'length_radius_ratio'."
        )
      }
    },
    polynomial_cylinder = {
      if (missing_arg("length_body")) {
        stop_missing("Polynomial cylinder requires 'length_body'.")
      }
      if (missing_arg("radius_body")) {
        stop_missing("Polynomial cylinder requires 'radius_body'.")
      }
      if (missing_arg("polynomial")) {
        stop_missing("Polynomial cylinder requires 'polynomial' coefficients.")
      }
    },
    arbitrary = {
      # Allow arbitrary inputs; validation handled downstream
    },
    NULL
  )
}

#' Filter shape arguments
#' @param shape_fun Shape-generating function.
#' @param arguments Input functions that will parameterize the shape-generating
#' function via 'shape_fun'.
#' @return Filtered argument list.
#'
#' @keywords internal
#' @noRd
.filter_shape_args <- function(shape_fun, arguments) {
  # ============================================================================
  # Get argument names
  arg_names <- names(formals(shape_fun))
  # ============================================================================
  # Get the overlapping names
  filtered_args <- arguments[arg_names]
  # ============================================================================
  # Prune the arguments
  Filter(Negate(is.null), filtered_args)
}
