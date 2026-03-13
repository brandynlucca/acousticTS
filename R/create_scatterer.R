################################################################################
# CREATE SCATTERER FUNCTIONS
################################################################################
################################################################################
# Create SBF-class object
################################################################################
#' Manually generate a SBF-class object.
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
#' @param theta_units Angular units.
#' @param length_units Length units.
#' @param ID Optional metadata identifier.
#' @param body_shape Optional pre-built Shape for the body; overrides
#'   x_body/w_body/zU_body/zL_body if supplied.
#' @param bladder_shape Optional pre-built Shape for the bladder; overrides
#'   x_bladder/w_bladder/zU_bladder/zL_bladder if supplied.
#' @details
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
  # Generate shape position matrix =============================================
  # Create body shape field ====================================================
  # Body shape can be a Shape or manual inputs --------------------------------
  body_shape_obj <- if (methods::is(body_shape, "Shape")) {
    body_shape
  } else {
    arbitrary(
      x_body = x_body, w_body = w_body,
      zU_body = zU_body, zL_body = zL_body
    )
  }
  bladder_shape_obj <- if (methods::is(bladder_shape, "Shape")) {
    bladder_shape
  } else {
    arbitrary(
      x_bladder = x_bladder, w_bladder = w_bladder,
      zU_bladder = zU_bladder, zL_bladder = zL_bladder
    )
  }
  # Normalize to x/w/zU/zL rows with suffix -----------------------------------
  normalize_rpos <- function(shape_obj, suffix) {
    pos <- acousticTS::extract(shape_obj, "position_matrix")
    col_or_default <- function(candidates, default = NULL) {
      for (nm in candidates) {
        if (!is.null(colnames(pos)) && nm %in% colnames(pos)) {
          return(pos[, nm])
        }
      }
      if (!is.null(default)) return(default)
      pos[, 1]
    }
    x <- col_or_default(c(paste0("x_", suffix), "x_body", "x_bladder", "x"))
    w <- col_or_default(
      c(paste0("w_", suffix), "w_body", "w_bladder", "w", "y_body", "y"),
      default = rep(0, length(x))
    )
    zU <- col_or_default(
      c(paste0("zU_", suffix), "zU_body", "zU_bladder", "zU"),
      default = rep(NA_real_, length(x))
    )
    zL <- col_or_default(
      c(paste0("zL_", suffix), "zL_body", "zL_bladder", "zL"),
      default = rep(NA_real_, length(x))
    )
    rpos <- rbind(x, w, zU, zL)
    rownames(rpos) <- if (suffix == "body") {
      c("x_body", "w_body", "zU_body", "zL_body")
    } else {
      c("x_bladder", "w_bladder", "zU_bladder", "zL_bladder")
    }
    rpos
  }
  # Define body shape ==========================================================
  body <- list(
    rpos = normalize_rpos(body_shape_obj, "body"),
    sound_speed = sound_speed_body,
    density = density_body,
    g = g_body,
    h = h_body,
    theta = theta_body
  )
  # Define bladder shape =======================================================
  bladder <- list(
    rpos = normalize_rpos(bladder_shape_obj, "bladder"),
    sound_speed = sound_speed_bladder,
    density = density_bladder,
    g = g_bladder,
    h = h_bladder,
    theta = theta_bladder
  )
  # Validate segment counts (need at least one interval) ----------------------
  if (ncol(body$rpos) < 2) {
    stop("Body shape must have at least two points (>=1 segment).",
         call. = FALSE)
  }
  if (ncol(bladder$rpos) < 2) {
    stop("Bladder shape must have at least two points (>=1 segment).",
         call. = FALSE)
  }
  # Define shape parameters ====================================================
  shape_parameters <- list(
    body = list(
      shape = paste0(class(body_shape_obj)),
      length = max(body$rpos[1, ]),
      n_segments = ncol(body$rpos)
    ),
    bladder = list(
      shape = paste0(class(bladder_shape_obj)),
      length = max(bladder$rpos[1, ], na.rm = TRUE) -
        min(bladder$rpos[1, ], na.rm = TRUE),
      n_segments = ncol(bladder$rpos)
    ),
    length_units = length_units,
    theta_units = theta_units
  )
  # Create metadata field ======================================================
  metadata <- list(ID = ifelse(!is.null(ID),
    ID,
    "UID"
  ))
  # Create FLS-class object ====================================================
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
#' @param diameter_units Units for diameter. Defaults to "m".
#' @param theta_units Units for direction. Defaults to "radians".
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
  # Define user input or default object ID =====================================
  metadata <- list(
    ID = ifelse(!is.null(ID),
      ID,
      "Calibration sphere"
    ),
    Material = material
  )
  # Create sphere object to define definitions =================================
  sphere_shape <- sphere(
    radius_body = diameter / 2,
    n_segments = n_segments,
    diameter_units = "m"
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
    diameter_units = diameter_units,
    theta_units = theta_units
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
#' Manually generate a FLS object.
#' @param shape Optional input argument that dictates shape-type, if desired,
#' for generalized shapes.
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
#' @param theta_units Units used for orientation. Defaults to "radians".
#' @param length_units Units used for position vector. Defaults to "m".
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
#' Material properties can be provided either as contrasts (`g_body`/`h_body`)
#' or as absolute density/sound speed (`density_body`/`sound_speed_body`).
#' Downstream models will derive contrasts automatically when only absolute
#' values are supplied.
#'
#' @seealso \code{\link{FLS}}
#'
#' @importFrom methods new
#' @keywords scatterer_type_generation
#' @export
fls_generate <- function(shape = "arbitrary",
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
  # Collect shape information if provided =====================================
  if (!methods::is(shape, "Shape") && shape == "arbitrary") {
    if (
      (is.null(x_body) && is.null(y_body) && is.null(z_body)) &&
        is.null(length_body)
    ) {
      stop("Body shape is not appropriately parameterized.")
    }
  }
  # Generate shape position matrix ============================================
  arg_pull <- lapply(as.list(match.call()), eval, parent.frame())
  shape_input <- .resolve_shape(shape, arg_pull)
  # Define shape parameters ==================================================
  pos_mat <- acousticTS::extract(shape_input, "position_matrix")
  shape_params <- acousticTS::extract(shape_input, "shape_parameters")
  shape_parameters <- list(
    length = max(pos_mat[, 1], na.rm = TRUE) - min(pos_mat[, 1], na.rm = TRUE),
    radius = max(shape_params$radius, na.rm = TRUE),
    n_segments = length(pos_mat[, 1]) - 1,
    length_units = length_units,
    theta_units = theta_units
  )
  # Add prolate spheroidal parameters, if necessary ++++++++++++++++++++++++++++
  if (methods::is(shape_input, "ProlateSpheroid")) {
    shape_parameters$semimajor_length <- shape_params$semimajor_length
    shape_parameters$semiminor_length <- shape_params$semiminor_length
  }
  # Add cylindrical parameters, if necessary +++++++++++++++++++++++++++++++++++
  if (methods::is(shape_input, "Cylinder")) {
    shape_parameters$radius_curvature_ratio <- (
      shape_params$radius_curvature_ratio
    )
  }
  if (is.character(shape) && shape == "arbitrary") {
    shape_parameters[["shape"]] <- "Arbitrary"
  } else if (methods::is(shape_input, "Shape")) {
    shape_parameters[["shape"]] <- paste0(class(shape_input))
  }

  if (methods::is(shape, "Sphere") ||
      (is.character(shape) && shape == "sphere")) {
    shape_parameters[["radius_shape"]] <- acousticTS::extract(
      shape_input,
      "shape_parameters"
    )$radius_shape
  }

  if ((is.character(shape) && shape == "cylinder") ||
      methods::is(shape, "Cylinder")) {
    shape_parameters$taper_order <- shape_input@shape_parameters$taper_order
  }
  # Check material properties length/availability ============================
  if (is.null(g_body) && is.null(density_body)) {
    stop("Supply either 'g_body' or 'density_body' for FLS objects.",
         call. = FALSE)
  }
  if (is.null(h_body) && is.null(sound_speed_body)) {
    stop("Supply either 'h_body' or 'sound_speed_body' for FLS objects.",
         call. = FALSE)
  }
  # g / density length checks -------------------------------------------------
  if (!is.null(g_body) && length(g_body) > 1 &&
      length(g_body) != length(shape_input@position_matrix[, 1]) - 1) {
    stop(
      paste0(
        "Vector input for 'g_body' with ", length(g_body),
        " elements does not match the expected number of segments (",
        length(shape_input@position_matrix[, 1]) - 1, ")"
      )
    )
  }
  if (!is.null(density_body) && length(density_body) > 1 &&
      length(density_body) != length(shape_input@position_matrix[, 1]) - 1) {
    stop(
      paste0(
        "Vector input for 'density_body' with ", length(density_body),
        " elements does not match the expected number of segments (",
        length(shape_input@position_matrix[, 1]) - 1, ")"
      )
    )
  }
  # h / sound speed length checks --------------------------------------------
  if (!is.null(h_body) && length(h_body) > 1 &&
      length(h_body) != length(shape_input@position_matrix[, 1]) - 1) {
    stop(
      paste0(
        "Vector input for 'h_body' with ", length(h_body),
        " elements does not match the expected number of segments (",
        length(shape_input@position_matrix[, 1]) - 1, ")"
      )
    )
  }
  if (!is.null(sound_speed_body) && length(sound_speed_body) > 1 &&
      length(sound_speed_body) != length(shape_input@position_matrix[, 1]) - 1) {
    stop(
      paste0(
        "Vector input for 'sound_speed_body' with ", length(sound_speed_body),
        " elements does not match the expected number of segments (",
        length(shape_input@position_matrix[, 1]) - 1, ")"
      )
    )
  }
  # Create body slot =========================================================
  body <- list(
    rpos = t(pos_mat),
    radius = shape_params$radius,
    radius_curvature_ratio = radius_curvature_ratio,
    theta = theta_body,
    g = g_body,
    h = h_body,
    density = density_body,
    sound_speed = sound_speed_body
  )
  # Create metadata field ======================================================
  metadata <- list(ID = ifelse(!is.null(ID),
    ID,
    "UID"
  ))
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
#' Create GAS object
#'
#' @inheritParams fls_generate
#' @param shape Optional pre-made shape input. Default is a sphere.
#' @param radius_body Radius (m). For non-canonical shapes, this would be the
#' maximum or mean radius at the scatterer midsection.
#' @param h_fluid Sound speed contrast of fluid relative to surrounding
#' medium (h).
#' @param g_fluid Density contrast of fluid relative to surrounding density (g).
#' @param sound_speed_fluid Optional fluid sound speed (m/s).
#' @param density_fluid Optional fluid density (m/s).
#' @param radius_units Diameter units. Defaults to "m".
#' @param n_segments Number of body segments.
#' @examples
#' shape_gas <- sphere(radius_body = 0.01, n_segments = 60)
#' gas_generate(shape = shape_gas, g_fluid = 0.0012, h_fluid = 0.22)
#' @return
#' GAS-class object
#'
#' @seealso \code{\link{GAS}}
#'
#' @importFrom methods new
#' @keywords scatterer_type_generation
#' @export
gas_generate <- function(shape = "sphere",
                         h_fluid = 0.2200,
                         g_fluid = 0.0012,
                         sound_speed_fluid = NULL,
                         density_fluid = NULL,
                         theta_body = pi / 2,
                         ID = NULL,
                         radius_units = "m",
                         theta_units = "radians",
                         n_segments = 100, ...) {
  # Collect shape information if provided ======================================
  arg_pull <- lapply(as.list(match.call()), eval, parent.frame())
  shape_input <- .resolve_shape(shape, arg_pull)
  # Create metadata field ======================================================
  metadata <- list(ID = ifelse(!is.null(ID),
    ID,
    "UID"
  ))
  # Create body shape field ====================================================
  body <- list(
    rpos = shape_input@position_matrix,
    radius = shape_input@shape_parameters$radius,
    theta = theta_body,
    g = g_fluid,
    h = h_fluid
  )
  # Define shape parameters ====================================================
  shape_parameters <- list(
    length = max(
      max(acousticTS::extract(shape_input, "position_matrix")[, 1])
    ),
    radius = shape_input@shape_parameters$radius,
    n_segments = n_segments,
    radius_units = radius_units,
    theta_units = theta_units,
    shape = class(shape_input)
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
#' Generate ESS shape
#' @inheritParams fls_generate
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
ess_generate <- function(shape = "sphere",
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
  # Create metadata field ======================================================
  metadata <- list(ID = ifelse(!is.null(ID), ID, "UID"))
  # Create container for different position matrices ===========================
  rpos <- list()
  # Create shape fields ========================================================
  if (methods::is(shape, "Shape")) {
    rpos[["shell"]] <- shape
  } else {
    arg_pull <- arg_pull <- lapply(as.list(match.call()), eval, parent.frame())
    # Build shell shape --------------------------------------------------------
    if (shape == "arbitrary") {
      rpos[["shell"]] <- arbitrary(
        x_body = x_body,
        y_body = y_body,
        z_body = z_body,
        radius_body = radius_shell
      )
    } else {
      # Reuse shape resolver for shell
      shell_call <- arg_pull
      shell_call$radius_body <- radius_shell
      rpos[["shell"]] <- .resolve_shape(shape, shell_call)
      # Optional inner fluid shape when thickness provided
      if (!is.null(shell_thickness)) {
        fluid_call <- arg_pull
        fluid_call$radius_body <- radius_shell - shell_thickness
        rpos[["fluid"]] <- .resolve_shape(shape, fluid_call)
      }
    }
  }
  # Assign elastic properties ==================================================
  elastic_params <- Filter(
    Negate(is.null),
    list(E = E, nu = nu, G = G, K = K)
  )
  # Iterate through to calculate the requisite parameters ++++++++++++++++++++++
  repeat {
    param_len <- length(elastic_params)
    if (is.null(elastic_params$nu)) {
      elastic_params$nu <- tryCatch(
        pois(elastic_params$K, elastic_params$E, elastic_params$G),
        error = function(e) NULL
      )
    }
    if (is.null(elastic_params$K)) {
      elastic_params$K <- tryCatch(
        bulk(elastic_params$E, elastic_params$G, elastic_params$nu),
        error = function(e) NULL
      )
    }
    if (is.null(elastic_params$E)) {
      elastic_params$E <- tryCatch(
        young(elastic_params$K, elastic_params$G, elastic_params$nu),
        error = function(e) NULL
      )
    }
    if (is.null(elastic_params$G)) {
      elastic_params$G <- tryCatch(
        shear(elastic_params$K, elastic_params$E, elastic_params$nu),
        error = function(e) NULL
      )
    }
    if (is.null(elastic_params$lambda)) {
      elastic_params$lambda <- tryCatch(
        lame(
          elastic_params$K, elastic_params$E, elastic_params$G,
          elastic_params$nu
        ),
        error = function(e) NULL
      )
    }
    if (length(elastic_params) == param_len) break
  }
  # Handle material properties =================================================
  material_properties <- list()
  # Validate parameter combinations ++++++++++++++++++++++++++++++++++++++++++++
  conflicts <- list(
    c("g_shell", "density_shell"),
    c("h_shell", "sound_speed_shell"),
    c("g_fluid", "density_fluid"),
    c("h_fluid", "sound_speed_fluid")
  )
  for (pair in conflicts) {
    if (!is.null(get(pair[1])) && !is.null(get(pair[2]))) {
      stop(paste("Cannot specify both", pair[1], "and", pair[2]))
    }
  }
  # Shell material properties ++++++++++++++++++++++++++++++++++++++++++++++++++
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
  # Fluid material properties ++++++++++++++++++++++++++++++++++++++++++++++++++
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
  # Shape parameters field =====================================================
  shape_parameters <- list(
    shell = list(
      diameter = shell[["radius"]] * 2,
      radius = shell[["radius"]],
      length_units = length_units
    ),
    fluid = list(
      diameter = ifelse(is.null(fluid[["radius"]]),
        NA,
        fluid[["radius"]] * 2
      ),
      radius = ifelse(is.null(fluid[["radius"]]),
        NA,
        fluid[["radius"]]
      ),
      length_units = length_units
    ),
    shape = class(rpos[["shell"]])
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
#' Filter shape arguments
#' @param shape_fun Shape-generating function.
#' @param arguments Input functions that will parameterize the shape-generating
#' function via 'shape_fun'.
#' @return Filtered argument list.
#'
#' @keywords internal
#' @noRd
.resolve_shape <- function(shape_fun, arguments) {
  if (methods::is(shape_fun, "Shape")) {
    return(shape_fun)
  }
  # Resolve function object if character supplied
  fun_obj <- if (is.character(shape_fun)) get(shape_fun, mode = "function") else shape_fun
  shape_name <- if (is.character(shape_fun)) shape_fun else deparse(substitute(shape_fun))
  .validate_shape_requirements(shape_name, arguments)
  if (shape_fun != "arbitrary") {
    true_args <- .filter_shape_args(fun_obj, arguments)
  } else {
    true_args <- arguments
  }
  do.call(fun_obj, true_args)
}

.validate_shape_requirements <- function(shape_name, arguments) {
  name <- tolower(shape_name)
  missing_arg <- function(key) {
    if (!(key %in% names(arguments))) return(TRUE)
    value <- arguments[[key]]
    is.null(value) || (is.atomic(value) && all(is.na(value)))
  }
  stop_missing <- function(msg) stop(msg, call. = FALSE)
  switch(
    name,
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
        stop_missing("Cylinder requires either 'radius_body' or 'length_radius_ratio'.")
      }
    },
    prolate_spheroid = {
      length_missing <- missing_arg("length_body") && missing_arg("semimajor_length")
      radius_missing <- missing_arg("radius_body") && missing_arg("semiminor_length")
      if (length_missing) {
        stop_missing("Prolate spheroid requires 'length_body' or 'semimajor_length'.")
      }
      if (radius_missing && missing_arg("length_radius_ratio")) {
        stop_missing("Prolate spheroid requires 'radius_body', 'semiminor_length', or 'length_radius_ratio'.")
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
