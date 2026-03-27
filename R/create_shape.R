################################################################################
################################################################################
# FUNCTIONS FOR GENERATIGN CANONICAL & PRE-DEFINED SHAPES
################################################################################
################################################################################
################################################################################
# Arbitrary (or pre-generated) body shape parameters
################################################################################
#' Creates arbitrary body shape from user inputs
#' @param ... Any coordinate or parameter vectors (e.g., x_body, y_body,
#'   z_body, radius_body, w_body, zU_body, zL_body, custom_*). Inputs are
#'   stored and padded to a common length; validation happens downstream in
#'   model-specific code.
#' @param length_units Units for body length. Defaults to meters: "m"
#' @examples
#' # Symmetric arbitrary shape using x/y/z + radius
#' arbitrary(
#'   x_body = c(0, 0.01), y_body = c(0, 0), z_body = c(0, 0),
#'   radius_body = c(0, 0.002)
#' )
#'
#' # Dorsal/ventral style inputs
#' arbitrary(
#'   x_body = c(0, 0.015),
#'   w_body = c(0.005, 0.0075),
#'   zU_body = c(0.001, 0.002),
#'   zL_body = c(-0.001, -0.002)
#' )
#' @return An \code{\link{Arbitrary}} shape object containing the padded
#'   position matrix and stored shape metadata.
#' @seealso \code{\link{Arbitrary}}
#'
#' @keywords shape_generation
#' @rdname arbitrary
#' @export
arbitrary <- function(...,
                      length_units = "m") {
  # Gather arguments ===========================================================
  args <- list(...)
  coord_names <- names(args)
  # If an rpos matrix is provided, use it directly =============================
  if (!is.null(args$rpos) && is.matrix(args$rpos)) {
    position_matrix <- args$rpos
  } else {
    # Collect all numeric vectors; pad to common length ========================
    numeric_vecs <- args[vapply(args, is.numeric, logical(1))]
    # Check coordinate aliases =================================================
    all_aliases <- Filter(Negate(is.null), list(
      x = .validate_coord_aliases(coord_names, "^x"),
      y = .validate_coord_aliases(coord_names, "^y"),
      w = .validate_coord_aliases(coord_names, "^w"),
      z = .validate_coord_aliases(coord_names, "^z(?![UL])"),
      zU = .validate_coord_aliases(coord_names, "^(zU)"),
      zL = .validate_coord_aliases(coord_names, "^(zL)"),
      a = .validate_coord_aliases(coord_names,
                                  "^radius(?!_curvature)\\w*")
    ))
    if (length(all_aliases) < 2) {
      stop(
        "'Arbitrary' shape-class requires at least two numeric vectors with ",
        "valid coordinates (e.g., 'x_body' and 'radius_body'). Current ",
        "detected coordinates from input: '",
        paste(unlist(all_aliases), collapse="', '"), "'.",
        call. = FALSE
      )
    }
    # Default column order hint: prioritize common names when present ==========
    position_matrix <- do.call(cbind, numeric_vecs[unlist(all_aliases)])
    colnames(position_matrix) <- names(all_aliases)
  }
  # Update shape parameters ====================================================
  shape_parameters <- list(
    n_segments = max(nrow(position_matrix) - 1, 1),
    length_units = length_units,
    diameter_units = length_units
  )
  # Check for radius +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  if ("a" %in% colnames(position_matrix)) {
    shape_parameters$radius <- position_matrix[, "a"]
  } else {
    shape_parameters$radius <- NA
  }
  # Back-fill zU and zL if not already present -- assuming uniformity ++++++++++
  is_zU <- "zU" %in% colnames(position_matrix)
  is_zL <- "zL" %in% colnames(position_matrix)
  is_a <- "a" %in% colnames(position_matrix)
  if (is_a & !is_zU & !is_zL) {
    position_matrix <- cbind(
      position_matrix,
      cbind(zU=position_matrix[,"a"], zL=-position_matrix[,"a"])
    )
  } else if(is_a & !is_zU & is_zL) {
    position_matrix <- cbind(
      position_matrix,
      cbind(zU=-position_matrix[,"zL"])
    )
  } else if(is_a & is_zU & !is_zL) {
    position_matrix <- cbind(
      position_matrix,
      cbind(zL=-position_matrix[,"zU"])
    )
  }
  shape_parameters <- c(
    shape_parameters,
    list(
      mean_radius = mean(shape_parameters$radius),
      max_radius = max(shape_parameters$radius)
    )
  )
  # Generate new shape objects =================================================
  methods::new(
    "Arbitrary",
    position_matrix = position_matrix,
    shape_parameters = shape_parameters
  )
}

.validate_coord_aliases <- function(x, pattern) {
  # Extract ====================================================================
  mask <- grepl(pattern, x, perl = TRUE)
  aliases <- x[mask]
  # Validate length ============================================================
  if (length(aliases) > 1) {
    stop(
      "The following arbitrary entries overlap and are ambiguous. There ",
      "should only be 1 entry corresponding to 'x', 'y', 'z', 'zU', and/or ",
      "'zL'. The following overlapped: '",
      paste(aliases, collapse = "', '"),
      "'."
    )
  } else if(length(aliases) == 0) {
    return(NULL)
  }
  aliases
}

################################################################################
# Sphere
################################################################################
#' Creates a sphere.
#' @param radius_body Object radius (m).
#' @param n_segments Number of segments to discretize object shape. Defaults to
#'    1e2 segments.
#' @param diameter_units Default is "m" for meters
#' @examples
#' sphere(radius_body = 0.01, n_segments = 50)
#' @return
#' Creates position vector for a spherical object of a defined radius.
#' @seealso \code{\link{Sphere}}
#'
#' @keywords shape_generation
#' @rdname sphere
#' @export
sphere <- function(radius_body,
                   n_segments = 1e2,
                   diameter_units = "m") {
  # Define semi-major or x-axis ================================================
  diameter <- radius_body * 2
  x_nodes <- seq(
    from = 0,
    to = diameter,
    length.out = n_segments + 1
  )
  # Along-semimajor radii ======================================================
  radius_output <- sqrt(pmax(
    radius_body * radius_body -
      (x_nodes - radius_body) * (x_nodes - radius_body),
    0
  ))
  # Maintain the historical trailing-to-leading x ordering =====================
  x_edges <- rev(x_nodes)
  radius_output <- rev(radius_output)
  # Generate position matrix ===================================================
  position_matrix <- cbind(
    x = x_edges,
    y = rep(0, length(x_edges)),
    z = rep(0, length(x_edges)),
    zU = radius_output,
    zL = -radius_output
  )
  # Generate shape parameters list =============================================
  shape_parameters <- list(
    radius_shape = radius_output,
    diameter_shape = diameter,
    radius = diameter / 2,
    n_segments = n_segments,
    diameter_units = diameter_units
  )
  # Generate new shape object ==================================================
  return(methods::new("Sphere",
    position_matrix = position_matrix,
    shape_parameters = shape_parameters
  ))
}
################################################################################
# Prolate spheroid
################################################################################
#' Creates a prolate spheroid.
#'
#' @param length_body Semi-major axis length (m).
#' @param semimajor_length Optional alias for semi-major axis length (m).
#' @param radius_body Semi-minor axis length (m). This can also be stylized as
#'    the "maximum radius" of the scattering object.
#' @param semiminor_length Optional alias for semi-minor axis length (m).
#' @param length_radius_ratio Optional ratio input when radius is not explicitly
#'    known.
#' @param n_segments Number of segments to discretize object shape. Defaults to
#'    18 segments.
#' @param length_units Units for body matrix (defaults to m).
#' @examples
#' prolate_spheroid(
#'   length_body = 0.04, radius_body = 0.004, n_segments = 60
#' )
#' prolate_spheroid(
#'   semimajor_length = 0.05, semiminor_length = 0.003
#' )
#' @return
#' Creates the position vector for a prolate spheroid object of defined
#'    semi-major and -minor axes.
#' @seealso \code{\link{ProlateSpheroid}}
#'
#' @keywords shape_generation
#' @rdname prolate_spheroid
#' @export
prolate_spheroid <- function(length_body = NULL,
                             radius_body = NULL,
                             length_radius_ratio = NULL,
                             semimajor_length = NULL,
                             semiminor_length = NULL,
                             n_segments = 100,
                             length_units = "m") {
  # Allow alternative argument names for clarity ===============================
  length_val <- if (!is.null(semimajor_length)) {
    semimajor_length * 2
  } else {
    length_body
  }
  radius_val <- if (!is.null(semiminor_length)) {
    semiminor_length
  } else {
    radius_body
  }
  # Define maximum radius ======================================================
  max_radius <- .calculate_max_radius(
    radius_val, length_val, length_radius_ratio
  )
  # Define semi-major or x-axis ================================================
  x_nodes <- seq(0, length_val, length.out = n_segments + 1)
  # Normalize the node positions for ellipse ==================================
  curved_x_nodes <- (x_nodes - length_val / 2) / (length_val / 2)
  # Compute the radius at each node ============================================
  radius_output <- max_radius * sqrt(pmax(1 - curved_x_nodes^2, 0))
  # Maintain the historical trailing-to-leading x ordering =====================
  x_edges <- rev(x_nodes)
  radius_output <- rev(radius_output)
  # Generate position matrix ===================================================
  position_matrix <- cbind(
    x = x_edges,
    y = rep(0, length(x_edges)),
    z = rep(0, length(x_edges)),
    zU = radius_output,
    zL = -radius_output
  )
  # Generate shape parameters list =============================================
  shape_parameters <- list(
    length = max(position_matrix[, 1]),
    radius = radius_output,
    semimajor_length = max(position_matrix[, 1]) / 2,
    semiminor_length = max_radius,
    length_radius_ratio = max(position_matrix[, 1]) / max_radius,
    n_segments = n_segments,
    length_units = length_units
  )
  # Generate new shape object ==================================================
  return(methods::new("ProlateSpheroid",
    position_matrix = position_matrix,
    shape_parameters = shape_parameters
  ))
}
################################################################################
# Oblate spheroid
################################################################################
#' Creates an oblate spheroid.
#'
#' @param length_body Body-axis length (m).
#' @param semiminor_length Optional alias for body-axis semi-length (m).
#' @param radius_body Maximum equatorial radius (m).
#' @param semimajor_length Optional alias for maximum equatorial radius (m).
#' @param length_radius_ratio Optional ratio input when radius is not explicitly
#'    known.
#' @param n_segments Number of segments to discretize object shape. Defaults to
#'    100 segments.
#' @param length_units Units for body matrix (defaults to m).
#' @examples
#' oblate_spheroid(
#'   length_body = 0.012, radius_body = 0.01, n_segments = 60
#' )
#' oblate_spheroid(
#'   semiminor_length = 0.006, semimajor_length = 0.01
#' )
#' @return
#' Creates the position vector for an oblate spheroid object of defined
#'    body-axis length and maximum equatorial radius.
#' @seealso \code{\link{OblateSpheroid}}
#'
#' @keywords shape_generation
#' @rdname oblate_spheroid
#' @export
oblate_spheroid <- function(length_body = NULL,
                            radius_body = NULL,
                            length_radius_ratio = NULL,
                            semimajor_length = NULL,
                            semiminor_length = NULL,
                            n_segments = 100,
                            length_units = "m") {
  # Allow alternative argument names for clarity ===============================
  length_val <- if (!is.null(semiminor_length)) {
    semiminor_length * 2
  } else {
    length_body
  }
  radius_val <- if (!is.null(semimajor_length)) {
    semimajor_length
  } else {
    radius_body
  }
  # Define maximum radius ======================================================
  max_radius <- .calculate_max_radius(
    radius_val, length_val, length_radius_ratio
  )
  # Validate oblate aspect ratio ===============================================
  if ((length_val / 2) > max_radius) {
    stop(
      "'OblateSpheroid' requires the body-axis semi-length to be <= the ",
      "equatorial radius. Use 'prolate_spheroid()' when the body axis is the ",
      "major axis.",
      call. = FALSE
    )
  }
  # Define along-body node positions ===========================================
  x_nodes <- seq(0, length_val, length.out = n_segments + 1)
  # Normalize the node positions for ellipse ==================================
  axial_half_length <- length_val / 2
  curved_x_nodes <- (x_nodes - axial_half_length) / axial_half_length
  # Compute the radius at each node ============================================
  radius_output <- max_radius * sqrt(pmax(1 - curved_x_nodes^2, 0))
  # Maintain the historical trailing-to-leading x ordering =====================
  x_edges <- rev(x_nodes)
  radius_output <- rev(radius_output)
  # Generate position matrix ===================================================
  position_matrix <- cbind(
    x = x_edges,
    y = rep(0, length(x_edges)),
    z = rep(0, length(x_edges)),
    zU = radius_output,
    zL = -radius_output
  )
  # Generate shape parameters list =============================================
  shape_parameters <- list(
    length = max(position_matrix[, 1]),
    radius = radius_output,
    semimajor_length = max_radius,
    semiminor_length = axial_half_length,
    length_radius_ratio = max(position_matrix[, 1]) / max_radius,
    n_segments = n_segments,
    length_units = length_units
  )
  # Generate new shape object ==================================================
  return(methods::new("OblateSpheroid",
    position_matrix = position_matrix,
    shape_parameters = shape_parameters
  ))
}
################################################################################
# Elongated cylinder
################################################################################
#' Creates a cylinder.
#'
#' @param length_body Length (m).
#' @param radius_body Maximum/uniform radius (m).
#' @param length_radius_ratio Optional ratio input when radius is not explicitly
#'    known.
#' @param taper Optional input that is the degree of taper to round ends of
#'    the cylinder.
#' @param radius_curvature_ratio Optional curvature ratio metadata for rounded
#'    cylinder-end workflows.
#' @param n_segments Number of segments to discretize object shape. Defaults to
#'    1e2 segments.
#' @param length_units Units (default is meters, "m").
#' @usage cylinder(length_body, radius_body = NULL, length_radius_ratio = NULL,
#'   taper = NULL, radius_curvature_ratio = NULL, n_segments = 100,
#'   length_units = "m")
#' @examples
#' cylinder(length_body = 0.05, radius_body = 0.003, n_segments = 80)
#' @return
#' Creates the position vector for a tapered or untapered cylinder.
#' @seealso \code{\link{Cylinder}}
#'
#' @keywords shape_generation
#' @rdname cylinder
#' @export
cylinder <- function(length_body,
                     radius_body = NULL,
                     length_radius_ratio = NULL,
                     taper = NULL,
                     radius_curvature_ratio = NULL,
                     n_segments = 1e2,
                     length_units = "m") {
  # Define maximum radius ======================================================
  max_radius <- .calculate_max_radius(
    radius_body, length_body, length_radius_ratio
  )
  # Define normalized x-axis ===================================================
  x_n_axis <- seq(-1, 1, length.out = n_segments + 1)
  # Define tapered radius vector, if applicable ================================
  if (!is.null(taper)) {
    tapering <- sqrt(pmax(1 - x_n_axis^taper, 0))
  } else {
    tapering <- rep(1, n_segments + 1)
  }
  radius_output <- max_radius * tapering
  x_nodes <- x_n_axis * length_body / 2 + length_body / 2
  # Maintain the historical trailing-to-leading x ordering =====================
  x_edges <- rev(x_nodes)
  radius_output <- rev(radius_output)
  # Generate position matrix ===================================================
  position_matrix <- cbind(
    x = x_edges,
    y = rep(0, length(x_edges)),
    z = rep(0, length(x_edges)),
    zU = radius_output,
    zL = -radius_output
  )
  # Generate shape parameters list =============================================
  shape_parameters <- list(
    length = max(position_matrix[, 1]),
    radius = radius_output,
    length_radius_ratio = max(position_matrix[, 1]) / max_radius,
    n_segments = n_segments,
    taper_order = ifelse(is.null(taper),
      NA,
      taper
    ),
    radius_curvature_ratio = ifelse(is.null(radius_curvature_ratio),
                                    NA,
                                    radius_curvature_ratio),
    length_units = length_units
  )
  # Generate new shape object ==================================================
  return(methods::new("Cylinder",
    position_matrix = position_matrix,
    shape_parameters = shape_parameters
  ))
}
################################################################################
# Polynomial cylinder
################################################################################
#' Creates a polynomial deformed cylinder
#'
#' @inheritParams cylinder
#' @param polynomial Polynomial coefficient vector.
#' @examples
#' # We can use the polynomial coefficients defined in Smith et al. (2013) to
#' # define the position vector of a sub-Arctic krill.
#' poly_vec <- c(0.83, 0.36, -2.10, -1.20, 0.63, 0.82, 0.64)
#' # Create the position vector
#' # This outputs a list containing "rpos" and "radius"
#' pos <- polynomial_cylinder(
#'   length_body = 15e-3, radius_body = 2e-3,
#'   polynomial = poly_vec
#' )
#' str(pos)
#' @return
#' Creates the position vector for a polynomial deformed cylinder.
#' @references
#' Smith, J.N., Ressler, P.H., and Warren, J.D. 2013. A distorted wave Born
#' approximation target strength model for Bering Sea euphausiids. ICES Journal
#' of Marine Science, 70(1): 204-214. https://doi.org/10.1093/icesjms/fss140
#'
#' @seealso \code{\link{PolynomialCylinder}}
#'
#' @keywords shape_generation
#' @rdname polynomial_cylinder
#' @export
polynomial_cylinder <- function(length_body,
                                radius_body,
                                n_segments = 1e2,
                                polynomial,
                                length_units = "m") {
  # Define normalized x-axis ===================================================
  x_n_axis <- seq(-1, 1, length.out = n_segments + 1)
  # Evaluate polynomial coefficients ===========================================
  n_order <- rev(seq_len(length(polynomial))) - 1
  poly_fun <- paste0(polynomial, paste0("*x_n_axis^", n_order), collapse = "+")
  # Define radius ==============================================================
  radius_output <- abs(eval(parse(text = poly_fun))) * radius_body
  # Define output x-axis =======================================================
  x_axis <- x_n_axis * length_body / 2 + length_body / 2
  # Generate position vector, 'rpos' ===========================================
  position_matrix <- cbind(
    x = x_axis,
    y = rep(0, length(x_axis)),
    z = rep(0, length(x_axis))
  )
  # Generate shape parameters list =============================================
  shape_parameters <- list(
    length = max(position_matrix[, 1]),
    radius = radius_output,
    length_radius_ratio = max(position_matrix[, 1]) /
      max(radius_output),
    n_segments = n_segments,
    length_units = length_units
  )
  # Generate new shape object ==================================================
  return(methods::new("PolynomialCylinder",
    position_matrix = position_matrix,
    shape_parameters = shape_parameters
  ))
}
################################################################################
# Wrapper function for generating shape
################################################################################
#' A wrapper function that automatically creates generalized and/or canonical
#' shapes for TS modeling.
#'
#' @param shape Shape. Details for shape specification are provided under
#'    'Details', including mandatory arguments.
#' @param ... Additional input arguments for subsequent shape generation
#' functions.
#' @examples
#' create_shape("sphere", radius_body = 0.01)
#' create_shape(
#'   "prolate_spheroid", length_body = 0.04, radius_body = 0.004
#' )
#' create_shape(
#'   "oblate_spheroid", length_body = 0.012, radius_body = 0.01
#' )
#' create_shape(
#'   "cylinder", length_body = 0.05, radius_body = 0.003
#' )
#'
#' @details
#' The \strong{shape} argument specifies what shape for the function to
#' generate the desired shape for TS modeling. Options currently include:
#' \tabular{rlllll}{
#'  \tab \strong{Object shape} \tab \strong{shape = ...} \tab  \tab
#'  \strong{Parameters} \tab \strong{Root function}\cr
#'  \tab \emph{Discrete/tapered cylinder} \tab "cylinder" \tab
#'  \tab length, radius \tab
#'  \code{\link[=cylinder]{cylinder(...)}}\cr
#'  \tab \emph{Polynomial cylinder} \tab "polynomial_cylinder" \tab
#'  \tab length, radius, polynomial \tab
#'  \code{\link[=polynomial_cylinder]{polynomial_cylinder(...)}}\cr
#'  \tab \emph{Oblate spheroid} \tab "oblate_spheroid" \tab
#'  \tab length, radius \tab
#'  \code{\link[=oblate_spheroid]{oblate_spheroid(...)}}\cr
#'  \tab \emph{Prolate spheroid} \tab "prolate_spheroid" \tab
#'  \tab length, radius \tab
#'  \code{\link[=prolate_spheroid]{prolate_spheroid(...)}}\cr
#'  \tab \emph{Sphere} \tab "sphere" \tab \tab radius \tab
#'  \code{\link[=sphere]{sphere(...)}}\cr
#' }
#'
#' \subsection{Model Parameter Definitions}{
#' \itemize{
#'  \item \strong{length}: the x-axis length of the shape.
#'  \item \strong{radius}: the radius of the shape when applicable.
#'  \item \strong{length_radius_ratio}: the length-to-radius ratio (L/A), which
#'  specifically refers to the radius at the mid-point of the cylinder and
#'  should be the maximum value. A typical L/A ratio in the literature is 16
#'  for krill.
#'  \item \strong{taper}: the taper order (n), which parameterizes the tapering
#'  function reported by Chu \emph{et al.} (1993) to create a tapered cylinder.
#'  The tapering order will converge on a prolate and oblate spheroid when
#'  L > 2a and L < 2a, respectively, and n = 2. A typical taper order in the
#'  literature is 10.
#'  \item \strong{polynomial}: the vector of arbitrary polynomial coefficients
#'  to generate a deformed cylinder as reported by Smith \emph{et al.} (2013).
#'  Although listed as a mandatory argument for the polynomial cylinder
#'  function, it has a default setting that uses the sixth-degree polynomial
#'  coefficients reported by Smith \emph{et al.} (2013).
#' }
#' }
#'
#' @return
#' A \code{\link{Shape}} object.
#'
#' @keywords shape_generation
#' @rdname create_shape
#' @export
create_shape <- function(shape, ...) {
  # Pull argument input names ==================================================
  arg_pull <- as.list(match.call())
  # Filter out inappropriate parameters ========================================
  true_args <- .filter_shape_args(shape, arg_pull)
  # Initialize =================================================================
  shape_out <- do.call(shape, true_args)
  # Return shape ===============================================================
  shape_out
}
