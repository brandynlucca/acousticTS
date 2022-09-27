################################################################################
# FUNCTIONS FOR GENERATIGN CANONICAL & PRE-DEFINED SHAPES
################################################################################
################################################################################
# Sphere
################################################################################
#' Creates a sphere.
#'
#' @param radius Object radius (m).
#' @param n_segments Number of segments to discretize object shape. Defaults to
#'    1e3 segments.
#' @usage
#' sphere(radius, n_segments)
#' @return
#' Creates position vector for a spherical object of a defined radius.
#' @export
sphere <- function(radius,
                   n_segments = 1e3) {
  # Define semi-major or x-axis ================================================
  diameter <- radius * 2
  x_axis <- seq(0, diameter, length.out = n_segments + 1)
  # Map radius along x-axis ====================================================
  radius_out <- sqrt((radius)^2 - (x_axis - radius)^2)
  # Generate position vector, 'rpos' ===========================================
  rpos <- rbind(x = x_axis - radius_out,
                radiusU = radius_out,
                radiusL = -rev(radius_out))
  return(rpos)
}
################################################################################
# Prolate spheroid
################################################################################
#' Creates a prolate spheroid.
#'
#' @param length Semi-major axis length (m).
#' @param radius Semi-minor axis length (m). This can also be stylized as the
#'    "maximum radius" of the scattering object.
#' @param length_radius_ratio Optional ratio input when radius is not explicitly
#'    known.
#' @param n_segments Number of segments to discretize object shape. Defaults to
#'    1e2 segments.
#' @usage
#' prolate_spheroid(length, radius, length_radius_ratio, n_segments)
#' @return
#' Creates the position vector for a prolate spheroid object of defined
#'    semi-major and -minor axes.
#' @export
prolate_spheroid <- function(length,
                             radius,
                             length_radius_ratio = NULL,
                             n_segments = 1e2) {
  # Define maximum radius ======================================================
  if(missing(radius) & !is.null(length_radius_ratio)) {
    max_radius <- length / length_radius_ratio
  }else if(!missing(radius)) {
    max_radius <- radius
  } else {
    stop("Radius/width and/or length-to-radius ratio are missing.")
  }
  # Define semi-major or x-axis ================================================
  x_axis <- seq(0, length, length.out = n_segments + 1)
  # Generate prolate spheroid shape ============================================
  radius_output <- max_radius * sqrt(1 - ((x_axis - length / 2) / (length / 2))^2)
  # Generate position vector, 'rpos' ===========================================
  rpos <- list(rpos = data.frame(x = x_axis,
                                 y = rep(0, length(x_axis)),
                                 z = rep(0, length(x_axis))),
               radius = radius_output)
  return(rpos)
}
################################################################################
# Elongated cylinder
################################################################################
#' Creates a cylinder.
#'
#' @param length Length (m).
#' @param radius Maximum/uniform radius (m).
#' @param length_radius_ratio Optional ratio input when radius is not explicitly
#'    known.
#' @param taper Optional input that is the degree of taper to round ends of
#'    the cylinder.
#' @param n_segments Number of segments to discretize object shape. Defaults to
#'    1e2 segments.
#' @usage
#' cylinder(length, radius, length_radius_ratio, taper, n_segments)
#' @return
#' Creates the position vector for a tapered or untapered cylinder.
#' @export
cylinder <- function(length,
                     radius,
                     length_radius_ratio = NULL,
                     taper = NULL,
                     n_segments = 1e2) {
  # Define maximum radius ======================================================
  if(missing(radius) & !is.null(length_radius_ratio)) {
    max_radius <- length / length_radius_ratio
  }else if(!missing(radius)) {
    max_radius <- radius
  } else {
    stop("Radius/width and/or length-to-radius ratio are missing.")
  }
  # Define normalized x-axis ===================================================
  x_n_axis <- seq(-1, 1, length.out = n_segments + 1)
  # Define tapered radius vector, if applicable
  tapered <- ifelse(!is.null(taper), sqrt(1 - x_n_axis^taper), 1)
  radius_output <- max_radius * tapered
  # Define output x-axis =======================================================
  x_axis <- x_n_axis * length / 2 + length / 2
  # Generate position vector, 'rpos' ===========================================
  rpos <- list(rpos = data.frame(x = x_axis,
                                 y = rep(0, length(x_axis)),
                                 z = rep(0, length(x_axis))),
               radius = radius_output)
  return(rpos)
}
################################################################################
# Polynomial cylinder
################################################################################
#' Creates a polynomial deformed cylinder.
#'
#' @param polynomial Polynomial coefficient vector.
#' @inheritParams cylinder
#' @usage
#' polynomial_cylinder(length, radius, n_segments, polynomial)
#' @examples
#' \dontrun{
#' # We can use the polynomial coefficients defined in Smith et al. (2013) to
#' # define the position vector of a sub-Arctic krill.
#' poly_vec <- c(0.83, 0.36, -2.10, -1.20, 0.63, 0.82, 0.64)
#' # Create the position vector
#' # This outputs a list containing "rpos" and "radius"
#' pos <- polynomial_cylinder(length = 15e-3, radius = 2e-3, polynomial = poly_vec)
#' str(pos)
#' }
#' @return
#' Creates the position vector for a polynomial deformed cylinder.
#' @references
#' Smith, J.N., Ressler, P.H., and Warren, J.D. 2013. A distorted wave Born
#' approximation target strength model for Bering Sea euphausiids. ICES Journal
#' of Marine Science, 70(1): 204-214. https://doi.org/10.1093/icesjms/fss140
#' @export
polynomial_cylinder <- function(length,
                                radius,
                                n_segments = 1e3,
                                polynomial) {
  # Define normalized x-axis ===================================================
  x_n_axis <- seq(-1, 1, length.out = n_segments + 1)
  # Evaluate polynomial coefficients ===========================================
  n_order <- rev(seq_len(length(polynomial))) - 1
  poly_fun <- paste0(polynomial, paste0("*x_n_axis^", n_order), collapse = "+")
  # Define radius ==============================================================
  radius_output <- abs(eval(parse(text = poly_fun))) * radius
  # Define output x-axis =======================================================
  x_axis <- x_n_axis * length / 2 + length / 2
  # Generate position vector, 'rpos' ===========================================
  rpos <- list(rpos = data.frame(x = x_axis,
                                 y = rep(0, length(x_axis)),
                                 z = rep(0, length(x_axis))),
               radius = radius_output)
  return(rpos)
}
################################################################################
# Wrapper function for generating shape
################################################################################
#' A wrapper function that automatically creates generalized and/or canonical
#' shapes for TS modeling.
#'
#' @param shape Shape. Details for shape specification are provided under
#'    'Details', including mandatory arguments.
#' @param ... Additional input arguments for subsequent shape generation functions.
#'
#' @details
#' The \strong{shape} argument specifies what shape for the function to generate into
#' the desired shape for TS modeling. Options currently include:
#' \tabular{rlllll}{
#'  \tab \strong{Object shape} \tab \strong{shape = ...} \tab  \tab
#'  \strong{Parameters} \tab \strong{Root function}\cr
#'  \tab \emph{Discrete/tapered cylinder} \tab "cylinder" \tab \tab length, radius \tab
#'  \code{\link[=cylinder]{cylinder(...)}}\cr
#'  \tab \emph{Polynomial cylinder} \tab "polynomial_cylinder" \tab \tab length, radius, polynomial \tab
#'  \code{\link[=polynomial_cylinder]{polynomial_cylinder(...)}}\cr
#'  \tab \emph{Prolate spheroid} \tab "prolate_speroid" \tab \tab length, radius \tab
#'  \code{\link[=prolate_spheroid]{prolate_spheroid(...)}}\cr
#'  \tab \emph{Sphere} \tab "sphere" \tab \tab radius \tab
#'  \code{\link[=sphere]{sphere(...)}}\cr
#' }
#'
#' \subsection{Model Parameter Definitions}{
#' \itemize{
#'  \item \strong{length}: the x-axis length of the shape.
#'  \item \strong{radius}: the radius of the shape when applicable.
#'  \item \strong{length_radius_ratio}: the length-to-radius ratio (L/A), which specifically
#'  refers to the radius at the mid-point of the cylinder and should be
#'  the maximum value. A typical L/A ratio in the literature is 16 for krill.
#'  \item \strong{taper}: the taper order (n), which parameterizes the tapering
#'  function reported by Chu \emph{et al.} (1993) to create a tapered cylinder.
#'  The tapering order will converge on a prolate and oblate spheroid when
#'  L > 2a and L < 2a, respectively, and n = 2. A typical taper order in the
#'  literature is 10.
#'  \item \strong{polynomial}: the vector of arbitrary polynomial coefficients to
#'  generate a deformed cylinder as reported by Smith \emph{et al.} (2013).
#'  Although listed as a mandatory argument for the polynomial cylinder
#'  function, it has a default setting that uses the sixth-degree polynomial
#'  coefficients reported by Smith \emph{et al.} (2013).
#' }
#' }
#'
#' @return
#' Chu, D., Foote, K.G., and Stanton, T.K. 1993. Further analysis of target
#' strength measurements of Antarctic krill at 38 and 120 kHz: Comparison and
#' deformed cylinder model and inference of orientation distribution. The
#' Journal of the Acoustical Society of America, 93(5): 2985-2988.
#' https://doi.org/10.1121/1.405818
#'
#' Smith, J.N., Ressler, P.H., and Warren, J.D. 2013. A distorted wave Born
#' approximation target strength model for Bering Sea euphausiids. ICES Journal
#' of Marine Science, 70(1): 204-214. https://doi.org/10.1093/icesjms/fss140
#' @export
create_shape <- function(shape,
                         ...) {
  # Define shape output ========================================================
  shape_out <- switch(shape,
                      sphere = sphere(...),
                      prolate_sphoeroid = prolate_spheroid(...),
                      cylinder = cylinder(...),
                      polynomial_cylinder = polynomial_cylinder(...))
  return(shape_out)
}
