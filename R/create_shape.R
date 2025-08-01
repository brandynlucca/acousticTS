################################################################################
################################################################################
# FUNCTIONS FOR GENERATIGN CANONICAL & PRE-DEFINED SHAPES
################################################################################
################################################################################
################################################################################
# Arbitrary (or pre-generated) body shape parameters
################################################################################
#' Creates arbitrary body shape from user inputs
#' @param x_body x-axis (m)
#' @param y_body y-axis (m)
#' @param z_body z-axis (m)
#' @param radius_body Radius (m)
#' @param length_units Units for body length. Defaults to meters: "m"
#' @seealso \code{\link{Arbitrary}}
#' 
#' @keywords shape_generation
#' @rdname arbitrary
#' @export
arbitrary <- function( x_body ,
                       y_body ,
                       z_body ,
                       radius_body ,
                       length_units = "m" ) {
  position_matrix <- cbind( x = x_body ,
                            y = y_body ,
                            z = z_body ,
                            zU = z_body + radius_body ,
                            zL = z_body - radius_body )
  # Generate shape parameters list =============================================
  shape_parameters <- list( radius = radius_body , 
                            n_segments = length( x_body ) - 1 ,
                            diameter_units = length_units )
  # Generate new shape object ==================================================
  return( new( "Arbitrary" ,
               position_matrix = position_matrix ,
               shape_parameters = shape_parameters ) )
}
################################################################################
# Sphere
################################################################################
#' Creates a sphere.
#' @param radius_body Object radius (m).
#' @param n_segments Number of segments to discretize object shape. Defaults to
#'    1e2 segments.
#' @param diameter_units Default is "m" for meters
#' @usage
#' sphere(radius_body, n_segments, diameter_units )
#' @return
#' Creates position vector for a spherical object of a defined radius.
#' @seealso \code{\link{Sphere}}
#' 
#' @keywords shape_generation
#' @rdname sphere
#' @importFrom utils head tail
#' @export
sphere <- function( radius_body ,
                    n_segments = 1e2 ,
                    diameter_units = "m" ) {
  # Define semi-major or x-axis ================================================
  diameter <- radius_body * 2
  semi_major <- seq( from = 0 ,
                     to = diameter ,
                     length.out = n_segments + 1 )
  # Along-semimajor radii ======================================================
  along_radius <- sqrt( radius_body * radius_body - ( semi_major - radius_body ) * ( semi_major - radius_body ) )
  # Segment midpoints ==========================================================
  x_mids <- ( head( semi_major , -1 ) + tail( semi_major , -1 ) ) / 2
  # Apply to the radii =========================================================
  radius_output <-( head( along_radius , -1 ) + tail( along_radius , -1 ) ) / 2
  # Assign "zeroth" radius =====================================================
  if ( which.max( semi_major ) == length( semi_major ) ) {
    x_edges <- rev( semi_major ) 
    radius_output <- rev( radius_output )
  }
  radius_output <- c( 0.0, radius_output )
  # Generate position matrix ===================================================
  position_matrix <- cbind( x = x_edges ,
                            y = rep( 0 , length( x_edges ) ) ,
                            z = rep( 0 , length( x_edges ) ) ,
                            zU = radius_output ,
                            zL = - radius_output )
  # Generate shape parameters list =============================================
  shape_parameters <- list(
    radius_shape = radius_output ,
    diameter_shape = diameter ,
    radius = diameter / 2 ,
    n_segments = n_segments ,
    diameter_units = diameter_units )
  # Generate new shape object ==================================================
  return( new( "Sphere" ,
               position_matrix = position_matrix ,
               shape_parameters = shape_parameters ) )
}
################################################################################
# Prolate spheroid
################################################################################
#' Creates a prolate spheroid.
#'
#' @param length_body Semi-major axis length (m).
#' @param radius_body Semi-minor axis length (m). This can also be stylized as the
#'    "maximum radius" of the scattering object.
#' @param length_radius_ratio Optional ratio input when radius is not explicitly
#'    known.
#' @param n_segments Number of segments to discretize object shape. Defaults to
#'    18 segments.
#' @param length_units Units for body matrix (defaults to m).
#' @param theta_units Units for body orientation (defaults to radians).
#' @return
#' Creates the position vector for a prolate spheroid object of defined
#'    semi-major and -minor axes.
#' @seealso \code{\link{ProlateSpheroid}}
#' 
#' @keywords shape_generation
#' @rdname prolate_spheroid
#' @importFrom utils head tail
#' @export
prolate_spheroid <- function( length_body ,
                              radius_body ,
                              length_radius_ratio = NULL ,
                              n_segments = 18 ,
                              length_units = "m" ,
                              theta_units = "radians" ) {
  # Define maximum radius ======================================================
  if ( missing( radius_body ) & ! is.null( length_radius_ratio ) ) {
    max_radius <- length_body / length_radius_ratio
  } else if ( ! missing( radius_body ) ) {
    max_radius <- radius_body
  } else {
    stop("Radius/width and/or length-to-radius ratio are missing.")
  }
  # Define semi-major or x-axis ================================================
  # x_axis <- seq( 0 , length_body , length.out = n_segments + 1 )
  x_edges <- seq( 0 , length_body , length.out = n_segments + 1 )
  # Get the segment midpoints ==================================================
  x_mids <- ( head( x_edges , -1 ) + tail( x_edges , -1 ) ) / 2
  # Normalize the midpoints for ellipse ========================================
  curved_x_mids <- ( x_mids - length_body / 2 ) / ( length_body / 2 )
  # Compute the radius at each segment midpoint ================================
  radius_output <- max_radius * sqrt( 1 - curved_x_mids ^ 2 )
  # Generate prolate spheroid shape ============================================
  # curved_x_axis <- ( ( x_axis - length_body / 2 ) / ( length_body / 2 ) ) 
  # radius_output <- max_radius * sqrt( 1 - curved_x_axis * curved_x_axis )
  # Assign "zeroth" radius =====================================================
  if ( which.max( x_edges ) == length( x_edges ) ) {
    x_edges <- rev( x_edges ) 
    radius_output <- rev( radius_output )
  }
  radius_output <- c( 0.0 , radius_output )
  # Generate position matrix ===================================================
  position_matrix <- cbind( x = x_edges , 
                            y = rep( 0 , length( x_edges ) ) ,
                            z = rep( 0 , length( x_edges ) ) ,
                            zU = radius_output ,
                            zL = - rev( radius_output ) )
  # Generate shape parameters list =============================================
  shape_parameters <- list(
    length = max( position_matrix[ , 1 ] ) ,
    radius = radius_output ,
    length_radius_ratio = max( position_matrix[ , 1 ] ) / max( radius_output ) ,
    n_segments = n_segments ,
    length_units = length_units
  )
  # Generate new shape object ==================================================
  return( new( "ProlateSpheroid" ,
               position_matrix = position_matrix ,
               shape_parameters = shape_parameters ) )
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
#' @param n_segments Number of segments to discretize object shape. Defaults to
#'    1e2 segments.
#' @param length_units Units (default is meters, "m").
#' @usage
#' cylinder(length_body, radius_body, length_radius_ratio, 
#' taper, n_segments, length_units)
#' @return
#' Creates the position vector for a tapered or untapered cylinder.
#' @seealso \code{\link{Cylinder}}
#' 
#' @keywords shape_generation
#' @rdname cylinder
#' @importFrom utils head tail
#' @export
cylinder <-  function( length_body ,
                       radius_body ,
                       length_radius_ratio ,
                       taper ,
                       n_segments = 1e2 ,
                       length_units = "m" ) {
  # Define maximum radius ======================================================
  if ( missing( radius_body ) ) {
    max_radius <- length_body / length_radius_ratio
  } else if ( ! missing( radius_body ) ) {
    max_radius <- radius_body
  } else {
    stop("Radius/width and/or length-to-radius ratio are missing.")
  }
  # Define normalized x-axis ===================================================
  x_n_axis <- seq( -1 , 1 , length.out = n_segments + 1 )
  # Define tapered radius vector, if applicable ================================
  if ( ! missing( taper ) ) {
    tapering <- sqrt( 1 - x_n_axis ^ taper )
  } else {
    tapering <- rep( 1 , n_segments + 1 )
  }
  radius_tapered <- max_radius * tapering
  # Segment midpoints ==========================================================
  x_edges <- x_n_axis * length_body / 2 + length_body / 2
  x_mids <- ( head( x_edges , -1 ) + tail( x_edges , -1 ) ) / 2
  # Apply to the radii =========================================================
  radius_mids <-( head( radius_tapered , -1 ) + tail( radius_tapered , -1 ) ) / 2
  # Assign "zeroth" radius =====================================================
  if ( which.max( x_edges ) == length( x_edges ) ) {
    x_edges <- rev( x_edges ) 
    radius_mids <- rev( radius_mids )
  }
  radius_output <- c( 0.0 , radius_mids )
  # Generate position matrix ===================================================
  position_matrix <- cbind( x = x_edges ,
                            y = rep( 0 , length( x_edges ) ) ,
                            z = rep( 0 , length( x_edges ) ) ,
                            zU = radius_output ,
                            zL = - radius_output )
  # Generate shape parameters list =============================================
  shape_parameters <- list(
    length = max( position_matrix[ , 1 ] ) ,
    radius = radius_output ,
    length_radius_ratio = max( position_matrix[ , 1 ] ) /
      max( radius_output ) ,
    n_segments = n_segments ,
    taper_order = ifelse( missing( taper ) ,
                          NA ,
                          taper ) ,
    length_units = length_units
  )
  # Generate new shape object ==================================================
  return( new( "Cylinder" ,
               position_matrix = position_matrix ,
               shape_parameters = shape_parameters ) )
}
################################################################################
# Polynomial cylinder
################################################################################
#' Creates a polynomial deformed cylinder
#'
#' @inheritParams cylinder
#' @param polynomial Polynomial coefficient vector.
#' @usage
#' polynomial_cylinder(length_body, radius_body, n_segments, polynomial, 
#' length_units)
#' @examples
#' \dontrun{
#' # We can use the polynomial coefficients defined in Smith et al. (2013) to
#' # define the position vector of a sub-Arctic krill.
#' poly_vec <- c(0.83, 0.36, -2.10, -1.20, 0.63, 0.82, 0.64)
#' # Create the position vector
#' # This outputs a list containing "rpos" and "radius"
#' pos <- polynomial_cylinder(length_body = 15e-3, radius_body = 2e-3, polynomial = poly_vec)
#' str(pos)
#' }
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
                                polynomial ,
                                length_units = "m" ) {
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
  position_matrix <- cbind( x = x_axis ,
                            y = rep( 0 , length( x_axis ) ) ,
                            z = rep( 0 , length( x_axis ) ) )
  # Generate shape parameters list =============================================
  shape_parameters <- list(
    length = max( position_matrix[ , 1 ] ) ,
    radius = radius_output ,
    length_radius_ratio = max( position_matrix[ , 1 ] ) /
      max( radius_output ) ,
    n_segments = n_segments ,
    length_units = length_units
  )
  # Generate new shape object ==================================================
  return( new( "PolynomialCylinder" ,
               position_matrix = position_matrix ,
               shape_parameters = shape_parameters ) )
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
#'  \tab \emph{Prolate spheroid} \tab "prolate_spheroid" \tab \tab length, radius \tab
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
#' 
#' @keywords shape_generation
#' @rdname create_shape
#' @export
create_shape <- function( shape , ... ) {
  # Pull argument input names ==================================================
  arg_pull <- as.list( match.call( ) )
  # Grab input arguments =======================================================
  arg_list <- names( formals( shape ) )
  # Filter out inappropriate parameters ========================================
  arg_full <- arg_pull[ arg_list ] 
  true_args <- Filter( Negate( is.null) , arg_full )
  # Initialize =================================================================
  shape_out <- do.call( shape , true_args )
  # Return shape ===============================================================
  return( shape_out )
}