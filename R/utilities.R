################################################################################
################################################################################
# UTILITY FUNCTIONS FOR VARIOUS DATA WRANGLING AND MANIPULATION OPERATIONS
################################################################################
################################################################################
# Accessor functions
################################################################################
################################################################################
#' Primary accessor function for dredging specific data from scatterer objects
#' @param object Scatterer-class object.
#' @param feature Feature of interest (e.g. body).
#' @export
extract <- function( object , feature ) {
  return( methods::slot( object , feature ) )
}
################################################################################
################################################################################
# Uniformly bend body shape
################################################################################
################################################################################
#' Support function for bending scatterer body shape and position matrix
#' @param input Dataframe or scatterer-class object
#' @param radius_curvature Radius of curvature that can be parameterized either 
#' as a ratio relative to body length or actual measurement
#' @param mode Either "ratio" or "measurement"
#' 
#' @rdname brake
#' @export
brake <- function( input , radius_curvature , mode = "ratio" ) {
  # Check object type ==========================================================
  class_type <- typeof( input )
  # Bend shape =================================================================
  output <- switch( class_type ,
                    list = brake_df( input , radius_curvature , mode ) ,
                    S4 = brake_scatterer( input , radius_curvature , mode ) )
  # Return output ==============================================================
  return( output )
}
################################################################################
#' Support function for bending scatterer position matrix dataframe
#' @param body_df Dataframe object containing body shape information
#' @param radius_curvature Radius of curvature that can be parameterized either 
#' as a ratio relative to body length or actual measurement
#' @param mode Either "ratio" or "measurement"
#' 
#' @keywords internal
#' @rdname brake_df
#' @export
brake_df <- function( body_df , radius_curvature , mode = "ratio" ) {
  # Extract object body shape ==================================================
  rpos <- body_df$rpos
  L <- max( rpos[ 1 , ] )
  n_segments <- ncol( rpos )
  # Pull in value for use ======================================================
  radius_curvature_new <- switch( mode ,
                                  ratio = radius_curvature , 
                                  measurement = radius_curvature / L )
  # Check angle-per-segment ====================================================
  arc_angle_per_segment <- L / ( radius_curvature_new * L * ( n_segments - 1 ) )
  # ---- Compare against threshold of 10 (arc) degrees in radians ++++++++++++++
  if ( arc_angle_per_segment > pi / 8 ) {
    warning(
      paste0(
        "Arc angle per segment [", round( arc_angle_per_segment , 5 ) , "] ",
        "exceeds pi/8; increase 'n_segments' [", n_segments, "] ",
        "or decrease the effective radius of curvature ('radius_curvature') ",
        "relative to body length [", radius_curvature_new, "]. Otherwise, ",
        "model predictions may be unstable and/or unreliable."
      )
    )
  }
  # Estimate the position angle of the bent cylinder ===========================
  # Geometry described in Stanton (1989) for total arc length ++++++++++++++++++
  gamma_max <- 0.5 / radius_curvature_new
  theta <- seq(-gamma_max, gamma_max, length.out = n_segments)
  # Rescale to the original body size ==========================================
  x_new <- (radius_curvature_new * L) * sin(theta) + (L / 2)
  z_new <- (radius_curvature_new * L) * (1 - cos(theta))
  # Determine position matrix direction ========================================
  rpos_direction <- ifelse( which.max( body_df$rpos[ 1 , ] ) == 1 ,
                            "REV" ,
                            "FWD" )
  x_direction <- switch( rpos_direction , 
                         REV = rev( x_new ) ,
                         FWD = x_new )
  z_direction <- switch( rpos_direction ,
                         REV = - rev( z_new ) ,
                         FWD = - z_new )
  # Calculate new arc length ===================================================
  arc_lengths <- sapply( 1:( length( x_direction ) - 1 ) , function( i ) {
    chord_length <- sqrt(
      ( x_direction[ i + 1 ] - x_direction[ i ] )^2 + 
        ( z_direction[ i + 1 ] - z_direction[ i ] )^2
    )
    cl_theta <- 2 * asin( chord_length / (2 * ( radius_curvature_new * L ) ) )
    ( radius_curvature_new * L ) * cl_theta
  } )
  total_arc_length <- sum(arc_lengths)
  # Check for segment thickness relative to curvature ==========================
  # ---- Get center coordinates of the osculating circle [Q] +++++++++++++++++++
  Q <- c( L / 2 , radius_curvature_new * L )
  # ---- Calculate the distances +++++++++++++++++++++++++++++++++++++++++++++++
  Q_d <- sqrt( colSums( ( rbind( x_direction , z_direction ) - Q )^2 ) )
  if( any( body_df$radius > Q_d ) ) {
    warning(
      "One or more body segments [n=", sum( body_df$radius > Q_d ), "] have a ",
      "thickness (radius) greater than their distance to the center of the ",
      "osculating circle used to curve the body shape. This means the ",
      "segment(s) will overlap the center of the arc, causing self-",
      "intersection or unrealistic geometry. Consider reducing the segment ",
      "thickness, increasing the radius of curvature to reduce the shape ",
      "bend, or decreasing the number of segments to avoid overlap."
    )
  }
  # Update object ==============================================================
  rpos[ c( 1 , 3 ) , ] <- rbind( x_direction , z_direction )
  body_df_new <- body_df 
  body_df_new$rpos <- rpos
  body_df_new$arc_length <- total_arc_length
  # Return object ==============================================================
  return( body_df_new )
}
################################################################################
#' Support function for bending scatterer body shape scatterer object
#' @param object Scatterer-class object
#' @param radius_curvature Radius of curvature that can be parameterized either 
#' as a ratio relative to body length or actual measurement
#' @param mode Either "ratio" or "measurement"
#' @rdname brake_scatterer
#' @importFrom methods slot<-
#' @export
brake_scatterer <- function( object , radius_curvature , mode = "ratio" ) {
  # Extract object body shape ==================================================
  body <- extract( object , "body" )
  # Pull in value for use ======================================================
  body_curved <- brake_df( body , radius_curvature , mode )
  # Update object ==============================================================
  slot( object , "body" ) <- body_curved
  # Return object ==============================================================
  return( object )
}
################################################################################
#' Support rotation function for KRM (swimbladder)
#' @inheritParams body_rotation
#' @keywords internal
#' @export
bladder_rotation <- function(sum_rpos, rpos, theta, k_length){
  v <- (sum_rpos[1, ]*cos(theta) + sum_rpos[3, ]*sin(theta)) / 2
  v <- matrix(data = rep(v, each = k_length),
              ncol = length(v),
              nrow = k_length)
  delta_u <- diff(rpos[1, ]) * sin(theta)
  return(list(v = v, delta_u = delta_u))
}
################################################################################
#' Support rotating function for KRM (body)
#' @param sum_rpos Summed position matrix
#' @param rpos Position matrix
#' @param theta Orientation angle
#' @param k_length Length of wavenumber vector
#' @keywords internal
#' @export
body_rotation <- function( sum_rpos , rpos , theta , k_length ) {
  dorsal_sum <- base::matrix( data = base::rep( sum_rpos[ 3 , ] , each = k_length ) ,
                              ncol = base::length( sum_rpos[ 3 , ] ) ,
                              nrow = k_length )
  ventral_sum <- base::matrix( data = base::rep( sum_rpos[ 4 , ] , each = k_length ) ,
                               ncol = base::length( sum_rpos[ 4 , ] ) ,
                               nrow = k_length )
  vbU <- ( dorsal_sum * base::cos( theta ) + dorsal_sum * base::sin( theta ) ) / 2
  vbL <- ( ventral_sum * base::cos( theta ) + ventral_sum * base::sin( theta ) ) / 2
  delta_u <- base::diff( rpos[ 1 , ] ) * base::sin( theta )
  return( base::list( vbU = vbU , vbL = vbL , delta_u = delta_u ) )
}
################################################################################
#' Format data for the modal series solution model into the appropriate matrix
#' @param v Vector input.
#' @param limit Modal series limit.
#' @keywords internal
#' @export
modal_matrix <- function( v , limit ) {
  base::return(
    base::matrix(
      data = base::rep( v , each = limit + 1 ) ,
      ncol = base::length( v ) ,
      nrow = limit + 1
    )
  )
}
################################################################################
#' Discretize vector into separate intervals of different length
#' @param x1 Desired vector/interval
#' @param x0 Original or initial vector/interval
#' @export
segmentize <- function(x1, x0) {
  
  # For each consecutive pair of original coordinates
  segments_to_update <- Map(function(start, end) {
    # Find points between this pair
    points_between <- which(x1 < start & x1 > end)
    
    if (length(points_between) > 0) {
      # Redistribute them evenly between start and end
      new_values <- seq(
        from = start, 
        to = end, 
        length.out = length(points_between) + 2)[2:(length(points_between) + 1)]
      list(indices = points_between, values = new_values)
    } else {
      NULL
    }
  }, x0[-length(x0)], x0[-1])
  
  # Apply all updates
  segments_to_update <- segments_to_update[!sapply(segments_to_update, is.null)]
  
  lapply(segments_to_update, function(update) {
    x1[update$indices] <<- update$values
  })
  
  # Return
  return(x1)
}
################################################################################
#' Resample shape for SDWBA model with piecewise constant radius
#' 
#' This function resamples the shape of a fluid-like scatterer (FLS) object for 
#' use in stochastic distorted wave Born approximation (SDWBA) calculations. 
#' The resampling preserves the overall shape of the scatterer while creating 
#' a new representation with the specified number of segments. The radius 
#' assignment uses a stepwise algorithm to maintain piecewise constant radius 
#' values across segments.
#' 
#' @param object FLS-class object to resample
#' @param n_segments Number of segments in the resampled shape
#' 
#' @return FLS object with resampled shape
#' @export
sdwba_resample <- function(object, n_segments) {
  # Check input
  if (!inherits(object, "FLS")) {
    stop("Object must be of class FLS")
  }
  
  # Extract shape data
  body <- extract( object , "body" )
  orig_rpos <- body$rpos
  orig_radius <- body$radius
  orig_x <- orig_rpos[ 1 , ]
  n_orig <- length( orig_x )
  
  # Create new x-coordinate positions (evenly spaced)
  x_new <- seq( body$rpos[1 , 1 ] , 
                body$rpos[1, dim( body$rpos )[ 2 ] ] , 
                length.out = n_segments + 1 )
  
  # Find the closest x-coordinate value
  x_nearest <- sapply( orig_x , function( x ) which.min( abs( x - x_new ) ) )
  
  # Update `x_new`
  x_new[ x_nearest[ 2 : ( n_orig - 1 ) ] ] <- orig_x[ 2 : ( n_orig - 1 ) ]
  
  # Segmentize the new coordinates
  x_new_seg <- segmentize( x_new , orig_x )
  # ---- Initialize the new position matrix
  new_rpos <- rbind( x = x_new_seg )
  
  # Interpolate coordinates using splines (vectorized)
  rpos_interp <- apply( body$rpos[ 2 : 3 , ] , 
                       1 ,
                       function( y ) stats::spline( x = body$rpos[ 1 , ] , 
                                             y = y , 
                                             xout = x_new_seg )$y )
  # ---- And transpose
  new_rpos <- rbind( new_rpos ,
                     t( rpos_interp ) )
  
  # Initialize radius array
  new_radius <- numeric(n_segments + 1)
  decreasing <- orig_x[ 1 ] > orig_x[ length( orig_x ) ]
  # ---- Update stepwise-assigned radius values [inplace operation]
  lapply( 1 : ( n_orig - 1 ), function( i ) {
    if( decreasing ) {
      indices <- which( x_new_seg <= orig_x[ i ] ) + 1
      indices <- indices[ indices <= length( new_radius ) ]
    } else {
      indices <- which( x_new_seg >= orig_x[ i ] & 
                          x_new_seg < orig_x[ i + 1] )
    }
    new_radius[ indices ] <<- body$radius[ i + 1 ]
  })
  
  # Add the last component to the position matrix
  new_rpos <- rbind(
    new_rpos,
    rbind( zL=new_rpos[ 3, ] - new_radius ,
           zU=new_rpos[ 3, ] + new_radius )
  )
  
  # Update object with new shape
  object@body$rpos <- new_rpos
  object@body$radius <- new_radius
  object@shape_parameters$n_segments <- n_segments
  object@shape_parameters$length <- abs( diff( range( new_rpos[ 1 , ] ) ) )
  
  return(object)
}
