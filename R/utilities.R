################################################################################
################################################################################
# Utility functions for data wrangling and organization
################################################################################
################################################################################
#' Support function for bending scatterer body shape and position matrix
#' @param input Dataframe or scatterer-class object
#' @param radius_curvature Radius of curvature that can be parameterized either 
#' as a ratio relative to body length or actual measurement
#' @param mode Either "ratio" or "measurement"
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
#' @rdname brake_df
#' @export
brake_df <- function( body_df , radius_curvature , mode = "ratio" ) {
  # Extract object body shape ==================================================
  rpos <- body_df$rpos
  L <- max( rpos[ 1 , ] )
  n_segments <- ncol( rpos )
  # Pull in value for use ======================================================
  radius_curvature <- switch( mode ,
                              ratio = radius_curvature , 
                              measurement = radius_curvature / L )
  # Estimate the position angle of the bent cylinder ===========================
  # Geometry described in Stanton (1989) for total arc length ++++++++++++++++++
  gamma_max <- 0.5 / radius_curvature
  # Simulate new x-axis over -1 to 1 ===========================================
  x_normal <- seq( from = -1 , to = 1 , length.out = n_segments )
  # Simulate arc positions =====================================================
  x_position <- sin( gamma_max ) * x_normal
  z_position <- 1 - sqrt( 1 - x_position * x_position )
  # Renormalize new x- and z-axis positions ====================================
  x_ratio <- x_position * 1 / gamma_max 
  z_ratio <- z_position * 1 / gamma_max
  # Rescale to the original body size ==========================================
  x_new <- x_ratio * L / 2 + L / 2
  z_new <- z_ratio * L
  # Calculate new arc length ===================================================
  arc_length <- sum( sqrt( diff( x_new ) * diff( x_new ) + diff( z_new ) * diff( z_new ) ) )
  # Determine position matrix direction ========================================
  rpos_direction <- ifelse( which.max( body$rpos[ 1 , ] ) == 1 ,
                            "REV" ,
                            "FWD" )
  x_direction <- switch( rpos_direction , 
                         REV = rev( x_new ) ,
                         FWD = x_new )
  z_direction <- switch( rpos_direction ,
                         REV = - rev( z_new ) ,
                         FWD = - z_new )
  # Update object ==============================================================
  rpos[ c( 1 , 3 ) , ] <- rbind( x_direction , z_direction )
  body_df_new <- body_df 
  body_df_new$rpos <- rpos
  body_df_new$arc_length <- arc_length
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
#' @export
brake_scatterer <- function( object , radius_curvature , mode = "ratio" ) {
  # Extract object body shape ==================================================
  body <- extract( object , "body" )
  shape <- extract( object , "shape_parameters" )
  # Pull in value for use ======================================================
  radius_curvature <- switch( mode ,
                              ratio = radius_curvature , 
                              measurement = radius_curvature / shape$length )
  # Estimate the position angle of the bent cylinder ===========================
  # Geometry described in Stanton (1989) for total arc length ++++++++++++++++++
  gamma_max <- 0.5 / radius_curvature
  # Simulate new x-axis over -1 to 1 ===========================================
  x_normal <- seq( from = -1 , to = 1 , length.out = shape$n_segments )
  # Simulate arc positions =====================================================
  x_position <- sin( gamma_max ) * x_normal
  z_position <- 1 - sqrt( 1 - x_position * x_position )
  # Renormalize new x- and z-axis positions ====================================
  x_ratio <- x_position * 1 / gamma_max 
  z_ratio <- z_position * 1 / gamma_max
  # Rescale to the original body size ==========================================
  x_new <- x_ratio * shape$length / 2 + shape$length / 2
  z_new <- z_ratio * shape$length
  # Calculate new arc length ===================================================
  arc_length <- sum( sqrt( diff( x_new ) * diff( x_new ) + diff( z_new ) * diff( z_new ) ) )
  # Determine position matrix direction ========================================
  rpos_direction <- ifelse( which.max( body$rpos[ 1 , ] ) == 1 ,
                            "REV" ,
                            "FWD" )
  x_direction <- switch( rpos_direction , 
                         REV = rev( x_new ) ,
                         FWD = x_new )
  z_direction <- switch( rpos_direction ,
                         REV = - rev( z_new ) ,
                         FWD = - z_new )
  # Update object ==============================================================
  body$rpos[ c( 1 , 3 ) , ] <- rbind( x_direction , z_direction )
  slot( object , "body" )$rpos <- body$rpos
  slot( object , "body" )$arc_length <- arc_length
  # Return object ==============================================================
  return( object )
}
################################################################################
#' Support rotation function for KRM (swimbladder)
#' @inheritParams body_rotation
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
#' Primary accessor function for dredging specific data from scatterer objects
#' @param object Scatterer-class object.
#' @param feature Feature of interest (e.g. body).
#' @export
extract <- function( object , feature ) {
  base::return( methods::slot( object , feature ) )
}
################################################################################
#' Format data for the modal series solution model into the appropriate matrix
#' @param v Vector input.
#' @param limit Modal series limit.
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