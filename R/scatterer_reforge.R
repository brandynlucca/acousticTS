################################################################################
# FORGE FUNCTIONS FOR MANIPULATING SCATTERER SHAPES 
################################################################################
################################################################################
# PRIMARY FORGE GENERATION FUNCTION 
################################################################################
# Forge scatterer object based on requested inputs 
#' @param x Along-body, or x-axis, vector
#' @param y Across-body, or y-axis, vector
#' @param z Dorsoventral, or z-axis, vector
#' @param a Radius vector 
#' @param scatterer This sets the desired scatterer class, such as FLS or GAS.
#' @param shape This will generate a canonical shape if something other than 
#' "arbitrary" is input. 
#' @export
forge <- function( x , y , z , a , scatterer ,  shape = "arbitrary" ) {
  
}

#' Resizing function for targets.
#' @param object Scatterer-class object.
#' @rdname reforge
#' @export
setGeneric( "reforge" ,
            function( object , ... )
              standardGeneric( "reforge" ) )
#' Resizing function for swimbladdered targets
#' @param object SBF-class object.
#' @param full_size new factored resize of body.
#' @param bladder_size new factored resize of bladder.
#' @export
setMethod("reforge",
          signature(object="SBF"),
          function(object,
                   full_size=1.0,
                   bladder_size=1.0) {

            if(full_size != 1.0){
              lscale <- full_size
              body_rpos <- extract(object, "body")$rpos
              bladder_rpos <- extract(object, "bladder")$rpos
              bladder_origin <- (bladder_rpos[1,1] - min(body_rpos[1, ])) /
                (max(body_rpos[1, ]) - min(body_rpos[1, ]))
              bladder_rpos[1, ] <- bladder_rpos[1, ] - min(bladder_rpos[1, ])
              mscale <- cbind(c(1,0,0,0),
                              c(0,1,0,0),
                              c(0,0,1,0),
                              c(0,0,0,1)) * lscale
              body_rpos <- t(t(body_rpos) %*% mscale)
              bladder_rpos <- t(t(bladder_rpos) %*% mscale)
              bladder_rpos[1, ] <- bladder_rpos[1, ] + bladder_origin*max(body_rpos[1, ])
              slot(object, "body")$rpos <- body_rpos
              slot(object, "bladder")$rpos <- bladder_rpos
            }

            if(bladder_size != 1.0){
              sbscale <- 1 - bladder_size
              body_rpos <- extract(object, "body")$rpos
              bladder_rpos <- extract(object, "bladder")$rpos
              bladder_origin <- (bladder_rpos[1,1] - min(body_rpos[1, ])) /
                (max(body_rpos[1, ]) - min(body_rpos[1, ]))

              vert_dist <- (bladder_rpos[3, ] - bladder_rpos[4, ]) * sbscale
              bladder_rpos[3, ] <- bladder_rpos[3, ] - 0.5*vert_dist
              bladder_rpos[4, ] <- bladder_rpos[4, ] + 0.5*vert_dist
              bladder_rpos[2, ] <- bladder_rpos[2, ] * sbscale

              #diagnostics


              # if(sum(body_rpos[4, body_rpos[1, ] >= min(bladder_rpos[1, ]) &
              #                  body_rpos[1, ] <= max(bladder_rpos[1, ])] -
              #        bladder_rpos[4, ] > 0) > 0 |
              #    sum(bladder_rpos[3, ] -
              #        body_rpos[3, body_rpos[1, ] >= min(bladder_rpos[1, ]) &
              #                  body_rpos[1, ] <= min(bladder_rpos[1, ])] > 0) > 0){
              #   stop("Swimbladder boundary exceeds boundary of scatterer body shape.
              #      Consider adjusting the swimbladder rescaling factor to amend.")
              # }else{
              #   return(object)
              # }

              slot(object, "bladder")$rpos <- bladder_rpos
            }

            return(object)
          })
#' Reforge FLS-class object.
#' @param object FLS-class object.
#' @param length  New body length resize.
#' @param radius New radius size
#' @param n_segments New number of segments
#' @param length_radius_ratio_constant Keep length-to-radius ratio based on new length
#' @export
setMethod( "reforge",
           signature( object = "FLS" ) ,
           function( object ,
                     length = NA , 
                     radius = NA , 
                     length_radius_ratio_constant = T ,
                     n_segments = NA ) {
             ###################################################################
             # Determine rescaling factors =====================================
             # Determine new number of cylinders +++++++++++++++++++++++++++++++
             if( ! missing( n_segments ) ) {
               # Parse shape ===================================================
               shape <- extract( object , "shape_parameters" )
               # Parse body ====================================================
               body <- extract( object , "body" )
               x_new <- seq( from = body$rpos[ 1 , 1 ] , 
                             to = body$rpos[ 1 , shape$n_segments ] , 
                             length.out = n_segments )
               rpos_new <- base::rbind(
                 x_new , 
                 base::t(
                   base::sapply( 2 : base::nrow( body$rpos ) , 
                                 FUN = function( i ) { 
                                   stats::approx(
                                     x = body$rpos[ 1 , ] ,
                                     y = body$rpos[ i , ] ,
                                     xout = x_new
                                   )
                                 }$y
                   )
                 )
               )
               radius_new <- stats::approx(
                 x = body$rpos[ 1 , ] ,
                 y = body$radius ,
                 xout = x_new
               )$y
               # Update metadata +++++++++++++++++++++++++++++++++++++++++++++++
               methods::slot( object , "body" )$rpos <- rpos_new
               methods::slot( object , "body" )$radius <- radius_new
               methods::slot( object , "shape_parameters" )$n_segments <- n_segments
             }
             # Determine new length ++++++++++++++++++++++++++++++++++++++++++
             if( ! base::missing( length ) ) {
               # Parse shape ===================================================
               shape <- acousticTS::extract( object , "shape_parameters" )
               # Parse body ====================================================
               body <- acousticTS::extract( object , "body" )
               new_scale <- length / shape$length
               matrix_rescale <- base::diag( x = 1 , 
                                             nrow = base::nrow( body$rpos ) ,
                                             ncol = base::nrow( body$rpos ) ) * new_scale
               rpos_new <- base::t( base::t( body$rpos ) %*% matrix_rescale )
               # New radius based on constant ratio or adjust ++++++++++++++++++
               if ( length_radius_ratio_constant ) {
                 if( base::missing( radius ) ) {
                   radius_new <- body$radius * new_scale
                 } else {
                   radius_rescale <- radius / shape$radius
                   radius_new <- radius * radius_rescale
                 }
               }
               # Update metadata +++++++++++++++++++++++++++++++++++++++++++++++
               methods::slot( object , "body" )$rpos <- rpos_new
               methods::slot( object , "body" )$radius <- radius_new
               methods::slot( object , "shape_parameters" )$length <- base::max( rpos_new[ 1 , ] )
               methods::slot( object , "shape_parameters" )$radius <- base::max( radius_new )
             }
             # Return object ===================================================
             return( object )
           } )