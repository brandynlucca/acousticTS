################################################################################
# FORGE FUNCTIONS FOR MANIPULATING SCATTERER SHAPES 
################################################################################
################################################################################
# PRIMARY FORGE GENERATION FUNCTION 
################################################################################
#' Resizing function for targets.
#' @param object Scatterer-class object.
#' @param length_body Updated body length when applicable
#' @param width_body Updated body width when applicable
#' @param height_body Update body height/depth when applicable
#' @param length_bladder Updated bladder length when applicable
#' @param width_bladder Updated bladder width when applicable
#' @param height_bladder Updated bladder height/depth when applicable
#' @param radius_body Updated body radius when applicable
#' @param radius_bladder Updated bladder radius when applicable 
#' @param bladder_inflation_factor Proportional bladder volume
#' @param body_bladder_ratio_constant Maintain body-to-bladder volume ratio when resizing
#' @param length_radius_ratio_constant Maintain body-to-radius ratio when resizing
#' @param n_segments_body Number of segments along the body
#' @param n_segments_bladder Number of segments along the bladder
#' @param ... Scatterer-specific arguments
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
setMethod( "reforge",
           signature( object = "SBF" ) ,
           function( object ,
                     length_body = NA ,
                     width_body = NA ,
                     height_body = NA ,
                     length_bladder = NA ,
                     width_bladder = NA ,
                     height_bladder = NA ,
                     radius_body = NA ,
                     bladder_inflation_factor = 1.0 ,
                     # body_bladder_ratio_constant = T ,
                     isometric_body = T ,
                     isometric_bladder = T ,
                     n_segments_body = NA , 
                     n_segments_bladder = NA ) {
             ###################################################################
             # Determine rescaling factors =====================================
             # Parse body ======================================================
             body <- extract( object , "body" )
             rpos_b <- body$rpos
             # Parse bladder ===================================================
             bladder <- extract( object , "bladder" )
             rpos_sb <- bladder$rpos
             # Parse shape =====================================================
             shape <- extract( object , "shape_parameters" )
             # Number of segments ==============================================
             # Body ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
             if ( ! is.na( n_segments_body ) ) {
               x_body_new <- seq( from = rpos_b[ 1 , 1 ] ,
                                  to = rpos_b[ 1 , shape$body$n_segments ] ,
                                  length.out = n_segments_body )
               rpos_body_new <- rbind( x_body_new ,
                                       t.default(
                                         sapply( 2 : nrow( rpos_b ) ,
                                                 FUN = function( i ) {
                                                   approx( x = rpos_b[ 1 , ] ,
                                                           y = rpos_b[ i , ] ,
                                                           xout = x_body_new ) }$y ) ) )
               rpos_b <- rpos_body_new
             }
             # Bladder ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
             if ( ! is.na( n_segments_bladder ) ) {
               x_bladder_new <- seq( from = rpos_sb[ 1 , 1 ] ,
                                     to = rpos_sb[ 1 , shape$bladder$n_segments ] ,
                                     length.out = n_segments_bladder )
               rpos_bladder_new <- rbind( x_bladder_new ,
                                          t.default(
                                            sapply( 2 : nrow( rpos_sb ) ,
                                                    FUN = function( i ) {
                                                      approx( x = rpos_sb[ 1 , ] ,
                                                              y = rpos_sb[ i , ] ,
                                                              xout = x_bladder_new ) }$y ) ) )
               rpos_sb <- rpos_bladder_new
             }
             # Determine new length & radius vectors +++++++++++++++++++++++++++
             # Body ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
             body_height <- max( rpos_b[ 3 , ] - rpos_b[ 4 , ] )
             body_dims <- c( shape$body$length , 
                             max( body$rpos[ 2 , ] ) , 
                             max( body_height ) ,
                             max( body_height ) )
             body_new <- c( length_body , 
                            width_body , 
                            height_body , 
                            height_body )
             body_dims_rat <- c( 1 , 1 , 1 , 1 )
             bladder_height <- rpos_sb[ 3 , ] - rpos_sb[ 4 , ]
             bladder_dims <- c( max( bladder$rpos[ 1 , ] ) - min( bladder$rpos[ 1 , ] ) , 
                                max( bladder$rpos[ 2 , ] ) , 
                                max( bladder_height ) ,
                                max( bladder_height ) )
             bladder_new <- c( length_bladder , 
                               width_bladder , 
                               height_bladder ,
                               height_bladder )
             bladder_dims_rat <- c( 1 , 1 , 1 , 1 )
             
             if ( any( is.na( body_new ) ) ) {
               body_idx <- which( ! is.na( body_new ) )
               
               if( isometric_body & ! is.na( body_new[ 1 ] ) ) {
                 molt <- body_new[ 1 ] / body_dims[ 1 ]
                 rpos_b <- t( t( rpos_b ) %*% diag( molt , 
                                                    nrow = nrow( rpos_b ) , 
                                                    ncol = nrow( rpos_b ) ) )
               } else {
                 molt <- body_dims_rat
                 molt[ body_idx ] <- body_new[ body_idx ] / body_dims[ body_idx ]
                 rpos_b <- t( t( rpos_b ) %*% diag( molt , 
                                                    nrow = nrow( rpos_b ) , 
                                                    ncol = nrow( rpos_b ) ) )
               }
             } else {
               molt <- body_dims_rat
               molt <- body_new / body_dims
               rpos_b <- t( t( rpos_b ) %*% diag( molt , 
                                                  nrow = nrow( rpos_b ) , 
                                                  ncol = nrow( rpos_b ) ) )
             }
             # Bladder ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
             if ( any( is.na( bladder_new ) ) ) {
               bladder_idx <- which( ! is.na( bladder_new ) )
               
               if( isometric_bladder & ! is.na( bladder_new[ 1 ] ) ) {
                 molt <- bladder_new[ 1 ] / bladder_dims[ 1 ]
                 rpos_sb <- t( t( rpos_sb ) %*% diag( molt , 
                                                      nrow = nrow( rpos_sb ) , 
                                                      ncol = nrow( rpos_sb ) ) )
               } else if( any( ! is.na( bladder_new ) ) ) {
                 molt <- bladder_dims_rat
                 molt[ bladder_idx ] <- bladder_new[ bladder_idx ] / bladder_dims[ bladder_idx ]
                 rpos_sb <- t( t( rpos_sb ) %*% diag( molt , 
                                                      nrow = nrow( rpos_sb ) , 
                                                      ncol = nrow( rpos_sb ) ) )
               } else if( any( ! is.na( body_new ) ) ) {
                 molt <- body_new[ 1 ] / body_dims[ 1 ] 
                 rpos_sb <- t( t( rpos_sb ) %*% diag( molt , 
                                                      nrow = nrow( rpos_sb ) , 
                                                      ncol = nrow( rpos_sb ) ) )
               }
             } else {
               molt <- bladder_dims_rat
               molt <- bladder_new / bladder_dims
               rpos_sb <- t( t( rpos_sb ) %*% diag( molt , 
                                                    nrow = nrow( rpos_sb ) , 
                                                    ncol = nrow( rpos_sb ) ) )
             }
             # Bladder inflation factor ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
             x_bladder_origin <- bladder$rpos[ 1 , 1 ] / max( body$rpos[ 1 , ] )
             xsb_start <- x_bladder_origin * max( rpos_b[ 1 , ] )
             xsb_offset <- rpos_sb[ 1 , 1 ] - xsb_start
             bladder_refill <- bladder_height * bladder_inflation_factor - bladder_height
             rpos_sb[ 1 , ] <- rpos_sb[ 1 , ] - xsb_offset
             rpos_sb[ 2 , ] <- rpos_sb[ 2 , ] * bladder_inflation_factor
             rpos_sb[ 3 , ] <- rpos_sb[ 3 , ] * bladder_inflation_factor
             rpos_sb[ 4 , ] <- rpos_sb[ 4 , ] * bladder_inflation_factor
             # Update object ===================================================
             slot( object , "body" )$rpos <- rpos_b
             slot( object , "bladder" )$rpos <- rpos_sb
             slot( object , "shape_parameters" )$body$length <- max( rpos_b[ 1, ] )
             slot( object , "shape_parameters" )$bladder$length <- max( rpos_sb[ 1 , ] )
             slot( object , "shape_parameters" )$body$n_segments <- length( rpos_b[ 1 , ] )
             slot( object , "shape_parameters" )$bladder$n_segments <- length( rpos_sb[ 1 , ] )
             # Return object ===================================================
             return( object )
           } )
################################################################################
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
                     length , 
                     radius , 
                     length_radius_ratio_constant = T ,
                     n_segments ) {
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
               rpos_new <- rbind(
                 x_new , 
                 t.default(
                   sapply( 2 : nrow( body$rpos ) , 
                           FUN = function( i ) { 
                             approx(x = body$rpos[ 1 , ] ,
                                    y = body$rpos[ i , ] ,
                                    xout = x_new ) }$y
                   )
                 )
               )
               radius_new <- approx( x = body$rpos[ 1 , ] ,
                                     y = body$radius ,
                                     xout = x_new )$y
               # Update metadata +++++++++++++++++++++++++++++++++++++++++++++++
               slot( object , "body" )$rpos <- rpos_new
               slot( object , "body" )$radius <- radius_new
               slot( object , "shape_parameters" )$n_segments <- n_segments
             }
             # Determine new length ++++++++++++++++++++++++++++++++++++++++++
             if( ! missing( length ) ) {
               # Parse shape ===================================================
               shape <- acousticTS::extract( object , "shape_parameters" )
               # Parse body ====================================================
               body <- acousticTS::extract( object , "body" )
               new_scale <- length / shape$length
               matrix_rescale <- diag( x = 1 , 
                                       nrow = nrow( body$rpos ) ,
                                       ncol = nrow( body$rpos ) ) * new_scale
               rpos_new <- t.default( t.default( body$rpos ) %*% matrix_rescale )
               # New radius based on constant ratio or adjust ++++++++++++++++++
               if ( length_radius_ratio_constant ) {
                 if( missing( radius ) ) {
                   radius_new <- body$radius * new_scale
                 } else {
                   radius_rescale <- radius / shape$radius
                   radius_new <- body$radius * radius_rescale
                 }
               }
               # Update metadata +++++++++++++++++++++++++++++++++++++++++++++++
               slot( object , "body" )$rpos <- rpos_new
               slot( object , "body" )$radius <- radius_new
               slot( object , "shape_parameters" )$length <- max( rpos_new[ 1 , ] )
               slot( object , "shape_parameters" )$radius <- max( radius_new )
             }
             # Return object ===================================================
             return( object )
           } )