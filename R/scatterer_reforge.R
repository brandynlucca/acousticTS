#' Resizing function for targets.
#' @param object Scatterer-class object.
#' @param ... Additional inputs.
#' @rdname reforge
#' @export
setGeneric("reforge", function(object, ...)
  standardGeneric("reforge"))

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
#' @param full_size New body length resize.
#' @export
setMethod( "reforge",
           signature( object = "FLS" ) ,
           function( object,
                     new_length ) {
             lscale <- new_length / acousticTS::extract( object , "shape_parameters" )$length
             body_rpos <- acousticTS::extract( object , "body" )$rpos[ 1 : 3 , ]
             radius <- acousticTS::extract( object , "body" )$radius
             mscale <- base::cbind( base::c( 1 , 0 , 0 ) ,
                                    base::c( 0 , 1 , 0 ) ,
                                    base::c( 0 , 0 , 1 ) ) * lscale
             body_rpos <- base::t( base::t( body_rpos ) %*% mscale )
             radius <- radius * lscale
             slot( object , "body")$rpos <- body_rpos
             slot( object , "body")$radius <- radius
             slot( object , "shape_parameters")$length <- new_length
             return( object )
           } )
