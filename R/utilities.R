################################################################################
################################################################################
# Utility functions for data wrangling and organization
################################################################################
################################################################################
#' Support rotatin function for KRM (swimbladder)
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
#'
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
