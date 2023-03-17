################################################################################
################################################################################
# Mathematical functions
################################################################################
################################################################################
#' Along-matrix summing function
#' @param rpos Position vector
#' @param iterations Number of iterations
#' @rdname along_sum
#' @export
along_sum <- function( rpos , iterations ) {
  output <- rpos[ , 1 : ( iterations - 1 ) ] + rpos[ , 2 : iterations ]
  base::return( output )
}
################################################################################
#' Convert between degrees and radians
#' @param x A real value in degrees or radians
#' for radians.
#' @usage
#' radians( x )
#' degrees( x )
#' @return
#' Converts degrees to radians or radians to degrees
#' @export
radians <- function( x ) x * 180.0 / pi
#' @rdname radians
degrees <- function( x ) x * pi / 180.0
################################################################################
#' Vectorized Euclidean norm function
#' @param x A matrix with numeric, real values.
#' @usage
#' vecnorm(x)
#' @examples
#' values <- matrix(c(1,2,3), ncol=3)
#' vecnorm(values)
#' # 3.741657
#' @return
#' Calculates the Euclidean norm of a vector.
#' @rdname vecnorm
#' @export
vecnorm <- function( x ) {
  base::return( base::sqrt( base::rowSums( x ^ 2 ) ) )
}
