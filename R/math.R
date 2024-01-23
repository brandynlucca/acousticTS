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
  return( output )
}
################################################################################
#' Numerical integration via adaptive quadrature of complex values
#' @param integral Input integration function that is indexed so that it can be
#' used within `base::apply`.
#' @param x Indexing argument for multi-row objects
#' @param y Indexing argument for multi-column objects
#' @rdname contour_integrate
#' @export
contour_integrate <- function( integral , x , y ) {
  complex( real = integrate( function( s ) Re( integral( s , x , y ) ) ,
                             lower = 0 , upper = 1 )$value ,
           imaginary = integrate( function( s ) Im( integral( s , x , y ) ) ,
                                  lower = 0 , upper = 1 )$value )
}
################################################################################
#' Wrapper function incorporating phase deviation into contour integration
#' @param integral Input integration function that is indexed so that it can be
#' used within `base::apply`
#' @param x Indexing argument for multi-row objects
#' @param y Indexing argument for multi-column objects
#' @param n_iterations Number of phase deviations to average and summarize
#' @param integral Integral function used for numerical integration via adaptive
#' quadrature
#' @param phase_sd Phase standard deviation 
#' @rdname phase_integrate
#' @export
phase_integrate <- function( x , y , n_iterations , integral , phase_sd ) {
  rng <- rnorm( n_iterations , 0 , 1 )
  phase <- exp( 1i * rng * phase_sd )
  contour_integrate( integral , x , y ) * phase
}
################################################################################
#' Convert angular measurements from radians to degrees
#' @param x A real value in radians
#' @usage
#' degrees( x )
#' @examples
#' orientation <- pi / 2 # radians
#' degrees( orientation ) # this should return a value equal to 90 degrees
#' @return
#' Angle in degrees
#' @rdname degrees
#' @export
degrees <- function( x ) x * 180.0 / pi
################################################################################
#' Convert angular measurements from degrees to radians.
#' @param x A real value in degrees
#' @usage
#' radians( x )
#' @examples
#' orientation <- 90 # degrees
#' radians( orientation ) # this should return a value equal to pi / 2 radians#' 
#' @return
#' Angle in radians.
#' @rdname radians
#' @export
radians <- function( x ) x * pi / 180.0
################################################################################
#' Calculates the Euclidean norm across each row of a given matrix.
#' @param x A matrix with numeric, real values.
#' @usage
#' vecnorm( x )
#' @examples
#' values <- matrix( c( 1 , 2 , 3 ) , ncol = 3 )
#' vecnorm( values ) # should yield 3.741657
#' @return
#' Calculates the Euclidean norm of a vector.
#' @rdname vecnorm
#' @export
vecnorm <- function( x ) sqrt( rowSums( x * x ) )