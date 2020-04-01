# All formulas for calculation of TS derived from:
#

#' Calculates the Euclidean norm of a vector.
#'
#' @param x A vector with numeric, real values.
#' @usage
#' vecnorm(x)
#' @examples
#' values <- c(1,2,3)
#' vecnorm(values)
#' [1] 3.741657
#' @return
#' Calculates the Euclidean norm of a vector.
#' @export

#Euclidean vector norm
vecnorm <- function(x){sqrt(sum(x**2))}

#' Toggle between radians and degrees.
#'
#' @param x A real value in degrees or radians
#' @param d The value input-type. Two input types: "deg" for degrees and "rad" for radians.
#' @usage
#' degrad(x,d)
#' @examples
#' x <- 180 #degrees
#' degrad(x, "deg")
#' [1] 3.141593
#' x <- pi #radians
#' degrad(x, "rad")
#' [1] 180
#' @return
#' Converts degrees to radians or radians to degrees
#' @export

#Toggle between radians and dgrees
degrad <- function(x,d){
  if(d == "deg"){
    value <- x*pi/180.0
  }else if(d == "rad"){
    value <- x*180.0/pi
  }
  return(value)
}

#' Calculate the acoustic wavenumber based on the sound speed of water.
#'
#' @param c Sound speed (m/s)
#' @param f Frequency (Hz)
#' @usage
#' kcalc(f,c)
#' @examples
#' c <- 1500 #m/s
#' f <- 120e3 #Hz
#' kcalc(f,c)
#' [1] 502.6547
#' @return
#' Calculates the acoustic wavenumber based on the sound speed of water
#' @export

#Calculate acoustic wavenumber based on the sound speed of water
kcalc <- function(f,c){2*pi*f/c}

#' Resize animal to maintain shape based on length.
#'
#' @param shape Desired object/animal shape.
#' @param length New length (m).
#' @usage
#' resize(shape, length)
#' @return
#' Rescales the shape of an animal based on a desired length.
#' @export

resize <- function(shape, length){
  lscale <- length/max(shape@rpos[1,]) #grab current length of shape and calculate scale ratio
  mscale <- cbind(c(1,0,0),c(0,1,0),c(0,0,1)) * lscale #calculate position matrix scale
  shape@rpos <- t(t(shape@rpos) %*% mscale) #rescale length of shape
  shape@a <- shape@a * lscale #scale radius based on same length ratio
  return(shape)
}


#' Manipulate object
#'
#' @param shape Desired object/animal shape.
#' @param curve Curve (boolean, T/F).
#' @param pc Radius of curvature (pc)
#' @param theta Orientation
#' @param length New Length (m)
#' @usage
#'
#' @return
#' Returns manipulated object.
#' @export
Shapely <- function(shape, curve=F, pc=0.0, theta=shape@theta, length=shape@L){
  if(curve == T){
    shape@curve <- T; shape@pc <- ifelse(pc == 0.0, 3.0, pc)
  }else if(curve == F & shape@curve == T){
    shape@curve <- F
  }
  if(shape@theta != theta){
    shape@theta <- theta
  }
  if(shape@L != length){
    lscale <- length/max(shape@rpos[1,])
    mscale <- cbind(c(1,0,0),c(0,1,0),c(0,0,1)) * lscale
    shape@rpos <- t(t(shape@rpos) %*% mscale)
    shape@a <- shape@a * lscale
  }
    return(shape)
}
