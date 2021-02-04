#' Calculates the Euclidean norm of a vector.
#'
#' @param x A vector with numeric, real values.
#' @usage
#' vecnorm(x)
#' @examples
#' values <- c(1,2,3)
#' vecnorm(values)
#' # 3.741657
#' @return
#' Calculates the Euclidean norm of a vector.
#' @export

#Euclidean vector norm
vecnorm <- function(x){sqrt(sum(x**2))}

#' Toggle between radians and degrees.
#'
#' @param x A real value in degrees or radians
#' @param d The value input-type. Two input types: "deg" for degrees and "rad" 
#' for radians.
#' @usage
#' degrad(x,d)
#' @examples
#' x <- 180 #degrees
#' degrad(x, "deg")
#' # 3.141593
#' x <- pi #radians
#' degrad(x, "rad")
#' # 180
#' @return
#' Converts degrees to radians or radians to degrees
#' @export

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
#' # 502.6547
#' @return
#' Calculates the acoustic wavenumber based on the sound speed of water
#' @export

#Calculate acoustic wavenumber based on the sound speed of water
kcalc <- function(f,c){2*pi*f/c}

#' Manipulate scatterer object.
#'
#' @param shape Input scatterer shape.
#' @param curve Curve (boolean, T/F).
#' @param pc Radius of curvature (pc)
#' @param theta Orientation
#' @param length New Length (m)
#' @usage
#' Shapely(shape, curve, pc, theta, length)
#' @return
#' Returns manipulated object.
#' @export
Shapely <- function(shape, curve=shape@curve, pc=0.0, theta=shape@theta, 
                    length=shape@L){
  if(curve == T){
    shape@curve <- T; shape@pc <- ifelse(pc == 0.0, 3.3, pc)
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
    shape@L <- max(shape@rpos[1,])
  }
    return(shape)
}

#' Toggle between log- and linear-domain for backscatter values.
#' @param x A real value in degrees or radians
#' @description
#' To convert from TS to \eqn{\sigma_bs}
#' TS2sigma(x)
#' To convert from \eqn{\sigma_bs} to TS
#' sigma2TS(x)
#' @examples
#' TS <- -85 #dB re: 1 m^2
#' (sigma <- TS2sigma(TS)) #m^2; convert TS to sigma_bs
#' # 3.162278e-09
#' sigma2TS(sigma) #dB re: 1 m^2; convert sigma_bs to TS
#' # -85
#' @return
#' Converts TS to sigma_bs and vice versa.
#' @rdname TS2sigma
#' @export
TS2sigma <- function(x){
  return(10^(x/10))
}

#' Convert sigma to TS.
#' @rdname TS2sigma
#' @export
sigma2TS <- function(x){
  return(10*log10(x))
}

#' Primary accessor function. 
#' @export
pull <- function(object, parameter){
  return(slot(object, parameter))
}

#' Extract position matrix from object for plotting and other purposes.
#' @export
pos_matrix <- function(object){
  out_df <- data.frame(t(pull(object, "rpos")))
  out_df$a <- pull(object, "a")
  return(out_df)
}

#' #' Generic for plotting scatterer-class objects.
#' #' @export
setGeneric("plot", function(x,y,...) standardGeneric("plot"))

#' Method for plotting scatterer-class objects.
#' @export
#' @import graphics
setMethod("plot", signature(x="scattering", y="missing"), 
          definition = function(x, nudge_y=1, ...) {
  coord <- pos_matrix(x)
  z_up <- coord$z - median(coord$z) #center on 0 vertically 
  a_hi <- z_up + coord$a; a_lo <- z_up - coord$a
  par(ask=F)
  plot(coord$x, z_up, type='l', 
       ylim=c(min(a_lo)*(1-(1-nudge_y)), max(a_hi)*nudge_y),
       xlab="Along-body axis (m)",, ylab="Height / Width (m)", ...)
  lines(coord$x, a_lo, lty=2)
  lines(coord$x, a_hi, lty=2)
  segments(y0=a_lo, y1=a_hi,
           x0=coord$x, x1=coord$x, lty=3, col="gray30")
})
