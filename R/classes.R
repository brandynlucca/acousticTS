#' Fluid-like scatterer (FLS) object/class.
#'
#' @description
#' A S4 class that provides slots to contain relevant animal metadata for parameterizing models for fluid-like scatterers (FLS) partitioned
#' into discretized cylinders. This, specifically, includes a position matrix, radius, material properties (g, h), orientation,
#' animal shape, and body curvature. This class is used within the DWBA and DFCM model functions. In the future, this will also allow for
#' converting one class of scatterer into another for seemless usage for model comparisons.
#'
#' @param rpos A 3xN matrix, where N is the number of columns, or segments, in our scatterer object. The rows represent the x-, y-, and z-axis
#' shape values in that order. Units in meters.
#' @param a An N-length radius vector. Units in meters.
#' @param g Density contrast. Units are dimensionless.
#' @param h Sound speed contrast. Units are dimensionless.
#' @param theta Animal orientation, where \eqn{\pi/2} is broadside incidence. Units in radians.
#' @param curve Body curvature/flexure. Boolean (T/F) which flags the appropriate model formulation.
#' @param pc Radius of curvature ratio. Units are dimensionless.
#' @param L Length of target. Units in meters.
#' @param ncyl Number of segments/cylinders comprising the shape. Units should be an integer value.
#'
#' @details
#' FLS objects can be created using \link[acousticTS]{FLSwrite} and \link[acousticTS]{FLSgenerate}. These objects can also be read using
#' built-in datasets (e.g., \link[acousticTS]{mcgehee}) or using \link[acousticTS]{FLSread}. The \link[acousticTS]{SDWBA}, \link[acousticTS]{SDWBA.sim},
#' \link[acousticTS]{DFCM}, and \link[acousticTS]{Shapely} functions all use FLS-objects for their respective inputs.
#'
#' @seealso
#' \link{SDWBA}
#' \code{\link{SDWBA}}
#' \link[acousticTS]{SDWBA}
#'
#'
#' @export

FLS <- setClass("FLS",slots=c(rpos="matrix",a="numeric",g="numeric",h="numeric",theta="numeric",curve="logical",pc="numeric",L="numeric",ncyl="numeric"))

#' Calls in a *.csv file as a FLS object
#'
#' @param file A *.csv file formatted with the following columns: x, y, z, a, g, h, theta [optional], pc [optional]
#' @usage
#' FLSread(file)
#' @return
#' Calls in an FLS-class object from a *.csv file
#' @export

FLSread <- function(file){
  animal <- read.csv(file, header=T) #Call in *.csv file; assumes headers are present
  return(new("FLS", rpos=as.matrix(rbind(animal$x,animal$y,animal$z)),a=animal$a,
             g=animal$g[1],h=animal$h[1],
             theta=ifelse(length(animal$theta) > 0, animal$theta, pi/2),
             curve=F,pc=0.0,L=max(animal$x),ncyl=length(animal$x)))}

#' Manually generate a FLS object.
#'
#' @param x Vector containing x-axis body (m) shape data.
#' @param y Vector containing y-axis body (m) shape data.
#' @param z Vector containing z-axis body (m) shape data.
#' @param a Vector containing radii (m).
#' @param g Density contrast.
#' @param h Soundspeed contrast
#' @param theta Orientation of the target relative to the transmit source (\eqn{\theta}). Broadside incidence is considered 90 degrees, or pi/2.
#' Default value is pi/2; input should be in radians.
#' @usage
#' FLSgenerate(x,y,z,a,g,h,theta)
#' @examples
#' #Manually parameterize shape
#' x <- seq(1,10,1)*1e-3; y <- rep(0,10); z <- c(seq(1,5,1),rev(seq(1,5,1)))*1e-4
#' a <- z/2
#' g <- 1.036
#' h <- 1.0279
#' new_target <- FLSgenerate(x=x,y=y,z=z,a=a,g=g,h=h)
#' #Let's model where sound speed (c) is 1500 m/s, frequency is 120 kHz, with no phase deviation
#' c <- 1500
#' freq <- 120e3
#' SDWBA(shape=new_target, c=c, frequency=freq)
#' [1] -114.0107
#' @return
#' Calls in an FLS-class object from a *.csv file
#' @export
#Manual shape creation
##Inputs required:
##x,y,z for position matrix, a for radius, g and h for material properties, theta for tilt angle (default = 0.0 deg)
FLSgenerate <- function(x,y,z,a,g,h,theta=pi/2,curve=F,pc=0.0){
  return(new("FLS", rpos=as.matrix(rbind(x,y,z)),a=a,g=g,h=h,theta=theta,curve=curve,pc=pc,L=max(x),ncyl=length(x)))}

#' @export
FLSwrite <- function(shape, filename=paste(getwd(),"/target_shape_",Sys.Date(),".csv",sep="")){
  object <- data.frame(x=shape@x,y=shape@y,z=shape@z,a=shape@a,h=rep(shape@h,shape@ncyl),g=rep(shape@g,shape@ncyl))
  write.csv(shape, file=filename)
}

#' Swimbladdered fish (SBF) object/class.
#'
#' @description
#' A S4 class that provides slots to contain relevant animal metadata for parameterizing models for swimbladdered fish/scatterers (SBF)
#' partitioend into two sets of discretized cylinders: the body and swimbladder. Both shapes comprise a position matrix, material properties
#' (sound speed, c, and density, rho), orientation (theta), the overall body length of the animal, and the number of cylinders. This class is
#' used within the KRM and GBZ model functions.
#' @export
SBF <- setClass("SBF", slots=c(body="matrix",bladder="matrix",theta="numeric",pb="numeric",
                               cb="numeric",psb="numeric",csb="numeric",L="numeric",
                               ncylb="numeric",ncylsb="numeric"))

#' Calls in a *.csv file as a SBF
#'
#' @param file A *.csv file formatted with the following columns: xb, wb, zbU, zbL, xsb, wsb, zsbU, zsbL, pb, cb, psb, csb, theta [optional].
#' @usage
#' SBFread(file)
#' @return
#' Calls in an SBF-class object from a *.csv file
#' @export
SBFread <- function(file){
  animal <- read.csv(file, header=T)
  return(new("SBF",
             body=as.matrix(rbind(animal$xb, animal$wb, animal$zbU, animal$zbL)),
             bladder=as.matrix(rbind(animal$xsb[!is.na(animal$xsb)], animal$wsb[!is.na(animal$xsb)],
                                     animal$zsbU[!is.na(animal$xsb)], animal$zsbL[!is.na(animal$xsb)])),
             pb=animal$pb[1], cb=animal$cb[1], psb=animal$psb[1], csb=animal$csb[1],
             theta=ifelse(length(animal$theta) > 0, animal$theta, pi/2), L=max(animal$xb),
             ncylb=length(animal$xb[!is.na(animal$xb)]),ncylsb=length(animal$xsb[!is.na(animal$xsb)])))}
