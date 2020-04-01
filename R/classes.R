#' Fluid-filled scatterer (FFS) object/class.
#'
#' @description
#' A S4 class that provides slots to contain relevant animal metadata for parameterizing models for fluid-filled scatterers (FFS) partitioned
#' into discretized cylinders. This, specifically, includes a position matrix, radius, material properties (g, h), orientation,
#' animal shape, and body curvature. This class is used within the DWBA and DFCM model functions. In the future, this will also allow for
#' converting one class of scatterer into another for seemless usage for model comparisons.
#' @export

FFS <- setClass("FFS",slots=c(rpos="matrix",a="numeric",g="numeric",h="numeric",theta="numeric",shape="character",pc="numeric",L="numeric",ncyl="numeric"))

#' Calls in a *.csv file as a FFS object
#'
#' @param file A *.csv file formatted with the following columns: x, y, z, a, g, h, theta [optional], pc [optional]
#' @usage
#' FFSread(file)
#' @return
#' Calls in an FFS-class object from a *.csv file
#' @export

FFSread <- function(file){
  animal <- read.csv(file, header=T) #Call in *.csv file; assumes headers are present
  return(new("FFS", rpos=as.matrix(rbind(animal$x,animal$y,animal$z)),a=animal$a,
             g=animal$g[1],h=animal$h[1],
             theta=ifelse(length(animal$theta) > 0, animal$theta, pi/2),
             shape="straight",pc=0.0,L=max(animal$x),ncyl=length(animal$x)))}

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
