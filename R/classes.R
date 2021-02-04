#' Scattering-class objects.
#' 
#' @description
#' The group of S4 class objects that represents different scatterer-types, such
#' as fluid-like (FLS) and gas-bearing (GAS) animals. In general, each of these 
#' classes contain relevant animal metadata for parameterizing models, such 
#' as orientation (\eqn{\theta}), material properties (\eqn{\rho} and g, 
#' c and h), and body morphology (e.g., \eqn{r_0}).
#' 
#' @details
#' There are currently several scattering-classes available for modeling using
#' the \eqn{acousticTS} package:
#' 
#' \itemize{
#'  \item \strong{FLS}: fluid-like 
#'  scatterers (\code{\link[acousticTS]{FLS}}) defined 
#'  by the FLS-class. These generally describe scatterers (or body structures) 
#'  that are fluid-like, meaning that the acoustic material properties (i.e., 
#'  sound speed and density) are similar to that of the surrounding medium. 
#'  Consequently, things like thin flesh and crustaceans are typically modeled 
#'  as FLS objects. Current models available to FLS objects includes the 
#'  \link[acousticTS]{SDWBA}, \link[acousticTS]{KRM}, and 
#'  \link[acousticTSs]{DCM}.
#'  \item \strong{SBF}: `swimbladdered fish (SBF)` defined by the SBF-class. 
#'  \item \strong{CAL}: `solid and calibration spheres (CAL)` defined by the 
#'  CAL-class. 
#'  \item \strong{GAS}: `gas-bearing scatterers (GAS)` defined by the GAS-class.
#'  \item \strong{ESS}: `elastic shelled scatterers (ESS)` defined by the 
#'  ESS-class.
#' }
#' 
#' @export
scattering <- setClass("scattering", representation("VIRTUAL"))

#' Fluid-like scatterer (FLS) object/class.
#'
#' @description
#' A S4 class that provides slots to contain relevant animal metadata for 
#' parameterizing models for fluid-like scatterers (FLS) partitioned into 
#' discretized cylinders. This, specifically, includes a position matrix, 
#' radius, material properties (g, h), orientation, animal shape, and body 
#' curvature. This class is used within the DWBA and DFCM model functions. In 
#' the future, this will also allow for converting one class of scatterer into 
#' another for seemless usage for model comparisons.
#'
#' @param rpos A 3xN matrix, where N is the number of columns, or segments, in 
#' our scatterer object. The rows represent the x-, y-, and z-axis
#' shape values in that order. Units in meters.
#' @param a An N-length radius vector. Units in meters.
#' @param g Density contrast. Units are dimensionless.
#' @param h Sound speed contrast. Units are dimensionless.
#' @param theta Animal orientation, where \eqn{\pi/2} is broadside incidence. 
#' Units in radians.
#' @param curve Body curvature/flexure. Boolean (T/F) which flags the a
#' ppropriate model formulation.
#' @param pc Radius of curvature ratio. Units are dimensionless.
#' @param L Length of target. Units in meters.
#' @param ncyl Number of segments/cylinders comprising the shape. Units should 
#' be an integer value.
#'
#' @details
#' FLS objects can be created using 
#' \link[acousticTS]{FLSwrite} and 
#' \link[acousticTS]{FLSgenerate}. These objects can also be 
#' read using
#' built-in datasets (e.g., \link[acousticTS]{mcgehee}) or using 
#' \link[acousticTS]{FLSread}. The 
#' \link[acousticTS]{SDWBA}, 
#' \link[acousticTS]{SDWBA.sim}, \link[acousticTS]{DCM}, and 
#' \link[acousticTS]{Shapely} functions all use FLS-objects for their 
#' respective 
#' inputs.
#'
#' @export
#' @import methods

FLS <- setClass("FLS",slots=c(rpos="matrix",a="numeric",
                              g="numeric",h="numeric",
                              theta="numeric",curve="logical",pc="numeric",
                              L="numeric",ncyl="numeric"),
                contains="scattering")

#' Swimbladdered fish (SBF) object/class.
#'
#' @description
#' A S4 class that provides slots to contain relevant animal metadata for 
#' parameterizing models for swimbladdered fish/scatterers (SBF)
#' partitioend into two sets of discretized cylinders: the body and 
#' swimbladder. Both shapes comprise a position matrix, material properties
#' (sound speed, c, and density, rho), orientation (theta), the overall body 
#' length of the animal, and the number of cylinders. This class is
#' used within the KRM and GBZ model functions.
#' @export
SBF <- setClass("SBF", slots=c(body="matrix",bladder="matrix",theta="numeric",
                               pb="numeric",
                               cb="numeric",psb="numeric",csb="numeric",
                               L="numeric",
                               ncylb="numeric",ncylsb="numeric"),
                contains="scattering")

#' Solid and calibration sphere (CAL) object/class.
#'
#' @description
#' A S4 class that provides slots to contain relevant metadata for solid sphere
#' objects belonging to the CAL-class scatterers.
#' @param material Material-type for the soldi sphere. See 'Details' built-in
#' material options. 
#' @param a Spherical radius (m).
#' @param c1 Longitudinal sound speed (m/s).
#' @param c2 Transversal sound speed (m/s).
#' @param rho1 Density (kg/m^3)
#'
#' @details
#' There are several options for the \strong{material} argument:
#' \tabular{rlllll}{
#'  \strong{Material} \tab \strong{Argumment} \tab \strong{c1} \tab \strong{c2} 
#'  \tab \strong{\eqn{\rho1}}\cr
#'  \emph{Tungsten carbide} \tab "WC" \tab 6853 \tab 4171 \tab 8360\cr
#'  \emph{Stainless steel} \tab "steel" \tab 5980 \tab 3297 \tab 7970\cr
#'  \emph{Brass} \tab "brass" \tab 4372 \tab 2100 \tab 8360\cr
#'  \emph{Copper} \tab "Cu" \tab 4760 \tab 2288.5 \tab 8947\cr
#'  \emph{Aluminum} \tab "Al" \tab 6260 \tab 3080 \tab 2700\cr
#' }
#' 
#' @export
CAL <- setClass("CAL", slots=c(material="character", 
                               a="numeric",
                               c1="numeric", c2="numeric", rho1="numeric"),
                contains="scattering")

#' Elastic shelled scatterer (ESS) object/class.
#'
#' @description
#' A S4 class that provides slots to contain relevant metadata for elastic 
#' shelled scatterers/objects belonging to the ESS-class scatterers.
#' 
#' @param a Radius (m).
#' @param g Density contrast. 
#' @param h Sound speed contrast.
#' @param theta Orientation (\eqn{\theta}, radians)
#' 
#' @param L Optional. Length of the longitudinal axis (m) for non-spherical 
#' shapes.
#' @param spherical Optional. Boolean for specifying whether scatterer is 
#' spherical or oblong. 
#' 
#' @export
ESS <- setClass("ESS", slots=c(L="numeric",
                               a="numeric", 
                               g="numeric", h="numeric", theta="numeric",
                               spherical="logical"), 
                contains="scattering")
