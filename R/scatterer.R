#' Scattering-class objects.
#'
#' @description
#'
#' The \eqn{acousticTS} package uses a variety of defined S4-class objects
#' comprising different types of scatterers, such as fish with gas-filled
#' swimbladders (\link[acousticTS]{SBF}) and fluid-like crustaceans
#' (\link[acousticTS]{FLS}).
#'
#' @section Data Organization:
#' \describe{
#'    \item{\code{metadata}:}{A \code{list} comprising any identifying
#'    information associated with the scatterer. The default metadata entry
#'    includes \code{ID} that uses a default \code{character} value of
#'    \code{"UID"} (i.e. \code{ID = "UID"}). This can otherwise be formatted
#'    in any manner for book keeping purposes.}
#'    \item{\code{body/bladder}:}{A \code{list} that includes information
#'    relevant to the scatterer's position vector, material properties,
#'    tilt/orientation, etc. For some scatterers, this may only include
#'    \code{body}, but other targets may have an additional parameter such as
#'    \code{bladder}. Generally, each entry includes:
#'    \itemize{
#'    \item{rpos}: the relevant position vector (\ifelse{html}{\out{r<sub>0</sub>}}{\eqn{r_{0}}})
#'    that includes axes such as \code{x}, \code{y}, \code{z}, etc., that depend
#'    on the type of scatterer being used.
#'    \item{\code{radius}:} in some cases this includes the radius measurements for each
#'    cylinders depending on the type of scatterer object.
#'    \item{\code{theta}:} the orientation of the scatterer relative to the transmitting
#'    transducer or sound source (\ifelse{html}{\out{&theta;<sub>animal</sub>}}{\eqn{\theta_{animal}}})
#'    that can be represented either by degrees or radians, although all functions
#'    require radians.
#'    \item{\code{g, h}:} material properties that represent the density and sound speed
#'    contrasts (g and h, respectively) relative to the ambient/surrounding fluid.
#'    Some targets may instead have standard sound speed
#'    (\ifelse{html}{\out{c<sub>animal</sub>}}{\eqn{c_{animal}}}, m \ifelse{html}{\out{s<sup>-1</sup>}}{\eqn{s^{-1}}})
#'    and density (\ifelse{html}{\out{&rho;<sub>animal</sub>}}{\eqn{\rho_{animal}}}, kg \ifelse{html}{\out{m<sup>3</sup>}}{\eqn{m^3}}).
#'    }
#'    }
#'    \item{\code{shape_parameters}:}{A \code{list} that includes metadata pertaining
#'    to the shape of the scatterer and any other features of interest (e.g. gas-filled
#'    swimbladder). Generally, each entry includes: overall body length, the
#'    number of discrete cylinders that make up the shape (if applicable),
#'    and the units related to both \ifelse{html}{\out{&theta;<sub>animal</sub>}}{\eqn{\theta_{animal}}}
#'    (e.g. rad, \ifelse{html}{\out{&deg;}}{\eqn{\degree}}) and length (e.g. mm, m).}
#'    \item{\code{model_parameters}:}{A \code{list} that contains relevant model
#'    parameterization once an object has been initialized for modeling
#'    \ifelse{html}{\out{&sigma;<sub>bs</sub>}}{\eqn{\sigma_{bs}}}. This is typically
#'    broken up into three categories:
#'    \itemize{
#'    \item{\code{parameters}:} A \code{list} that includes information such as
#'    frequency (Hz), acoustic wavenumber (i.e. k), etc.
#'    \item{\code{medium}:} A \code{data.frame} including information such as
#'    the material properties of the ambient medium.
#'    \item{\code{scatterer}:} A \code{list} containing summarized information
#'    used to parameterize certain scattering models.}}
#'    \item{\code{model}:}{A \code{list} that collects model results from one
#'    or more models in the linear domain (i.e. \ifelse{html}{\out{&sigma;<sub>bs</sub>}}{\eqn{\sigma_{bs}}}).
#'    }
#'    }
#' @section Supported Scatterers:
#' \describe{
#'    \item{\code{Calibration spheres (CAL)}}{\link[acousticTS]{CAL}}
#'    \item{\code{Fluid-like scatterers (FLS)}}{\link[acousticTS]{FLS}}
#'    \item{\code{Swimbladdered fish (SBF)}}{\link[acousticTS]{SBF}}
#'    \item{\code{Elastic shell scatterers (ESS)}}{\link[acousticTS]{ESS}}
#' }
#' @rdname scatterer
#' @export
setClass("scatterer",
         slots=c(metadata="list",
                 model_parameters="list"),
         contains="VIRTUAL")

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
SBF <- setClass("SBF",
                slots=c(metadata="list",
                        model_parameters="list",
                        model="list",
                        body="list",
                        bladder="list",
                        shape_parameters="list"),
                contains="scatterer")

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
#' @rdname FLS
#' @export
FLS <- setClass("FLS",
                slots=c(metadata="list",
                        model_parameters="list",
                        model="list",
                        body="list",
                        shape_parameters="list"),
                contains="scatterer")

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
CAL <- setClass("CAL",
                slots = c(metadata = "list",
                          model_parameters = "list",
                          model = "list",
                          body = "list",
                          shape_parameters = "list"),
                contains = "scatterer")

#' Elastic shelled scatterer (ESS) object/class.
#'
#' @description
#' A S4 class that provides slots to contain relevant metadata for elastic
#' shelled scatterers/objects belonging to the ESS-class scatterers.
#'
#' @export
ESS <- setClass("ESS",
                slots=c(metadata = "list",
                        model_parameters = "list",
                        model = "list",
                        shell = "list",
                        body = "list",
                        shape_parameters = "list"),
                contains = "scatterer")

#' GAS scatterer.
#' @export
GAS <- setClass("GAS",
                slots = c(metadata = "list",
                          model_parameters = "list",
                          model = "list",
                          body = "list",
                          shape_parameters = "list"),
                contains = "scatterer")
