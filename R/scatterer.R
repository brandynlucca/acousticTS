################################################################################
################################################################################
# OVERALL SCATTERING OBJECTS
################################################################################
################################################################################
#' Scatterer-class object for target strength estimation
#'
#' @description
#' The \eqn{acousticTS} package uses a variety of defined S4-class objects
#' comprising different types of scatterers, such as fish with gas-filled
#' swimbladders (\link{SBF}) and fluid-like crustaceans
#' (\link{FLS}).
#' @slot metadata List containing relevant metadata
#' @slot model_parameters Model parameters necessary for predicting TS
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
#'    \item{rpos}: the relevant position vector
#'    (\ifelse{html}{\out{r<sub>0</sub>}}{\eqn{r_{0}}})
#'    that includes axes such as \code{x}, \code{y}, \code{z}, etc., that depend
#'    on the type of scatterer being used.
#'    \item{\code{radius}:} in some cases this includes the radius measurements
#'    for each cylinders depending on the type of scatterer object.
#'    \item{\code{theta}:} the orientation of the scatterer relative to the
#'    transmitting transducer or sound source
#'    (\ifelse{html}{\out{&theta;<sub>animal</sub>}}{\eqn{\theta_{animal}}})
#'    that can be represented either by degrees or radians, although all
#'    functions require radians.
#'    \item{\code{g, h}:} material properties that represent the density and
#'    sound speed contrasts (g and h, respectively) relative to the
#'    ambient/surrounding fluid. Some targets may instead have standard
#'    sound speed
#'    (\ifelse{html}{\out{c<sub>animal</sub>}}{\eqn{c_{animal}}},
#'    m \ifelse{html}{\out{s<sup>-1</sup>}}{\eqn{s^{-1}}})
#'    and density
#'    (\ifelse{html}{\out{&rho;<sub>animal</sub>}}{\eqn{\rho_{animal}}},
#'    kg \ifelse{html}{\out{m<sup>3</sup>}}{\eqn{m^3}}).
#'    }
#'    }
#'    \item{\code{shape_parameters}:}{A \code{list} that includes metadata
#'    pertaining to the shape of the scatterer and any other features of
#'    interest (e.g. gas-filled swimbladder). Generally, each entry includes:
#'    overall body length, the number of discrete cylinders that make up the
#'    shape (if applicable), and the units related to both
#'    \ifelse{html}{\out{&theta;<sub>animal</sub>}}{\eqn{\theta_{animal}}}
#'    (e.g. rad, \ifelse{html}{\out{&deg;}}{\eqn{^\circ}}) and length
#'    (e.g. mm, m).}
#'    \item{\code{model_parameters}:}{A \code{list} that contains relevant model
#'    parameterization once an object has been initialized for modeling
#'    \ifelse{html}{\out{&sigma;<sub>bs</sub>}}{\eqn{\sigma_{bs}}}. This is
#'    typically broken up into three categories:
#'    \itemize{
#'    \item{\code{parameters}:} A \code{list} that includes information such as
#'    frequency (Hz), acoustic wavenumber (i.e. k), etc.
#'    \item{\code{medium}:} A \code{data.frame} including information such as
#'    the material properties of the ambient medium.
#'    \item{\code{scatterer}:} A \code{list} containing summarized information
#'    used to parameterize certain scattering models.}}
#'    \item{\code{model}:}{A \code{list} that collects model results from one
#'    or more models in the linear domain (i.e.
#'    \ifelse{html}{\out{&sigma;<sub>bs</sub>}}{\eqn{\sigma_{bs}}}).
#'    }
#'    }
#' @section Supported Scatterers:
#' \itemize{
#'    \item \code{Calibration spheres} (\link{CAL})
#'    \item \code{Elastic-shelled scatterers} (\link{ESS})
#'    \item \code{Fluid-like scatterers} (\link{FLS})
#'    \item \code{Gas-filled scatterers} (\link{GAS})
#'    \item \code{Swimbladdered fish} (\link{SBF})
#' }
#'
#' @keywords scatterer_type
#'
#' @aliases Scatterer
#' @rdname Scatterer-class
#' @export
setClass("Scatterer",
  slots = c(
    metadata = "list",
    model_parameters = "list"
  )
)
################################################################################
################################################################################
# GAS/FLUID-FILLED SCATTERERS
################################################################################
################################################################################
#' Generic gas-filled scatterer (GAS) object/class.
#' @description
#' A S4 class that provides slots to contain relevant metadata for gas-bearing
#' scatterers belonging to the GAS-class. This object can include simple
#' gas-filled bubbles to other scatterers with gas occlusions, swimbladders,
#' and other internal features, if applicable. The default behavior for this
#' type of object is to only reference the gaseous/fluid feature with
#' exceptions that are model-dependent. See \link{Scatterer} for a more
#' detailed description on how this S4 object is organized.
#'
#' @seealso \code{\link{Scatterer}}
#'
#' @keywords scatterer_type
#'
#' @rdname GAS-class
#' @export
GAS <- setClass("GAS",
  slots = c(
    metadata = "list",
    model_parameters = "list",
    model = "list",
    body = "list",
    shape_parameters = "list"
  ),
  contains = "Scatterer"
)
################################################################################
#' Swimbladdered fish (SBF) object/class.
#' @description
#' A S4 class that provides slots to contain relevant animal metadata for
#' parameterizing models for swimbladdered fish (SBF) that are partitioned into
#' two sets of discretized cylinders: the body and the swimbladder. Both shapes
#' comprise independent position matrices, material properties, orientations,
#' and other relevant shape-related data and metadata. See
#' \link{Scatterer} for a more detailed description on how
#' this S4 object is organized.
#'
#' @seealso \link{Scatterer}
#'
#' @keywords scatterer_type
#'
#' @rdname SBF-class
#' @export
SBF <- setClass("SBF",
  slots = c(
    metadata = "list",
    model_parameters = "list",
    model = "list",
    body = "list",
    bladder = "list",
    shape_parameters = "list"
  ),
  contains = "GAS"
)
################################################################################
################################################################################
# ELASTIC-SHELLED SCATTERERS
################################################################################
################################################################################
#' Elastic shelled scatterer (ESS) object/class.
#' @description
#' A S4 class that provides slots to contain relevant metadata for elastic
#' shelled scatterers/objects belonging to the ESS-class. This object can be
#' created using values for both an outer shell and internal tissues, if
#' applicable. The default behavior for this type of this object is to only
#' reference the outer shell with few exceptions that are model-dependent. See
#' \link{Scatterer} for a more detailed description on how this S4 object is
#' organized.
#'
#' @seealso \link{Scatterer}
#'
#' @keywords scatterer_type
#'
#' @rdname ESS-class
#' @export
ESS <- setClass("ESS",
  slots = c(
    metadata = "list",
    model_parameters = "list",
    model = "list",
    shell = "list",
    fluid = "list",
    shape_parameters = "list"
  ),
  contains = "Scatterer"
)
################################################################################
#' Solid and calibration sphere (CAL) object/class.
#' @description
#' A S4 class that provides slots to contain relevant metadata for solid sphere
#' objects belonging to the CAL-class scatterers. This object is created using
#' parameters specific to the outer shell. The default behavior of this object
#' is to only reference these outer elastic shell properties with few
#' exceptions that are model-dependent. See \link{Scatterer} for a more
#' detailed description on how this S4 object is organized.
#'
#' @seealso \link{Scatterer}
#'
#' @keywords scatterer_type
#'
#' @rdname CAL-class
#' @export
CAL <- setClass("CAL",
  slots = c(
    metadata = "list",
    model_parameters = "list",
    model = "list",
    body = "list",
    shape_parameters = "list"
  ),
  contains = "ESS"
)
################################################################################
################################################################################
# FLUID-LIKE SCATTERERS
################################################################################
################################################################################
#' Fluid-like scatterer (FLS) object/class.
#' @description
#' A S4 class that provides slots to contain relevant metadata for scatterers
#' similar to the surrounding fluid medium (i.e fluid-like) belonging to
#' FLS-class scatterers. See \link{Scatterer} for a more
#' detailed description on how this S4 object is organized.
#'
#' @seealso \link{Scatterer}
#'
#' @keywords scatterer_type
#'
#' @rdname FLS-class
#' @export
FLS <- setClass("FLS",
  slots = c(
    metadata = "list",
    model_parameters = "list",
    model = "list",
    body = "list",
    shape_parameters = "list"
  ),
  contains = "Scatterer"
)
