################################################################################
################################################################################
# CLASSES FOR CANONICAL AND PRE-DEFINED SHAPE
################################################################################
################################################################################
#' Generic scattering shape object used throughout this package.
#' @description
#' A S4 class that provides slots to contain relevant shape data and metadata
#' for a variety of arbitrary and canonical shapes and geometries. See
#' \link{Scatterer} for a more detailed description on how
#' this S4 object interacts with generic \link{Scatterer}
#' objects.
#' @slot position_matrix Position matrix that provides the 2D representation of
#' the body shape
#' @slot shape_parameters A list of additional shape specifications
#'
#' @keywords shapes internal
#'
#' @name Shape
#' @aliases Shape-class
#' @rdname Shape-class
#' @exportClass Shape
Shape <- setClass("Shape",
  slots = c(
    position_matrix = "matrix",
    shape_parameters = "list"
  )
)
################################################################################
# Arbitrary (or pre-generated) body shape
################################################################################
#' Arbitrary body shape
#'
#' @seealso \code{\link{Shape}}
#'
#' @keywords shapes internal
#'
#' @name Arbitrary
#' @aliases Arbitrary-class
#' @rdname Arbitrary-class
#' @exportClass Arbitrary
Arbitrary <- setClass("Arbitrary",
  contains = "Shape"
)
################################################################################
# Sphere
################################################################################
#' Spherical body shape
#'
#' @seealso \code{\link{Shape}}
#'
#' @keywords shapes internal
#'
#' @name Sphere
#' @aliases Sphere-class
#' @rdname Sphere-class
#' @exportClass Sphere
Sphere <- setClass("Sphere",
  contains = "Shape"
)
################################################################################
# Prolate spheroid
################################################################################
#' Prolate spheroidal body shape
#'
#' @seealso \code{\link{Shape}}
#'
#' @keywords shapes internal
#'
#' @name ProlateSpheroid
#' @aliases ProlateSpheroid-class
#' @rdname ProlateSpheroid-class
#' @exportClass ProlateSpheroid
ProlateSpheroid <- setClass("ProlateSpheroid",
  contains = "Shape"
)
################################################################################
# Oblate spheroid
################################################################################
#' Oblate spheroidal body shape
#'
#' @seealso \code{\link{Shape}}
#'
#' @keywords shapes internal
#'
#' @name OblateSpheroid
#' @aliases OblateSpheroid-class
#' @rdname OblateSpheroid-class
#' @exportClass OblateSpheroid
OblateSpheroid <- setClass("OblateSpheroid",
  contains = "Shape"
)
################################################################################
# Cylinder
################################################################################
#' Cylindrical body shape
#'
#' @seealso \code{\link{Shape}}
#'
#' @keywords shapes internal
#'
#' @name Cylinder
#' @aliases Cylinder-class
#' @rdname Cylinder-class
#' @exportClass Cylinder
Cylinder <- setClass("Cylinder",
  contains = "Shape"
)
################################################################################
# Polynomial cylinder
################################################################################
#' Cylindrical body shape deformed using a polynomial
#'
#' @seealso \code{\link{Shape}}
#'
#' @keywords shapes internal
#'
#' @name PolynomialCylinder
#' @aliases PolynomialCylinder-class
#' @rdname PolynomialCylinder-class
#' @exportClass PolynomialCylinder
PolynomialCylinder <- setClass("PolynomialCylinder",
  contains = "Cylinder"
)
