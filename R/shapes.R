
################################################################################
################################################################################
# CLASSES FOR CANONICAL AND PRE-DEFINED SHAPE
################################################################################
################################################################################
#' Generic scattering shape object used throughout this package.
#' @description
#' A S4 class that provides slots to contain relevant shape data and metadata for 
#' a variety of arbitrary and canonical shapes and geometries. See 
#' \link{Scatterer} for a more detailed description on how 
#' this S4 object interacts with generic \link{Scatterer} 
#' objects.
#' @slot position_matrix Position matrix that provides the 2D representation of 
#' the body shape
#' @slot shape_parameters A list of additional shape specifications 
#' 
#' @keywords shapes, internal
#' 
#' @rdname Shape-class
#' @aliases Shape
#' 
#' @export
Shape <- setClass( "Shape" ,
                   slots = c(position_matrix = "matrix" ,
                             shape_parameters = "list" ) )
################################################################################
# Arbitrary (or pre-generated) body shape
################################################################################
#' Arbitrary body shape
#' 
#' @seealso \code{\link{Shape}}
#' 
#' @keywords shapes, internal
#'
#' @rdname Arbitrary-class
#' @aliases Arbitrary
#' 
#' @export
Arbitrary <- setClass( "Arbitrary" ,
                       contains = "Shape" )
################################################################################
# Sphere
################################################################################
#' Spherical body shape
#' 
#' @seealso \code{\link{Shape}}
#' 
#' @keywords shapes, internal
#' 
#' @aliases Sphere
#' @rdname Sphere-class
#' @export
Sphere <- setClass( "Sphere" ,
                    contains = "Shape" )
################################################################################
# Prolate spheroid
################################################################################
#' Prolate spheroidal body shape
#' 
#' @seealso \code{\link{Shape}}
#' 
#' @keywords shapes, internal
#' 
#' @aliases ProlateSpheroid
#' @rdname ProlateSpheroid-class
#' @export
ProlateSpheroid <- setClass( "ProlateSpheroid" ,
                             contains = "Shape" )
################################################################################
# Cylinder
################################################################################
#' Cylindrical body shape
#' 
#' @seealso \code{\link{Shape}}
#' 
#' @keywords shapes, internal
#'
#' @aliases Cylinder
#' @rdname Cylinder-class
#' @export
Cylinder <- setClass( "Cylinder" ,
                      contains = "Shape" )
################################################################################
# Polynomial cylinder
################################################################################
#' Cylindrical body shape deformed using a polynomial 
#' 
#' @seealso \code{\link{Shape}}
#' 
#' @keywords shapes, internal
#'
#' @aliases PolynomialCylinder
#' @rdname PolynomialCylinder-class
#' @export
PolynomialCylinder <- setClass( "PolynomialCylinder" ,
                                contains = "Cylinder" )
