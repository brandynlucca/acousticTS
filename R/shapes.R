
################################################################################
################################################################################
# CLASSES FOR CANONICAL AND PRE-DEFINED SHAPE
################################################################################
################################################################################
#' Generic scattering shape object used throughout this package.
#' @description
#' A S4 class that provides slots to contain relevant shape data and metadata for 
#' a variety of arbitrary and canonical shapes and geometries. See 
#' \link[acousticTS]{scatterer-class} for a more detailed description on how 
#' this S4 object interacts with generic \link[acousticTS]{scatterer-class} 
#' objects.
#' @slot position_matrix Position matrix that provides the 2D representation of 
#' the body shape
#' @slot shape_parameters A list of additional shape specifications 
#' @rdname shape
#' @export
setClass( "shape" ,
          slots = c(position_matrix = "matrix" ,
                    shape_parameters = "list" ) )
################################################################################
# Arbitrary (or pre-generated) body shape
################################################################################
#' Arbitrary body shape
#' @rdname Arbitrary-class
#' @export
Arbitrary <- setClass( "Arbitrary" ,
                       contains = "shape" )
################################################################################
# Sphere
################################################################################
#' Spherical body shape
#' @rdname Sphere-class
#' @export
Sphere <- setClass( "Sphere" ,
                    contains = "shape" )
################################################################################
# Prolate spheroid
################################################################################
#' Prolate spheroidal body shape
#' @rdname ProlateSpheroid-class
#' @export
ProlateSpheroid <- setClass( "ProlateSpheroid" ,
                             contains = "shape" )
################################################################################
# Cylinder
################################################################################
#' Cylindrical body shape
#' @rdname Cylinder-class
#' @export
Cylinder <- setClass( "Cylinder" ,
                      contains = "shape" )
################################################################################
# Polynomial cylinder
################################################################################
#' Cylindrical body shape deformed using a polynomial 
#' @rdname PolynomialCylinder-class
#' @export
PolynomialCylinder <- setClass( "PolynomialCylinder" ,
                                contains = "Cylinder" )
