
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
#' @rdname Shape
#' @export
setClass( "Shape" ,
          slots = c(position_matrix = "matrix" ,
                    shape_parameters = "list" ) )
################################################################################
# Arbitrary (or pre-generated) body shape
################################################################################
#' Arbitrary body shape
#' @name Arbitrary-class
#' @rdname Arbitrary-class
#' @export
Arbitrary <- setClass( "Arbitrary" ,
                       contains = "Shape" )
################################################################################
# Sphere
################################################################################
#' Spherical body shape
#' @name Sphere-class
#' @rdname Sphere-class
#' @export
Sphere <- setClass( "Sphere" ,
                    contains = "Shape" )
################################################################################
# Prolate spheroid
################################################################################
#' Prolate spheroidal body shape
#' @name ProlateSpheroid-class
#' @rdname ProlateSpheroid-class
#' @export
ProlateSpheroid <- setClass( "ProlateSpheroid" ,
                             contains = "Shape" )
################################################################################
# Cylinder
################################################################################
#' Cylindrical body shape
#' @name Cylinder-class
#' @rdname Cylinder-class
#' @export
Cylinder <- setClass( "Cylinder" ,
                      contains = "Shape" )
################################################################################
# Polynomial cylinder
################################################################################
#' Cylindrical body shape deformed using a polynomial 
#' @name PolynomialCylinder-class
#' @rdname PolynomialCylinder-class
#' @export
PolynomialCylinder <- setClass( "PolynomialCylinder" ,
                                contains = "Cylinder" )
