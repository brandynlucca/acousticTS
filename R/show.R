################################################################################
# Show functions
################################################################################
################################################################################
# Methods for "show(...)" for each scattering class object
################################################################################
################################################################################
#' Generic function for show(...) for different scatterers.
#' @param object Scattering object.
#' @import graphics
#' @import stats
#' @import grDevices
#' @export
setMethod( f = "show" ,
           signature = "scatterer" ,
           definition = function( object ) {
# Detect scatterer type ========================================================
             sc_type <- class( object )
# Toggle through scatterer types ===============================================
             switch( sc_type ,
                     FLS = fls_show( object ),
                     SBF = sbf_show( object ) ,
                     CAL = cal_show( object ) ,
                     GAS = gas_show( object ) ,
                     ESS = ess_show( object ) )
           } )
################################################################################
#' show(...) for FLS-class objects.
#' @param object FLS-class object.
#' @export
fls_show <- function( object ) {
  # Print out informational text ===============================================
  # Parse metadata =============================================================
  meta <- acousticTS::extract( object ,
                               "metadata" )
  # Parse shape ================================================================
  shape <- acousticTS::extract( object ,
                                "shape_parameters" )
  # Parse body =================================================================
  body <- acousticTS::extract( object ,
                               "body" )
  # Print object summary information ===========================================
  base::cat(
    base::paste0( methods::is( object )[[ 1 ]], "-object" ), "\n" ,
    " Fluid-like scatterer \n " ,
    " ID:" ,
    base::paste0( meta$ID ) , "\n" ,
    "Body dimensions:\n" ,
    " Length:" , base::paste0( base::round( shape$length , 3 ) ,
                               " " ,
                               shape$length_units ) ,
    base::paste0( "(n = " , shape$n_segments , " cylinders)" ) , "\n" ,
    " Mean radius:", base::paste0( base::round( base::mean( body$radius ) , 4 ) ,
                                   " " ,
                                   shape$length_units ) , "\n" ,
    " Max radius:" , base::paste0( base::round( base::max( body$radius ) , 4 ) ,
                                   " " ,
                                   shape$length_units ) , "\n" ,
    "Shape parameters:\n" ,
    " Defined shape:" ,
    paste0( shape$shape ) , "\n" ,
    " L/a ratio:" ,
    base::paste0( base::round( shape$length / shape$radius , 1 ) ) , "\n" ,
    " Taper order:" ,
    base::paste0(  shape$taper_order ) , "\n" ,
    "Material properties:\n" ,
    base::paste0( " g: " , base::round( base::mean( body$g ), 4 ) ) , "\n" ,
    base::paste0( " h: " , base::round( base::mean( body$h ) , 4 ) ) , "\n" ,
    "Body orientation (relative to transducer face/axis):" ,
    base::paste0( base::round( body$theta , 3 ) ,
                  " " ,
                  shape$theta_units ) )
}
#' show(...) for GAS_class objects
#' @param object GAS-class object
#' @export
gas_show <- function ( object ) {
  # Print out informational text ===============================================
  # Parse metadata =============================================================
  meta <- acousticTS::extract( object ,
                               "metadata" )
  # Parse shape =============================== =================================
  shape <- acousticTS::extract( object ,
                                "shape_parameters" )
  # Parse body =================================================================
  body <- acousticTS::extract( object ,
                               "body" )
  cat(
    paste0( methods::is( object )[[ 1 ]], "-object" ) , "\n" ,
    " Gas- and fluid-filled scatterer \n " ,
    " ID:" , paste0( meta$ID ) , "\n" ,
    "Body dimensions:\n" ,
    " Diameter:" ,
    paste0( shape$radius * 2 ,
            " " ,
            shape$radius_units ) , "\n" ,
    " Radius:" ,
    paste0( shape$radius ,
            " " ,
            shape$radius_units ) , "\n" ,
    "Material properties:\n" ,
    paste0( " g: " , round( mean( body$g ) , 4 ) ) , "\n" ,
    paste0( " h: " , round( mean( body$h ) , 4 ) ) , "\n" 
  )
}
#' show(...) for SBF-class objects.
#' @param object SBF_class object.
#' @export
sbf_show <- function(object) {
  # Print out informational text ===============================================
  # Parse metadata =============================================================
  meta <- acousticTS::extract( object ,
                               "metadata")
  # Parse shape ================================================================
  shape <- acousticTS::extract( object ,
                                "shape_parameters" )
  # Parse body =================================================================
  body <- acousticTS::extract( object ,
                               "body" )
  # Parse bladder ==============================================================
  bladder <- acousticTS::extract( object ,
                                  "bladder" )
  # Print object summary information ===========================================
  base::cat(
    base::paste0( methods::is( object )[[ 1 ]], "-object" ), "\n" ,
    " Fluid-like scatterer \n " ,
    " ID:" ,
    base::paste0( meta$ID ) , "\n" ,
    "Body dimensions:\n" ,
    " Length:" , base::paste0( base::round( shape$body$length , 3 ) ,
                               " " ,
                               shape$length_units ) ,
    base::paste0( "(n = " , shape$body$n_segments , " cylinders)" ) , "\n" ,
    " Mean radius:", base::paste0( base::round( base::mean( body$rpos[ 2 , ] / 2 ) , 4 ) ,
                                   " " ,
                                   shape$length_units ) , "|" ,
    "Max radius:" , base::paste0( base::round( base::max( body$rpos[ 2 , ] / 2 ) , 4 ) ,
                                  " " ,
                                  shape$length_units ) , "\n" ,
    "Bladder dimensions:\n" ,
    " Length:" , base::paste0( base::round( shape$bladder$length , 3 ) ,
                               " " ,
                               shape$length_units ) ,
    base::paste0( "(n = " , shape$bladder$n_segments , " cylinders)" ) , "\n" ,
    " Mean radius:", base::paste0( base::round( base::mean( bladder$rpos[ 2 , ] / 2 ) , 4 ) ,
                                   " " ,
                                   shape$length_units ) , "|" ,
    "Max radius:" , base::paste0( base::round( base::max( bladder$rpos[ 2 , ] / 2 ) , 4 ) ,
                                  " " ,
                                  shape$length_units ) , "\n" ,
    "Body material properties:\n" ,
    base::paste0( " Density: " , base::round( base::mean( body$density ), 4 ) ) ,
    "kg m^-3" , "|" ,
    base::paste0( "Sound speed: " , base::round( base::mean( body$sound_speed ) , 4 ) ) ,
    "m s^-1" , "\n" ,
    "Bladder fluid material properties:\n" ,
    base::paste0( " Density: " , base::round( base::mean( bladder$density ), 4 ) ) ,
    "kg m^-3" , "|" ,
    base::paste0( "Sound speed: " , base::round( base::mean( bladder$sound_speed ) , 4 ) ) ,
    "m s^-1" , "\n" ,
    "Body orientation (relative to transducer face/axis):" ,
    base::paste0( base::round( body$theta , 3 ) ,
                  " " ,
                  shape$theta_units ) )
}
#' show(...) for CAL-class objects.
#' @param object CAL-class object.
#' @export
cal_show <- function( object ) {
  # Parse metadata =============================================================
  meta <- acousticTS::extract( object ,
                               "metadata")
  # Parse shape ================================================================
  shape <- acousticTS::extract( object ,
                                "shape_parameters" )
  # Parse body =================================================================
  body <- acousticTS::extract( object ,
                               "body" )
  # Print object summary information ===========================================
  cat( paste0( is( object )[[ 1 ]] , "-object" ) , "\n" ,
       "Calibration sphere" , "\n" ,
       " ID:" ,
       paste0( meta$ID ) , "\n" ,
       "Material:" ,
       paste0( meta$Material ) , "\n" ,
       " Sphere longitudinal sound speed:" ,
       paste0( body$sound_speed_longitudinal ) , "m/s" ,  "\n" ,
       " Sphere transversal sound speed:" ,
       paste0( body$sound_speed_transversal ) , "m/s" , "\n" ,
       " Sphere density:" ,
       paste0( body$density ) , "kg/m^3" , "\n" ,
       "Diameter:" ,
       paste0( shape$diameter ,
               " " ,
               shape$diameter_units ) , "\n" ,
       " Radius:" ,
       paste0( shape$radius ,
               " " ,
               shape$diameter_units) , "\n" ,
       "Propagation direction of the incident sound wave:" ,
       paste0( round( body$theta , 3 ) ,
               " " ,
               shape$theta_units ) )
}
#' show(...) for ESS-class objects.
#' @param object ESS-class object.
#' @export
ess_show <- function( object ) {
  # Parse metadata =============================================================
  meta <- acousticTS::extract( object ,
                               "metadata")
  # Parse shape ================================================================
  shape <- acousticTS::extract( object ,
                                "shape_parameters" )
  # Parse shell ================================================================
  shell <- acousticTS::extract( object ,
                                "shell" )
  # Parse fluid ================================================================
  fluid <- acousticTS::extract( object ,
                                "fluid" )
  # Create material properties string ==========================================
  # Shell ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  shell_material_props <- shell[ 
    names( shell ) %in% c( "sound_speed" , "density" , "g" , "h" , "K" , 
                           "E" , "G" , "nu" ) 
  ]
  shell_material_text <- if( length( shell_material_props ) > 0 ) {
    prop_strings <- mapply( function( name , value ) {
      clean_name <- gsub( "_" , " " , name )  # Remove underscores
      clean_name <- switch( clean_name ,
                            "density" = "Density" ,
                            "sound speed" = "Sound speed" ,
                            "K" = "Bulk modulus (K)" ,
                            "E" = "Young's modulus (E)" ,
                            "G" = "Shear modulus (G)" ,
                            "nu" = "Poisson's ratio" ,
                            clean_name )
      units <- switch( name ,
                       "density" = " kg m^-3" ,
                       "sound_speed" = " m s^-1" ,
                       "K" = " Pa" ,
                       "E" = " Pa" ,
                       "G" = " Pa" ,
                       "" )
      paste0( clean_name , ": " , round( value , 4 ) , units )
    } , names( shell_material_props ) , shell_material_props , SIMPLIFY = F )
    paste( "   " , prop_strings , collapse = "\n " )
  } else {
    "   None specified"
  }
  # Internal fluids ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  fluid_material_props <- fluid[
    names( fluid ) %in% c( "sound_speed" , "density" , "g" , "h" )     
  ]
  fluid_material_text <- if( length( shell_material_props ) > 0 ) {
    prop_strings <- mapply( function( name , value ) {
      clean_name <- gsub( "_" , " " , name )  # Remove underscores
      clean_name <- switch( clean_name ,
                            "density" = "Density" ,
                            "sound speed" = "Sound speed" ,
                            clean_name )
      units <- switch( name ,
                       "density" = " kg m^-3" ,
                       "sound_speed" = " m s^-1" ,
                       "" )
      paste0( clean_name , ": " , round( value , 4 ) , units )
    } , names( fluid_material_props ) , fluid_material_props , SIMPLIFY = F )
    paste( "   " , prop_strings , collapse = "\n " )
  } else {
    "   None specified"
  }
  # Print object summary information ===========================================
  cat( paste0( is( object )[[ 1 ]] , "-object" ) , "\n" ,
       "Elastic-shelled scatterer" , "\n" ,
       " ID:" ,
       paste0( meta$ID ) , "\n" ,
       "Material:" ,
       paste0( meta$Material ) , "\n" ,
       "  Shell: \n" ,
       shell_material_text , " \n" ,
       "  Internal fluid-like body: \n" ,
       fluid_material_text , " \n" ,
       "Shape: \n" ,
       "  Shell: \n" ,
       "    Radius:" , paste0( shape$shell$radius , 
                               " " , 
                               shape$shell$length_units ) , " \n" ,
       "    Diameter:" , paste0( shape$shell$diameter , 
                                 " " , 
                                 shape$shell$length_units ) , " \n" ,
       "    Outer thickness:" , paste0( shell$shell_thickness , 
                                        " " , 
                                        shape$shell$length_units ) , "\n" ,
       "  Internal fluid-like body: \n" ,
       "    Radius:" , paste0( shape$fluid$radius , 
                               " " , 
                               shape$fluid$length_units ) , " \n" ,
       "    Diameter:" , paste0( shape$fluid$diameter , 
                                 " " , 
                                 shape$shell$length_units ) , " \n" ,
       "Propagation direction of the incident sound wave:" ,
       paste0( round(shell$theta , 3 ) , " radians" )
  )
}
