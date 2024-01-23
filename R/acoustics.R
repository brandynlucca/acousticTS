################################################################################
# ACOUSTIC AND SIGNAL PROCESSING VARIABLE CALCULATIONS
################################################################################
################################################################################
# GENERIC ACOUSTIC VARIABLES
################################################################################
################################################################################
#' Calculate the acoustic wavenumber (k) based on the sound speed of water.
#' @param sound_speed Sound speed (c, m/s)
#' @param frequency Frequency (f, Hz)
#' @return
#' Calculates the acoustic wavenumber (k) based on the sound speed of water.
#' @rdname k
#' @export
k <- function( frequency , sound_speed ) 2 * pi * frequency / sound_speed
################################################################################
#' Calculates the linear backscattering coefficient (sigma_bs) from the linear
#' scattering length, f_bs.
#' @param f_bs Linear scattering length (m), or related expression
#' @return
#' Returns the linear backscattering coefficient that can be then converted into 
#' TS.
#' @rdname sigma_bs
#' @export
sigma_bs <- function( f_bs ) abs( f_bs ) * abs( f_bs )
################################################################################
#' Convert backscatter values from log- to linear-domain.
#' @description
#' The `linear(...)` function converts a given value into the linear domain, while
#' the `db(...)` function converts inputs into the log domain.
#' @param value Logarithmic (e.g. TS) or linear (\eqn{\sigma_bs}) value
#' @param coefficient Optional. Numeric coefficient preceding the logarithm. 
#' Default is 10.
#' @return
#' Transforms the backscattering response into either the log (`db`) or 
#' linear (`linear`) domains.
#' @rdname linear
#' @export
linear <- function( value , coefficient = 10 ) coefficient ^ ( value / coefficient )
#' @rdname linear
#' @export
db <- function( value , coefficient = 10 ) coefficient * log( value , base = coefficient ) 
################################################################################
################################################################################
# REFLECTIVITY VARIABLES, COEFFICIENTS, AND EQUATIONS
################################################################################
################################################################################
#' Plane wave/plane interface reflection coefficient
#' @param interface1 Dataframe object containing density (kg/m^3) and sound
#' speed (m/s) values for a boundary/interface (1)
#' @param interface2 Dataframe object containing density (kg/m^3) and sound
#' speed (m/s) values for a boundary/interface (2)
#' @param mode Two options: coefficient calculation for "DWBA" and "KRM"
#' @export
reflection_coefficient <- function( interface1 , interface2 , mode = "DWBA" ) {
  # Calculate acoustic impedance of the first interface ========================
  Z1 <- interface1$density * interface1$sound_speed
  # Calculate acoustic impedance of the second interface =======================
  Z2 <- interface2$density * interface2$sound_speed
  # Calculate the reflection coefficient =======================================
  R <- ( Z2 - Z1 ) / ( Z2 + Z1 )
  # Output =====================================================================
  return( R )
}
################################################################################
#' Transmission coefficient for transmission between two mediums
#' @param interface1 Dataframe object containing density (kg/m^3) and sound
#' speed (m/s) values for a boundary/interface (1)
#' @param interface2 Dataframe object containing density (kg/m^3) and sound
#' speed (m/s) values for a boundary/interface (2)
#' @export
transmission_coefficient <- function( interface1 , interface2 ) {
  # Calculate acoustic impedance of the first interface ========================
  Z1 <- interface1$density * interface1$sound_speed
  # Calculate acoustic impedance of the second interface =======================
  Z2 <- interface2$density * interface2$sound_speed
  # Calculate the reflection coefficient =======================================
  T12 <- ( 2 * ( Z2 / Z1 ) ) / ( 1 + ( Z2 / Z1 ) )
  # Output =====================================================================
  return( T12 )
}
################################################################################
#' Calculate the compressibility material properties of a scatterer's tissue or
#' shell (kappa)
#' @param interface1 Dataframe object containing density (kg/m^3) and sound
#' speed (m/s) values for a boundary/interface (1)
#' @param interface2 Dataframe object containing density (kg/m^3) and sound
#' speed (m/s) values for a boundary/interface (2)
#' @export
kappa <- function( interface1 , interface2 ) {
  # Calculate acoustic compressibility of the first interface ==================
  K1 <- ( interface1$density * interface1$sound_speed ^ 2 ) ^ ( - 1 )
  # Calculate acoustic compressibility of the second interface =================
  K2 <- ( interface2$density * interface2$sound_speed ^ 2 ) ^ ( - 1 )
  # Calculate the reflection coefficient =======================================
  gamma_kappa <- ( K2 - K1 ) / K1
  # Output =====================================================================
  return( gamma_kappa )
}
################################################################################
################################################################################
# ELASTICITY CALCULATIONS AND EQUATIONS
################################################################################
################################################################################
#' Calculate the Poisson's ratio (\eqn{\nu}).
#' @description
#' Calculate the Poisson's ratio (\eqn{\nu}) from two of the three other elastic
#' moduli to calculate the Lamé's parameter. When more than two values are input, 
#' the function will default to using the bulk (K) and Young's (E) moduli. This
#' assumes that the input values represent 3D material properties.
#' @param K Bulk modulus (K, GPa).
#' @param E Young's modulus (E, GPa).
#' @param G Shear modulus (GPa).
#' @return
#' Returns a unitless ratio known as Poisson's ratio (\eqn{\nu}).
#' @rdname pois
#' @export
pois <- function( K , E , G ) {
  if ( missing( K ) & ( ! missing( E ) & ! missing( G ) ) ) E / ( 2 * G ) - 1 
  else if ( missing( E ) & ( ! missing( K ) & ! missing( G ) ) ) ( 3 * K - 2 * G ) / ( 2 * ( 3 * K + G ) )
  else if ( missing( G ) & ( ! missing( K ) & ! missing( E ) ) ) ( 3 * K - E ) / ( 6 * K )
  else if ( ! missing( K ) & ! missing( E ) & ! missing( G ) ) ( 3 * K - E ) / ( 6 * K )
  else stop( "At least two elasticity moduli values are required to calculate Poisson's ratio." )
}
################################################################################
#' Calculate the bulk modulus (K).
#' @description
#' Calculate the bulk modulus (K) from two of the three other elastic
#' moduli to calculate the Lamé's parameter. When more than two values are input, 
#' the function will default to using Young's (E) and shear (G) moduli. This
#' assumes that the input values represent 3D material properties.
#' @param E Young's modulus (GPa).
#' @param G Shear modulus (GPa).
#' @param nu Poisson's ratio (unitless).
#' @return
#' Returns an estimate for the bulk modulus (K).
#' @rdname bulk
#' @export
bulk <- function( E , G , nu ) {
  if ( missing( nu ) & ( ! missing( E ) & ! missing( G ) ) ) E * G / ( 3 * ( 3 * G - E ) )
  else if ( missing( E ) & ( ! missing( nu ) & ! missing( G ) ) ) 2 * G * ( 1 + nu ) / ( 3 * ( 1 - 2 * nu ) )
  else if ( missing( G ) & ( ! missing( nu ) & ! missing( E ) ) ) E / ( 3 * ( 1 - 2 * nu ) )
  else if ( ! missing( K ) & ! missing( E ) & ! missing( G ) ) E * G / ( 3 * ( 3 * G - E ) )
  else stop( "At least two elasticity moduli values are required to calculate the bulk modulus." )
}
################################################################################
#' Calculate Young's modulus (E).
#' @description
#' Calculate the Young's modulus (E) from two of the three other elastic
#' moduli to calculate the Lamé's parameter. When more than two values are input, 
#' the function will default to using the bulk (K) and shear (G) moduli. This
#' assumes that the input values represent 3D material properties.
#' @param K Bulk modulus (GPa).
#' @param G Shear modulus (GPa).
#' @param nu Poisson's ratio (unitless).
#' @return
#' Returns an estimate for the Young's modulus (E).
#' @rdname young
#' @export
young <- function( K , G , nu ) {
  if ( missing( nu ) & ( ! missing( K ) & ! missing( G ) ) ) 9 * K * G / ( 3 * K + G )
  else if ( missing( G ) & ( ! missing( K ) & ! missing( nu ) ) ) 3 * K * ( 1 - 2 * nu )
  else if ( missing( K ) & ( ! missing( G ) & ! missing( nu ) ) ) 2 * G * ( 1 + nu )
  else if ( ! missing( K ) & ! missing( G ) & ! missing( nu ) ) 9 * K * G / ( 3 * K + G )
  else stop( "At least two elasticity moduli values are required to calculate Young's modulus." )
}
################################################################################
#' Calculate the shear modulus (G).
#' @description
#' Calculate the shear modulus (G) from two of the three other elastic
#' moduli to calculate the Lamé's parameter. When more than two values are input, 
#' the function will default to using the bulk (K) and Young's (E) moduli. This
#' assumes that the input values represent 3D material properties.
#' @param K Bulk modulus (GPa).
#' @param E Young's modulus (GPa).
#' @param nu Poisson's ratio (unitless).
#' @return
#' Returns an estimate for the shear modulus (G).
#' @rdname shear
#' @export
shear <- function( K , E , nu ) {
  if ( missing( nu ) & ( ! missing( K ) & ! missing( E ) ) ) 3 * K * E / ( 9 * K - E )
  else if ( missing( E ) & ( ! missing( K ) & ! missing( nu ) ) ) 3 * K * ( 1 - 2 * nu ) / ( 2 * ( 1 + nu ) )
  else if ( missing( K ) & ( ! missing( E ) & ! missing( nu ) ) ) E / ( 2 * ( 1 + nu ) )
  else if ( ! missing( K ) & ! missing( E ) & ! missing( nu ) ) 3 * K * E / ( 9 * K - E )
  else stop( "At least two elasticity moduli values are required to calculate the shear modulus." )
}
################################################################################
#' Calculate Lamé's first parameter (\eqn{\lambda}).
#' @description
#' Calculate Lamé's first parameter (\eqn{\lambda}) from two of the four other 
#' elastic moduli. When more than two values are input, the function will default 
#' to using the bulk (K) and Young's (E) moduli. This assumes that the input 
#' values represent 3D material properties.
#' @param K Bulk modulus (GPa).
#' @param E Young's modulus (GPa).
#' @param G Shear modulus (GPa).
#' @param nu Poisson's ratio (unitless).
#' @return
#' Returns an Lamé's first parameter (\eqn{\lambda}).
#' @rdname lame
#' @export
lame <- function( K , E , G , nu ) {
  if ( ! missing( K ) & ! missing( E ) & missing( G ) & missing( nu ) ) {
    3 * K * ( 3 * K - E ) / ( 9 * K - E )
  } else if ( ! missing( K ) & missing( E ) & ! missing( G ) & missing ( nu ) ) {
    K - 2 * G / 3 
  } else if ( ! missing( K ) & missing( E ) & missing( G ) & ! missing( nu ) ) {
    3 * K * nu / ( 1 + nu )
  } else if ( missing( K ) & ! missing( E ) & ! missing( G ) & missing( nu ) ) {
    G * ( E - 2 * G ) / ( 3 * G - E )
  } else if ( missing( K ) & ! missing( E ) & missing( G ) & ! missing( nu ) ) {
    E * nu / ( ( 1 + nu ) * ( 1 - 2 * nu ) )
  } else if ( missing( K ) & missing( E ) & ! missing( G ) & ! missing( nu ) ) {
    2 * G * nu / ( 1 - 2 * nu )
  } else if ( ! missing( K ) & ! missing( E ) & ! missing( G ) & missing( nu ) ) {
    3 * K * ( 3 * K - E ) / ( 9 * K - E )
  } else if ( missing( K ) & ! missing( E ) & ! missing( G ) & ! missing( nu ) ) {
    G * ( E - 2 * G ) / ( 3 * G - E )
  } else if ( ! missing( K ) & missing( E ) & ! missing( G ) & ! missing ( nu ) ) {
    K - 2 * G / 3 
  } else if ( ! missing( K ) & ! missing( E ) & missing( G ) & ! missing ( nu ) ) {
    3 * K * ( 3 * K - E ) / ( 9 * K - E )
  } else if ( ! missing( K ) & ! missing( E ) & ! missing( G ) & missing ( nu ) ) {
    3 * K * ( 3 * K - E ) / ( 9 * K - E )
  } else {
    stop( "At least two elasticity moduli values are required to calculate the the Lame' parameter." )
  }
}
################################################################################
################################################################################
#' Calculate the mass density material properties of a scatterer's tissue or
#' shell (kappa)
#' @param interface1 Dataframe object containing density (kg/m^3) and sound
#' speed (m/s) values for a boundary/interface (1)
#' @param interface2 Dataframe object containing density (kg/m^3) and sound
#' speed (m/s) values for a boundary/interface (2)
#' @export
rho <- function( interface1 , interface2 ) ( interface2$density - interface1$density ) / interface2$density
################################################################################
#' Wrapper function to model acoustic target strength
#' @param object Scatterer-class object.
#' @param frequency Frequency (Hz).
#' @param model Model name.
#' @param verbose Prints current procedural step occurring from model initialization
#' to calculating TS. Defaults to F.
#' @param ... Additional optional model inputs/parameters.
#' @export
target_strength <- function( object ,
                             frequency ,
                             model , 
                             verbose = F , ... ) {
  # Ignore model name case =====================================================
  model <- tolower( model )
  ts_model <- gsub( "(_.*)" , "\\L\\1" , paste0( toupper( model ) ) , perl = T )
  ts_model <- ifelse( ts_model %in% c( "CALIBRATION" ,
                                      "HIGH_pass_stanton" ) , 
                      tolower( ts_model ) , 
                      ts_model )
  # Pull argument input names ==================================================
  arg_pull <- as.list( match.call( ) )
  if( length( model ) == 1 ) {
    # Initialize objects to input model parameters =============================
    # Grab input arguments =====================================================
    init_name <- paste0( model , "_initialize" )
    arg_list <- names( formals( init_name ) )
    # Filter out inappropriate parameters ======================================
    arg_full <- arg_pull[ arg_list ] 
    true_args <- Filter( Negate( is.null) , arg_full )
    # Initialize ===============================================================
    object_copy <- do.call( init_name , true_args )
    slot( object ,
          "model_parameters" )[ ts_model ] <- extract( object_copy , "model_parameters" )[ ts_model ]
    slot( object ,
          "model" )[ ts_model ] <- extract( object_copy , "model" )[ ts_model ]
    
    # object <- switch( model ,
    #                   mss_anderson = mss_anderson_initialize( object,
    #                                                           frequency , ... ) ,
    #                   calibration = calibration_initialize( object,
    #                                                         frequency , ... ) ,
    #                   dcm = dcm_initialize( object,
    #                                         frequency , ... ) ,
    #                   dwba = dwba_initialize( object ,
    #                                           frequency , ... ) ,
    #                   dwba_curved = dwba_curved_initialize( object ,
    #                                                         frequency ) ,
    #                   sdwba = sdwba_initialize( object ,
    #                                             frequency , ... ) ,
    #                   sdwba_curved = sdwba_curved_initialize( object ,
    #                                                           frequency , ... ) ,
    #                   krm = krm_initialize( object ,
    #                                         frequency , ... ) ,
    #                   stanton_high_pass = stanton_high_pass_initialize( object ,
    #                                                                     frequency , ... ) )
    # cat( toupper( model ) , "model for" , paste0( class(object) , "-object: " ,
    #                                                      extract( object , "metadata" )$ID ) , "initialized.\n" )
    # Determine which model to use =============================================
    # object <- switch( model,
    #                   mss_anderson = MSS_anderson( object ) ,
    #                   calibration = calibration( object ) ,
    #                   dcm = DCM( object ) ,
    #                   dwba = DWBA( object ) ,
    #                   dwba_curved = DWBA_curved( object ) ,
    #                   sdwba = SDWBA( object ) ,
    #                   sdwba_curved = SDWBA_curved( object ) ,
    #                   krm = KRM( object ) ,
    #                   stanton_high_pass = stanton_high_pass( object ) )
    object <- eval( parse( text = paste0( ts_model , "(object)" ) ) )
  } else {
    # Initialize objects to input model parameters =============================
    idx <- 1
    repeat {
      if( idx > length( model ) ) {
        break
      }
      # Pull correct formal arguments ==========================================
      model_name <- paste0( model[ idx ] , "_initialize" )
      arg_list <- names( formals( model_name ) )
      # Filter out inappropriate parameters ====================================
      arg_full <- arg_pull[ arg_list ] 
      true_args <- Filter( Negate( is.null) , arg_full )
      # Initialize =============================================================
      object_copy <- do.call( model_name , true_args )
      slot( object ,
            "model_parameters" )[ ts_model[ idx ] ] <- extract( object_copy , "model_parameters" )[ ts_model[ idx ] ]
      slot( object ,
            "model" )[ ts_model[ idx ] ] <- extract( object_copy , "model" )[ ts_model[ idx ] ]
      if ( verbose ) {
        cat( toupper( model[ idx ] ) , "model for" , paste0( class(object) , "-object: " ,
                                                             extract( object , "metadata" )$ID ) ,
             "initialized.\n\n" )
      }
      idx <- idx + 1
    }
    idx <- 1
    repeat {
      if( idx > length( model ) ) {
        break
      }
      if ( verbose ) {
        cat( "Beginning TS modeling via" , toupper( model[ idx ] ) ,
             "model for" , paste0( class(object) , "-object: " ,
                                   extract( object , "metadata" )$ID ) , "\n" )
      }
      # Parse correct model name ===============================================
     
      # Calculate modeled TS ===================================================
      object <- eval( parse( text = paste0( ts_model[ idx ] , "(object)" ) ) )
      if ( verbose ) {
        object_copy <- do.call( ts_model , list( object ) )
        cat( toupper( model[ idx ] ) , "TS model predictions for" , paste0( class(object) , "-object: " ,
                                                                            extract( object , "metadata" )$ID ) , "complete.\n\n" )
      }
      idx <- idx + 1
    }
  }
  # Output object ==============================================================
  return( object )
}
