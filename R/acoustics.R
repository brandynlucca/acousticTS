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
#' moduli to calculate the \enc{Lamé}{Lame}'s parameter. When more than two 
#' values are input, the function will default to using the bulk (K) and 
#' Young's (E) moduli. This assumes that the input values represent 3D material 
#' properties.
#' @param K Bulk modulus (K, Pa).
#' @param E Young's modulus (E, Pa).
#' @param G Shear modulus (Pa).
#' @return
#' Returns a dimensionless ratio known as Poisson's ratio (\eqn{\nu}).
#' @rdname pois
#' @export
pois <- function( K = NULL , E = NULL , G = NULL ) {
  # Check inputs ===============================================================
  inputs <- list( K = K , E = E , G = G )
  provided <- !sapply( inputs , is.null )
  # Checksum ===================================================================
  if ( sum( provided ) < 2 ) {
    stop(
      paste0( "At least two elasticity moduli values are required to " ,
              "calculate Poisson's ratio." )
    )
  }
  # Use first available combination ============================================
  if ( provided[ "E" ] && provided[ "G" ] ) {
    return( E / ( 2 * G ) - 1  ) 
  } else if ( provided[ "K" ] && provided[ "G" ] ) {
    return( 3 * K - 2 * G ) / ( 2 * ( 3 * K + G ) )
  } else if ( provided[ "K" ] && provided[ "E" ] ) {
    return( ( 3 * K - E ) / ( 6 * K ) )
  }
}
################################################################################
#' Calculate the bulk modulus (K).
#' @description
#' Calculate the bulk modulus (K) from two of the three other elastic
#' moduli to calculate the \enc{Lamé}{Lame}'s parameter. When more than two 
#' values are input, the function will default to using Young's (E) and shear 
#' (G) moduli. This assumes that the input values represent 3D material 
#' properties.
#' @param E Young's modulus (Pa).
#' @param G Shear modulus (Pa).
#' @param nu Poisson's ratio (Dimensionless).
#' @return
#' Returns an estimate for the bulk modulus (K).
#' @rdname bulk
#' @export
bulk <- function( E = NULL , G = NULL , nu = NULL ) {
  # Check inputs ===============================================================
  inputs <- list( E = E , G = G , nu = nu )
  provided <- !sapply( inputs , is.null )
  # Checksum ===================================================================
  if ( sum( provided ) < 2 ) {
    stop(
      paste0( "At least two elasticity moduli values are required to " ,
              "calculate the bulk modulus." )
    )
  }
  # Use first available combination ============================================
  if ( provided[ "E" ] && provided[ "G" ] ) {
    return( E * G / ( 3 * ( 3 * G - E ) ) )
  } else if ( provided[ "G" ] && provided[ "nu" ] ) {
    return( 2 * G * ( 1 + nu ) / ( 3 * ( 1 - 2 * nu ) ) )
  } else if ( provided[ "E" ] && provided[ "nu" ] ) {
    return( E / ( 3 * ( 1 - 2 * nu ) ) )
  }
}
################################################################################
#' Calculate Young's modulus (E).
#' @description
#' Calculate the Young's modulus (E) from two of the three other elastic
#' moduli to calculate the \enc{Lamé}{Lame}'s parameter. When more than two 
#' values are input, the function will default to using the bulk (K) and shear 
#' (G) moduli. This assumes that the input values represent 3D material 
#' properties.
#' @param K Bulk modulus (Pa).
#' @param G Shear modulus (Pa).
#' @param nu Poisson's ratio (Dimensionless).
#' @return
#' Returns an estimate for the Young's modulus (E).
#' @rdname young
#' @export
young <- function( K = NULL , G = NULL , nu = NULL ) {
  # Check inputs ===============================================================
  inputs <- list( K = K , G = G , nu = nu )
  provided <- !sapply( inputs , is.null )
  # Checksum ===================================================================
  if ( sum( provided ) < 2 ) {
    stop(
      paste0( "At least two elasticity moduli values are required to " ,
              "calculate Young's modulus." )
    )
  }
  # Use first available combination ============================================
  if ( provided[ "K" ] && provided[ "G" ] ) {
    return( 9 * K * G / ( 3 * K + G ) )
  } else if ( provided[ "K" ] && provided[ "nu" ] ) {
    return( 3 * K * ( 1 - 2 * nu ) )
  } else if ( provided[ "G" ] && provided[ "nu" ] ) {
    return( 2 * G * ( 1 + nu ) )
  }
}
################################################################################
#' Calculate the shear modulus (G).
#' @description
#' Calculate the shear modulus (G) from two of the three other elastic
#' moduli to calculate the \enc{Lamé}{Lame}'s parameter. When more than two 
#' values are input, the function will default to using the bulk (K) and 
#' Young's (E) moduli. This assumes that the input values represent 3D material 
#' properties.
#' @param K Bulk modulus (Pa).
#' @param E Young's modulus (Pa).
#' @param nu Poisson's ratio (Dimensionless).
#' @return
#' Returns an estimate for the shear modulus (G).
#' @rdname shear
#' @export
shear <- function( K = NULL , E = NULL , nu = NULL ) {
  # Check inputs ===============================================================
  inputs <- list( K = K , E = E , nu = nu )
  provided <- !sapply( inputs , is.null )
  # Checksum ===================================================================
  if ( sum( provided ) < 2 ) {
    stop(
      paste0( "At least two elasticity moduli values are required to " ,
              "calculate the shear modulus." )
    )
  }
  # Use first available combination ============================================
  if ( provided[ "K" ] && provided[ "E" ] ) {
    return( 3 * K * E / ( 9 * K - E ) )
  } else if ( provided[ "K" ] && provided[ "nu" ] ) {
    return( 3 * K * ( 1 - 2 * nu ) / ( 2 * ( 1 + nu ) ) )
  } else if ( provided[ "E" ] && provided[ "nu" ] ) {
    return( E / ( 2 * ( 1 + nu ) ) )
  }
}
################################################################################
#' Calculate \enc{Lamé}{Lame}'s first parameter (\eqn{\lambda}).
#' @description
#' Calculate \enc{Lamé}{Lame}'s first parameter (\eqn{\lambda}) from two of the 
#' four other elastic moduli. When more than two values are input, the function 
#' will default to using the bulk (K) and Young's (E) moduli. This assumes that 
#' the input values represent 3D material properties.
#' @param K Bulk modulus (Pa).
#' @param E Young's modulus (Pa).
#' @param G Shear modulus (Pa).
#' @param nu Poisson's ratio (Dimensionless).
#' @return
#' Returns an \enc{Lamé}{Lame}'s first parameter (\eqn{\lambda}).
#' @rdname lame
#' @export
lame <- function( K , E , G , nu ) {
  # Check inputs ===============================================================
  inputs <- list( K = K , E = E , G = G , nu = nu )
  provided <- !sapply( inputs , is.null )
  # Checksum ===================================================================
  if ( sum( provided ) < 2 ) {
    stop(
      paste0( "At least two elasticity moduli values are required to " ,
              "calculate the Lamé parameter." )
    )
  }
  # Use first available combination ============================================
  if ( provided[ "K" ] && provided[ "G" ] ) {
    return( K - 2 * G / 3 )
  } else if ( provided[ "E" ] && provided[ "nu" ] ) {
    return( E * nu / ( ( 1 + nu ) * ( 1 - 2 * nu ) ) )
  } else if ( provided[ "G" ] && provided[ "nu" ] ) {
    return( 2 * G * nu / ( 1 - 2 * nu ) )
  } else if ( provided[ "K" ] && provided[ "nu" ] ) {
    return( 3 * K * nu / ( 1 + nu ) )
  } else if ( provided[ "K" ] && provided[ "E" ] ) {
    return( 3 * K * ( 3 * K - E ) / ( 9 * K - E ) )
  } else if ( provided[ "E" ] && provided[ "G" ] ) {
    return( G * ( E - 2 * G ) / ( 3 * G - E ) )
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
#' to calculating TS. Defaults to FALSE.
#' @param ... Additional optional model inputs/parameters.
#' @export
target_strength <- function( object , frequency , model , verbose = FALSE , ... ) {
  # Validate inputs ============================================================
  if ( missing( object ) ) stop( "object is required" )
  if ( missing( frequency ) ) stop( "frequency is required" ) 
  if ( missing( model ) ) stop( "model is required" )
  
  # Store the object in a variable that won't conflict with model internals ===
  target_object <- object
  
  # Capture all arguments including ... =======================================
  arg_pull <- list( object = target_object , frequency = frequency , ... )
  
  # Handle model names (convert to uppercase for consistency) ==================
  model <- tolower( model )
  ts_model <- gsub( "(_.*)" , "\\L\\1" , paste0( toupper( model ) ) , perl = T )
  ts_model <- ifelse( ts_model %in% c( "CALIBRATION" ,
                                       "HIGH_pass_stanton" ) , 
                      tolower( ts_model ) , 
                      ts_model )
  
  # Initialize objects to input model parameters ==============================
  idx <- 1
  repeat {
    if( idx > length( model ) ) {
      break
    }
    
    # Pull correct formal arguments ==========================================
    model_name <- paste0( model[ idx ] , "_initialize" )
    
    # Check if initialization function exists ================================
    if ( !exists( model_name ) ) {
      stop( "Initialization function " , model_name , " not found for model " , model[ idx ] )
    }
    
    arg_list <- names( formals( model_name ) )
    
    # Filter out inappropriate parameters ====================================
    arg_full <- arg_pull[ arg_list ] 
    true_args <- Filter( Negate( is.null ) , arg_full )
    
    # Initialize ==============================================================
    object_copy <- do.call( model_name , true_args )
    
    # Store model parameters and results =====================================
    slot( target_object , "model_parameters" )[ ts_model[ idx ] ] <- extract( object_copy , "model_parameters" )[ ts_model[ idx ] ]
    slot( target_object , "model" )[ ts_model[ idx ] ] <- extract( object_copy , "model" )[ ts_model[ idx ] ]
    
    if( verbose ) {
      cat( toupper( model[ idx ] ) , "model for" , paste0( class( target_object ) , "-object: " ,
                                                           extract( target_object , "metadata" )$ID ) ,
           "initialized.\n\n" )
    }
    
    idx <- idx + 1
  }
  
  # Run the models =============================================================
  idx <- 1
  repeat {
    if( idx > length( model ) ) {
      break
    }
    
    if( verbose ) {
      cat( "Beginning TS modeling via" , toupper( model[ idx ] ) ,
           "model for" , paste0( class( target_object ) , "-object: " ,
                                 extract( target_object , "metadata" )$ID ) , "\n" )
    }
    
    # Check if model function exists =========================================
    if ( !exists( ts_model[ idx ] ) ) {
      stop( "Model function " , ts_model[ idx ] , " not found" )
    }
    
    # Calculate modeled TS using do.call instead of eval(parse(...)) =========
    target_object <- do.call( ts_model[ idx ] , list( object = target_object ) )
    
    if( verbose ) {
      cat( toupper( model[ idx ] ) , "TS model predictions for" , paste0( class( target_object ) , "-object: " ,
                                                                          extract( target_object , "metadata" )$ID ) , "complete.\n\n" )
    }
    
    idx <- idx + 1
  }
  
  # Output object ==============================================================
  return( target_object )
}

#' Calculate ka matrix for Goodman and Stern (1962) model
#' @param frequency Frequency vector
#' @param sound_speed_sw Seawater sound speed
#' @param sound_speed_fluid Fluid sound speed
#' @param sound_speed_longitudinal Longitudinal sound speed in shell
#' @param sound_speed_transversal Transversal sound speed in shell
#' @param radius_shell Shell radius
#' @param radius_fluid Fluid radius
#' @return Matrix of ka values
calculate_ka_matrix <- function( frequency , sound_speed_sw , 
                                 sound_speed_fluid , sound_speed_longitudinal , 
                                 sound_speed_transversal ,
                                 radius_shell , radius_fluid ) {
  
  k1 <- k( frequency , sound_speed_sw )
  k3 <- k( frequency , sound_speed_fluid )
  kL <- k( frequency , sound_speed_longitudinal )
  kT <- k( frequency , sound_speed_transversal )
  
  ka_matrix <- rbind(
    k1a_shell = k1 * radius_shell ,
    kLa_shell = kL * radius_shell ,
    kTa_shell = kT * radius_shell ,
    k1a_fluid = k1 * radius_fluid ,
    kTa_fluid = kT * radius_fluid ,
    kLa_fluid = kL * radius_fluid ,
    k3a_fluid = k3 * radius_fluid
  )
  
  return( ka_matrix )
}

#' Calculate alpha coefficients for Goodman-Stern model
#' @param bessel_cache Cached Bessel function values
#' @param ka_matrix_m Modal ka matrix
#' @param m Modal vector
#' @param lambda Lamé first parameter
#' @param mu Shear modulus
#' @param density_sw Seawater density
#' @param density_shell Shell density
#' @param density_fluid Fluid density
#' @return List of alpha coefficients
calculate_goodman_stern_alpha <- function( bessel_cache , ka_matrix_m , m , 
                                           lambda , mu , density_sw , 
                                           density_shell , density_fluid ) {
  # Compute the alpha values needed for the boundary matrices ==================
  a1 <- bessel_cache$k1a_shell$js * density_sw / density_shell
  a2 <- ka_matrix_m$k1a_shell * bessel_cache$k1a_shell$jsd
  a11 <- ka_matrix_m$k1a_shell * bessel_cache$k1a_shell$hs * 
    density_sw / density_shell
  a21 <-  ka_matrix_m$k1a_shell * bessel_cache$k1a_shell$hsd
  a12 <- ( lambda * bessel_cache$kLa_shell$js - 
             2 * mu * bessel_cache$kLa_shell$jsdd ) / ( lambda + 2 * mu )
  a22 <- ka_matrix_m$kLa_shell * bessel_cache$kLa_shell$jsd
  a32 <- 2 * ( ka_matrix_m$kLa_shell * 
                 bessel_cache$kLa_shell$jsd - bessel_cache$kLa_shell$js )
  a42 <- (lambda * bessel_cache$kLa_fluid$js - 
            2 * mu * bessel_cache$kLa_fluid$jsdd) / (lambda + 2*mu)
  a52 <- ka_matrix_m$kLa_fluid * bessel_cache$kLa_fluid$jsd
  a62 <- 2 * ( ka_matrix_m$kLa_fluid * 
                 bessel_cache$kLa_fluid$jsd - bessel_cache$kLa_fluid$js )
  a13 <- -2 * m * ( m + 1 ) * ka_matrix_m$kTa_shell^( -2 ) * 
    ( ka_matrix_m$kTa_shell * 
        bessel_cache$kTa_shell$jsd - bessel_cache$kTa_shell$js )
  a23 <- m * ( m + 1 ) * bessel_cache$kTa_shell$js
  a33 <- ka_matrix_m$kTa_shell^2 * 
    bessel_cache$kTa_shell$jsdd + (m + 2) * (m - 1) * bessel_cache$kTa_shell$js
  a43 <- -2 * m * (m + 1) * ka_matrix_m$kTa_fluid^( -2 ) * 
    ( ka_matrix_m$kTa_fluid * 
        bessel_cache$kTa_fluid$jsd - bessel_cache$kTa_fluid$js )
  a53 <- m * ( m + 1 ) * bessel_cache$kTa_fluid$js
  a63 <- ka_matrix_m$kTa_fluid ^ 2 * 
    bessel_cache$kTa_fluid$jsdd + ( m + 2 ) * ( m - 1 ) * 
    bessel_cache$kTa_fluid$js
  a14 <- ( lambda * bessel_cache$kLa_shell$ys - 
             2 * mu * bessel_cache$kLa_shell$ysdd ) / ( lambda + 2*mu )
  a24 <- ka_matrix_m$kLa_shell * bessel_cache$kLa_shell$ysd
  a34 <- 2 * ( ka_matrix_m$kLa_shell * 
                 bessel_cache$kLa_shell$ysd - bessel_cache$kLa_shell$ys)
  a44 <- ( lambda * bessel_cache$kLa_fluid$ys - 2 * mu * 
             bessel_cache$kLa_fluid$ysdd ) / (lambda + 2*mu)
  a54 <- ka_matrix_m$kLa_fluid * bessel_cache$kLa_fluid$ysd
  a64 <- 2 * ( ka_matrix_m$kLa_fluid * 
                 bessel_cache$kLa_fluid$ysd - bessel_cache$kLa_fluid$ys )
  a15 <- -2 * m * ( m + 1 ) * ka_matrix_m$kTa_shell^( -2 ) * 
    ( ka_matrix_m$kTa_shell * 
        bessel_cache$kTa_shell$ysd - bessel_cache$kTa_shell$ys )
  a25 <- m * ( m + 1 ) * bessel_cache$kTa_shell$ys
  a35 <- ka_matrix_m$kTa_shell^2 * bessel_cache$kTa_shell$ysdd + 
    ( m + 2 ) * ( m - 1 ) * bessel_cache$kTa_shell$ys
  a45 <- -2 * m * ( m + 1 ) * ka_matrix_m$kTa_fluid^( -2 ) * 
    ( ka_matrix_m$kTa_fluid * bessel_cache$kTa_fluid$ysd - 
        bessel_cache$kTa_fluid$ys )
  a55 <- m * ( m + 1 ) * bessel_cache$kTa_fluid$ys
  a65 <- ka_matrix_m$kTa_fluid^2 * 
    bessel_cache$kTa_fluid$ysdd + ( m + 2 ) * ( m - 1 ) * 
    bessel_cache$kTa_fluid$ys
  a46 <- bessel_cache$k3a_fluid$js * density_fluid / density_shell
  a56 <- ka_matrix_m$k3a_fluid * bessel_cache$k3a_fluid$jsd
  # Format and return list =====================================================
  return( list(
    a1 = a1 , a2 = a2 , a11 = a11 , a21 = a21 ,
    a12 = a12 , a22 = a22 , a32 = a32, a42 = a42 , a52 = a52 , a62 = a62 ,
    a13 = a13 , a23 = a23 , a33 = a33, a43 = a43 , a53 = a53 , a63 = a63 ,
    a14 = a14 , a24 = a24 , a34 = a34, a44 = a44 , a54 = a54 , a64 = a64 ,
    a15 = a15 , a25 = a25 , a35 = a35, a45 = a45 , a55 = a55 , a65 = a65 ,
    a46 = a46 , a56 = a56
  ) )
}
#' Calculate boundary condition matrices for Goodman and Stern (1962) model
#' @param alpha Alpha coefficient list
#' @param ka_matrix ka matrix
#' @param m Modal vector
#' @return List of boundary matrices for each frequency and modal order
calculate_goodman_stern_boundary_matrices <- function( alpha , ka_matrix , m ) {
  # Create template matrices ==================================================
  A_template_0 <- matrix( 0 , nrow = 4 , ncol = 4 )
  A_template_m <- matrix( 0 , nrow = 6 , ncol = 6 )
  # Create the boundary condition matrices for each frequency and modal order ==
  boundary_matrices <- lapply( 1 : ncol( ka_matrix ) , function( freq_idx ) {
    lapply( 1 : length( m ) , function( m_idx ) {
      
      # Case: m == 0 ============================================================
      if ( m[ m_idx ] == 0 ) {
        # Apply reduced template ================================================
        A_numerator <- A_template_0
        
        # Numerator =============================================================
        A_numerator[ 1 , 1 : 3 ] <- c(
          alpha$a1[ m_idx , freq_idx ] , alpha$a12[ m_idx , freq_idx ] , 
          alpha$a14[ m_idx , freq_idx ] 
        )
        A_numerator[ 2 , 1 : 3 ] <- c(
          alpha$a2[ m_idx , freq_idx ] , alpha$a22[ m_idx , freq_idx ] , 
          alpha$a24[ m_idx , freq_idx ]        
        )
        A_numerator[ 3 , 2 : 4 ] <- c(
          alpha$a42[ m_idx , freq_idx ] , alpha$a44[ m_idx , freq_idx ] , 
          alpha$a46[ m_idx , freq_idx ]        
        )
        A_numerator[ 4 , 2 : 4 ] <- c(
          alpha$a52[ m_idx , freq_idx ] , alpha$a54[ m_idx , freq_idx ] , 
          alpha$a56[ m_idx , freq_idx ]        
        )
        
        # Denominator ===========================================================
        A_denominator <- A_numerator
        A_denominator[ 1 : 2 , 1 ] <- c( 
          alpha$a11[ m_idx , freq_idx ] , alpha$a21[ m_idx , freq_idx ]
        )
        
        return( 
          list( A_numerator = A_numerator , A_denominator = A_denominator ) 
        )
      } else {
        # Apply reduced template ================================================
        A_numerator <- A_template_m
        
        # Numerator =============================================================
        A_numerator[ 1 , 1 : 5 ] <- c(
          alpha$a1[ m_idx , freq_idx ] , alpha$a12[ m_idx , freq_idx ] ,
          alpha$a13[ m_idx , freq_idx ] , alpha$a14[ m_idx , freq_idx ] ,
          alpha$a15[ m_idx , freq_idx ]
        )
        A_numerator[ 2 , 1 : 5 ] <- c(
          alpha$a2[ m_idx , freq_idx ] , alpha$a22[ m_idx , freq_idx ] ,
          alpha$a23[ m_idx , freq_idx ] , alpha$a24[ m_idx , freq_idx ] ,
          alpha$a25[ m_idx , freq_idx ]
        )
        A_numerator[ 3 , 2 : 5 ] <- c(
          alpha$a32[ m_idx , freq_idx ] , alpha$a33[ m_idx , freq_idx ] ,
          alpha$a34[ m_idx , freq_idx ] , alpha$a35[ m_idx , freq_idx ]
        )
        A_numerator[ 4 , 2 : 6 ] <- c(
          alpha$a42[ m_idx , freq_idx ] , alpha$a43[ m_idx , freq_idx ] ,
          alpha$a44[ m_idx , freq_idx ] , alpha$a45[ m_idx , freq_idx ] ,
          alpha$a46[ m_idx , freq_idx ]
        )
        A_numerator[ 5 , 2 : 6 ] <- c(
          alpha$a52[ m_idx , freq_idx ] , alpha$a53[ m_idx , freq_idx ] ,
          alpha$a54[ m_idx , freq_idx ] , alpha$a55[ m_idx , freq_idx ] ,
          alpha$a56[ m_idx , freq_idx ]
        )
        A_numerator[ 6 , 2 : 5 ] <- c(
          alpha$a62[ m_idx , freq_idx ] , alpha$a63[ m_idx , freq_idx ] ,
          alpha$a64[ m_idx , freq_idx ] , alpha$a65[ m_idx , freq_idx ]
        )
        
        # Denominator ===========================================================
        A_denominator <- A_numerator
        A_denominator[ 1 : 2 , 1 ] <- c( 
          alpha$a11[ m_idx , freq_idx ] , alpha$a21[ m_idx , freq_idx ]
        )
        
        return( 
          list( A_numerator = A_numerator , A_denominator = A_denominator ) 
        )
      }
    } )
  } )
  
  return( boundary_matrices )
}
