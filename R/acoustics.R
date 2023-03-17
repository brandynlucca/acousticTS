################################################################################
# ACOUSTIC AND SIGNAL PROCESSING VARIABLE CALCULATIONS
################################################################################
################################################################################
# Generic acoustic variables
################################################################################
#' Calculate the acoustic wavenumber (k) based on the sound speed of water.
#' @param sound_speed Sound speed (c, m/s)
#' @param frequency Frequency (f, Hz)
#' @usage
#' kcalc( frequency , sound_speed )
#' @return
#' Calculates the acoustic wavenumber (k) based on the sound speed of water
#' @rdname kcalc
#' @export
kcalc <- function( frequency , sound_speed ) {
  base::return( 2 * pi * frequency / sound_speed )
}
################################################################################
#' Calculates the linear backscattering coefficient (sigma_bs)
#' @param f_bs Linear scattering length (m), or related expression
#' @usage
#' sigma_bs( f_bs )
#' @return
#' Returns the linear backscattering coefficient that can be converted to TS
#' @rdname sigma_bs
#' @export
sigma_bs <- function( f_bs ) {
  base::return( base::abs( f_bs ) ^ 2 )
}
################################################################################
#' Convert backscatter values from log- to linear-domain
#' @param value Logarithmic (e.g. TS) or linear (\eqn{\sigma_bs}) value
#' @param coefficient Optional. Numeric coefficient preceding the logarithm. Default is 10.
#' @return
#' Calculates the acoustic wavenumber based on the sound speed of water
#' @export
linear <- function( value , coefficient = 10 ) {
  base::return( coefficient ^ ( value / coefficient ) )
}
#' @rdname linear
dB <- function( value , coefficient = 10 ){
  base::return( coefficient * base::log10( value ) )
}
################################################################################
################################################################################
# Scattering function variable and parameter calculations
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
  base::return( R )
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
  base::return( T12 )
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
  base::return( gamma_kappa )
}
################################################################################
#' Calculate the mass density material properties of a scatterer's tissue or
#' shell (kappa)
#' @param interface1 Dataframe object containing density (kg/m^3) and sound
#' speed (m/s) values for a boundary/interface (1)
#' @param interface2 Dataframe object containing density (kg/m^3) and sound
#' speed (m/s) values for a boundary/interface (2)
#' @export
rho <- function( interface1 , interface2 ) {
  # Calculate acoustic compressibility of the first interface ==================
  gamma_rho <- ( interface2$density - interface1$density ) / interface2$density
  # Output =====================================================================
  base::return( gamma_rho )
}
################################################################################
#' Wrapper function to model acoustic target strength
#' @param object Scatterer-class object.
#' @param frequency Frequency (Hz).
#' @param model Model name.
#' @param ... Additional optional model inputs/parameters.
#' @export
target_strength <- function( object ,
                             frequency ,
                             model , ... ) {
  # Initialize objects to input model parameters ===============================
  object <- switch( model ,
                    anderson = anderson_initialize( object,
                                                    frequency , ... ) ,
                    calibration = calibration_initialize( object,
                                                          frequency , ... ) ,
                    DCM = dcm_initialize( object,
                                          frequency ) ,
                    DWBA = dwba_initialize( object ,
                                            frequency ) ,
                    KRM = krm_initialize( object ,
                                          frequency ) ,
                    stanton_high_pass = stanton_high_pass_initialize( object ,
                                                                      frequency , ... ) )
  # Determine which model to use ===============================================
  object <- switch( model,
                    anderson = anderson_model( object ) ,
                    calibration = calibration( object ) ,
                    DCM = DCM( object ) ,
                    DWBA = DWBA( object ) ,
                    KRM = KRM( object ) ,
                    stanton_high_pass = stanton_high_pass( object ) )
  # Output object ==============================================================
  return( object )
}