################################################################################
################################################################################
# Model initialization functions
################################################################################
################################################################################
#' Initialize GAS-object for modal series solution.
#' @param object GAS-class object.
#' @param radius Radius of sphere (m).
#' @param g_body Density contrast for gas.
#' @param h_body Sound speed contrast for gas.
#' @param ka_limit Modal series limit (i.e. max "m"). The default is the maximum
#'    ka + 10.
#' @inheritParams krm_initialize
#' @export
mss_anderson_initialize <- function(object ,
                                    # Required input
                                    frequency ,
                                    # Optional body parameters
                                    radius = NULL ,
                                    g_body = NULL ,
                                    h_body = NULL ,
                                    # Optional medium inputs
                                    sound_speed_sw = 1500 ,
                                    density_sw = 1026 ,
                                    # Optional modal series limit
                                    ka_limit = NULL ) {
  # Detect object class ========================================================
  scatterer_type <- class( object )
  # Parse metadata =============================================================
  meta <- extract( object , "metadata")
  # Parse shape ================================================================
  shape <- extract( object , "shape_parameters" )
  # Parse body =================================================================
  body <- extract( object , "body" )
  # Define medium parameters ===================================================
  medium_params <- data.frame( sound_speed = sound_speed_sw ,
                               density = density_sw )
  # Define body radii compatible w/ model ======================================
  # bladder <- base::switch( scatterer_type ,
  #                          SBF = acousticTS::extract( object ,
  #                                                     "bladder" ) ,
  #                          GAS = body )
  # # Equivalent spherical radius ================================================
  # if ( scatterer_type == "SBF" ) {
  #   # Equivalent volume ========================================================
  #   acousticTS::along_sum( bladder$rpos , shape$bladder$n_segments )
  #   bladder$rpos[ 2 , ] / 2
  # }
  # aesr <- base::switch( scatterer_type ,
  #                       )
  # Define model parameters recipe =============================================
  model_params <- base::list(
    acoustics = base::data.frame(
      frequency = frequency ,
      # Wavenumber (medium) ====================================================
      k_sw = acousticTS::k( frequency ,
                                sound_speed_sw ) ,
      # Wavenumber (fluid) =====================================================
      k_f = acousticTS::k( frequency ,
                               body$h * sound_speed_sw ) ) ,
    ncyl_b = shape$n_segments )
  # Define ka limit ============================================================
  ka_limit <- base::ifelse( ! is.null( ka_limit ) ,
                            ka_limit ,
                            base::round(
                              base::max( model_params$acoustics$k_sw ) * body$radius ) + 20 )
  model_params$ka_limit <- ka_limit
  # Define body parameters recipe ==============================================
  body_params <- base::data.frame(
    radius = body$radius ,
    h = body$h ,
    g = body$g ,
    theta = body$theta )
  # Add model parameters slot to scattering object =============================
  methods::slot( object ,
                 "model_parameters" )$MSS_anderson <- base::list(
                   parameters = model_params ,
                   medium = medium_params ,
                   body = body_params )
  # Add model results slot to scattering object ================================
  methods::slot( object ,
                 "model" )$MSS_anderson <- base::data.frame(
                   frequency = frequency ,
                   sigma_bs = base::rep( NA ,
                                         length( frequency ) ) )
  # Output =====================================================================
  return( object )
}
################################################################################
################################################################################
# Model initialization functions
################################################################################
################################################################################
#' Initialize ESS-object for modal series solution.
#' @param object ESS-class object.
#' @param frequency Frequency vector (Hz).
#' @param sound_speed_sw Seawater sound speed (m/s). Default: 1500.
#' @param density_sw Seawater density (kg/m³). Default: 1026.
#' @param m_limit Modal series limit (i.e. max "m"). The default is the maximum
#'    ka + 10.
#' @export
mss_goodman_stern_initialize <- function( object ,
                                          frequency ,
                                          sound_speed_sw = 1500 ,
                                          density_sw = 1026 ,
                                          m_limit = NULL ) {
  # Detect object class ========================================================
  scatterer_type <- class( object )
  # Parse metadata =============================================================
  meta <- extract( object , "metadata")
  # Parse shape ================================================================
  shape <- extract( object , "shape_parameters" )
  # Parse shell ================================================================
  shell <- extract( object , "shell" )
  # Parse fluid ================================================================
  fluid <- extract( object , "fluid" )
  # Calculate wavenumbers ======================================================
  k1 <- k( frequency , sound_speed_sw )
  # Define the maximum Modal series iterator limit =============================
  m_limit <- ifelse( is.null( m_limit ) , 
                     max( round( k1 * shell$radius ) ) + 10 , 
                     m_limit )
  # Define medium parameters ===================================================
  medium_params <- data.frame( sound_speed = sound_speed_sw ,
                               density = density_sw )
  # Define model parameters recipe =============================================
  model_params <- base::list(
    acoustics = base::data.frame(
      frequency = frequency ,
      # Wavenumber (medium) ====================================================
      k_sw = k1 
    ) ,
    m_limit = m_limit ,
    m = 0 : m_limit
  )
  # Initialize lists for fluid and shell =======================================
  shell_params <- list()
  fluid_params <- list()
  # Establish the correct material properties ==================================
  # Shell ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  shell_params$sound_speed <- ifelse( "sound_speed" %in% names( shell ) ,
                                      shell$sound_speed ,
                                      ifelse( "h" %in% names( shell ) ,
                                              shell$h * sound_speed_sw ,
                                              NA ) )
  shell_params$density <- ifelse( "density" %in% names( shell ) ,
                                  shell$density ,
                                  ifelse( "g" %in% names( shell ) ,
                                          shell$g * density_sw ,
                                          NA ) )
  shell_params$nu <- ifelse( "nu" %in% names( shell ) , shell$nu , NA )
  shell_params$K <- ifelse( "K" %in% names( shell ) , shell$K , NA )
  shell_params$E <- ifelse( "E" %in% names( shell ) , shell$E , NA )
  shell_params$G <- ifelse( "G" %in% names( shell ) , shell$G , NA )
  shell_params$lambda <- ifelse( "lambda" %in% names( shell ) , 
                                 shell$lambda , NA )
  # Fluid ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  fluid_params$sound_speed <- ifelse( "sound_speed" %in% names( fluid ) ,
                                      fluid$sound_speed ,
                                      ifelse( "h" %in% names( fluid ) ,
                                              fluid$h * sound_speed_sw ,
                                              NA ) )
  fluid_params$density <- ifelse( "density" %in% names( fluid ) ,
                                  fluid$density ,
                                  ifelse( "g" %in% names( fluid ) ,
                                          fluid$g * density_sw ,
                                          NA ) )
  # Add morphometrics ==========================================================
  # Shell ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  shell_params$radius <- shell$radius 
  shell_params$shell_thickness <- shell$shell_thickness
  # Fluid ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  fluid_params$radius <- fluid$radius
  # Add model parameters slot to scattering object =============================
  methods::slot( object ,
                 "model_parameters" )$MSS_goodman_stern <- base::list(
                   parameters = model_params ,
                   shell = shell_params ,
                   fluid = fluid_params ,
                   medium = medium_params )
  # Add model results slot to scattering object ================================
  methods::slot( object ,
                 "model" )$MSS_goodman_stern <- base::data.frame(
                   frequency = frequency ,
                   sigma_bs = base::rep( NA ,
                                         length( frequency ) ) )
  # Output =====================================================================
  return( object )
}
################################################################################
#' Initialize CAL-class object for modeling.
#' @param object CAL-class object.
#' @inheritParams krm_initialize
#' @export
calibration_initialize <- function( object ,
                                    frequency ,
                                    sound_speed_sw = 1500 ,
                                    density_sw = 1026 ) {
  # Parse metadata =============================================================
  meta <- acousticTS::extract( object ,
                               "metadata")
  # Parse shape ================================================================
  shape <- acousticTS::extract( object ,
                                "shape_parameters" )
  # Parse body =================================================================
  body <- acousticTS::extract( object ,
                               "body" )
  # Define medium parameters ===================================================
  medium_params <- base::data.frame( sound_speed = sound_speed_sw ,
                                     density = density_sw )
  # Define model parameters recipe =============================================
  model_params <- base::list(
    acoustics = base::data.frame(
      frequency = frequency ,
      # Wavenumber (medium) ====================================================
      k_sw = acousticTS::k( frequency ,
                                sound_speed_sw ) ,
      # Wavenumber (longitudinal axis) =========================================
      k_l = acousticTS::k( frequency ,
                               body$sound_speed_longitudinal ) ,
      # Wavenumber (tranversal axis) ===========================================
      k_t = acousticTS::k( frequency ,
                               body$sound_speed_transversal ) ) ,
    parameters = base::data.frame(
      ncyl_b = shape$n_segments )
  )
  # Define body parameters recipe ==============================================
  body_params <- base::data.frame(
    diameter = body$diameter ,
    radius = body$radius ,
    sound_speed_longitudinal = body$sound_speed_longitudinal ,
    sound_speed_transversal = body$sound_speed_transversal ,
    density = body$density ,
    theta = body$theta )
  # Add model parameters slot to scattering object =============================
  methods::slot( object ,
                 "model_parameters" )$calibration <- base::list(
                   parameters = model_params ,
                   medium = medium_params ,
                   body = body_params )
  # Add model results slot to scattering object ================================
  methods::slot( object ,
                 "model" )$calibration <- base::data.frame(
                   frequency = frequency ,
                   sigma_bs = base::rep( NA ,
                                         length( frequency ) ) )
  # Output =====================================================================
  return( object )
}
################################################################################
#' Initialize object for modeling using the DCM.
#'
#' @param object FLS-class object.
#' @param frequency Transmit frequency (kHz)
#' @param radius_cylinder Optional input to override current shape radius.
#' @param radius_curvature Numeric input for the radius of curvature
#' @param radius_curvature_ratio Ratio between body length and the radius of
#'    curvature. Defaults to 3.0.
#' @param radius_cylinder_fun Defines which radius value will be used from the
#'    radius vector. Defaults to "max", but also accepts "mean" and "median".
#' @param length Body length (m).
#' @param g Density contrast.
#' @param h Sound speed contrast.
#' @param theta Body orientation relative to the incident sound wave.
#' @param sound_speed_sw Seawater sound speed
#' (\ifelse{html}{\out{c<sub>body</sub>}}{\eqn{c_{body}}},
#' m \ifelse{html}{\out{s<sup>-1</sup>}}{\eqn{s^{-1}}}).
#' @param density_sw Seawater density
#' (\ifelse{html}{\out{&rho;<sub>body</sub>}}{\eqn{\rho_{body}}},
#' kg \ifelse{html}{\out{m<sup>3</sup>}}{\eqn{m^3}})
#' @param alpha_B Numerical coefficient (\ifelse{html}{\out{&alpha;<sub>B</sub>}}{\eqn{\alpha_{B}}}).
#' @export
dcm_initialize <- function(object,
                           # Required inputs
                           frequency,
                           # Optional inputs
                           # Shape (radius)
                           radius_cylinder = NULL ,
                           radius_curvature = NULL ,
                           radius_curvature_ratio = 3.0 ,
                           radius_cylinder_fun = "max" ,
                           # Body parameters
                           length = NULL ,
                           g = NULL ,
                           h = NULL ,
                           theta = NULL ,
                           # Medium
                           sound_speed_sw = 1500 ,
                           density_sw = 1026 ,
                           # Numerical coefficient
                           alpha_B = 0.8 ) {
  # Parse metadata =============================================================
  meta <- acousticTS::extract( object , "metadata" )
  # Parse shape ================================================================
  shape <- acousticTS::extract( object , "shape_parameters" )
  # Parse body =================================================================
  body <- acousticTS::extract( object , "body" )
  # Determine radius to be used for DCM ========================================
  radius_uniform <- base::ifelse( !base::is.null( radius_cylinder ) ,
                                  radius_cylinder,
                                  base::switch( radius_cylinder_fun ,
                                                max = max( body$radius, na.rm = T ) ,
                                                mean = mean( body$radius, na.rm = T ) ,
                                                median = median( body$radius, na.rm = T ) ) )
  # Calculate radius of curvature either based on user input or ratio ==========
  radius_curvature <- base::ifelse( is.null( radius_curvature ) ,
                                    ifelse( ! is.null( body$radius_curvature_ratio ) ,
                                            body$radius_curvature_ratio * shape$length , 
                                            shape$length * radius_curvature_ratio ) ,
                                    radius_curvature )
  # Fill in remaining model inputs required by DCM =============================
  # Orientation ================================================================
  theta <- base::ifelse( !base::is.null( theta ) , theta , body$theta )
  # Density contrast, g ========================================================
  g <- base::ifelse( !base::is.null( g ) , g , body$g )
  # Sound speed contrast, h ====================================================
  h <- base::ifelse( !base::is.null( h ) , h , body$h )
  # Define model dataframe for body parameterization ===========================
  body_params <- base::data.frame( length = shape$length ,
                                   radius = radius_uniform ,
                                   radius_function = base::ifelse( !base::is.null( radius_cylinder ) ,
                                                                   "none",
                                                                   radius_cylinder_fun ) ,
                                   radius_curvature = radius_curvature ,
                                   theta = theta ,
                                   g = g ,
                                   h = h ,
                                   alpha_B = alpha_B )
  # Define model dataframe for medium parameterization =========================
  medium_params <-data.frame( sound_speed = sound_speed_sw ,
                                     density = density_sw )
  # Define model dataframe for acoustic parameterization =======================
  model_params <- list(
    acoustics = data.frame(
      frequency = frequency ,
      k_sw = k( frequency , sound_speed_sw ) ,
      k_b = k( frequency , sound_speed_sw * h ) ) ,
    n_s = shape$n_segments )
  # Tidy up model parameters to insert into object =============================
  slot( object ,
        "model_parameters" )$DCM <- list( parameters = model_params ,
                                          medium = medium_params ,
                                          body = body_params )
  # Add model results slot to scattering object ================================
  slot( object ,
        "model" )$DCM <- data.frame( frequency = frequency ,
                                     sigma_bs = rep( NA ,
                                                     length( frequency ) ) )
  # Output =====================================================================
  return(object)
}
################################################################################
#' Initialize FLS-class object for TS modeling.
#' @param object FLS-class object.
#' @inheritParams krm_initialize
#' @param theta Angle of incident soundwave (pi / 2 is broadside).
#' @export
dwba_initialize <- function( object ,
                             frequency ,
                             sound_speed_sw = 1500 ,
                             density_sw = 1026 ,
                             theta = pi / 2 ) {
  # Parse metadata =============================================================
  meta <- acousticTS::extract( object , "metadata" )
  # Parse shape ================================================================
  shape <- acousticTS::extract( object , "shape_parameters" )
  # Parse body =================================================================
  body <- acousticTS::extract( object , "body" )
  # Define medium parameters ===================================================
  medium_params <- base::data.frame( sound_speed = sound_speed_sw ,
                                     density = density_sw )
  # Define model parameters recipe =============================================
  model_params <- base::list(
    acoustics = base::data.frame(
      frequency = frequency ,
      # Wavenumber (medium) ====================================================
      k_sw = acousticTS::k( frequency ,
                            sound_speed_sw ) ,
      # Wavenumber (fluid) =====================================================
      k_f = acousticTS::k( frequency ,
                           body$h * sound_speed_sw ) ) ,
    ncyl_b = shape$n_segments )
  # Define body parameters recipe ==============================================
  body_params <- base::data.frame(
    radius = body$radius ,
    h = body$h ,
    g = body$g ,
    theta = body$theta )
  # Add model parameters slot to scattering object =============================
  methods::slot( object ,
                 "model_parameters" )$DWBA <- base::list( parameters = model_params ,
                                                          medium = medium_params ,
                                                          body = body_params )
  # Add model results slot to scattering object ================================
  methods::slot( object ,
                 "model" )$DWBA <- base::data.frame( frequency = frequency ,
                                                     sigma_bs = base::rep( NA ,
                                                                           base::length( frequency ) ) )
  # Output =====================================================================
  return( object )
}
################################################################################
#' Initialize FLS-class object for TS modeling.
#' @param object FLS-class object.
#' @param radius_curvature_ratio Radius of curvature ratio (length-to-curvature 
#' ratio).
#' @inheritParams dwba_initialize
#' @export
dwba_curved_initialize <- function( object ,
                                    frequency ,
                                    sound_speed_sw = 1500 ,
                                    density_sw = 1026 ,
                                    radius_curvature_ratio = NULL ,
                                    theta = pi / 2 ) {
  # Parse metadata =============================================================
  meta <- extract( object , "metadata" )
  # Parse shape ================================================================
  shape <- extract( object , "shape_parameters" )
  # Parse body =================================================================
  body <- extract( object , "body" )
  # Bend body ==================================================================
  body <- brake( body , 
                 radius_curvature = ifelse( ! is.null( radius_curvature_ratio ) ,
                                            radius_curvature_ratio , 
                                            body$radius_curvature_ratio ) , 
                 mode = "ratio" )
  # Define medium parameters ===================================================
  medium_params <- base::data.frame( sound_speed = sound_speed_sw ,
                                     density = density_sw )
  # Define model parameters recipe =============================================
  model_params <- base::list(
    acoustics = base::data.frame(
      frequency = frequency ,
      # Wavenumber (medium) ====================================================
      k_sw = acousticTS::k( frequency ,
                            sound_speed_sw ) ,
      # Wavenumber (fluid) =====================================================
      k_f = acousticTS::k( frequency ,
                           body$h * sound_speed_sw ) ) ,
    ncyl_b = shape$n_segments )
  # Define body parameters recipe ==============================================
  body_params <- list(
    rpos = body$rpos ,
    radius = body$radius ,
    h = body$h ,
    g = body$g ,
    theta = body$theta ,
    radius_curvature_ratio = body$radius_curvature_ratio , 
    arc_length = body$arc_length )
  # Update object to reflect curvature =========================================
  
  # Add model parameters slot to scattering object =============================
  methods::slot( object ,
                 "model_parameters" )$DWBA_curved <- base::list( parameters = model_params ,
                                                                 medium = medium_params ,
                                                                 body = body_params )
  # Add model results slot to scattering object ================================
  methods::slot( object ,
                 "model" )$DWBA_curved <- base::data.frame( frequency = frequency ,
                                                            sigma_bs = base::rep( NA ,
                                                                                  base::length( frequency ) ) )
  # Output =====================================================================
  return( object )
}
################################################################################
#' Initialize FLS-class object for SDWBA modeling
#' @param object FLS-class object.
#' @param n_iterations Number of iterations to repeat SDWBA
#' @param n_segments_init Reference number of body segments
#' @param phase_sd_init Reference phase deviation
#' @param length_init Reference body length
#' @param frequency_init Reference frequency
#' @inheritParams dwba_initialize
#' @export
sdwba_initialize <- function( object , 
                              frequency , 
                              sound_speed_sw = 1500 , 
                              density_sw = 1026 , 
                              n_iterations = 100 , 
                              n_segments_init = 14 , 
                              phase_sd_init = sqrt( 2 ) / 2 , 
                              length_init = 38.35e-3 , 
                              frequency_init = 120e3 ) {
  # Parse metadata =============================================================
  meta <- acousticTS::extract( object , "metadata" )
  # Parse shape ================================================================
  shape <- acousticTS::extract( object , "shape_parameters" )
  # Parse body =================================================================
  body <- acousticTS::extract( object , "body" )
  # Define medium parameters ===================================================
  medium_params <- base::data.frame( sound_speed = sound_speed_sw ,
                                     density = density_sw )
  # Define model parameters recipe =============================================
  model_params <- base::list(
    acoustics = base::data.frame(
      frequency = frequency ,
      # Wavenumber (medium) ====================================================
      k_sw = acousticTS::k( frequency ,
                                sound_speed_sw ) ,
      # Wavenumber (fluid) =====================================================
      k_f = acousticTS::k( frequency ,
                               body$h * sound_speed_sw ) ) ,
    n_segments = shape$n_segments 
  )
  # Define stochastic recipe ===================================================
  # First calculate new resampled body shape resolution ++++++++++++++++++++++++
  N_f <- base::ceiling( n_segments_init * ( frequency / frequency_init ) * 
                          ( shape$length / length_init ) )
  N_f_vec <- base::ifelse( N_f > n_segments_init , 
                           N_f , 
                           n_segments_init )
  N_f_idx <- base::unique(
    N_f_vec
  )
  # Now we calculate the new phase standard deviation value ++++++++++++++++++++
  phase_sd <- phase_sd_init * ( n_segments_init / N_f_vec ) * ( shape$length / length_init )
  # Create stochastic recipe +++++++++++++++++++++++++++++++++++++++++++++++++++
  stochastic_params <- lapply( 1 : length( N_f_idx ) ,
                               FUN = function( i ) {
                                 idx <- which( N_f_vec == N_f_idx[ i ] )
                                 object_new <- sdwba_resample( object , 
                                                               n_segments = N_f_idx[ i ] )
                                 body <- acousticTS::extract( object_new , "body" )
                                 n_segments <- N_f_idx[ i ]
                                 phase_sd <- phase_sd[ i ] 
                                 acoustics <- model_params$acoustics[ idx , ]
                                 return( list(
                                   meta_params = data.frame(
                                     n_iterations = n_iterations , 
                                     N0 = n_segments_init , 
                                     f0 = frequency_init , 
                                     L0 = length_init , 
                                     p0 = phase_sd_init
                                   ) ,
                                   body_params = body , 
                                   n_segments = n_segments ,
                                   acoustics = acoustics
                                 ) )
                               }
  )
  # Add model parameters slot to scattering object =============================
  methods::slot( object ,
                 "model_parameters" )$SDWBA <- base::list( parameters = stochastic_params ,
                                                           medium = medium_params ,
                                                           body = body )
  # Add model results slot to scattering object ================================
  methods::slot( object ,
                 "model" )$SDWBA <- base::data.frame( frequency = frequency ,
                                                      f_bs = base::rep( NA ,
                                                                        base::length( frequency ) ) ,
                                                      sigma_bs = base::rep( NA ,
                                                                            base::length( frequency ) ) ,
                                                      TS = base::rep( NA , 
                                                                      base::length( frequency ) ) )
  # Output =====================================================================
  return( object )
}
################################################################################
#' Initialize FLS-class object for SDWBA modeling
#' @param object FLS-class object.
#' @param n_iterations Number of iterations to repeat SDWBA
#' @param n_segments_init Reference number of body segments
#' @param phase_sd_init Reference phase deviation
#' @param length_init Reference body length
#' @param frequency_init Reference frequency
#' @inheritParams dwba_initialize
#' @export
sdwba_curved_initialize <- function( object , 
                                     frequency , 
                                     sound_speed_sw = 1500 , 
                                     density_sw = 1026 , 
                                     n_iterations = 100 , 
                                     n_segments_init = 14 , 
                                     phase_sd_init = sqrt( 2 ) / 2 ,
                                     length_init = 38.35e-3 , 
                                     frequency_init = 120e3 ) {
  # Parse metadata =============================================================
  meta <- acousticTS::extract( object , "metadata" )
  # Parse shape ================================================================
  shape <- acousticTS::extract( object , "shape_parameters" )
  # Parse body =================================================================
  body <- acousticTS::extract( object , "body" )
  # Bend body ==================================================================
  object_copy <- object
  object_copy <- brake( object_copy ,
                        radius_curvature = body$radius_curvature_ratio , 
                        mode = "ratio" )
  # Define medium parameters ===================================================
  medium_params <- base::data.frame( sound_speed = sound_speed_sw ,
                                     density = density_sw )
  # Define model parameters recipe =============================================
  model_params <- base::list(
    acoustics = base::data.frame(
      frequency = frequency ,
      # Wavenumber (medium) ====================================================
      k_sw = acousticTS::k( frequency ,
                                sound_speed_sw ) ,
      # Wavenumber (fluid) =====================================================
      k_f = acousticTS::k( frequency ,
                               body$h * sound_speed_sw ) ) ,
    n_segments = shape$n_segments 
  )
  # Define stochastic recipe ===================================================
  # First calculate new resampled body shape resolution ++++++++++++++++++++++++
  N_f <- base::ceiling( n_segments_init * ( frequency / frequency_init ) * 
                          ( shape$length / length_init ) )
  N_f_vec <- base::ifelse( N_f > n_segments_init , 
                           N_f , 
                           n_segments_init )
  N_f_idx <- base::unique(
    N_f_vec
  )
  # Now we calculate the new phase standard deviation value ++++++++++++++++++++
  phase_sd <- phase_sd_init * ( n_segments_init / N_f_vec ) * ( shape$length / length_init )
  # Create stochastic recipe +++++++++++++++++++++++++++++++++++++++++++++++++++
  stochastic_params <- lapply( 1 : length( N_f_idx ) ,
                               FUN = function( i ) {
                                 idx <- which( N_f_vec == N_f_idx[ i ] )
                                 object_new <- sdwba_resample( object_copy , 
                                                               n_segments = N_f_idx[ i ] )
                                 body <- extract( object_new , "body" )
                                 n_segments <- N_f_idx[ i ]
                                 phase_sd <- phase_sd[ i ] 
                                 acoustics <- model_params$acoustics[ idx , ]
                                 return( list(
                                   meta_params = data.frame(
                                     n_iterations = n_iterations , 
                                     N0 = n_segments_init , 
                                     f0 = frequency_init , 
                                     L0 = length_init , 
                                     p0 = phase_sd_init
                                   ) ,
                                   body_params = body , 
                                   n_segments = n_segments ,
                                   acoustics = acoustics
                                 ) )
                               }
  )
  # Add model parameters slot to scattering object =============================
  methods::slot( object ,
                 "model_parameters" )$SDWBA_curved <- base::list( parameters = stochastic_params ,
                                                                  medium = medium_params ,
                                                                  body = body )
  # Add model results slot to scattering object ================================
  methods::slot( object ,
                 "model" )$SDWBA_curved <- base::data.frame( frequency = frequency ,
                                                             f_bs = base::rep( NA ,
                                                                               base::length( frequency ) ) ,
                                                             sigma_bs = base::rep( NA ,
                                                                                   base::length( frequency ) ) ,
                                                             TS = base::rep( NA , 
                                                                             base::length( frequency ) ) )
  # Output =====================================================================
  return( object )
}
################################################################################
#' Initialize SBF-class object for KRM calculations.
#'
#' @param object SBF-class object
#' @param frequency Frequency (Hz).
#' @param sound_speed_sw Seawater sound speed.
#' @param density_sw Seawater density.
#' @param density_body Optional flesh density input.
#' @param density_swimbladder Optional gas density input.
#' @param sound_speed_body Optional flesh sound speed input.
#' @param sound_speed_swimbladder Optional gas sound speed input.
#' @param theta_body Optional orientation input (relative to incident sound wave).
#' @param theta_swimbladder Optional orientation input (relative to incident sound wave).
#' @export
krm_initialize <- function( object ,
                            frequency ,
                            sound_speed_sw = 1500 ,
                            density_sw = 1026 ,
                            density_body = NULL ,
                            density_swimbladder = NULL ,
                            sound_speed_body = NULL ,
                            sound_speed_swimbladder = NULL ,
                            theta_body = NULL ,
                            theta_swimbladder = NULL ) {
  # Detect scatterer object class ==============================================
  scatterer_type <- base::class( object )
  # Parse metadata =============================================================
  meta <- acousticTS::extract( object , "metadata" )
  # Parse body =================================================================
  body <- acousticTS::extract( object , "body" )
  # Parse shape ================================================================
  shape <- acousticTS::extract( object , "shape_parameters" )
  # Define medium parameters ===================================================
  medium_params <- base::data.frame( sound_speed = sound_speed_sw ,
                                     density = density_sw )
  # Define model parameters recipe =============================================
  model_params <- base::list(
    acoustics = base::data.frame(
      # Frequency (Hz) =========================================================
      frequency = frequency ,
      # Wavenumber (medium) ====================================================
      k_sw = acousticTS::k( frequency ,
                                sound_speed_sw )
    )
  )
  # Class-specific model parameters ============================================
  if ( scatterer_type == "FLS" ) {
    # Define body parameters recipe ============================================
    body_params <- base::data.frame( length = shape$length ,
                                     theta = body$theta ,
                                     density = density_sw * body$g ,
                                     sound_speed = sound_speed_sw * body$h )
    # Body wave number =========================================================
    model_params$acoustics$k_b <- acousticTS::k( frequency ,
                                                     sound_speed_sw * body$h )
    # Define body segment count ================================================
    model_params$ns_b <- base::ncol( body$rpos )
    # Slot new model parameters input for KRM ==================================
    methods::slot( object , "model_parameters" )$KRM <- base::list( parameters = model_params ,
                                                                    medium = medium_params ,
                                                                    body = body_params )
  } else if ( scatterer_type == "SBF" ) {
    # Define body parameters recipe ============================================
    body_params <- base::data.frame( length = shape$body$length ,
                                     theta = body$theta ,
                                     density = body$density ,
                                     sound_speed = body$sound_speed )
    # Pull bladder information =================================================
    bladder <- acousticTS::extract( object , "bladder" )
    # Define bladder parameters recipe =========================================
    bladder_params <- base::data.frame( length = shape$bladder$length ,
                                        theta = bladder$theta ,
                                        density = bladder$density ,
                                        sound_speed = bladder$sound_speed )
    # Body wave number =========================================================
    model_params$acoustics$k_b <- acousticTS::k( frequency ,
                                                     body$sound_speed )
    # Define body segment count ================================================
    model_params$ns_b <- shape$body$n_segments
    # Define bladder segment count =============================================
    model_params$ns_sb <- shape$bladder$n_segments
    # Slot new model parameters input for KRM ==================================
    methods::slot( object , "model_parameters" )$KRM <- base::list( parameters = model_params ,
                                                                    medium = medium_params ,
                                                                    body = body_params ,
                                                                    bladder = bladder_params )
  }
  # Slot new model output for KRM ==============================================
  methods::slot( object , "model" )$KRM <- base::data.frame( frequency ,
                                                             sigma_bs = base::rep( NA ,
                                                                                   base::length( frequency ) ) )
  # Output =====================================================================
  return( object )
}
################################################################################
#' Initialize object for Stanton high-pass approximation
#' @param object ESS-class object.
#' @param radius_shell Radius of shell.
#' @param g_shell Optional shell density contrast.
#' @param h_shell Optional shell sound speed contrast.
#' @inheritParams krm_initialize
#' @export
high_pass_stanton_initialize <- function( object ,
                                          frequency ,
                                          radius_shell ,
                                          g_shell ,
                                          h_shell ,
                                          sound_speed_sw = 1500 ,
                                          density_sw = 1026 ) {
  # Extract model body parameters ==============================================
  shell <- extract( object , "shell" )
  # Define medium parameters ===================================================
  medium_params <- data.frame( sound_speed = sound_speed_sw ,
                               density = density_sw )
  # Define acoustic parameters =================================================
  acoustics <- data.frame( frequency = frequency,
                           k_sw = k( frequency , sound_speed_sw ) )
  # Fill out material properties ===============================================
  shell$h <- ifelse( "h" %in% names( shell ) ,
                     shell$h ,
                     ifelse( "sound_speed" %in% names( shell ) ,
                             shell$sound_speed / medium_params$sound_speed ,
                             NA )
  )
  shell$g <- ifelse( "g" %in% names( shell ) ,
                     shell$g ,
                     ifelse( "density" %in% names( shell ) ,
                             shell$density / medium_params$density ,
                             NA )
  )
  # Calculate wavenumber =======================================================
  acoustics$k_b <- acoustics$k_sw * shell$h
  # Define model dataframe for body parameterization ===========================
  shell_params <- data.frame( radius = shell$radius ,
                              diameter = shell$radius * 2 ,
                              g = shell$g ,
                              h = shell$h )
  # Tidy up model parameters to insert into object =============================
  model_params <- list( shell = shell_params ,
                        medium = medium_params ,
                        acoustics = acoustics )
  # Define object and model parameters =========================================
  slot( object , "model_parameters" )$high_pass_stanton <- model_params
  slot( object , "model" )$high_pass_stanton <- data.frame( frequency ,
                                                            sigma_bs = rep( NA ,
                                                                            length( frequency ) ) )
  return( object )
}