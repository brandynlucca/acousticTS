################################################################################
################################################################################
# Model initialization functions
################################################################################
################################################################################
#' Initialize GAS-object for modal series solution.
#' @param object GAS-class object.
#' @param radius_body Radius of sphere (m).
#' @param g_body Density contrast for gas.
#' @param h_body Sound speed contrast for gas.
#' @param ka_limit Modal series limit (i.e. max "m"). The default is the maximum
#'    ka + 10.
#' @inheritParams krm_initialize
#' @export
anderson_initialize <- function(object ,
                                        # Required input
                                        frequency ,
                                        # Optional body parameters
                                        radius_body = NULL ,
                                        g_body = NULL ,
                                        h_body = NULL ,
                                        # Optional medium inputs
                                        sound_speed_sw = 1500 ,
                                        density_sw = 1026 ,
                                        # Optional modal series limit
                                        ka_limit = NULL ) {
  # Detect object class ========================================================
  scatterer_type <- base::class( object )
  # Parse shape ================================================================
  shape <- acousticTS::extract( object ,
                                "shape_parameters" )
  # Parse body =================================================================
  body <- acousticTS::extract( object ,
                               "body" )
  # Define medium parameters ===================================================
  medium_params <- base::data.frame( sound_speed = sound_speed_sw ,
                                     density = density_sw )
  # Define body radii compatible w/ model ======================================
  bladder <- base::switch( scatterer_type ,
                           SBF = acousticTS::extract( object ,
                                                      "bladder" ) ,
                           GAS = body )
  # Equivalent spherical radius ================================================
  if ( scatterer_type == "SBF" ) {
    # Equivalent volume ========================================================
    acousticTS::along_sum( bladder$rpos , shape$bladder$n_segments )
    bladder$rpos[ 2 , ] / 2
  }
  aesr <- base::switch( scatterer_type ,
                        )



  # Define model parameters recipe =============================================
  model_params <- base::list(
    acoustics = base::data.frame(
      frequency = frequency ,
      # Wavenumber (medium) ====================================================
      k_sw = acousticTS::kcalc( frequency ,
                                sound_speed_sw ) ,
      # Wavenumber (fluid) =====================================================
      k_f = acousticTS::kcalc( frequency ,
                               body$h * sound_speed_sw ) ) ,
    ncyl_b = shape$n_segments )
  # Define ka limit ============================================================
  ka_limit <- base::ifelse( !base::is.null( ka_limit ) ,
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
                 "model_parameters" )$anderson <- base::list(
                   parameters = model_params ,
                   medium = medium_params ,
                   body = body_params )
  # Add model results slot to scattering object ================================
  methods::slot( object ,
                 "model" )$anderson <- base::data.frame(
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
      k_sw = acousticTS::kcalc( frequency ,
                                sound_speed_sw ) ,
      # Wavenumber (longitudinal axis) =========================================
      k_l = acousticTS::kcalc( frequency ,
                               body$sound_speed_longitudinal ) ,
      # Wavenumber (tranversal axis) ===========================================
      k_t = acousticTS::kcalc( frequency ,
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
#' @param radius_curvature_ratio Ratio between body length and the radius of
#'    curvature. Defaults to 3.0.
#' @param radius_cylinder_fun Defines which radius value will be used from the
#'    radius vector. Defaults to "max", but also accepts "mean" and "median".
#' @param length Body length (m).
#' @param g Density contrast.
#' @param h Sound speed contrast.
#' @param theta Body orientation relative to the indicent sound wave.
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
                                                max = base::max( body$radius, na.rm = T ) ,
                                                mean = base::mean( body$radius, na.rm = T ) ,
                                                median = base::median( body$radius, na.rm = T ) ) )
  # Calculate radius of curvature either based on user input or ratio ==========
  radius_curvature <- base::ifelse( !base::is.null( radius_curvature ) ,
                                    radius_curvature,
                                    shape$length * radius_curvature_ratio )
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
  medium_params <- base::data.frame( sound_speed = sound_speed_sw ,
                                     density = density_sw )
  # Define model dataframe for acoustic parameterization =======================
  model_params <- base::list(
    acoustics = base::data.frame(
      frequency = frequency ,
      k_sw = acousticTS::kcalc( frequency , sound_speed_sw ) ,
      k_b = acousticTS::kcalc( frequency , sound_speed_sw * h ) ) ,
    n_s = shape$n_segments )
  # Tidy up model parameters to insert into object =============================
  methods::slot( object ,
                 "model_parameters" )$DCM <- base::list( parameters = model_params ,
                                                         medium = medium_params ,
                                                         body = body_params )
  # Add model results slot to scattering object ================================
  methods::slot( object ,
                 "model" )$DCM <- base::data.frame( frequency = frequency ,
                                                    sigma_bs = base::rep( NA ,
                                                                          base::length( frequency ) ) )
  # Output =====================================================================
  return(object)
}
################################################################################
#' Initialize FLS-class object for TS modeling.
#' @param object FLS-class object.
#' @inheritParams krm_initialize
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
      k_sw = acousticTS::kcalc( frequency ,
                                sound_speed_sw ) ,
      # Wavenumber (fluid) =====================================================
      k_f = acousticTS::kcalc( frequency ,
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
#' @param theta Optional orientation input (relative to incident sound wave).
#' @param body_resize Factor resize for body.
#' @param swimbladder_resize Factor resize for swimbladder.
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
      k_sw = acousticTS::kcalc( frequency ,
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
    model_params$acoustics$k_b <- acousticTS::kcalc( frequency ,
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
    model_params$acoustics$k_b <- acousticTS::kcalc( frequency ,
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
stanton_high_pass_initialize <- function(object,
                                         frequency,
                                         radius_shell = NULL,
                                         g_shell = NULL,
                                         h_shell = NULL,
                                         sound_speed_sw = 1500,
                                         density_sw = 1026) {
  # Extract model body parameters ==============================================
  shell <- extract(object, "shell")
  # Define medium parameters ===================================================
  medium_params <- data.frame(sound_speed = sound_speed_sw,
                              density = density_sw)
  # Define acoustic parameters =================================================
  acoustics <- data.frame(frequency = frequency,
                          k_sw = kcalc(frequency, sound_speed_sw))
  acoustics$k_b <- acoustics$k_sw * extract(object, "shell")$h
  # Fill in remaining model inputs required by DCM =============================
  radius <- ifelse(!is.null(radius_shell), radius_shell, shell$radius)
  # Density contrast, g ========================================================
  g <- ifelse(!is.null(g_shell), g_shell, shell$g)
  # Sound speed contrast, h ====================================================
  h <- ifelse(!is.null(h_shell), h_shell, shell$h)
  # Define model dataframe for body parameterization ===========================
  shell_params <- data.frame(radius = radius,
                             diameter = radius * 2,
                             g = g,
                             h = h)
  # Tidy up model parameters to insert into object =============================
  model_params <- list(shell = shell_params,
                       medium = medium_params,
                       acoustics = acoustics)
  # Define object and model parameters =========================================
  slot(object, "model_parameters")$stanton_high_pass <- model_params
  return(object)
}
