################################################################################
# CREATE SCATTERER FUNCTIONS
################################################################################
################################################################################
# Create SBF-class object
################################################################################
#' Manually generate a SBF-class object.
#'
#' @param x_body Vector containing along-body axis (m).
#' @param w_body Vector containing across-body axis (m).
#' @param zU_body Vector containing dorsal-body axis (m).
#' @param zL_body Vector containing ventral-body axis (m).
#' @param x_bladder Vector containing along-bladder axis (m).
#' @param w_bladder Vector containing across-bladder axis (m).
#' @param zU_bladder Vector containing dorsal-bladder axis (m).
#' @param zL_bladder Vector containing ventral-bladder axis (m).
#' @param density_body Flesh density (\ifelse{html}{\out{&rho;<sub>body</sub>}}{\eqn{\rho_{body}}}, kg \ifelse{html}{\out{m<sup>3</sup>}}{\eqn{m^3}}).
#' @param sound_speed_body Flesh sound speed (\ifelse{html}{\out{c;<sub>body</sub>}}{\eqn{c_{body}}}, m \ifelse{html}{\out{s<sup>-1</sup>}}{\eqn{s^{-1}}}).
#' @param density_bladder Bladder density (\eqn{\rho}, kg \ifelse{html}{\out{m<sup>3</sup>}}{\eqn{m^3}}).
#' @param sound_speed_bladder Bladder sound speed (c, m \eqn{s^-1}.
#' @param theta_body Angle of body relative to wavefront (\eqn{\theta_body}, radians).
#' @param theta_bladder Angle of body relative to wavefront (\eqn{\theta_bladder}, radians).
#' @param theta_units Angular units.
#' @param length_units Angular units.
#' @param ID Angular units.
#'
#' @return
#' Generates a SBF-class object.
#' 
#' @seealso \code{\link{SBF}}
#' @export
sbf_generate <- function( x_body ,
                          w_body ,
                          zU_body ,
                          zL_body ,
                          x_bladder ,
                          w_bladder ,
                          zU_bladder ,
                          zL_bladder ,
                          sound_speed_body ,
                          sound_speed_bladder ,
                          density_body ,
                          density_bladder ,
                          theta_body = pi / 2 ,
                          theta_bladder = pi / 2 ,
                          theta_units = "radians" ,
                          length_units = "m" ,
                          ID = NULL ) {
  # Generate shape position matrix =============================================
  # Create body shape field ====================================================
  shape_body <- "Arbitrary"
  shape_bladder <- "Arbitrary"
  # Define body shape ==========================================================
  body <- base::list( rpos = base::rbind( x = x_body[ !base::is.na( x_body ) ] ,
                                          w = w_body[ !base::is.na( w_body ) ] ,
                                          zU = zU_body[ !base::is.na( zU_body ) ] ,
                                          zL = zL_body[ !base::is.na( zL_body ) ] ) ,
                      sound_speed = sound_speed_body[ !base::is.na( sound_speed_body ) ] ,
                      density = density_body[ !base::is.na( density_body ) ] ,
                      theta = theta_body[ !base::is.na( theta_body ) ] )
  # Define bladder shape =========================================================
  bladder <- base::list( rpos = base::rbind( x = x_bladder[ !base::is.na( x_bladder ) ] ,
                                             w = w_bladder[ !base::is.na( w_bladder ) ] ,
                                             zU = zU_bladder[ !base::is.na( zU_bladder ) ] ,
                                             zL = zL_bladder[ !base::is.na( zL_bladder ) ] ) ,
                         sound_speed = sound_speed_bladder[ !base::is.na( sound_speed_bladder ) ] ,
                         density = density_bladder[ !base::is.na( density_bladder ) ] ,
                         theta = theta_bladder[ !base::is.na( theta_bladder ) ] )
  # Define shape parameters ====================================================
  shape_parameters <- base::list(
    body = base::list(
      shape = base::ifelse( base::class( shape_body ) == "character" ,
                            shape_body ,
                            "Arbitrary" ) ,
      length = base::max( body$rpos[ 1 , ] ) ,
      n_segments = base::length( body$rpos[ 1 , ] )
    ) ,
    bladder = base::list(
      shape = base::ifelse( base::class( shape_bladder ) == "character" ,
                            shape_bladder ,
                            "Arbitrary" ) ,
      length = base::max( bladder$rpos[ 1 , ] ) - base::min( bladder$rpos[ 1 , ] ) ,
      n_segments = base::length( bladder$rpos[ 1 , ] )
    ) ,
    length_units = length_units ,
    theta_units = theta_units
  )
  # Create metadata field ======================================================
  metadata <- base::list( ID = base::ifelse ( !base::is.null( ID ) ,
                                              ID ,
                                              "UID" ) )
  # Create FLS-class object ====================================================
  return( methods::new( "SBF" ,
                        metadata = metadata ,
                        model_parameters = base::list( ) ,
                        model = base::list( ) ,
                        body = body ,
                        bladder = bladder ,
                        shape_parameters = shape_parameters ) )
}
################################################################################
# Create CAL-class object
################################################################################
#' Generate a CAL-class object.
#' @param material Material-type for a solid sphere. See `Details` for available
#' options. Default is tungsten carbide (WC).
#' @param diameter Spherical diameter (m).
#' @param n_segments Number of segments to discretize object shape.
#' @param sound_speed_longitudinal Longitudinal sound speed (m/s).
#' @param sound_speed_transversal Transversal sound speed (m/s).
#' @param density_sphere Density (kg/m^3).
#' @param theta_sphere Backscattering direction (Default: pi radians).
#' @param ID Optional metadata ID input.
#' @param diameter_units Units for diameter. Defaults to "m".
#' @param theta_units Units for direction. Defaults to "radians".
#' @param material Material-type for the soldi sphere. See 'Details' built-in
#' material options.
#'
#' @details
#' There are several options for the \strong{material} argument:
#' \tabular{rlllll}{
#'  \strong{Material} \tab \strong{Argument} \tab \strong{c1} \tab \strong{c2}
#'  \tab \strong{\eqn{\rho1}}\cr
#'  \emph{Tungsten carbide} \tab "WC" \tab 6853 \tab 4171 \tab 14900\cr
#'  \emph{Stainless steel} \tab "steel" \tab 5980 \tab 3297 \tab 7970\cr
#'  \emph{Brass} \tab "brass" \tab 4372 \tab 2100 \tab 8360\cr
#'  \emph{Copper} \tab "Cu" \tab 4760 \tab 2288.5 \tab 8947\cr
#'  \emph{Aluminum} \tab "Al" \tab 6260 \tab 3080 \tab 2700\cr
#' }
#' @return
#' Generates a CAL-class object.
#' 
#' @seealso \code{\link{CAL}}
#' @export
cal_generate <- function( material = "WC" ,
                          diameter = 38.1e-3 ,
                          sound_speed_longitudinal = NULL ,
                          sound_speed_transversal = NULL ,
                          density_sphere = NULL ,
                          theta_sphere = pi ,
                          ID = NULL ,
                          diameter_units = "m" ,
                          theta_units = "radians" ,
                          n_segments = 1e2 ) {
  # Define user input or default object ID =====================================
  metadata <- list(
    ID = ifelse( !is.null( ID ) ,
                 ID ,
                 "Calibration sphere" ),
    Material = material )
  # Create sphere object to define definitions =================================
  sphere_shape <- sphere( radius_body = diameter / 2 ,
                          n_segments = n_segments ,
                          diameter_units = "m" )
  # Define calibration sphere body shape =======================================
  body <- list( rpos = sphere_shape@position_matrix ,
                diameter = diameter ,
                radius = diameter / 2 ,
                theta = theta_sphere )
  # Define material properties =================================================
  material_properties <- base::switch(
    material ,
    Cu = list(sound_speed_longitudinal = 4760 ,
              sound_speed_transversal = 2288.5 ,
              density = 8947 ) ,
    WC = list(sound_speed_longitudinal = 6853 ,
              sound_speed_transversal = 4171 ,
              density = 14900 ) ,
    Al = list(sound_speed_longitudinal = 6260 ,
              sound_speed_transversal = 3080 ,
              density = 2700 ) ,
    steel = list(sound_speed_longitudinal = 5610 ,
                 sound_speed_transversal = 3120 ,
                 density = 7800 ) ,
    brass = list(sound_speed_longitudinal = 4372 ,
                 sound_speed_transversal = 2100 ,
                 density = 8360 )
  )
  if( !is.null( sound_speed_longitudinal ) ) {
    material_properties$sound_speed_longitudinal <- sound_speed_longitudinal
  }
  if( !is.null( sound_speed_transversal ) ) {
    material_properties$sound_speed_transversal <- sound_speed_transversal
  }
  if( !is.null( density_sphere ) ) {
    material_properties$density <- density_sphere
  }
  # Append material properties to the shape body ===============================
  body <- base::append(
    body ,
    material_properties
  )
  # Define shape parameters ====================================================
  shape_parameters <- base::list(
    diameter = diameter ,
    radius_body = diameter / 2 ,
    n_segments = n_segments ,
    diameter_units = diameter_units ,
    theta_units = theta_units
  )
  # Generate calibration sphere object =========================================
  return( new( "CAL" ,
               metadata = metadata ,
               model_parameters = base::list( ) ,
               model = base::list( ) ,
               body = body ,
               shape_parameters = shape_parameters ) )
}
################################################################################
# Create FLS-class object
################################################################################
#' Manually generate a FLS object.
#' @param shape Optional input argument that dictates shape-type, if desired, for
#' generalized shapes.
#' @param x_body Vector containing x-axis body (m) shape data.
#' @param y_body Vector containing y-axis body (m) shape data.
#' @param z_body Vector containing z-axis body (m) shape data.
#' @param length_body Optional input for a generic length value input.
#' @param radius_body Vector containing radii (m).
#' @param n_segments Number of body segments.
#' @param radius_curvature_ratio Length-to-curvature ratio (pc/L).
#' @param g_body Density contrast. This can either be a single value (i.e. 
#' homogenous) or a vector of values (i.e. inhomogenous).
#' @param h_body Soundspeed contrast. This can either be a single value (i.e. 
#' homogenous) or a vector of values (i.e. inhomogenous).
#' @param theta_body Orientation of the target relative to the transmit source
#' (\eqn{\theta}). Broadside incidence is considered 90 degrees, or pi/2.
#' Default value is pi/2; input should be in radians.
#' @param theta_units Units used for orientation. Defaults to "radians".
#' @param length_units Units used for position vector. Defaults to "m".
#' @param ID Optional metadata entry.
#' @param ... Additional parameters.
#' @return
#' FLS-class object
#' 
#' @seealso \code{\link{FLS}}
#' @import methods
#' @export
fls_generate <- function( shape = "arbitrary" ,
                          x_body = NULL ,
                          y_body = NULL,
                          z_body = NULL ,
                          length_body = NULL ,
                          radius_body = NULL ,
                          radius_curvature_ratio = NULL ,
                          n_segments = 18 ,
                          g_body ,
                          h_body ,
                          theta_body = pi / 2 ,
                          ID = NULL ,
                          length_units = "m" ,
                          theta_units = "radians" , ... ) {
  # Collect shape information if provided ======================================
  if( !is( shape , "Shape" ) ) {
    if ( shape != "arbitrary" ) {
      if ( base::is.null( length_body ) )
        base::stop( "Body shape is not appropriately parameterized." )
    } else if ( base::is.null( x_body ) ) {
      base::stop( "Body shape is not appropriately parameterized." )
    }
  }
  # Generate shape position matrix =============================================
  # Create body shape field ====================================================
  if ( is( shape , "Shape" ) ) {
    shape_input <- shape
  } else {
    if ( shape == "arbitrary" ) {
      shape_input <- arbitrary( x_body = x_body ,
                                y_body = y_body ,
                                z_body = z_body ,
                                radius_body = radius_body )
    } else {
      # Pull argument input names ==============================================
      arg_pull <- as.list( match.call( ) )
      # Grab input arguments ===================================================
      arg_list <- names( formals( shape ) )
      # Filter out inappropriate parameters ====================================
      arg_full <- arg_pull[ arg_list ] 
      true_args <- Filter( Negate( is.null) , arg_full )
      # Initialize =============================================================
      shape_input <- do.call( shape , true_args )
    }
  }
    # Define shape parameters ==================================================
    shape_parameters <- base::list(
      length = base::max( shape_input@position_matrix[ , 1 ] , na.rm = T ) -
        base::min( shape_input@position_matrix[ , 1 ] , na.rm = T ) ,
      radius = base::max( shape_input@shape_parameters$radius , na.rm = T ) ,
      n_segments = length( shape_input@position_matrix[ , 1 ]  ) - 1 ,
      length_units = length_units ,
      theta_units = theta_units ,
      shape = base::ifelse( is.character( shape ) | is( shape, "Arbitrary" ) ,
                            "Arbitrary" ,
                            class( shape ) ) )
    if ( is( shape , "Sphere" ) || ( is.character( shape ) && shape == "Sphere" ) ) {
      shape_parameters[[ "radius_shape" ]] <- extract( shape_input ,
                                                       "shape_parameters" )$radius_shape
    }
    
    if( is.character( shape ) ){
      if ( shape == "cylinder" ) {
        shape_parameters$taper_order <- shape_input@shape_parameters$taper_order
      }
    }
    # Check material properties length =========================================
    # g ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    if ( length( g_body ) > 1 && 
         length( g_body ) != length( shape_input@position_matrix[ , 1 ] ) - 1 ){
      stop(
        paste0("Vector input for 'g_body' with ", length( g_body ) ,
               " elements does not match the expected number of segments (",
               length( shape_input@position_matrix[ , 1 ] ) - 1 , ")"
        )
      )
    }
    # h ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    if ( length( h_body ) > 1 && 
         length( h_body ) != length( shape_input@position_matrix[ , 1 ] ) - 1 ){
      stop(
        paste0("Vector input for 'h_body' with ", length( h_body ) ,
               " elements does not match the expected number of segments (",
               length( shape_input@position_matrix[ , 1 ] ) - 1 , ")"
        )
      )
    }
    # Create body slot =========================================================
    body <- list( rpos = t.default( shape_input@position_matrix ) ,
                  radius = shape_input@shape_parameters$radius ,
                  radius_curvature_ratio = radius_curvature_ratio ,
                  theta = theta_body ,
                  g = g_body ,
                  h = h_body )
  # Create metadata field ======================================================
  metadata <- base::list( ID = base::ifelse ( !base::is.null( ID ) ,
                                              ID ,
                                              "UID" ) )
  # Create FLS-class object ====================================================
  return( new( "FLS" ,
               metadata = metadata ,
               model_parameters = base::list( ) ,
               model = base::list( ) ,
               body = body ,
               shape_parameters = shape_parameters ) )
}
################################################################################
# Create GAS-class object
################################################################################
#' Create GAS object
#'
#' @inheritParams fls_generate
#' @param shape Optional pre-made shape input. Default is a sphere.
#' @param radius_body Radius (m). For non-canonical shapes, this would be the 
#' maximum or mean radius at the scatterer midsection.
#' @param h_fluid Sound speed contrast of fluid relative to surrounding
#' medium (h).
#' @param g_fluid Density contrast of fluid relative to surrounding density (g).
#' @param sound_speed_fluid Optional fluid sound speed (m/s).
#' @param density_fluid Optional fluid density (m/s).
#' @param radius_units Diameter units. Defaults to "m".
#' @param n_segments Number of body segments.
#' @return
#' GAS-class object
#' 
#' @seealso \code{\link{GAS}}
#' 
#' @import methods
#' @export
gas_generate <- function( shape = "sphere" ,
                          radius_body ,
                          h_fluid = 0.2200 ,
                          g_fluid = 0.0012 ,
                          sound_speed_fluid = NULL ,
                          density_fluid = NULL ,
                          theta_body = pi / 2 ,
                          ID = NULL ,
                          radius_units = "m" ,
                          theta_units = "radians" ,
                          n_segments = 100 ) {
  # Collect shape information if provided ======================================
  if ( base::is.null( radius_body ) & base::class( shape ) == "character" ) {
    stop( "Canonical shape generation requires 'double' input for radius_body argument." )
  }
  # Pull argument input names ==================================================
  arg_pull <- as.list( match.call( ) )
  # Grab input arguments =======================================================
  arg_list <- names( formals( shape ) )
  # Filter out inappropriate parameters ========================================
  arg_full <- arg_pull[ arg_list ] 
  true_args <- Filter( Negate( is.null ) , arg_full )
  # Initialize =================================================================
  shape_input <- do.call( shape , true_args )
  # Create metadata field ======================================================
  metadata <- base::list( ID = base::ifelse ( !base::is.null( ID ) ,
                                              ID ,
                                              "UID" ) )
  # # Create body shape field ==================================================
  body <- base::list( rpos = shape_input@position_matrix ,
                      radius = base::ifelse( shape %in%
                                               base::c( "sphere" ,
                                                        "prolate_spheroid" ) ,
                                             radius_body ,
                                             base::mean(
                                               base::abs(
                                                 base::diff(
                                                   base::t(
                                                     shape_input@position_matrix[ , base::c( 2, 3 ) ]
                                                     ) ) ) ) ) ,
                      theta = theta_body ,
                      g = g_fluid ,
                      h = h_fluid )
  # Define shape parameters ====================================================
  shape_parameters <- base::list(
    radius_body = body$radius ,
    n_segments = n_segments ,
    radius_units = radius_units ,
    theta_units = theta_units ,
    shape = base::ifelse( base::class( shape ) == "character" ,
                          shape ,
                          "Arbitrary" )
  )
  # Create GAS-class object ====================================================
  return( methods::new( "GAS" ,
                        metadata = metadata ,
                        model_parameters = base::list( ) ,
                        model = base::list( ) ,
                        body = body ,
                        shape_parameters = shape_parameters ) )
}
################################################################################
# Create GAS-class object
################################################################################
#' Generate ESS shape
#' @inheritParams fls_generate
#' @param radius_shell Radius of shell (m).
#' @param shell_thickness Optional shell thickness (m).
#' @param g_fluid Optional density contrast for fluid-like body.
#' @param density_fluid Optional density for fluid-like body (kg/m³).
#' @param h_fluid Optional sound speed contrast for fluid-like body.
#' @param sound_speed_fluid Optional sound speed for fluid-like body (m/s).
#' @param g_shell Density contrast for the shell.
#' @param density_shell Optional density for the shell (kg/m³).
#' @param h_shell Sound speed contrast for the shell.
#' @param sound_speed_shell Optional sound speed for the shell (m/s).
#' @param E Young's modulus (Pa) of the shell material.
#' @param nu Poisson's ratio (Dimensionless) of the shell material. 
#' @param G Shear modulus (Pa) of the shell material.
#' @param K Bulk modulus (Pa) of the shell material. 
#' @param theta_shell Object orientation relative to incident sound wave.
#' @return ESS-class object
#' 
#' @seealso \code{\link{ESS}}
#' 
#' @export
ess_generate <- function( shape = "sphere" ,
                          x_body = NULL ,
                          y_body = NULL ,
                          z_body = NULL ,
                          radius_shell ,
                          shell_thickness = NULL ,
                          g_fluid = NULL ,
                          density_fluid = NULL ,
                          h_fluid = NULL ,
                          sound_speed_fluid = NULL ,
                          g_shell = NULL ,
                          density_shell = NULL ,
                          h_shell = NULL ,
                          sound_speed_shell = NULL ,
                          E = NULL ,
                          G = NULL ,
                          K = NULL ,
                          nu = NULL ,
                          theta_shell = pi / 2 ,
                          ID = NULL ,
                          theta_units = "radians" ,
                          length_units = "m" ) {
  # Create metadata field ======================================================
  metadata <- list( ID = ifelse( !is.null( ID ) , ID , "UID" ) )
  # Create shape fields ========================================================
  if ( is( shape , "Shape" ) ) {
    shell_rpos <- shape
  } else {
    if ( shape == "Arbitrary" ) {
      shell_rpos <- arbitrary( x_body = x_body ,
                               y_body = y_body ,
                               z_body = z_body ,
                               radius_body = radius_shell )
    } else {
      # Pull argument input names ==============================================
      arg_pull <- as.list( match.call( ) )
      # Create container for different position matrices =======================
      rpos <- list()
      # Rename radius_shell to radius for shape functionality ==================
      if ( "radius_shell" %in% names( arg_pull ) ) {
        # Create copy ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        subarg_pull <- arg_pull
        # Create shape +++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        subarg_pull$radius_body <- arg_pull$radius_shell
        # Grab input arguments +++++++++++++++++++++++++++++++++++++++++++++++++
        arg_list <- names( formals( shape ) )
        # Filter out inappropriate parameters ++++++++++++++++++++++++++++++++++
        arg_full <- subarg_pull[ arg_list ] 
        true_args <- Filter( Negate( is.null) , arg_full )
        # Initialize +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        rpos[[ "shell" ]] <- do.call( shape , true_args )        
      }
      # Create radius_fluid, if present for shape functionality ================
      if ( "shell_thickness" %in% names( arg_pull ) && 
           !is.null( arg_pull$shell_thickness ) ) {
        # Create copy ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        subarg_pull <- arg_pull
        # Create shape +++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        subarg_pull$radius_body <- arg_pull$radius_shell - arg_pull$shell_thickness
        # Grab input arguments +++++++++++++++++++++++++++++++++++++++++++++++++
        arg_list <- names( formals( shape ) )
        # Filter out inappropriate parameters ++++++++++++++++++++++++++++++++++
        arg_full <- subarg_pull[ arg_list ] 
        true_args <- Filter( Negate( is.null) , arg_full )
        # Initialize +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        rpos[[ "fluid" ]] <- do.call( shape , true_args )             
      }
    }
  }
  # Assign elastic properties ==================================================
  elastic_params <- Filter( Negate( is.null ) , 
                            list( E = E, nu = nu, G = G, K = K ) )
  # Iterate through to calculate the requisite parameters ++++++++++++++++++++++
  repeat {
    param_len <- length( elastic_params )
    if ( is.null( elastic_params$nu ) ) elastic_params$nu <- tryCatch(
      pois( elastic_params$K , elastic_params$E , elastic_params$G ) , 
      error = function( e ) NULL )
    if ( is.null( elastic_params$K ) ) elastic_params$K <- tryCatch(
      bulk( elastic_params$E , elastic_params$G , elastic_params$nu ) , 
      error = function( e ) NULL )
    if ( is.null( elastic_params$E ) ) elastic_params$E <- tryCatch(
      young( elastic_params$K , elastic_params$G , elastic_params$nu ) , 
      error = function( e ) NULL )
    if ( is.null( elastic_params$G ) ) elastic_params$G <- tryCatch(
      shear( elastic_params$K , elastic_params$E , elastic_params$nu ) , 
      error = function( e ) NULL )
    if ( is.null( elastic_params$lambda ) ) elastic_params$lambda <- tryCatch(
      lame( elastic_params$K , elastic_params$E , elastic_params$G , 
            elastic_params$nu ) , 
      error = function( e ) NULL )
    if ( length( elastic_params ) == param_len ) break
  }
  # Handle material properties =================================================
  material_properties <- list()
  # Validate parameter combinations ++++++++++++++++++++++++++++++++++++++++++++
  conflicts <- list(
    c( "g_shell" , "density_shell" ),
    c( "h_shell" , "sound_speed_shell" ) , 
    c( "g_fluid" , "density_fluid" ) ,
    c( "h_fluid" , "sound_speed_fluid" )
  )
  for ( pair in conflicts ) {
    if ( !is.null( get( pair[ 1 ] ) ) && !is.null( get( pair[ 2 ] ) ) ) {
      stop( paste( "Cannot specify both" , pair[ 1 ] , "and", pair[ 2 ] ) )
    }
  }
  # Shell material properties ++++++++++++++++++++++++++++++++++++++++++++++++++
  material_properties[[ "shell" ]] <- c(
    Filter( Negate( is.null ) , 
            list( g = g_shell , density = density_shell , 
                  h = h_shell , sound_speed = sound_speed_shell ) ) , 
    elastic_params
  )
  # Fluid material properties ++++++++++++++++++++++++++++++++++++++++++++++++++
  material_properties[[ "fluid" ]] <- Filter( 
    Negate( is.null ) , 
    list( g = g_fluid , density = density_fluid ,
          h = h_fluid , sound_speed = sound_speed_fluid ) )
  # Finalize shell slot ========================================================
  shell <- c(
    list(
      rpos = rpos[[ "shell" ]]@position_matrix ,
      radius = rpos[[ "shell" ]]@shape_parameters$radius ,
      shell_thickness = ifelse( !is.null( shell_thickness ) ,
                                shell_thickness ,
                                NA ) ,
      theta = theta_shell 
    ) ,
    material_properties[[ "shell" ]]
  )
  # Finalize fluid slot ========================================================
  fluid <- c(
    list(
      rpos = if ( "fluid" %in% names( rpos ) ) {
        rpos[[ "fluid" ]]@position_matrix
      } else { NULL } ,
      radius = if ( "fluid" %in% names( rpos ) ) {
        rpos[[ "fluid" ]]@shape_parameters$radius
      } else { NULL } ,
      theta = theta_shell
    ) ,
    material_properties[[ "fluid" ]]
  )
  # Shape parameters field =====================================================
  shape_parameters <- list(
    shell = list(
      diameter = radius_shell * 2 ,
      radius = radius_shell ,
      length_units = length_units 
    ) ,
    fluid = list(
      diameter = ifelse( is.null( fluid[[ "radius" ]] ) ,
                         NA ,
                         fluid[[ "radius" ]] * 2 ) ,
      radius = ifelse( is.null( fluid[[ "radius" ]] ) ,
                       NA ,
                       fluid[[ "radius" ]] ) ,
      length_units = length_units 
    )
  )
  # Create ESS-class object ====================================================
  return( new( "ESS" ,
               metadata = metadata ,
               model_parameters = list( ) ,
               model = list( ) ,
               shell = shell ,
               fluid = fluid ,
               shape_parameters = shape_parameters ) )
}
