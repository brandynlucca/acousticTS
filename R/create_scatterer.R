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
#' @param rho_body Flesh density (\ifelse{html}{\out{&rho;<sub>body</sub>}}{\eqn{\rho_{body}}}, kg \ifelse{html}{\out{m<sup>3</sup>}}{\eqn{m^3}}).
#' @param c_body Flesh sound speed (\ifelse{html}{\out{c;<sub>body</sub>}}{\eqn{c_{body}}}, m \ifelse{html}{\out{s<sup>-1</sup>}}{\eqn{s^{-1}}}).
#' @param rho_bladder Bladder density (\eqn{\rho}, kg \ifelse{html}{\out{m<sup>3</sup>}}{\eqn{m^3}}).
#' @param c_bladder Bladder sound speed (c, m \eqn{s^-1}.
#' @param theta_body Angle of body relative to wavefront (\eqn{\theta_body}, radians).
#' @param theta_bladder Angle of body relative to wavefront (\eqn{\theta_bladder}, radians).
#' @param theta_units Angular units.
#' @param length_units Angular units.
#' @param ID Angular units.
#'
#' @return
#' Generates a SBF-class object.
#' @export
sbf_generate <- function(x_body, w_body, zU_body, zL_body,
                         x_bladder, w_bladder, zU_bladder, zL_bladder,
                         c_body, c_bladder,
                         rho_body, rho_bladder,
                         theta_body=pi/2, theta_bladder=pi/2,
                         theta_units="radians",
                         length_units="m",
                         ID=NULL){
  metadata <- list(ID=ifelse(!is.null(ID),
                             ID,
                             "UID"))
  body <- list(rpos=rbind(x=x_body[!is.na(x_body)],
                          w=w_body[!is.na(w_body)],
                          zU=zU_body[!is.na(zU_body)],
                          zL=zL_body[!is.na(zL_body)]),
               sound_speed=c_body[!is.na(c_body)],
               density=rho_body[!is.na(rho_body)],
               theta=theta_body[!is.na(theta_body)])
  bladder <- list(rpos=rbind(x=x_bladder[!is.na(x_bladder)],
                             w=w_bladder[!is.na(w_bladder)],
                             zU=zU_bladder[!is.na(zU_bladder)],
                             zL=zL_bladder[!is.na(zL_bladder)]),
                  sound_speed=c_bladder[!is.na(c_bladder)],
                  density=rho_bladder[!is.na(rho_bladder)],
                  theta=theta_bladder[!is.na(theta_bladder)])
  shape_parameters <- list(
    body=list(
      length=abs(max(x_body, na.rm=T)-min(x_body, na.rm=T)),
      ncyl=length(x_body[!is.na(x_body)]-1),
      theta_units=theta_units,
      length_units=length_units
    ),
    bladder=list(
      length=abs(max(x_bladder, na.rm=T)-min(x_bladder, na.rm=T)),
      ncyl=length(x_bladder[!is.na(x_bladder)]-1),
      theta_units=theta_units,
      length_units=length_units
    )
  )

  return(new("SBF",
             metadata=metadata,
             body=body,
             bladder=bladder,
             model_parameters=list(),
             shape_parameters=shape_parameters))
}
################################################################################
# Create CAL-class object
################################################################################
#' Generate a CAL-class object.
#' @param material Material-type for a solid sphere. See `Details` for available
#' options. Default is tungsten carbide (WC).
#' @param diameter Spherical diameter (m).
#' @param sound_speed_longitudinal Longitudinal sound speed (m/s).
#' @param sound_speed_transversal Transversal sound speed (m/s).
#' @param density_sphere Density (kg/m^3).
#' @param ID Optional metadata ID input.
#' @param diameter_units Units for diameter. Defaults to "m".
#' @return
#' Generates a CAL-class object.
#' @export
cal_generate <- function(material = "WC",
                         diameter = 38.1e-3,
                         sound_speed_longitudinal = NULL,
                         sound_speed_transversal = NULL,
                         density_sphere = NULL,
                         ID = NULL,
                         diameter_units = "m") {
  metadata <- list(ID = ifelse(!is.null(ID), ID, "Calibration sphere"),
                   Material = material)
  sphere_rpos <- sphere(diameter)

  body <- list(rpos = sphere_rpos,
               diameter = diameter,
               radius = diameter / 2)

  if(is.null(sound_speed_longitudinal) & is.null(sound_speed_transversal) &
     is.null(density_sphere)){
    material_properties <- switch(material,
                                  Cu = list(sound_speed_longitudinal = 4760,
                                            sound_speed_transversal = 2288.5,
                                            density = 8947),
                                  WC = list(sound_speed_longitudinal = 6853,
                                            sound_speed_transversal = 4171,
                                            density = 14900),
                                  Al = list(sound_speed_longitudinal = 6260,
                                            sound_speed_transversal = 3080,
                                            density = 2700),
                                  steel = list(sound_speed_longitudinal = 5610,
                                               sound_speed_transversal = 3120,
                                               density = 7800),
                                  brass = list(sound_speed_longitudinal = 4372,
                                               sound_speed_transversal = 2100,
                                               density = 8360))
  } else {
    material_properties <- list(sound_speed_longitudinal = sound_speed_longitudinal,
                                sound_speed_transversal = sound_speed_transversal,
                                density = density_sphere)
  }

  body <- append(body, material_properties)

  shape_parameters <- list(body=list(diameter = diameter,
                                     radius = diameter / 2,
                                     ncyl = length(body$rpos[1, ]) - 1,
                                     diameter_units = diameter_units))

  return(new("CAL",
             metadata = metadata,
             model_parameters = list(),
             model = list(),
             body = body,
             shape_parameters = shape_parameters))
}
################################################################################
# Create FLS-class object
################################################################################
#' Manually generate a FLS object.
#'
#' @param x_body Vector containing x-axis body (m) shape data.
#' @param y_body Vector containing y-axis body (m) shape data.
#' @param z_body Vector containing z-axis body (m) shape data.
#' @param radius_body Vector containing radii (m).
#' @param g_body Density contrast.
#' @param h_body Soundspeed contrast
#' @param theta_body Orientation of the target relative to the transmit source
#' (\eqn{\theta}). Broadside incidence is considered 90 degrees, or pi/2.
#' Default value is pi/2; input should be in radians.
#' @param theta_units Units used for orientation. Defaults to "radians".
#' @param length_units Units used for position vector. Defaults to "m".
#' @param ID Optional metadata entry.
#' @return
#' Calls in an FLS-class object from a *.csv file
#' @import methods
#' @export
fls_generate <- function(x_body,
                         y_body,
                         z_body,
                         radius_body,
                         g_body,
                         h_body,
                         theta_body=pi/2,
                         theta_units="radians",
                         length_units="m",
                         radius_curvature=NULL,
                         ID=NULL){

  metadata <- list(ID=ifelse(!is.null(ID), ID, "UID"))
  body <- list(rpos=rbind(x=x_body,
                          y=y_body,
                          z=z_body),
               radius=radius_body,
               theta=theta_body,
               g=g_body,
               h=h_body)
  shape_parameters <- list(body=list(length=abs(max(x_body, na.rm=T)-min(x_body, na.rm=T)),
                                     ncyl=length(x_body)-1,
                                     theta_units=theta_units,
                                     length_units=length_units))
  return(new("FLS",
             metadata=metadata,
             model_parameters = list(),
             model = list(),
             body=body,
             shape_parameters=shape_parameters))
}
################################################################################
# Create GAS-class object
################################################################################
#' Create GAS object
#'
#' @inheritParams fls_generate
#' @param diameter_units Diameter units. Defaults to "m".
#' @return
#' Creates a FLS-class object from a *.csv file
#' @import methods
#' @export
gas_generate <- function(radius_body,
                         g_body = 0.0012,
                         h_body = 0.22,
                         ID = NULL,
                         diameter_units = "m") {
  # Create metadata field =================================================
  metadata <- list(ID = ifelse(!is.null(ID), ID, "UID"))
  # Create body shape field ===============================================
  body <- list(rpos = sphere(radius_body),
               radius = radius_body,
               diameter = radius_body * 2,
               g = g_body,
               h = h_body)
  # Shape parameters field ================================================
  shape_parameters <- list(body = list(diameter = radius_body * 2,
                                       radius = radius_body,
                                       ncyl = length(body$rpos[1, ]) - 1,
                                       diameter_units = diameter_units))
  # Create GAS-class object ===============================================
  return(new("GAS",
             metadata = metadata,
             model_parameters = list(),
             model = list(),
             body = body,
             shape_parameters = shape_parameters))
}

################################################################################
# Create GAS-class object
################################################################################
#' Generate ESS shape
#'
#' @inheritParams fls_generate
#' @param radius_shell Radius of shell (m).
#' @param shell_thickness Optional shell thickness (m).
#' @param g_fluid Optional density contrast for fluid-like body.
#' @param h_fluid Optional sound speed contrast for fluid-like body.
#' @param g_shell Density contrast for the shell.
#' @param h_shell Sound speed contrast for the shell.
#' @param theta_shell Object orientation relative to incident sound wave.
#' @export
ess_generate <- function(x_body = NULL,
                         y_body = NULL,
                         z_body = NULL,
                         radius_shell,
                         shell_thickness = NULL,
                         g_fluid = NULL,
                         h_fluid = NULL,
                         g_shell,
                         h_shell,
                         theta_shell = pi / 2,
                         ID = NULL,
                         theta_units = "radians",
                         length_units = "m") {
  # Create metadata field ======================================================
  metadata <- list(ID = ifelse(!is.null(ID), ID, "UID"))
  # Create shell shape field ===================================================
  if(is.null(x_body) & is.null(y_body) & is.null(z_body)) {
    shell_rpos <- sphere(radius_shell)
  } else {
    shell_rpos <- rbind(x = x_body,
                        y = y_body,
                        z = z_body)
  }
  shell <- list(rpos = shell_rpos,
                radius = radius_shell,
                g = g_shell,
                h = h_shell,
                theta = theta_shell,
                shell_thickness = ifelse(!is.null(shell_thickness),
                                         shell_thickness,
                                         NA))
  # Create body shape field ====================================================
  body <- list(g = ifelse(!is.null(g_fluid), g_fluid, NA),
               h = ifelse(!is.null(h_fluid), h_fluid, NA))
  # Shape parameters field =====================================================
  shape_parameters <- list(body = list(diameter = radius_shell * 2,
                                       radius = radius_shell,
                                       ncyl = length(body$rpos[1, ]) - 1,
                                       length_units = length_units))
  # Create ESS-class object ====================================================
  return(new("ESS",
             metadata = metadata,
             model_parameters = list(),
             model = list(),
             shell = shell,
             body = body,
             shape_parameters = shape_parameters))
}
