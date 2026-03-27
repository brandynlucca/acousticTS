# SCATTERING CLASS FIXTURES
# CONSTANTS

# ----> SPHERE
fixture_sphere <- function(boundary_type) {

  # Medium
  density_sw <- 1026.8
  sound_speed_sw <- 1477.3

  # Shape
  shell_radius <- 0.01
  shell_thickness <- 0.001
  radius_body <- shell_radius
  radius_shell <- shell_radius
  sphere_shape <- sphere(radius_body = shell_radius, n_segments = 80)

  # Shape
  switch(boundary_type,
         fixed_rigid = ess_generate(shape = sphere_shape,
                                    radius_shell = shell_radius,
                                    shell_thickness = shell_thickness,
                                    density_fluid = 1028.9,
                                    sound_speed_fluid = 1480.3),
         pressure_release = ess_generate(shape = sphere_shape,
                                         radius_shell = shell_radius,
                                         shell_thickness = shell_thickness,
                                         density_fluid = 1028.9,
                                         sound_speed_fluid = 1480.3),
         gas_filled = gas_generate(shape = sphere_shape,
                                   density_fluid = 1.24,
                                   sound_speed_fluid = 345.0),
         liquid_filled = fls_generate(
           shape = sphere_shape,
           g_body = 1028.9 / density_sw,
           h_body = 1480.3 / sound_speed_sw
           ),
         shelled_pressure_release = ess_generate(
           shape = sphere_shape,
           radius_shell = shell_radius,
           shell_thickness = shell_thickness,
           density_shell = 1028.9,
           sound_speed_shell = 1480.3
         ),
         shelled_gas = ess_generate(shape = sphere_shape,
                                    radius_shell = shell_radius,
                                    shell_thickness = shell_thickness,
                                    density_shell = 1070.0,
                                    sound_speed_shell = 1570.0,
                                    density_fluid = 1.24,
                                    sound_speed_fluid = 345.0),
         shelled_liquid = ess_generate(shape = sphere_shape,
                                       radius_shell = shell_radius,
                                       shell_thickness = shell_thickness,
                                       density_shell = 1028.9,
                                       sound_speed_shell = 1480.3,
                                       density_fluid = 1031.0,
                                       sound_speed_fluid = 1483.3)
        )
}

# ----> CYLINDER
fixture_cylinder <- function(boundary_type) {

  # Shape
  length_body <- 0.07
  radius_body <- 0.01
  cylinder_shape <- cylinder(
    length_body = length_body,
    radius_body = radius_body,
    n_segments = 80
  )

  # Medium
  density_sw <- 1026.8
  sound_speed_sw <- 1477.3

  # Shape
  switch(boundary_type,
         fixed_rigid = fls_generate(shape = cylinder_shape,
                                    g_body = NA,
                                    h_body = NA),
         pressure_release = fls_generate(shape = cylinder_shape,
                                         g_body = NA,
                                         h_body = NA),
         gas_filled = gas_generate(shape = cylinder_shape,
                                   density_fluid = 1.24,
                                   sound_speed_fluid = 345.0),
         liquid_filled = fls_generate(
            shape = cylinder_shape,
            g_body = 1028.9 / density_sw,
            h_body = 1480.3 / sound_speed_sw
         )
  )
}

# ----> PROLATE SPHEROID
fixture_ps <- function(boundary_type) {

  # Shape
  length_body <- 0.14
  radius_body <- 0.01
  prolate_shape <- prolate_spheroid(
    length_body = length_body,
    radius_body = radius_body,
    n_segments = 100
  )

  # Medium
  density_sw <- 1026.8
  sound_speed_sw <- 1477.3

  # Shape
  switch(boundary_type,
         fixed_rigid = fls_generate(shape = prolate_shape,
                                    g_body = NA,
                                    h_body = NA),
         pressure_release = fls_generate(shape = prolate_shape,
                                         g_body = NA,
                                         h_body = NA),
         gas_filled = gas_generate(
           shape = prolate_shape,
           density_fluid = 1.24,
           sound_speed_fluid = 345.0
         ),
         liquid_filled = fls_generate(
           shape = prolate_shape,
           g_body = 1028.9 / density_sw,
            h_body = 1480.3 / sound_speed_sw
         )
  )
}

# SPOOF MODEL INITIALIZATION
spoof_initialize <- function(object) {
  object
}
