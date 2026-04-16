library(acousticTS)

compute_essms_modal_coeffs <- function(density_fluid, sound_speed_fluid) {
  object <- ess_generate(
    shape = sphere(radius_body = 0.01, n_segments = 80),
    radius_shell = 0.01,
    shell_thickness = 0.0008,
    density_shell = 2565,
    sound_speed_shell = 3750,
    density_fluid = density_fluid,
    sound_speed_fluid = sound_speed_fluid,
    E = 7.0e10,
    nu = 0.32
  )

  object <- essms_initialize(
    object,
    frequency = 38e3,
    sound_speed_sw = 1477.3,
    density_sw = 1026.8,
    m_limit = 8L
  )

  model <- extract(object, "model_parameters")$ESSMS
  shell <- model$shell
  fluid <- model$fluid
  acoustics <- model$parameters$acoustics

  sound_speed_longitudinal <- sqrt((shell$lambda + 2 * shell$G) / shell$density)
  sound_speed_transversal <- sqrt(shell$G / shell$density)

  ka_matrix <- .calculate_ka_matrix(
    frequency = acoustics$frequency,
    sound_speed_sw = model$medium$sound_speed,
    sound_speed_fluid = fluid$sound_speed,
    sound_speed_longitudinal = sound_speed_longitudinal,
    sound_speed_transversal = sound_speed_transversal,
    radius_shell = shell$radius,
    radius_fluid = fluid$radius
  )

  as.vector(
    elastic_shell_boundary_conditions(
      ka_matrix = ka_matrix,
      m_limit = acoustics$m_limit,
      lambda = shell$lambda,
      mu = shell$G,
      rho_ratio_sw = model$medium$density / shell$density,
      rho_ratio_fl = fluid$density / shell$density
    )[1L, ]
  )
}

test_that("ESSMS shell coefficients respond to the Stanton fluid generalization", {
  coeff_equal_fluids <- compute_essms_modal_coeffs(
    density_fluid = 1026.8,
    sound_speed_fluid = 1477.3
  )
  coeff_contrasted_fluids <- compute_essms_modal_coeffs(
    density_fluid = 1077.3,
    sound_speed_fluid = 1575
  )

  expect_false(
    isTRUE(all.equal(
      coeff_equal_fluids,
      coeff_contrasted_fluids,
      tolerance = 1e-8
    ))
  )

  expect_true(all(is.finite(Mod(coeff_equal_fluids))))
  expect_true(all(is.finite(Mod(coeff_contrasted_fluids))))
})
