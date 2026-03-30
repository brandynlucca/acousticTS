library(acousticTS)

test_that("hpa_initialize validates method inputs and supported shapes", {
  sphere_obj <- gas_generate(
    shape = sphere(radius_body = 0.01, n_segments = 40),
    density_fluid = 1.24,
    sound_speed_fluid = 345
  )
  cylinder_obj <- fls_generate(
    shape = cylinder(length_body = 0.02, radius_body = 0.001, n_segments = 20),
    g_body = 1.02,
    h_body = 1.03,
    theta_body = 1.1
  )
  oblate_obj <- gas_generate(
    shape = oblate_spheroid(
      length_body = 0.02,
      radius_body = 0.01,
      n_segments = 40
    ),
    density_fluid = 1.24,
    sound_speed_fluid = 345
  )

  expect_error(
    acousticTS:::hpa_initialize(sphere_obj, frequency = 38000, method = "bad"),
    "'method' must be one of the following"
  )
  expect_error(
    acousticTS:::hpa_initialize(
      sphere_obj,
      frequency = 38000,
      deviation_fun = function(ka, extra) ka + extra
    ),
    "must only comprise a single argument"
  )
  expect_error(
    acousticTS:::hpa_initialize(
      sphere_obj,
      frequency = 38000,
      deviation_fun = c(1, 2)
    ),
    "must be scalar"
  )
  expect_error(
    acousticTS:::hpa_initialize(
      sphere_obj,
      frequency = 38000,
      null_fun = function(ka, extra) ka + extra
    ),
    "must only comprise a single argument"
  )
  expect_error(
    acousticTS:::hpa_initialize(
      sphere_obj,
      frequency = 38000,
      null_fun = c(1, 2)
    ),
    "must be scalar"
  )
  expect_error(
    acousticTS:::hpa_initialize(cylinder_obj, frequency = 38000, method = "johnson"),
    "requires scatterer to one of shape-type 'Sphere'"
  )
  expect_error(
    acousticTS:::hpa_initialize(oblate_obj, frequency = 38000, method = "stanton"),
    "requires scatterer to one of the following shape-types"
  )
})

test_that("hpa helper formulas match their analytical expressions", {
  acoustics <- data.frame(frequency = 1000, k_sw = 2)

  sphere_body <- list(
    shape = "Sphere",
    radius = c(0.08, 0.1),
    g = 1.2,
    h = 0.9,
    theta = 1.1,
    length = 0.2
  )
  cylinder_body <- list(
    shape = "Cylinder",
    radius = 0.03,
    g = 1.05,
    h = 0.98,
    theta = 1.1,
    length = 0.3
  )
  prolate_body <- list(
    shape = "ProlateSpheroid",
    radius = 0.03,
    g = 1.05,
    h = 0.98,
    theta = 1.1,
    length = 0.3
  )

  alpha_sphere <- (1 - sphere_body$g * sphere_body$h^2) /
    (3 * sphere_body$g * sphere_body$h^2) +
    (1 - sphere_body$g) / (1 + 2 * sphere_body$g)
  alpha_cylinder <- (1 - cylinder_body$g * cylinder_body$h^2) /
    (2 * cylinder_body$g * cylinder_body$h^2) +
    (1 - cylinder_body$g) / (1 + cylinder_body$g)

  expect_equal(acousticTS:::.alpha_pi(sphere_body), alpha_sphere)
  expect_equal(acousticTS:::.alpha_pi(cylinder_body), alpha_cylinder)

  ka_sphere <- acoustics$k_sw * max(sphere_body$radius)
  johnson_expected <- (sphere_body$radius^2 * ka_sphere^4 * alpha_sphere^2) /
    (1 + 1.5 * ka_sphere^4)
  expect_equal(
    acousticTS:::.johnson_hp(acoustics, sphere_body, alpha_sphere),
    johnson_expected
  )

  sphere_parameters <- list(
    null_fun = function(ka) ka + 1,
    deviation_fun = function(ka) ka^2 + 1
  )
  sphere_r <- (sphere_body$g * sphere_body$h - 1) /
    (sphere_body$g * sphere_body$h + 1)
  sphere_expected <- (
    max(sphere_body$radius)^2 * ka_sphere^4 * alpha_sphere^2 *
      sphere_parameters$null_fun(ka_sphere)
  ) / (
    1 + (4 * ka_sphere^4 * alpha_sphere^2) /
      (sphere_r^2 * sphere_parameters$deviation_fun(ka_sphere))
  )
  expect_equal(
    acousticTS:::.stanton_hp(
      acoustics = acoustics,
      body = sphere_body,
      parameters = sphere_parameters,
      alpha = alpha_sphere
    ),
    sphere_expected
  )

  ka_prolate <- acoustics$k_sw * prolate_body$radius
  prolate_r <- (prolate_body$g * prolate_body$h - 1) /
    (prolate_body$g * prolate_body$h + 1)
  prolate_expected <- ((1 / 9) * prolate_body$length^2 * ka_prolate^4 *
    alpha_cylinder^2 * 0.8) / (
    1 + ((16 / 9) * ka_prolate^4 * alpha_cylinder^2) /
      (prolate_r^2 * 1.4)
  )
  expect_equal(
    acousticTS:::.stanton_hp(
      acoustics = acoustics,
      body = prolate_body,
      parameters = list(null_fun = 0.8, deviation_fun = 1.4),
      alpha = alpha_cylinder
    ),
    prolate_expected
  )

  straight_expected <- {
    ka <- acoustics$k_sw * cylinder_body$radius
    s <- sin(acoustics$k_sw * cylinder_body$length * cos(cylinder_body$theta)) /
      (acoustics$k_sw * cylinder_body$length * cos(cylinder_body$theta))
    Ka <- ka * sin(cylinder_body$theta)
    r <- (cylinder_body$g * cylinder_body$h - 1) /
      (cylinder_body$g * cylinder_body$h + 1)

    (0.25 * cylinder_body$length^2 * Ka^4 * alpha_cylinder^2 * s^2 * 0.75) /
      (1 + (pi * Ka^3 * alpha_cylinder^2) / (r^2 * 1.25))
  }
  expect_equal(
    acousticTS:::.stanton_cylinder(
      k = acoustics$k_sw,
      l = cylinder_body$length,
      a = cylinder_body$radius,
      theta = cylinder_body$theta,
      rho_c = numeric(),
      r = (cylinder_body$g * cylinder_body$h - 1) /
        (cylinder_body$g * cylinder_body$h + 1),
      alpha = alpha_cylinder,
      gnull = 0.75,
      fdevs = 1.25
    ),
    straight_expected
  )

  bent_expected <- {
    ka <- acoustics$k_sw * cylinder_body$radius
    rho_c <- 3
    r_c <- rho_c * cylinder_body$length
    H <- 0.5 + 0.5 * rho_c * sin(1 / rho_c)
    r <- (cylinder_body$g * cylinder_body$h - 1) /
      (cylinder_body$g * cylinder_body$h + 1)

    (0.25 * cylinder_body$length^2 * ka^4 * alpha_cylinder^2 * H^2 * 0.75) /
      (1 + (cylinder_body$length^2 * ka^4 * alpha_cylinder^2 * H^2) /
        (r_c * cylinder_body$radius * r^2 * 1.25))
  }
  expect_equal(
    acousticTS:::.stanton_cylinder(
      k = acoustics$k_sw,
      l = cylinder_body$length,
      a = cylinder_body$radius,
      theta = cylinder_body$theta,
      rho_c = 3,
      r = (cylinder_body$g * cylinder_body$h - 1) /
        (cylinder_body$g * cylinder_body$h + 1),
      alpha = alpha_cylinder,
      gnull = 0.75,
      fdevs = 1.25
    ),
    bent_expected
  )
})

test_that("hpa initialization and evaluation wire sphere, cylinder, and ESS objects", {
  gas_obj <- gas_generate(
    shape = sphere(radius_body = 0.01, n_segments = 40),
    density_fluid = 1.24,
    sound_speed_fluid = 345
  )
  johnson_initialized <- acousticTS:::hpa_initialize(
    gas_obj,
    frequency = c(38000, 120000),
    method = "johnson"
  )
  johnson_out <- acousticTS:::HPA(johnson_initialized)
  johnson_model <- johnson_out@model_parameters$HPA
  johnson_alpha <- acousticTS:::.alpha_pi(johnson_model$body)
  johnson_expected <- acousticTS:::.johnson_hp(
    acoustics = johnson_model$parameters$acoustics,
    body = johnson_model$body,
    alpha = johnson_alpha
  )

  expect_s4_class(johnson_out, "GAS")
  expect_equal(johnson_out@model$HPA$sigma_bs, johnson_expected)
  expect_equal(johnson_out@model$HPA$TS, 10 * log10(johnson_expected))

  fls_obj <- fls_generate(
    shape = cylinder(length_body = 0.02, radius_body = 0.001, n_segments = 20),
    g_body = 1.02,
    h_body = 1.03,
    theta_body = 1.1
  )
  stanton_initialized <- acousticTS:::hpa_initialize(
    fls_obj,
    frequency = 38000,
    method = "stanton",
    deviation_fun = function(ka) ka + 1,
    null_fun = function(ka) ka^2 + 1
  )
  stanton_out <- acousticTS:::HPA(stanton_initialized)
  stanton_model <- stanton_out@model_parameters$HPA
  stanton_alpha <- acousticTS:::.alpha_pi(stanton_model$body)
  stanton_expected <- acousticTS:::.stanton_hp(
    acoustics = stanton_model$parameters$acoustics,
    body = stanton_model$body,
    parameters = stanton_model$parameters,
    alpha = stanton_alpha
  )

  expect_s4_class(stanton_out, "FLS")
  expect_equal(stanton_out@model$HPA$sigma_bs, stanton_expected)
  expect_equal(stanton_out@model$HPA$TS, 10 * log10(stanton_expected))

  ess_obj <- ess_generate(
    shape = sphere(radius_body = 0.01, n_segments = 40),
    shell_thickness = 1e-3,
    sound_speed_shell = 3750,
    sound_speed_fluid = 1575,
    density_shell = 2565,
    density_fluid = 1077.3,
    K = 70e9,
    nu = 0.32
  )
  ess_initialized <- acousticTS:::hpa_initialize(
    ess_obj,
    frequency = 38000,
    method = "stanton"
  )

  expect_s4_class(ess_initialized, "ESS")
  expect_true("HPA" %in% names(ess_initialized@model_parameters))
  expect_identical(
    as.character(ess_initialized@model_parameters$HPA$body$shape),
    "Sphere"
  )
})
