library(acousticTS)

make_vesms_reference_object <- function() {
  radius_gas <- 1e-3
  radius_shell <- radius_gas + 0.02e-3
  shear_shell <- 0.2e6
  lambda_shell <- 2.4e9
  bulk_shell <- lambda_shell + 2 * shear_shell / 3

  ess_generate(
    shape = sphere(radius_body = radius_shell, n_segments = 80),
    radius_shell = radius_shell,
    shell_thickness = radius_shell - radius_gas,
    density_shell = 1040,
    density_fluid = 80,
    sound_speed_fluid = 325,
    G = shear_shell,
    K = bulk_shell
  )
}

test_that("vesms_initialize stores the viscous-layer radius and scaffolding", {
  radius_gas <- 1e-3
  radius_shell <- radius_gas + 0.02e-3
  shear_shell <- 0.2e6
  lambda_shell <- 2.4e9
  bulk_shell <- lambda_shell + 2 * shear_shell / 3

  obj <- ess_generate(
    shape = sphere(radius_body = radius_shell, n_segments = 80),
    radius_shell = radius_shell,
    shell_thickness = radius_shell - radius_gas,
    density_shell = 1040,
    density_fluid = 80,
    sound_speed_fluid = 325,
    G = shear_shell,
    K = bulk_shell
  )

  out <- vesms_initialize(
    object = obj,
    frequency = c(1000, 38000, 120000),
    sound_speed_sw = 1500,
    density_sw = 1027,
    sound_speed_viscous = 1510,
    density_viscous = 1040,
    shear_viscosity_viscous = 3,
    bulk_viscosity_viscous = 3,
    m_limit = 2
  )

  expect_s4_class(out, "ESS")
  expect_true("VESMS" %in% names(out@model_parameters))
  expect_true("VESMS" %in% names(out@model))
  expect_equal(
    out@model_parameters$VESMS$viscous$radius,
    0.00419542498515355,
    tolerance = 1e-15
  )
})

test_that("VESMS reproduces the local reference implementation case", {
  radius_gas <- 1e-3
  radius_shell <- radius_gas + 0.02e-3
  shear_shell <- 0.2e6
  lambda_shell <- 2.4e9
  bulk_shell <- lambda_shell + 2 * shear_shell / 3

  obj <- ess_generate(
    shape = sphere(radius_body = radius_shell, n_segments = 80),
    radius_shell = radius_shell,
    shell_thickness = radius_shell - radius_gas,
    density_shell = 1040,
    density_fluid = 80,
    sound_speed_fluid = 325,
    G = shear_shell,
    K = bulk_shell
  )

  out <- target_strength(
    object = obj,
    frequency = c(1000, 38000, 120000),
    model = "vesms",
    sound_speed_sw = 1500,
    density_sw = 1027,
    sound_speed_viscous = 1510,
    density_viscous = 1040,
    shear_viscosity_viscous = 3,
    bulk_viscosity_viscous = 3,
    m_limit = 2
  )

  expect_equal(
    out@model$VESMS$TS,
    c(-116.009657062216, -55.8090704683511, -65.2551347115176),
    tolerance = 0.03
  )
})

test_that("vesms_initialize respects explicit viscous-layer controls", {
  radius_gas <- 1e-3
  radius_shell <- radius_gas + 0.02e-3
  shear_shell <- 0.2e6
  lambda_shell <- 2.4e9
  bulk_shell <- lambda_shell + 2 * shear_shell / 3

  obj <- ess_generate(
    shape = sphere(radius_body = radius_shell, n_segments = 80),
    radius_shell = radius_shell,
    shell_thickness = radius_shell - radius_gas,
    density_shell = 1040,
    density_fluid = 80,
    sound_speed_fluid = 325,
    G = shear_shell,
    K = bulk_shell
  )

  out_radius <- vesms_initialize(
    object = obj,
    frequency = c(1000, 38000, 120000),
    sound_speed_sw = 1500,
    density_sw = 1027,
    sound_speed_viscous = 1510,
    density_viscous = 1040,
    shear_viscosity_viscous = 3,
    radius_viscous = 0.0015,
    m_limit = 4
  )

  expect_equal(out_radius@model_parameters$VESMS$viscous$radius, 0.0015)
  expect_equal(out_radius@model_parameters$VESMS$viscous$bulk_viscosity, 3)
  expect_equal(out_radius@model_parameters$VESMS$parameters$acoustics$m_limit, c(4L, 4L, 4L))

  out_thickness <- vesms_initialize(
    object = obj,
    frequency = c(1000, 38000, 120000),
    sound_speed_sw = 1500,
    density_sw = 1027,
    sound_speed_viscous = 1510,
    density_viscous = 1040,
    shear_viscosity_viscous = 3,
    viscous_thickness = 0.0015 - radius_shell,
    bulk_viscosity_viscous = 5,
    m_limit = c(2, 3, 4)
  )

  expect_equal(out_thickness@model_parameters$VESMS$viscous$radius, 0.0015)
  expect_equal(out_thickness@model_parameters$VESMS$viscous$bulk_viscosity, 5)
  expect_equal(out_thickness@model_parameters$VESMS$parameters$acoustics$m_limit, c(2L, 3L, 4L))
})

test_that("vesms_initialize rejects inconsistent geometry and viscous-layer inputs", {
  radius_gas <- 1e-3
  radius_shell <- radius_gas + 0.02e-3
  shear_shell <- 0.2e6
  lambda_shell <- 2.4e9
  bulk_shell <- lambda_shell + 2 * shear_shell / 3

  sphere_obj <- ess_generate(
    shape = sphere(radius_body = radius_shell, n_segments = 80),
    radius_shell = radius_shell,
    shell_thickness = radius_shell - radius_gas,
    density_shell = 1040,
    density_fluid = 80,
    sound_speed_fluid = 325,
    G = shear_shell,
    K = bulk_shell
  )

  expect_error(
    vesms_initialize(
      object = sphere_obj,
      frequency = c(1000, 38000),
      sound_speed_sw = 1500,
      density_sw = 1027,
      sound_speed_viscous = 1510,
      density_viscous = 1040,
      shear_viscosity_viscous = 3,
      radius_viscous = 0.0015,
      viscous_thickness = 0.0004
    ),
    "Specify at most one"
  )

  expect_error(
    vesms_initialize(
      object = sphere_obj,
      frequency = c(1000, 38000),
      sound_speed_sw = 1500,
      density_sw = 1027,
      sound_speed_viscous = 1510,
      density_viscous = 1027,
      shear_viscosity_viscous = 3
    ),
    "Neutral-buoyancy estimation"
  )

  nonspherical_obj <- sphere_obj
  nonspherical_obj@shape_parameters$shape <- "Cylinder"

  expect_error(
    vesms_initialize(
      object = nonspherical_obj,
      frequency = c(1000, 38000),
      sound_speed_sw = 1500,
      density_sw = 1027,
      sound_speed_viscous = 1510,
      density_viscous = 1040,
      shear_viscosity_viscous = 3,
      radius_viscous = 0.0015
    ),
    "spherical ESS geometry"
  )
})

test_that("VESMS enforces supported object type and layered-property requirements", {
  sphere_obj <- make_vesms_reference_object()

  expect_error(
    target_strength(
      fls_generate(
        shape = sphere(radius_body = 1e-3, n_segments = 80),
        density_body = 1040,
        sound_speed_body = 1500
      ),
      frequency = c(1000, 38000),
      model = "vesms",
      sound_speed_sw = 1500,
      density_sw = 1027,
      sound_speed_viscous = 1510,
      density_viscous = 1040,
      shear_viscosity_viscous = 3,
      radius_viscous = 0.0015
    ),
    "requires an elastic-shelled scatterer"
  )

  missing_shell_props <- sphere_obj
  missing_shell_props@shell$G <- NULL
  expect_error(
    vesms_initialize(
      object = missing_shell_props,
      frequency = c(1000, 38000),
      sound_speed_sw = 1500,
      density_sw = 1027,
      sound_speed_viscous = 1510,
      density_viscous = 1040,
      shear_viscosity_viscous = 3,
      radius_viscous = 0.0015
    ),
    "shell to have finite density and shear modulus"
  )

  missing_lambda <- sphere_obj
  missing_lambda@shell$K <- NULL
  missing_lambda@shell$E <- NULL
  missing_lambda@shell$nu <- NULL
  missing_lambda@shell$lambda <- NULL
  expect_error(
    vesms_initialize(
      object = missing_lambda,
      frequency = c(1000, 38000),
      sound_speed_sw = 1500,
      density_sw = 1027,
      sound_speed_viscous = 1510,
      density_viscous = 1040,
      shear_viscosity_viscous = 3,
      radius_viscous = 0.0015
    ),
    "requires Lam"
  )

  missing_gas_props <- sphere_obj
  missing_gas_props@fluid$sound_speed <- NULL
  missing_gas_props@fluid$h <- NULL
  expect_error(
    vesms_initialize(
      object = missing_gas_props,
      frequency = c(1000, 38000),
      sound_speed_sw = 1500,
      density_sw = 1027,
      sound_speed_viscous = 1510,
      density_viscous = 1040,
      shear_viscosity_viscous = 3,
      radius_viscous = 0.0015
    ),
    "inner fluid slot to represent a gas core"
  )
})

test_that("VESMS validates viscous-layer arguments and modal caps", {
  sphere_obj <- make_vesms_reference_object()

  expect_error(
    vesms_initialize(
      object = sphere_obj,
      frequency = c(1000, 38000),
      sound_speed_sw = 1500,
      density_sw = 1027,
      sound_speed_viscous = 1510,
      density_viscous = 1040,
      shear_viscosity_viscous = 3,
      radius_viscous = 0.0010
    ),
    "radius to be a finite scalar larger than the shell radius"
  )

  expect_error(
    vesms_initialize(
      object = sphere_obj,
      frequency = c(1000, 38000),
      sound_speed_sw = 1500,
      density_sw = 1027,
      sound_speed_viscous = -1510,
      density_viscous = 1040,
      shear_viscosity_viscous = 3,
      radius_viscous = 0.0015
    ),
    "requires positive finite viscous-layer sound speed"
  )

  expect_error(
    vesms_initialize(
      object = sphere_obj,
      frequency = c(1000, 38000),
      sound_speed_sw = 1500,
      density_sw = 1027,
      sound_speed_viscous = 1510,
      density_viscous = 1040,
      shear_viscosity_viscous = 3,
      bulk_viscosity_viscous = -1,
      radius_viscous = 0.0015
    ),
    "finite non-negative bulk viscosity"
  )

  expect_error(
    vesms_initialize(
      object = sphere_obj,
      frequency = c(1000, 38000, 120000),
      sound_speed_sw = 1500,
      density_sw = 1027,
      sound_speed_viscous = 1510,
      density_viscous = 1040,
      shear_viscosity_viscous = 3,
      radius_viscous = 0.0015,
      m_limit = c(2, 3)
    ),
    "non-negative scalar or a vector with one value per frequency"
  )
})

test_that("VESMS gives consistent end-to-end results for equivalent outer-radius inputs", {
  sphere_obj <- make_vesms_reference_object()
  radius_shell <- sphere_obj@shell$radius
  radius_viscous <- 0.0015
  frequency <- c(1000, 38000, 120000)

  out_radius <- target_strength(
    object = sphere_obj,
    frequency = frequency,
    model = "vesms",
    sound_speed_sw = 1500,
    density_sw = 1027,
    sound_speed_viscous = 1510,
    density_viscous = 1040,
    shear_viscosity_viscous = 3,
    bulk_viscosity_viscous = 5,
    radius_viscous = radius_viscous,
    m_limit = c(2, 3, 4)
  )

  out_thickness <- target_strength(
    object = sphere_obj,
    frequency = frequency,
    model = "vesms",
    sound_speed_sw = 1500,
    density_sw = 1027,
    sound_speed_viscous = 1510,
    density_viscous = 1040,
    shear_viscosity_viscous = 3,
    bulk_viscosity_viscous = 5,
    viscous_thickness = radius_viscous - radius_shell,
    m_limit = c(2, 3, 4)
  )

  expect_true(all(c(
    "ka_viscous", "ka_shell", "ka_gas", "f_bs", "sigma_bs", "TS"
  ) %in% names(out_radius@model$VESMS)))
  expect_equal(
    out_radius@model_parameters$VESMS$viscous$radius,
    out_thickness@model_parameters$VESMS$viscous$radius,
    tolerance = 1e-12
  )
  expect_equal(out_radius@model$VESMS$f_bs, out_thickness@model$VESMS$f_bs, tolerance = 1e-12)
  expect_equal(out_radius@model$VESMS$TS, out_thickness@model$VESMS$TS, tolerance = 1e-12)
})
