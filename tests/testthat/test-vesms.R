library(acousticTS)

test_that("vesms_initialize stores the viscous-layer radius and scaffolding", {
  radius_gas <- 1e-3
  radius_shell <- radius_gas + 0.02e-3
  shear_shell <- 0.2e6
  lambda_shell <- 2.4e9
  bulk_shell <- lambda_shell + 2 * shear_shell / 3

  obj <- ess_generate(
    shape = "sphere",
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
    shape = "sphere",
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
