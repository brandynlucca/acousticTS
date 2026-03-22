library(acousticTS)

test_that("pcdwba_initialize stores the expected model scaffolding", {
  obj <- fls_generate(
    shape = cylinder(
      length_body = 0.015,
      radius_body = 0.001,
      taper = 10,
      radius_curvature_ratio = 3,
      n_segments = 50
    ),
    g_body = 1.02,
    h_body = 1.02,
    theta_body = pi / 2
  )

  out <- pcdwba_initialize(
    object = obj,
    frequency = c(38000, 120000),
    sound_speed_sw = 1500,
    density_sw = 1026
  )

  expect_s4_class(out, "FLS")
  expect_true("PCDWBA" %in% names(out@model_parameters))
  expect_true("PCDWBA" %in% names(out@model))
  expect_equal(out@model$PCDWBA$frequency, c(38000, 120000))
})

test_that("pcdwba matches the bent-cylinder reference case", {
  obj <- fls_generate(
    shape = cylinder(
      length_body = 0.015,
      radius_body = 0.001,
      taper = 10,
      radius_curvature_ratio = 3,
      n_segments = 50
    ),
    g_body = 1.02,
    h_body = 1.02,
    theta_body = pi / 2
  )

  out_regular <- target_strength(
    object = obj,
    frequency = c(38000, 120000),
    model = "pcdwba",
    sound_speed_sw = 1500,
    density_sw = 1026
  )

  out_regular_df <- out_regular@model$PCDWBA

  expect_equal(
    Re(out_regular_df$f_bs),
    c(-2.63183454623927e-06, -3.72290353412917e-06),
    tolerance = 1e-14
  )
  expect_equal(
    Im(out_regular_df$f_bs),
    c(1.50552093122769e-09, -4.51576556149518e-09),
    tolerance = 1e-14
  )
  expect_equal(
    out_regular_df$TS,
    c(-111.59482691279, -108.582357947069),
    tolerance = 1e-11
  )

  obj_override <- fls_generate(
    shape = cylinder(
      length_body = 0.015,
      radius_body = 0.001,
      taper = 10,
      n_segments = 50
    ),
    g_body = 1.02,
    h_body = 1.02,
    theta_body = pi / 2
  )

  out_override <- target_strength(
    object = obj_override,
    frequency = c(38000, 120000),
    model = "pcdwba",
    sound_speed_sw = 1500,
    density_sw = 1026,
    radius_curvature_ratio = 3
  )

  expect_equal(
    out_override@model$PCDWBA$f_bs,
    out_regular_df$f_bs,
    tolerance = 1e-14
  )
  expect_equal(
    out_override@model$PCDWBA$TS,
    out_regular_df$TS,
    tolerance = 1e-12
  )
})
