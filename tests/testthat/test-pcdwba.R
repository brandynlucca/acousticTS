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
    tolerance = 1e-12
  )
  expect_equal(
    Im(out_regular_df$f_bs),
    c(1.50552093122769e-09, -4.51576556149518e-09),
    tolerance = 1e-12
  )
  expect_equal(
    out_regular_df$TS,
    c(-111.59482691279, -108.582357947069),
    tolerance = 1e-10
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
    tolerance = 1e-12
  )
  expect_equal(
    out_override@model$PCDWBA$TS,
    out_regular_df$TS,
    tolerance = 1e-10
  )
})

test_that("pcdwba helper functions normalize geometry inputs cleanly", {
  row_major <- rbind(
    x = c(0, 1, 2),
    z = c(0, 0, 0),
    a = c(0.2, 0.3, 0.2)
  )
  curved_body <- list(
    rpos = rbind(
      x = c(0, 1, 2),
      z = c(0, 0.5, 0),
      a = c(0.2, 0.3, 0.2)
    ),
    radius = c(0.2, 0.3, 0.2)
  )

  expect_equal(acousticTS:::.pcdwba_axis_row(row_major, c("z", "z_body")), 2L)
  expect_equal(
    acousticTS:::.pcdwba_axis_row(row_major, c("radius"), default = 1L),
    1L
  )
  expect_true(acousticTS:::.pcdwba_is_straight(list(rpos = row_major)))
  expect_false(acousticTS:::.pcdwba_is_straight(curved_body))

  expect_equal(
    acousticTS:::.pcdwba_node_property(1.2, n_nodes = 3, label = "h"),
    rep(1.2, 3)
  )
  expect_equal(
    acousticTS:::.pcdwba_node_property(c(1.1, 1.2, 1.3), 3, "g"),
    c(1.1, 1.2, 1.3)
  )
  expect_equal(
    acousticTS:::.pcdwba_node_property(c(1.2, 1.4), 3, "g"),
    c(1.2, 1.2, 1.4)
  )
  expect_error(
    acousticTS:::.pcdwba_node_property(1:4, 3, "g"),
    "requires 'g' to be scalar, nodewise, or segmentwise"
  )

  regular_geometry <- acousticTS:::.pcdwba_regular_geometry(
    n_nodes = 5,
    radius_curvature_ratio = 3,
    taper_order = 2
  )
  expect_equal(length(regular_geometry$taper), 5)
  expect_equal(regular_geometry$taper[c(1, 5)], c(0, 0))
  expect_equal(length(regular_geometry$theta_tilt), 5)

  profile_geometry <- acousticTS:::.pcdwba_profile_geometry(curved_body)
  expected_length <- 2 * sqrt(1^2 + 0.5^2)
  expect_equal(profile_geometry$length_body, expected_length)
  expect_equal(profile_geometry$taper, c(2 / 3, 1, 2 / 3))
  expect_equal(length(profile_geometry$theta_tilt), 3)
  expect_equal(length(profile_geometry$gamma_t), 3)
  expect_equal(length(profile_geometry$dr_norm), 3)
  expect_equal(length(profile_geometry$r_pos), 3)

  expect_error(
    acousticTS:::.pcdwba_profile_geometry(
      list(
        rpos = rbind(
          x = c(1, 1, 1),
          z = c(0, 0, 0),
          a = c(0.2, 0.3, 0.2)
        ),
        radius = c(0.2, 0.3, 0.2)
      )
    ),
    "could not resolve a positive body length"
  )
})

test_that("pcdwba preparation rejects incompatible inputs and curves straight arbitrary profiles", {
  expect_error(
    acousticTS:::.pcdwba_prepare_geometry(
      gas_generate(shape = sphere(radius_body = 0.01, n_segments = 40)),
      radius_curvature_ratio = 3
    ),
    "requires a fluid-like scatterer"
  )

  obj <- fls_generate(
    shape = arbitrary(
      x_body = c(0, 0.01, 0.02),
      z_body = c(0, 0, 0),
      radius_body = c(0.001, 0.0015, 0.001)
    ),
    density_body = 1050,
    sound_speed_body = 1505
  )

  expect_error(
    acousticTS:::.pcdwba_prepare_geometry(
      obj,
      radius_curvature = 0.03,
      radius_curvature_ratio = 3
    ),
    "Specify at most one"
  )

  prepared <- acousticTS:::.pcdwba_prepare_geometry(
    obj,
    radius_curvature_ratio = 4
  )
  bent_body <- acousticTS::extract(prepared$object, "body")

  expect_s4_class(prepared$object, "FLS")
  expect_true(prepared$geometry$length_body > 0)
  expect_equal(prepared$max_radius, 0.0015)
  expect_equal(bent_body$radius_curvature_ratio, 4)
  expect_true(!is.null(bent_body$arc_length))
})
