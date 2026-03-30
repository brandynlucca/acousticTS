library(acousticTS)

.with_temp_pdf_device <- function(expr) {
  path <- tempfile(fileext = ".pdf")
  grDevices::pdf(path)
  on.exit({
    grDevices::dev.off()
    unlink(path)
  }, add = TRUE)
  force(expr)
}

test_that("plotting helpers cover empty, missing-column, fallback-label, and overlay branches", {
  local({
    testthat::local_mocked_bindings(
      extract = function(object, what) list(),
      .package = "acousticTS"
    )
    expect_error(
      acousticTS:::.collect_model_plot_data(list()),
      "no model results detected"
    )
  })

  local({
    testthat::local_mocked_bindings(
      extract = function(object, what) {
        list(
          SPHMS = data.frame(frequency = 38e3, TS = -40),
          TMM = data.frame(frequency = 38e3, TS = -41)
        )
      },
      .package = "acousticTS"
    )
    expect_error(
      acousticTS:::.collect_model_plot_data(list(), x_units = "k_sw"),
      "does not contain the requested plotting column"
    )
  })

  local({
    testthat::local_mocked_bindings(
      extract = function(object, what) {
        list(
          SPHMS = data.frame(frequency = 38e3, custom_x = 1, custom_y = 2),
          TMM = data.frame(frequency = 70e3, custom_x = 3, custom_y = 4)
        )
      },
      .package = "acousticTS"
    )
    plot_data <- acousticTS:::.collect_model_plot_data(
      list(),
      x_units = "custom_x",
      y_units = "custom_y"
    )
    expect_equal(plot_data$x_lab, "custom_x")
    expect_equal(plot_data$y_lab, "custom_y")

    .with_temp_pdf_device({
      local({
        testthat::local_mocked_bindings(
          .has_active_multi_panel_layout = function() TRUE,
          .package = "acousticTS"
        )
        expect_no_error(acousticTS:::.plot_model_results(plot_data))
      })
    })
  })

  position_matrix <- cbind(x = c(0, 1, 2), radius = c(0.1, 0.2, 0.1))
  .with_temp_pdf_device({
    outline <- acousticTS:::.plot_axisymmetric_outline(position_matrix, init = TRUE)
    expect_equal(outline$radius, c(0.1, 0.2, 0.1), tolerance = 1e-12)
    expect_no_error(acousticTS:::.plot_axisymmetric_outline(position_matrix, init = FALSE))
  })

  row_major_z_bounds <- rbind(
    x = c(0, 1, 2),
    zU = c(0.2, 0.3, 0.2),
    zL = c(-0.2, -0.3, -0.2)
  )
  row_major_flat <- rbind(
    x = c(0, 1, 2),
    y = c(0, 0, 0)
  )
  .with_temp_pdf_device({
    local({
      testthat::local_mocked_bindings(
        .has_active_multi_panel_layout = function() TRUE,
        .package = "acousticTS"
      )
      expect_no_error(acousticTS:::.plot_axisymmetric_outline(position_matrix, init = TRUE))
      expect_no_error(acousticTS:::.plot_row_major_segmented_body(
        row_major_z_bounds,
        aspect_ratio = "equal"
      ))
      expect_no_error(acousticTS:::.plot_row_major_segmented_body(row_major_flat))
    })
  })
})

test_that("SPHMS helper branches cover default boundaries, validation, and parameter fallbacks", {
  expect_no_error(acousticTS:::.sphms_validate_shape(list(shape = "Sphere")))
  expect_error(
    acousticTS:::.sphms_validate_shape(list(shape = "Cylinder")),
    "requires scatterer to be shape-type 'Sphere'"
  )

  cal_obj <- cal_generate()
  ess_obj <- fixture_sphere("shelled_liquid")
  fls_obj <- fixture_sphere("liquid_filled")
  gas_obj <- gas_generate(
    shape = sphere(radius_body = 0.01, n_segments = 20),
    density_fluid = 1.24,
    sound_speed_fluid = 345
  )

  expect_error(
    acousticTS:::.sphms_default_boundary(cal_obj, NULL),
    "Use 'model = \"calibration\"'"
  )
  expect_equal(acousticTS:::.sphms_default_boundary(ess_obj, NULL), "shelled_liquid")
  expect_equal(acousticTS:::.sphms_default_boundary(fls_obj, NULL), "liquid_filled")
  expect_equal(acousticTS:::.sphms_default_boundary(gas_obj, NULL), "gas_filled")
  expect_error(
    acousticTS:::.sphms_default_boundary(structure(list(), class = "mystery"), NULL),
    "Could not infer a default SPHMS boundary"
  )

  expect_equal(acousticTS:::.sphms_validate_boundary(fls_obj, "liquid_filled"), "liquid_filled")
  expect_error(
    acousticTS:::.sphms_validate_boundary(ess_obj, "bad_boundary"),
    "'ESS'-class"
  )
  expect_error(
    acousticTS:::.sphms_validate_boundary(structure(list(), class = "mystery"), "mystery"),
    "Only the following values for 'boundary' are available"
  )

  exterior <- list(density = 1100, sound_speed = 1500, radius = 0.01)
  body_params_fls <- acousticTS:::.sphms_body_parameters(fls_obj, exterior, 1480, 1026)
  expect_equal(body_params_fls$g31, 1100 / 1026, tolerance = 1e-12)
  expect_equal(body_params_fls$h31, 1500 / 1480, tolerance = 1e-12)

  body_params_ess <- acousticTS:::.sphms_body_parameters(
    ess_obj,
    extract(ess_obj, "shell"),
    1480,
    1026
  )
  expect_true(is.finite(body_params_ess$g21))
  expect_true(is.finite(body_params_ess$h21))
  expect_true(is.finite(body_params_ess$g31))
  expect_true(is.finite(body_params_ess$h31))
  local({
    testthat::local_mocked_bindings(
      extract = function(object, what) {
        if (what == "fluid") {
          list(g = 0.9, h = 0.8, radius = 0.005)
        } else {
          stop("unexpected extract target")
        }
      },
      .package = "acousticTS"
    )
    body_params_ess_fallback <- acousticTS:::.sphms_body_parameters(
      ess_obj,
      list(g = 1.1, h = 1.2, radius = 0.01),
      1480,
      1026
    )
    expect_equal(body_params_ess_fallback$h32, 0.8 / 1.2, tolerance = 1e-12)
  })

  acoustics <- acousticTS:::.sphms_acoustics(38e3, 1480, body_params_fls)
  expect_true(all(c("k_sw", "k_shell", "k_fluid") %in% names(acoustics)))
  expect_equal(acousticTS:::.sphms_m_limit(11L, acoustics, body_params_fls), 11L)
  expect_true(acousticTS:::.sphms_m_limit(NULL, acoustics, body_params_ess) > 0)
  expect_true(acousticTS:::.sphms_m_limit(NULL, acoustics, body_params_fls) > 0)
})

test_that("PCDWBA helper branches cover axis resolution, profile handling, and validation", {
  straight_body <- list(
    rpos = rbind(
      x = c(0, 1, 2),
      y = c(0, 0, 0)
    )
  )
  curved_body <- list(
    rpos = rbind(
      x = c(0, 1, 2),
      z = c(0, 0.1, 0.2)
    ),
    radius = c(0.1, 0.2, 0.1)
  )

  expect_equal(acousticTS:::.pcdwba_axis_row(straight_body$rpos, c("x", "x_body")), 1L)
  expect_null(acousticTS:::.pcdwba_axis_row(straight_body$rpos, c("z", "z_body")))
  expect_equal(acousticTS:::.pcdwba_axis_row(straight_body$rpos, c("z", "z_body"), default = 9L), 9L)

  expect_true(acousticTS:::.pcdwba_is_straight(straight_body))
  expect_false(acousticTS:::.pcdwba_is_straight(curved_body))

  expect_equal(acousticTS:::.pcdwba_node_property(2, 3, "density"), c(2, 2, 2))
  expect_equal(acousticTS:::.pcdwba_node_property(c(1, 2, 3), 3, "density"), c(1, 2, 3))
  expect_equal(acousticTS:::.pcdwba_node_property(c(2, 3), 3, "density"), c(2, 2, 3))
  expect_error(
    acousticTS:::.pcdwba_node_property(1:4, 3, "density"),
    "requires 'density' to be scalar, nodewise, or segmentwise"
  )

  regular_default_taper <- acousticTS:::.pcdwba_regular_geometry(5, radius_curvature_ratio = 2)
  expect_equal(regular_default_taper$taper, rep(1, 5), tolerance = 1e-12)
  regular <- acousticTS:::.pcdwba_regular_geometry(5, radius_curvature_ratio = 2, taper_order = 2)
  expect_equal(length(regular$taper), 5)

  profile_straight <- acousticTS:::.pcdwba_profile_geometry(straight_body)
  expect_true(profile_straight$length_body > 0)
  profile <- acousticTS:::.pcdwba_profile_geometry(curved_body)
  expect_true(profile$length_body > 0)
  expect_error(
    acousticTS:::.pcdwba_profile_geometry(list(rpos = rbind(x = c(1, 1, 1), z = c(0, 0, 0)), radius = c(0.1, 0.1, 0.1))),
    "could not resolve a positive body length"
  )

  expect_error(
    acousticTS:::.pcdwba_prepare_geometry(gas_generate(shape = sphere(radius_body = 0.01, n_segments = 20), density_fluid = 1.24, sound_speed_fluid = 345)),
    "requires a fluid-like scatterer"
  )
  expect_error(
    acousticTS:::.pcdwba_prepare_geometry(fixture_sphere("liquid_filled"), radius_curvature = 1, radius_curvature_ratio = 2),
    "Specify at most one of 'radius_curvature' or 'radius_curvature_ratio'"
  )

  local({
    testthat::local_mocked_bindings(
      .pcdwba_prepare_geometry = function(object, radius_curvature = NULL, radius_curvature_ratio = NULL) {
        list(
          object = fixture_sphere("liquid_filled"),
          geometry = list()
        )
      },
      .derive_contrasts = function(body, sound_speed_sw, density_sw) {
        list(h = NA_real_, g = NA_real_)
      },
      .package = "acousticTS"
    )
    expect_error(
      acousticTS:::pcdwba_initialize(fixture_sphere("liquid_filled"), 38e3),
      "PCDWBA requires finite density and sound-speed contrasts."
    )
  })
})
