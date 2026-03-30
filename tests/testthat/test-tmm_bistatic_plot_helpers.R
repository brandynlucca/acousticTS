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

test_that("TMM geometry helpers cover coordinate conversions and lobe widths", {
  expect_equal(
    acousticTS:::.tmm_spherical_to_cartesian(pi / 2, 0),
    c(1, 0, 0),
    tolerance = 1e-12
  )

  roundtrip <- acousticTS:::.tmm_cartesian_to_spherical(c(0, -2, 0))
  expect_equal(roundtrip[["theta"]], pi / 2, tolerance = 1e-12)
  expect_equal(roundtrip[["phi"]], 3 * pi / 2, tolerance = 1e-12)

  expect_equal(
    acousticTS:::.tmm_cross_product(c(1, 0, 0), c(0, 1, 0)),
    c(0, 0, 1),
    tolerance = 1e-12
  )

  basis <- acousticTS:::.tmm_forward_basis(theta_body = 0, phi_body = 0)
  expect_equal(sum(basis$forward * basis$e1), 0, tolerance = 1e-12)
  expect_equal(sum(basis$forward * basis$e2), 0, tolerance = 1e-12)
  expect_equal(sum(basis$e1 * basis$e2), 0, tolerance = 1e-12)

  slice_dirs <- acousticTS:::.tmm_local_slice_directions(
    theta_body = pi / 2,
    phi_body = 0,
    psi_scatter = c(0, pi),
    alpha = 0
  )
  expect_equal(slice_dirs$theta_scatter[1], pi / 2, tolerance = 1e-12)
  expect_equal(slice_dirs$phi_scatter[1], 0, tolerance = 1e-12)
  expect_equal(slice_dirs$theta_scatter[2], pi / 2, tolerance = 1e-12)
  expect_equal(slice_dirs$phi_scatter[2], pi, tolerance = 1e-12)

  psi_grid <- acousticTS:::.tmm_forward_separation_matrix(
    theta_body = pi / 2,
    phi_body = 0,
    theta_scatter = c(pi / 2, pi / 2),
    phi_scatter = c(0, pi)
  )
  expect_equal(psi_grid[1, 1], 0, tolerance = 1e-12)
  expect_equal(psi_grid[1, 2], pi, tolerance = 1e-12)

  expect_error(
    acousticTS:::.tmm_lobe_width(
      psi_scatter = c(0, 1, 2),
      sigma_scat_dB = c(-10, -5, 0),
      drop_dB = 0
    ),
    "'drop_dB' must be a single positive numeric value."
  )
  expect_true(is.na(acousticTS:::.tmm_lobe_width(
    psi_scatter = c(0, 1, 2),
    sigma_scat_dB = c(-10, NA_real_, 0),
    center_psi = 1
  )))
  expect_equal(
    acousticTS:::.tmm_lobe_width(
      psi_scatter = c(0, 1, 2, 3),
      sigma_scat_dB = c(-10, -1, 0, -1),
      center_psi = 2,
      drop_dB = 3
    ),
    2,
    tolerance = 1e-12
  )
})

test_that("TMM sector and summary-input helpers validate their contracts", {
  defaults <- acousticTS:::.tmm_default_sectors()
  expect_equal(nrow(defaults), 3)
  expect_equal(
    acousticTS:::.tmm_validate_sectors(NULL)$sector,
    defaults$sector
  )

  expect_error(
    acousticTS:::.tmm_validate_sectors(data.frame(sector = "x", psi_min = 0)),
    "Missing: psi_max"
  )
  expect_error(
    acousticTS:::.tmm_validate_sectors(
      data.frame(sector = "", psi_min = 0, psi_max = 1)
    ),
    "'sectors\\$sector' must contain non-empty names."
  )
  expect_error(
    acousticTS:::.tmm_validate_sectors(
      data.frame(sector = "bad", psi_min = 1, psi_max = 1)
    ),
    "psi_min < psi_max <= pi"
  )

  custom_sectors <- data.frame(
    sector = c("forward", "backward"),
    psi_min = c(0, pi / 2),
    psi_max = c(pi / 2, pi),
    stringsAsFactors = FALSE
  )
  expect_equal(
    acousticTS:::.tmm_validate_sectors(custom_sectors),
    custom_sectors
  )

  solid_angle <- acousticTS:::.tmm_grid_solid_angle(
    theta_scatter = c(pi / 4, 3 * pi / 4),
    phi_scatter = c(pi / 2, 3 * pi / 2)
  )
  expect_equal(dim(solid_angle), c(2L, 2L))
  expect_true(all(solid_angle > 0))

  spherical_params <- list(
    parameters = list(
      coordinate_system = "spherical",
      acoustics = list(frequency = c(38e3, 70e3))
    ),
    body = list(theta_body = pi / 3, phi_body = pi / 4)
  )
  summary_inputs <- acousticTS:::.tmm_bistatic_summary_inputs(
    model_params = spherical_params,
    frequency = 70e3,
    theta_body = NULL,
    phi_body = NULL,
    n_psi = 7,
    sectors = custom_sectors
  )
  expect_equal(summary_inputs$idx, 2)
  expect_equal(summary_inputs$theta_body, pi / 3, tolerance = 1e-12)
  expect_equal(summary_inputs$phi_body, pi / 4, tolerance = 1e-12)
  expect_equal(summary_inputs$sectors, custom_sectors)

  expect_error(
    acousticTS:::.tmm_bistatic_summary_inputs(
      model_params = list(parameters = list(coordinate_system = "cylindrical")),
      frequency = NULL,
      theta_body = NULL,
      phi_body = NULL,
      n_psi = 7,
      sectors = NULL
    ),
    "Stored cylindrical TMM bistatic summaries are not available yet"
  )
  expect_error(
    acousticTS:::.tmm_bistatic_summary_inputs(
      model_params = spherical_params,
      frequency = 70e3,
      theta_body = NULL,
      phi_body = NULL,
      n_psi = 2,
      sectors = NULL
    ),
    "'n_psi' must be a single integer >= 3."
  )

  mask_open <- acousticTS:::.tmm_bistatic_sector_mask(
    psi_grid = matrix(c(0.1, 0.5, 2.0, 2.5), nrow = 2),
    sectors = custom_sectors,
    i = 1
  )
  mask_closed <- acousticTS:::.tmm_bistatic_sector_mask(
    psi_grid = matrix(c(0.1, 0.5, 2.0, 2.5), nrow = 2),
    sectors = custom_sectors,
    i = 2
  )
  expect_true(mask_open[1, 1])
  expect_true(mask_closed[2, 2])

  sector_integrals <- acousticTS:::.tmm_bistatic_sector_integrals(
    sectors = custom_sectors,
    psi_grid = matrix(c(0.1, 0.5, 2.0, 2.5), nrow = 2),
    sigma_scat = matrix(c(1, 2, 3, 4), nrow = 2),
    solid_angle = matrix(1, nrow = 2, ncol = 2)
  )
  expect_equal(sector_integrals$integrated_sigma_scat, c(3, 7))
})

test_that("TMM bistatic summary and product helpers assemble mocked sub-results", {
  object <- fls_generate(
    shape = sphere(radius_body = 0.01, n_segments = 20),
    density_body = 1028.9,
    sound_speed_body = 1480.3
  )
  sectors <- data.frame(
    sector = "all",
    psi_min = 0,
    psi_max = pi,
    stringsAsFactors = FALSE
  )
  grid <- list(
    frequency = 70e3,
    theta_scatter = c(0, pi / 2),
    phi_scatter = c(0, pi),
    sigma_scat = matrix(c(1, 2, 3, 4), nrow = 2)
  )

  local({
    testthat::local_mocked_bindings(
      .tmm_require_stored_blocks = function(object) {
        list(
          parameters = list(
            acoustics = list(frequency = c(38e3, 70e3)),
            coordinate_system = "spherical"
          ),
          body = list(theta_body = pi / 2, phi_body = pi)
        )
      },
      .tmm_warn_exploratory_cylinder_blocks = function(...) NULL,
      extract = function(object, what) list(shape = "mock"),
      .tmm_bistatic_summary_inputs = function(...) {
        list(idx = 2, theta_body = pi / 2, phi_body = pi, sectors = sectors)
      },
      .tmm_bistatic_summary_grid = function(...) {
        list(
          grid = grid,
          psi_grid = matrix(c(0, 1, 2, 3), nrow = 2),
          solid_angle = matrix(1, nrow = 2, ncol = 2)
        )
      },
      .tmm_bistatic_summary_slices = function(...) {
        list(
          psi_scatter = c(0, pi),
          forward_slice = data.frame(
            psi_scatter = c(0, pi),
            sigma_scat_dB = c(-8, -2)
          ),
          dorsal_slice = data.frame(
            psi_scatter = c(0, pi),
            sigma_scat_dB = c(-6, -3)
          )
        )
      },
      .tmm_bistatic_summary_points = function(...) {
        list(
          peak = list(theta = 0.4, phi = 0.5, sigma = 4, psi = 3),
          forward_point = list(
            sigma_scat = c(1, 2),
            sigma_scat_dB = c(0, 3)
          ),
          backscatter_point = list(sigma_scat = c(5, 6))
        )
      },
      .tmm_bistatic_sector_integrals = function(...) {
        data.frame(
          sector = "all",
          psi_min = 0,
          psi_max = pi,
          integrated_sigma_scat = 9
        )
      },
      .tmm_lobe_width = function(...) 0.25,
      .package = "acousticTS"
    )

    summary <- tmm_bistatic_summary(
      object,
      frequency = 70e3,
      n_theta = 5,
      n_phi = 7,
      n_psi = 9,
      include_grid = TRUE
    )

    expect_equal(summary$metrics$frequency, 70e3)
    expect_equal(summary$metrics$forward_sigma_scat, 2)
    expect_equal(summary$metrics$sigma_bs, 6)
    expect_equal(summary$metrics$peak_sigma_scat, 4)
    expect_equal(summary$metrics$backscatter_lobe_width, 0.25)
    expect_equal(summary$sector_integrals$integrated_sigma_scat, 9)
    expect_equal(summary$grid, grid)
  })

  local({
    testthat::local_mocked_bindings(
      .tmm_require_stored_blocks = function(object) {
        list(body = list(theta_body = pi / 4, phi_body = pi / 2))
      },
      .tmm_warn_exploratory_cylinder_blocks = function(...) NULL,
      .tmm_scalar_angle = function(value, default, name) {
        if (is.null(value)) default else value
      },
      tmm_scattering = function(...) list(f_scat = 1 + 1i, sigma_scat = 2),
      tmm_average_orientation = function(...) list(sigma_bs = 3),
      tmm_bistatic_summary = function(...) list(metrics = data.frame(flag = 1)),
      .package = "acousticTS"
    )

    products <- tmm_products(
      object,
      frequency = 70e3,
      orientation = structure(list(), class = "orientation_distribution"),
      bistatic_summary = TRUE,
      include_grid = TRUE,
      n_theta = 5,
      n_phi = 7,
      n_psi = 9
    )

    expect_true(all(c(
      "monostatic",
      "orientation_average",
      "bistatic_summary"
    ) %in% names(products)))
    expect_equal(products$monostatic$f_scat, 1 + 1i)
    expect_equal(products$orientation_average$sigma_bs, 3)
    expect_equal(products$bistatic_summary$metrics$flag, 1)
  })
})

test_that("TMM plotting helpers normalize quantities and build slice/grid data", {
  expect_equal(acousticTS:::.tmm_plot_quantity(), "sigma_scat_dB")
  expect_equal(acousticTS:::.tmm_plot_quantity(NULL), "sigma_scat_dB")
  expect_equal(acousticTS:::.tmm_plot_quantity("level_dB"), "sigma_scat_dB")
  expect_error(
    acousticTS:::.tmm_plot_quantity("bad"),
    "'arg' should be one of"
  )

  grid <- list(
    sigma_scat_dB = matrix(c(-10, -9, -8, -7), nrow = 2),
    sigma_scat = matrix(c(1, 2, 3, 4), nrow = 2),
    f_scat = matrix(complex(real = 1:4, imaginary = 4:1), nrow = 2)
  )
  expect_equal(acousticTS:::.tmm_grid_plot_data(grid)$quantity, "sigma_scat_dB")
  expect_equal(
    acousticTS:::.tmm_grid_plot_data(grid, "sigma_scat")$z,
    grid$sigma_scat
  )
  expect_equal(
    acousticTS:::.tmm_grid_plot_data(grid, "mod_f")$z,
    Mod(grid$f_scat)
  )
  expect_equal(
    acousticTS:::.tmm_grid_plot_data(grid, "phase")$z,
    Arg(grid$f_scat)
  )

  object <- fls_generate(
    shape = sphere(radius_body = 0.01, n_segments = 20),
    density_body = 1028.9,
    sound_speed_body = 1480.3
  )
  local({
    testthat::local_mocked_bindings(
      .tmm_require_stored_blocks = function(object) {
        list(
          parameters = list(acoustics = list(frequency = c(38e3, 70e3))),
          body = list(theta_body = pi / 3, phi_body = pi / 4)
        )
      },
      extract = function(object, what) list(shape = "mock"),
      .tmm_plot_frequency_index = function(...) 2,
      .tmm_scalar_angle = function(value, default, name) {
        if (is.null(value)) default else value
      },
      .tmm_scattering_points = function(model_params,
                                        frequency_idx,
                                        shape_parameters,
                                        theta_body,
                                        phi_body,
                                        theta_scatter,
                                        phi_scatter) {
        complex(
          real = seq_along(theta_body),
          imaginary = seq_along(theta_body) / 10
        )
      },
      .package = "acousticTS"
    )

    expect_error(
      acousticTS:::.tmm_scattering_slice_data(object, n_points = 1),
      "'n_points' must be a single integer >= 2."
    )

    theta_body_slice <- acousticTS:::.tmm_scattering_slice_data(
      object,
      frequency = 70e3,
      vary = "theta_body",
      quantity = "mod_f",
      n_points = 5
    )
    expect_equal(theta_body_slice$frequency, 70e3)
    expect_equal(theta_body_slice$vary, "theta_body")
    expect_equal(length(theta_body_slice$x), 5)

    phi_body_slice <- acousticTS:::.tmm_scattering_slice_data(
      object,
      frequency = 70e3,
      vary = "phi_body",
      quantity = "phase",
      n_points = 4,
      theta_scatter = pi / 3
    )
    expect_equal(phi_body_slice$vary, "phi_body")
    expect_equal(length(phi_body_slice$y), 4)

    theta_scatter_slice <- acousticTS:::.tmm_scattering_slice_data(
      object,
      frequency = 70e3,
      vary = "theta_scatter",
      quantity = "sigma_scat",
      n_points = 6,
      phi_scatter = 1.5 * pi
    )
    expect_equal(theta_scatter_slice$quantity, "sigma_scat")
    expect_equal(length(theta_scatter_slice$y), 6)

    phi_scatter_slice <- acousticTS:::.tmm_scattering_slice_data(
      object,
      frequency = 70e3,
      vary = "phi_scatter",
      quantity = "level_dB",
      n_points = 3
    )
    expect_equal(phi_scatter_slice$quantity, "sigma_scat_dB")
    expect_equal(length(phi_scatter_slice$y), 3)
  })
})

test_that("TMM plotting renderers and dispatcher cover polar, heatmap, and slice branches", {
  grid <- list(
    frequency = 70e3,
    theta_scatter = c(pi / 4, 3 * pi / 4),
    phi_scatter = c(pi / 2, 3 * pi / 2),
    sigma_scat_dB = matrix(c(-10, -9, -8, -7), nrow = 2),
    sigma_scat = matrix(c(1, 2, 3, 4), nrow = 2),
    f_scat = matrix(complex(real = 1:4, imaginary = 4:1), nrow = 2)
  )

  expect_error(
    acousticTS:::.plot_tmm_scattering_polar(
      modifyList(grid, list(sigma_scat_dB = matrix(NA_real_, nrow = 2, ncol = 2)))
    ),
    "does not contain finite plotted values"
  )

  .with_temp_pdf_device({
    expect_identical(
      acousticTS:::.plot_tmm_scattering_heatmap(grid),
      grid
    )
    expect_identical(
      acousticTS:::.plot_tmm_scattering_polar(grid, quantity = "phase"),
      grid
    )
  })

  expect_error(
    acousticTS:::.plot_tmm_scattering_slice(
      object = structure(list(), class = "Scatterer"),
      polar = TRUE,
      heatmap = TRUE
    ),
    "Only one of 'polar' or 'heatmap' can be TRUE"
  )

  local({
    testthat::local_mocked_bindings(
      tmm_scattering_grid = function(...) grid,
      .plot_tmm_scattering_polar = function(grid, quantity) {
        attr(grid, "quantity") <- quantity
        grid
      },
      .plot_tmm_scattering_heatmap = function(grid, quantity) {
        attr(grid, "quantity") <- quantity
        grid
      },
      .package = "acousticTS"
    )

    polar_grid <- acousticTS:::.plot_tmm_scattering_slice(
      object = structure(list(), class = "Scatterer"),
      polar = TRUE,
      quantity = "phase"
    )
    heatmap_grid <- acousticTS:::.plot_tmm_scattering_slice(
      object = structure(list(), class = "Scatterer"),
      heatmap = TRUE
    )
    expect_equal(attr(polar_grid, "quantity"), "phase")
    expect_equal(attr(heatmap_grid, "quantity"), "sigma_scat_dB")
  })

  local({
    testthat::local_mocked_bindings(
      .tmm_scattering_slice_data = function(...) {
        list(
          frequency = 70e3,
          x = c(0, 1, 2),
          y = c(3, 3, 3),
          x_lab = "x",
          y_lab = "y"
        )
      },
      .package = "acousticTS"
    )

    .with_temp_pdf_device({
      slice <- acousticTS:::.plot_tmm_scattering_slice(
        object = structure(list(), class = "Scatterer")
      )
      expect_equal(slice$y, c(3, 3, 3))
    })
  })

  local({
    testthat::local_mocked_bindings(
      .tmm_scattering_slice_data = function(...) {
        list(
          frequency = 70e3,
          x = c(0, 1, 2),
          y = c(NA_real_, NA_real_, NA_real_),
          x_lab = "x",
          y_lab = "y"
        )
      },
      .package = "acousticTS"
    )

    expect_error(
      acousticTS:::.plot_tmm_scattering_slice(
        object = structure(list(), class = "Scatterer")
      ),
      "does not contain any finite plotted values"
    )
  })
})
