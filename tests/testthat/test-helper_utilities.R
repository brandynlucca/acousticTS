library(acousticTS)

test_that("low-level geometry and formatting helpers return stable outputs", {
  rpos <- cbind(
    x = c(0, 1, 2),
    radius = c(0, 1, 0),
    z = c(0, 2, 4)
  )

  resampled <- acousticTS:::.resample_rpos(rpos, n_new = 5)
  expect_equal(resampled[, 1], seq(0, 2, length.out = 5))
  expect_equal(resampled[, 2], c(0, 0.5, 1, 0.5, 0))
  expect_equal(resampled[, 3], c(0, 1, 2, 3, 4))
  expect_equal(colnames(resampled), colnames(rpos))

  lines <- acousticTS:::.format_property_lines(
    properties = c(sound_speed = 1500.126, density = 1026),
    label_map = c(sound_speed = "Sound speed"),
    unit_map = c(sound_speed = " m/s", density = " kg/m^3"),
    digits = 2,
    indent = "  "
  )
  expect_equal(
    unname(lines),
    c("  Sound speed: 1500.13 m/s", "  density: 1026 kg/m^3")
  )
  expect_equal(acousticTS:::.format_property_lines(numeric()), character())
})

test_that("shape-component extraction falls back in the expected order", {
  position_matrix <- cbind(
    x_body = c(0, 1, 2),
    z_body = c(3, 4, 5)
  )

  expect_equal(
    acousticTS:::.extract_shape_component_column(
      position_matrix,
      candidates = c("z_body", "radius")
    ),
    c(3, 4, 5)
  )
  expect_equal(
    acousticTS:::.extract_shape_component_column(
      position_matrix,
      candidates = c("radius"),
      default = c(9, 9, 9)
    ),
    c(9, 9, 9)
  )
  expect_equal(
    acousticTS:::.extract_shape_component_column(
      position_matrix,
      candidates = c("radius")
    ),
    c(0, 1, 2)
  )
})

test_that("pkgdown current-page detection recognizes each vignette page type", {
  skip_if_not_installed("knitr")

  local({
    testthat::local_mocked_bindings(
      current_input = function() "index.Rmd",
      .package = "knitr"
    )
    expect_equal(acousticTS:::.model_family_current_page(), "Overview")
  })

  local({
    testthat::local_mocked_bindings(
      current_input = function() "dwba-theory.Rmd",
      .package = "knitr"
    )
    expect_equal(acousticTS:::.model_family_current_page(), "Theory")
  })

  local({
    testthat::local_mocked_bindings(
      current_input = function() "dwba-implementation.qmd",
      .package = "knitr"
    )
    expect_equal(acousticTS:::.model_family_current_page(), "Implementation")
  })

  local({
    testthat::local_mocked_bindings(
      current_input = function() stop("no current input"),
      .package = "knitr"
    )
    expect_null(acousticTS:::.model_family_current_page())
  })
})
