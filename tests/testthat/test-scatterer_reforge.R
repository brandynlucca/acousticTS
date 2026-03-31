library(acousticTS)

test_that("`.discover_reforge_params` returns expected keys.", {
  # Test FLS parameters
  expect_equal(
    .discover_reforge_params("FLS"),
    c(
      "body_scale", "body_target", "isometric_body", "n_segments_body",
      "length", "radius", "length_radius_ratio_constant", "n_segments"
    )
  )

  # Test BBF parameters
  expect_equal(
    .discover_reforge_params("BBF"),
    c(
      "body_scale", "body_target", "backbone_scale",
      "backbone_target", "isometric_body", "isometric_backbone",
      "maintain_ratio", "n_segments_body", "n_segments_backbone",
      "containment"
    )
  )

  # Test SBF parameters
  expect_equal(
    .discover_reforge_params("SBF"),
    c(
      "body_scale", "body_target", "swimbladder_scale",
      "swimbladder_target", "isometric_body", "isometric_swimbladder",
      "maintain_ratio", "swimbladder_inflation_factor", "n_segments_body",
      "n_segments_swimbladder", "containment"
    )
  )

  # Test class vectors and fallback discovery
  expect_equal(
    .discover_reforge_params(c("OTHER", "FLS")),
    .discover_reforge_params("FLS")
  )

  # Test erroneous input
  expect_equal(
    .discover_reforge_params("OTHER"),
    character(0)
  )
})

test_that("`reforge('BBF')` correctly updates composite body/backbone geometry", {
  fish <- bbf_generate(
    body_shape = arbitrary(
      x_body = c(0, 0.04, 0.08),
      zU_body = c(0.001, 0.004, 0.001),
      zL_body = c(-0.001, -0.004, -0.001)
    ),
    backbone_shape = cylinder(
      length_body = 0.05,
      radius_body = 0.0008,
      n_segments = 40
    ),
    density_body = 1070,
    sound_speed_body = 1570,
    density_backbone = 1900,
    sound_speed_longitudinal_backbone = 3500,
    sound_speed_transversal_backbone = 1700,
    x_offset_backbone = 0.015,
    z_offset_backbone = 0.0002
  )

  orig_shape <- extract(fish, "shape_parameters")
  orig_body <- extract(fish, "body")
  orig_backbone <- extract(fish, "backbone")

  # Explicit body/backbone resampling follows the multi-component reforge path.
  fish_hires <- reforge(
    fish,
    n_segments_body = 80,
    n_segments_backbone = 30
  )
  hires_shape <- extract(fish_hires, "shape_parameters")
  hires_body <- extract(fish_hires, "body")
  hires_backbone <- extract(fish_hires, "backbone")
  expect_equal(hires_shape$body$n_segments, 80)
  expect_equal(hires_shape$backbone$n_segments, 30)
  expect_equal(ncol(hires_body$rpos), 80)
  expect_equal(ncol(hires_backbone$rpos), 30)

  # Default maintain_ratio should resize the backbone with the body.
  scale_factor <- 1.5
  fish_scaled <- reforge(fish, body_scale = scale_factor)
  scaled_shape <- extract(fish_scaled, "shape_parameters")
  scaled_body <- extract(fish_scaled, "body")
  scaled_backbone <- extract(fish_scaled, "backbone")

  expect_equal(
    scaled_shape$body$length,
    orig_shape$body$length * scale_factor
  )
  expect_equal(
    scaled_shape$backbone$length,
    orig_shape$backbone$length * scale_factor
  )
  expect_equal(
    max(scaled_body$radius, na.rm = TRUE),
    max(orig_body$radius, na.rm = TRUE) * scale_factor
  )
  expect_equal(
    max(scaled_backbone$radius, na.rm = TRUE),
    max(orig_backbone$radius, na.rm = TRUE) * scale_factor
  )

  # The backbone's relative axial start should stay tied to the body.
  orig_start_ratio <- (
    min(orig_backbone$rpos["x", ], na.rm = TRUE) -
      min(orig_body$rpos["x", ], na.rm = TRUE)
  ) / orig_shape$body$length
  new_start_ratio <- (
    min(scaled_backbone$rpos["x", ], na.rm = TRUE) -
      min(scaled_body$rpos["x", ], na.rm = TRUE)
  ) / scaled_shape$body$length
  expect_equal(new_start_ratio, orig_start_ratio, tolerance = 1e-10)
  expect_equal(
    .reforge_relative_vertical_offset(orig_backbone$rpos, orig_body$rpos),
    .reforge_relative_vertical_offset(scaled_backbone$rpos, scaled_body$rpos),
    tolerance = 1e-10
  )

  # Body-only scaling should leave the backbone unchanged when ratio lock is off.
  fish_body_only <- reforge(
    fish,
    body_scale = c(length = 2, width = 1.25, height = 1.5),
    isometric_body = FALSE,
    maintain_ratio = FALSE
  )
  expect_equal(
    extract(fish_body_only, c("shape_parameters", "body", "length")),
    orig_shape$body$length * 2
  )
  expect_equal(
    extract(fish_body_only, c("shape_parameters", "backbone", "length")),
    orig_shape$backbone$length
  )
  expect_equal(
    max(extract(fish_body_only, "backbone")$radius, na.rm = TRUE),
    max(orig_backbone$radius, na.rm = TRUE)
  )

  # Guardrails and containment warnings should mirror the other composite paths.
  expect_error(reforge(fish), "Must specify at least one")
  expect_error(
    reforge(fish, body_scale = 2, body_target = c(length = 0.12)),
    "Specify only one of body_scale or body_target, not both."
  )
  expect_error(
    reforge(
      fish,
      backbone_scale = c(width = 1.5, height = 2),
      isometric_backbone = FALSE,
      maintain_ratio = FALSE
    ),
    "Backbone reforge must preserve a circular cross-section"
  )
  expect_warning(
    reforge(fish, backbone_scale = 20, maintain_ratio = FALSE),
    "Backbone exceeds body bounds at some positions."
  )
  expect_error(
    reforge(
      fish,
      backbone_scale = 20,
      maintain_ratio = FALSE,
      containment = "error"
    ),
    "Backbone exceeds body bounds at some positions."
  )
  expect_silent(
    reforge(
      fish,
      backbone_scale = 20,
      maintain_ratio = FALSE,
      containment = "ignore"
    )
  )
})

test_that("`reforge('FLS')` correctly updates body shape", {
  # Call in example FLS object
  data(krill, package = "acousticTS")

  # Get the original body and shape parameter parameters
  orig_shape <- extract(krill, "shape_parameters")
  orig_body <- extract(krill, "body")

  # Test 'n_segments' change
  krill_hires <- reforge(krill, n_segments = 100)
  expect_equal(extract(krill_hires, "shape_parameters")$n_segments, 100)
  expect_equal(length(extract(krill_hires, "body")$rpos[1, ]), 101)
  krill_hires_new <- reforge(krill, n_segments_body = 120)
  expect_equal(extract(krill_hires_new, "shape_parameters")$n_segments, 120)
  expect_equal(length(extract(krill_hires_new, "body")$rpos[1, ]), 121)

  # Test the new shared body-scale pathway
  krill_scaled <- reforge(krill, body_scale = 2)
  expect_equal(
    extract(krill_scaled, "shape_parameters")$length,
    orig_shape$length * 2
  )
  expect_equal(
    extract(krill_scaled, "shape_parameters")$radius,
    orig_shape$radius * 2
  )
  expect_equal(
    max(extract(krill_scaled, "body")$radius),
    max(orig_body$radius) * 2
  )

  target_radius <- 30e-3
  krill_targeted <- reforge(
    krill,
    body_target = c(radius = target_radius),
    isometric_body = FALSE
  )
  expect_equal(
    extract(krill_targeted, "shape_parameters")$length,
    orig_shape$length
  )
  expect_equal(
    extract(krill_targeted, "shape_parameters")$radius,
    target_radius
  )
  expect_equal(
    max(extract(krill_targeted, "body")$radius),
    target_radius
  )
  # Bent-body resizing should honor the true centerline arc length.
  bent_krill <- brake(krill, radius_curvature = 5)
  orig_arc <- .shape_arc_length(body = extract(bent_krill, "body"))
  target_arc <- orig_arc * 1.5
  bent_rescaled <- reforge(bent_krill, body_target = c(length = target_arc))
  bent_body <- extract(bent_rescaled, "body")
  projected_length <- .shape_length(
    position_matrix = bent_body$rpos,
    row_major = TRUE
  )
  expect_equal(
    extract(bent_rescaled, "shape_parameters")$length,
    target_arc,
    tolerance = 1e-3
  )
  expect_equal(bent_body$arc_length, target_arc, tolerance = 1e-3)
  expect_equal(.shape_arc_length(body = bent_body), target_arc, tolerance = 1e-3)
  expect_true(projected_length < target_arc)
  expect_error(
    reforge(krill, length = 0.04, body_scale = 2),
    "Use either the legacy length/radius arguments or the new body_scale/body_target arguments"
  )

  # Test 'length' scaling
  # ---- Case: Constant L/a
  new_length <- orig_shape$length * 2.0
  krill_stretched <- reforge(krill, length = new_length)
  expect_equal(
    extract(krill_stretched, "shape_parameters")$length,
    new_length
  )
  expect_equal(
    extract(krill_stretched, "shape_parameters")$length /
      extract(krill_stretched, "shape_parameters")$radius,
    orig_shape$length / orig_shape$radius
  )
  # ---- Case: Varied L/a
  krill_stretched <- reforge(krill,
    length = new_length,
    length_radius_ratio_constant = FALSE
  )
  expect_equal(
    extract(krill_stretched, "shape_parameters")$length,
    new_length
  )
  expect_equal(
    extract(krill_stretched, "shape_parameters")$length /
      extract(krill_stretched, "shape_parameters")$radius,
    orig_shape$length / orig_shape$radius * 2
  )

  # Test 'radius' scaling
  new_radius <- 30e-3
  mega_krill <- reforge(krill, radius = new_radius)
  expect_equal(
    extract(mega_krill, "shape_parameters")$radius,
    new_radius
  )
  expect_equal(
    max(extract(mega_krill, "body")$radius),
    30e-3
  )
  expect_equal(
    mean(extract(mega_krill, "body")$radius),
    mean(orig_body$radius) * (new_radius / orig_shape$radius)
  )
})

test_that("`reforge('SBF')` correctly updates body shape", {
  # Call in example SBF object
  data(sardine, package = "acousticTS")

  # Get the original body and shape parameter parameters
  orig_shape <- extract(sardine, "shape_parameters")
  orig_body <- extract(sardine, "body")
  orig_bladder <- extract(sardine, "bladder")

  # Test 'n_segments' change for body and bladder
  sardine_hires <- reforge(sardine,
    n_segments_body = 800,
    n_segments_swimbladder = 200
  )
  hires_shape <- extract(sardine_hires, "shape_parameters")
  hires_body <- extract(sardine_hires, "body")
  hires_bladder <- extract(sardine_hires, "bladder")
  expect_equal(hires_shape$body$n_segments, 800)
  expect_equal(hires_shape$bladder$n_segments, 200)
  expect_equal(length(hires_body$rpos[1, ]), 800)
  expect_equal(length(hires_bladder$rpos[1, ]), 200)

  ##############################################################################
  # Scaling dimensions
  # Test 'isometric_body' scaling [body]
  # ---- Case: Uniform across all dimensions
  scale_factor <- 2
  body_expand <- reforge(sardine, body_scale = scale_factor)

  orig_length <- orig_shape$body$length
  orig_width <- max(orig_body$rpos[2, ])
  orig_height <- max(orig_body$rpos[3, ] - orig_body$rpos[4, ])
  new_length <- extract(body_expand, "shape_parameters")$body$length
  new_width <- max(extract(body_expand, "body")$rpos[2, ])
  z_pos <- extract(body_expand, "body")$rpos[3:4, ]
  new_height <- max(head(z_pos, 1) - tail(z_pos, 1))

  expect_equal(new_length, orig_length * scale_factor)
  expect_equal(new_width, orig_width * scale_factor)
  expect_equal(new_height, orig_height * scale_factor)
  expect_equal(
    extract(body_expand, "shape_parameters")$body$length,
    new_length
  )

  orig_sblength <- orig_shape$bladder$length
  orig_sbwidth <- max(orig_bladder$rpos[2, ])
  orig_sbheight <- max(orig_bladder$rpos[3, ] - orig_bladder$rpos[4, ])
  new_sblength <- extract(body_expand, "shape_parameters")$bladder$length
  new_sbwidth <- max(extract(body_expand, "bladder")$rpos[2, ])
  z_sbpos <- extract(body_expand, "bladder")$rpos[3:4, ]
  new_sbheight <- max(head(z_sbpos, 1) - tail(z_sbpos, 1))

  expect_equal(new_sblength, orig_sblength * scale_factor)
  expect_equal(new_sbwidth, orig_sbwidth * scale_factor)
  expect_equal(new_sbheight, orig_sbheight * scale_factor)
  expect_equal(
    extract(body_expand, "shape_parameters")$bladder$length,
    new_sblength
  )
  expect_equal(
    .reforge_relative_vertical_offset(orig_bladder$rpos, orig_body$rpos),
    .reforge_relative_vertical_offset(
      extract(body_expand, "bladder")$rpos,
      extract(body_expand, "body")$rpos
    ),
    tolerance = 1e-10
  )
  # ---- Case: Defining just length -> applies uniformly across dimensions
  body_expand <- reforge(sardine, body_scale = c(length = scale_factor))

  new_length <- extract(body_expand, "shape_parameters")$body$length
  new_width <- max(extract(body_expand, "body")$rpos[2, ])
  z_pos <- extract(body_expand, "body")$rpos[3:4, ]
  new_height <- max(head(z_pos, 1) - tail(z_pos, 1))
  expect_equal(new_length, orig_length * scale_factor)
  expect_equal(new_width, orig_width * scale_factor)
  expect_equal(new_height, orig_height * scale_factor)
  expect_equal(
    extract(body_expand, "shape_parameters")$body$length,
    new_length
  )

  new_sblength <- extract(body_expand, "shape_parameters")$bladder$length
  new_sbwidth <- max(extract(body_expand, "bladder")$rpos[2, ])
  z_sbpos <- extract(body_expand, "bladder")$rpos[3:4, ]
  new_sbheight <- max(head(z_sbpos, 1) - tail(z_sbpos, 1))

  expect_equal(new_sblength, orig_sblength * scale_factor)
  expect_equal(new_sbwidth, orig_sbwidth * scale_factor)
  expect_equal(new_sbheight, orig_sbheight * scale_factor)
  expect_equal(
    extract(body_expand, "shape_parameters")$bladder$length,
    new_sblength
  )
  # ---- Case: Defining length and width -> overrides `isometric_body`
  expect_error(
    reforge(sardine, body_scale = c(length = 4, width = 1.5)),
    paste0(
      "'body_scale' contains more than 1 dimension while 'isometric_body' ",
      "is TRUE. When TRUE 'body_scale' must be a scalar value."
    )
  )
  # ---- Case: Defining length, width, and height -> overrides `isometric_body`
  expect_error(
    reforge(sardine, body_scale = c(length = 4, width = 1.5, height = 2)),
    paste0(
      "'body_scale' contains more than 1 dimension while 'isometric_body' ",
      "is TRUE. When TRUE 'body_scale' must be a scalar value."
    )
  )
  # ---- Case: Repeat previous but drop 'isometric_body'
  body_expand <- reforge(sardine,
    body_scale = c(length = 4, width = 1.5, height = 2.0),
    isometric_body = FALSE
  )

  new_length <- extract(body_expand, "shape_parameters")$body$length
  new_width <- max(extract(body_expand, "body")$rpos[2, ])
  z_pos <- extract(body_expand, "body")$rpos[3:4, ]
  new_height <- max(head(z_pos, 1) - tail(z_pos, 1))
  expect_equal(new_length, orig_length * 4)
  expect_equal(new_width, orig_width * 1.5)
  expect_equal(new_height, orig_height * 2)
  expect_equal(
    extract(body_expand, "shape_parameters")$body$length,
    new_length
  )

  new_sblength <- extract(body_expand, "shape_parameters")$bladder$length
  new_sbwidth <- max(extract(body_expand, "bladder")$rpos[2, ])
  z_sbpos <- extract(body_expand, "bladder")$rpos[3:4, ]
  new_sbheight <- max(head(z_sbpos, 1) - tail(z_sbpos, 1))

  expect_equal(new_sblength, orig_sblength * 4)
  expect_equal(new_sbwidth, orig_sbwidth * 1.5)
  expect_equal(new_sbheight, orig_sbheight * 2)
  expect_equal(
    extract(body_expand, "shape_parameters")$bladder$length,
    new_sblength
  )
  # ---- Case: Repeat previous but drop 'maintain_ratio'
  body_expand <- reforge(
    sardine,
    body_scale = c(length = 4, width = 1.5, height = 2.0),
    isometric_body = FALSE,
    maintain_ratio = FALSE
  )

  new_length <- extract(body_expand, "shape_parameters")$body$length
  new_width <- max(extract(body_expand, "body")$rpos[2, ])
  z_pos <- extract(body_expand, "body")$rpos[3:4, ]
  new_height <- max(head(z_pos, 1) - tail(z_pos, 1))
  expect_equal(new_length, orig_length * 4)
  expect_equal(new_width, orig_width * 1.5)
  expect_equal(new_height, orig_height * 2)
  expect_equal(
    extract(body_expand, "shape_parameters")$body$length,
    new_length
  )

  new_sblength <- extract(body_expand, "shape_parameters")$bladder$length
  new_sbwidth <- max(extract(body_expand, "bladder")$rpos[2, ])
  z_sbpos <- extract(body_expand, "bladder")$rpos[3:4, ]
  new_sbheight <- max(head(z_sbpos, 1) - tail(z_sbpos, 1))

  expect_equal(new_sblength, orig_sblength)
  expect_equal(new_sbwidth, orig_sbwidth)
  expect_equal(new_sbheight, orig_sbheight)
  expect_equal(
    extract(body_expand, "shape_parameters")$bladder$length,
    new_sblength
  )

  # Test 'isometric_bladder' scaling [bladder]
  # ---- Case: Uniform across all dimensions
  bladder_expand <- reforge(sardine,
    swimbladder_scale = 5
  )

  new_length <- extract(bladder_expand, "shape_parameters")$body$length
  new_width <- max(extract(bladder_expand, "body")$rpos[2, ])
  z_pos <- extract(bladder_expand, "body")$rpos[3:4, ]
  new_height <- max(head(z_pos, 1) - tail(z_pos, 1))
  expect_equal(new_length, orig_length * 5)
  expect_equal(new_width, orig_width * 5)
  expect_equal(new_height, orig_height * 5)
  expect_equal(
    extract(bladder_expand, "shape_parameters")$body$length,
    new_length
  )

  new_sblength <- extract(bladder_expand, "shape_parameters")$bladder$length
  new_sbwidth <- max(extract(bladder_expand, "bladder")$rpos[2, ])
  z_sbpos <- extract(bladder_expand, "bladder")$rpos[3:4, ]
  new_sbheight <- max(head(z_sbpos, 1) - tail(z_sbpos, 1))

  expect_equal(new_sblength, orig_sblength * 5)
  expect_equal(new_sbwidth, orig_sbwidth * 5)
  expect_equal(new_sbheight, orig_sbheight * 5)
  expect_equal(
    extract(bladder_expand, "shape_parameters")$bladder$length,
    new_sblength
  )
  # ---- Case: Defining just length -> applies uniformly across dimensions
  bladder_expand <- reforge(sardine,
    swimbladder_scale = c(length = 5)
  )

  new_length <- extract(bladder_expand, "shape_parameters")$body$length
  new_width <- max(extract(bladder_expand, "body")$rpos[2, ])
  z_pos <- extract(bladder_expand, "body")$rpos[3:4, ]
  new_height <- max(head(z_pos, 1) - tail(z_pos, 1))
  expect_equal(new_length, orig_length * 5)
  expect_equal(new_width, orig_width * 5)
  expect_equal(new_height, orig_height * 5)
  expect_equal(
    extract(bladder_expand, "shape_parameters")$body$length,
    new_length
  )

  new_sblength <- extract(bladder_expand, "shape_parameters")$bladder$length
  new_sbwidth <- max(extract(bladder_expand, "bladder")$rpos[2, ])
  z_sbpos <- extract(bladder_expand, "bladder")$rpos[3:4, ]
  new_sbheight <- max(head(z_sbpos, 1) - tail(z_sbpos, 1))

  expect_equal(new_sblength, orig_sblength * 5)
  expect_equal(new_sbwidth, orig_sbwidth * 5)
  expect_equal(new_sbheight, orig_sbheight * 5)
  expect_equal(
    extract(bladder_expand, "shape_parameters")$bladder$length,
    new_sblength
  )
  # ---- Case: Defining length and width -> overrides `isometric_body`
  expect_error(
    reforge(
      sardine,
      swimbladder_scale = c(length = 4, width = 1.5)
    ),
    paste0(
      "'swimbladder_scale' contains more than 1 dimension while ",
      "'isometric_swimbladder' is TRUE. When TRUE 'swimbladder_scale' must ",
      "be a scalar value."
    )
  )
  # ---- Case: Defining length, width, and height -> overrides `isometric_body`
  expect_error(
    reforge(sardine,
      swimbladder_scale = c(length = 4, width = 1.5, height = 2.0)
    ),
    paste0(
      "'swimbladder_scale' contains more than 1 dimension while ",
      "'isometric_swimbladder' is TRUE. When TRUE 'swimbladder_scale' must ",
      "be a scalar value."
    )
  )
  # ---- Case: Repeat previous but drop 'isometric_swimbladder'
  bladder_expand <- reforge(
    sardine,
    swimbladder_scale = c(
      length = 4,
      width = 1.5,
      height = 2.0
    ),
    isometric_swimbladder = FALSE
  )

  new_length <- extract(bladder_expand, "shape_parameters")$body$length
  new_width <- max(extract(bladder_expand, "body")$rpos[2, ])
  z_pos <- extract(bladder_expand, "body")$rpos[3:4, ]
  new_height <- max(head(z_pos, 1) - tail(z_pos, 1))
  expect_equal(new_length, orig_length * 4)
  expect_equal(new_width, orig_width * 1.5)
  expect_equal(new_height, orig_height * 2)
  expect_equal(
    extract(bladder_expand, "shape_parameters")$body$length,
    new_length
  )

  new_sblength <- extract(bladder_expand, "shape_parameters")$bladder$length
  new_sbwidth <- max(extract(bladder_expand, "bladder")$rpos[2, ])
  z_sbpos <- extract(bladder_expand, "bladder")$rpos[3:4, ]
  new_sbheight <- max(head(z_sbpos, 1) - tail(z_sbpos, 1))

  expect_equal(new_sblength, orig_sblength * 4)
  expect_equal(new_sbwidth, orig_sbwidth * 1.5)
  expect_equal(new_sbheight, orig_sbheight * 2)
  expect_equal(
    extract(bladder_expand, "shape_parameters")$bladder$length,
    new_sblength
  )
  # ---- Case: Repeat previous but drop 'maintain_ratio'
  suppressWarnings({
    bladder_expand <- reforge(
      sardine,
      swimbladder_scale = c(
        length = 4,
        width = 1.5,
        height = 2.0
      ),
      isometric_swimbladder = FALSE,
      maintain_ratio = FALSE
    )
  })
  expect_warning(
    reforge(sardine,
      swimbladder_scale = c(length = 4, width = 1.5, height = 2.0),
      isometric_swimbladder = FALSE,
      maintain_ratio = FALSE
    ),
    "Swimbladder exceeds body bounds at some positions."
  )
  expect_error(
    reforge(
      sardine,
      swimbladder_scale = c(length = 4, width = 1.5, height = 2.0),
      isometric_swimbladder = FALSE,
      maintain_ratio = FALSE,
      containment = "error"
    ),
    "Swimbladder exceeds body bounds at some positions."
  )
  expect_silent(
    reforge(
      sardine,
      swimbladder_scale = c(length = 4, width = 1.5, height = 2.0),
      isometric_swimbladder = FALSE,
      maintain_ratio = FALSE,
      containment = "ignore"
    )
  )

  new_length <- extract(bladder_expand, "shape_parameters")$body$length
  new_width <- max(extract(bladder_expand, "body")$rpos[2, ])
  z_pos <- extract(bladder_expand, "body")$rpos[3:4, ]
  new_height <- max(head(z_pos, 1) - tail(z_pos, 1))
  expect_equal(new_length, orig_length)
  expect_equal(new_width, orig_width)
  expect_equal(new_height, orig_height)
  expect_equal(
    extract(bladder_expand, "shape_parameters")$body$length,
    new_length
  )

  new_sblength <- extract(bladder_expand, "shape_parameters")$bladder$length
  new_sbwidth <- max(extract(bladder_expand, "bladder")$rpos[2, ])
  z_sbpos <- extract(bladder_expand, "bladder")$rpos[3:4, ]
  new_sbheight <- max(head(z_sbpos, 1) - tail(z_sbpos, 1))

  expect_equal(new_sblength, orig_sblength * 4)
  expect_equal(new_sbwidth, orig_sbwidth * 1.5)
  expect_equal(new_sbheight, orig_sbheight * 2.0)
  expect_equal(
    extract(bladder_expand, "shape_parameters")$bladder$length,
    new_sblength
  )

  ##############################################################################
  # Target dimensions
  # Test manual target dimension entry
  # ---- Case: unnamed numeric -> Error
  expect_error(
    reforge(sardine, body_target = 1),
    paste0(
      "'body_target' must either be a scalar or a named vector with ",
      "dimensions: 'length', 'width', 'height'."
    )
  )
  expect_error(
    reforge(sardine, swimbladder_target = 1),
    paste0(
      "'swimbladder_target' must either be a scalar or a named vector with ",
      "dimensions: 'length', 'width', 'height'."
    )
  )
  # ---- Case: Defining just length -> applies uniformly across dimensions
  body_expand <- reforge(sardine, body_target = c(length = 0.8))

  new_length <- extract(body_expand, "shape_parameters")$body$length
  scaling_length_ratio <- new_length / orig_shape$body$length
  new_width <- max(extract(body_expand, "body")$rpos[2, ])
  z_pos <- extract(body_expand, "body")$rpos[3:4, ]
  new_height <- max(head(z_pos, 1) - tail(z_pos, 1))
  expect_equal(new_length, orig_length * scaling_length_ratio)
  expect_equal(new_width, orig_width * scaling_length_ratio)
  expect_equal(new_height, orig_height * scaling_length_ratio)
  expect_equal(
    extract(body_expand, "shape_parameters")$body$length,
    new_length
  )

  new_sblength <- extract(body_expand, "shape_parameters")$bladder$length
  new_sbwidth <- max(extract(body_expand, "bladder")$rpos[2, ])
  z_sbpos <- extract(body_expand, "bladder")$rpos[3:4, ]
  new_sbheight <- max(head(z_sbpos, 1) - tail(z_sbpos, 1))

  expect_equal(new_sblength, orig_sblength * scaling_length_ratio)
  expect_equal(new_sbwidth, orig_sbwidth * scaling_length_ratio)
  expect_equal(new_sbheight, orig_sbheight * scaling_length_ratio)
  expect_equal(
    extract(body_expand, "shape_parameters")$bladder$length,
    new_sblength
  )
  # # ---- Case: Defining length and width -> overrides `isometric_body`
  expect_error(
    reforge(sardine, body_target = c(length = 0.8, width = 0.4)),
    paste0(
      "'body_target' contains more than 1 dimension while 'isometric_body' ",
      "is TRUE. When TRUE 'body_target' must be a scalar value."
    )
  )

  # ---- Case: Defining length/width/height -> overrides `isometric_body`
  expect_error(
    reforge(sardine, body_target = c(length = 0.8, width = 0.4, height = 0.1)),
    paste0(
      "'body_target' contains more than 1 dimension while 'isometric_body' ",
      "is TRUE. When TRUE 'body_target' must be a scalar value."
    )
  )
  # ---- Case: Repeat previous but drop 'isometric_body'
  body_expand <- reforge(sardine,
    body_target = c(length = 0.8, width = 0.4, height = 0.1),
    isometric_body = FALSE
  )

  new_length <- extract(body_expand, "shape_parameters")$body$length
  new_width <- max(extract(body_expand, "body")$rpos[2, ])
  z_pos <- extract(body_expand, "body")$rpos[3:4, ]
  new_height <- max(head(z_pos, 1) - tail(z_pos, 1))
  scaling_width_ratio <- new_width / orig_width
  scaling_height_ratio <- new_height / orig_height
  expect_equal(new_length, orig_length * scaling_length_ratio)
  expect_equal(new_width, orig_width * scaling_width_ratio)
  expect_equal(new_height, orig_height * scaling_height_ratio)
  expect_equal(
    extract(body_expand, "shape_parameters")$body$length,
    new_length
  )

  new_sblength <- extract(body_expand, "shape_parameters")$bladder$length
  new_sbwidth <- max(extract(body_expand, "bladder")$rpos[2, ])
  z_sbpos <- extract(body_expand, "bladder")$rpos[3:4, ]
  new_sbheight <- max(head(z_sbpos, 1) - tail(z_sbpos, 1))

  expect_equal(new_sblength, orig_sblength * scaling_length_ratio)
  expect_equal(new_sbwidth, orig_sbwidth * scaling_width_ratio)
  expect_equal(new_sbheight, orig_sbheight * scaling_height_ratio)
  expect_equal(
    extract(body_expand, "shape_parameters")$bladder$length,
    new_sblength
  )
  # ---- Case: Repeat previous but drop 'maintain_ratio'
  body_expand <- reforge(sardine,
    body_target = c(length = 0.8, width = 0.4, height = 0.1),
    isometric_body = FALSE,
    maintain_ratio = FALSE
  )

  new_length <- extract(body_expand, "shape_parameters")$body$length
  new_width <- max(extract(body_expand, "body")$rpos[2, ])
  z_pos <- extract(body_expand, "body")$rpos[3:4, ]
  new_height <- max(head(z_pos, 1) - tail(z_pos, 1))
  expect_equal(new_length, orig_length * scaling_length_ratio)
  expect_equal(new_width, orig_width * scaling_width_ratio)
  expect_equal(new_height, orig_height * scaling_height_ratio)
  expect_equal(
    extract(body_expand, "shape_parameters")$body$length,
    new_length
  )

  new_sblength <- extract(body_expand, "shape_parameters")$bladder$length
  new_sbwidth <- max(extract(body_expand, "bladder")$rpos[2, ])
  z_sbpos <- extract(body_expand, "bladder")$rpos[3:4, ]
  new_sbheight <- max(head(z_sbpos, 1) - tail(z_sbpos, 1))

  expect_equal(new_sblength, orig_sblength)
  expect_equal(new_sbwidth, orig_sbwidth)
  expect_equal(new_sbheight, orig_sbheight)
  expect_equal(
    extract(body_expand, "shape_parameters")$bladder$length,
    new_sblength
  )

  # Test 'isometric_bladder' scaling [bladder]
  # ---- Case: Defining just length -> applies uniformly across dimensions
  scaling_length_ratio <- 0.2 / orig_sblength
  scaling_width_ratio <- 0.3 / orig_sbheight
  scaling_height_ratio <- 0.1 / orig_sbheight
  bladder_expand <- reforge(sardine,
    swimbladder_target = c(length = 0.2)
  )

  new_length <- extract(bladder_expand, "shape_parameters")$body$length
  new_width <- max(extract(bladder_expand, "body")$rpos[2, ])
  z_pos <- extract(bladder_expand, "body")$rpos[3:4, ]
  new_height <- max(head(z_pos, 1) - tail(z_pos, 1))
  expect_equal(new_length, orig_length * scaling_length_ratio)
  expect_equal(new_width, orig_width * scaling_length_ratio)
  expect_equal(new_height, orig_height * scaling_length_ratio)
  expect_equal(
    extract(bladder_expand, "shape_parameters")$body$length,
    new_length
  )

  new_sblength <- extract(bladder_expand, "shape_parameters")$bladder$length
  new_sbwidth <- max(extract(bladder_expand, "bladder")$rpos[2, ])
  z_sbpos <- extract(bladder_expand, "bladder")$rpos[3:4, ]
  new_sbheight <- max(head(z_sbpos, 1) - tail(z_sbpos, 1))

  expect_equal(new_sblength, orig_sblength * scaling_length_ratio)
  expect_equal(new_sbwidth, orig_sbwidth * scaling_length_ratio)
  expect_equal(new_sbheight, orig_sbheight * scaling_length_ratio)
  expect_equal(
    extract(bladder_expand, "shape_parameters")$bladder$length,
    new_sblength
  )
  # ---- Case: Defining length and width -> overrides `isometric_body`
  expect_error(
    reforge(sardine,
      swimbladder_target = c(length = 0.2, width = 0.3)
    ),
    paste0(
      "'swimbladder_target' contains more than 1 dimension while ",
      "'isometric_swimbladder' is TRUE. When TRUE 'swimbladder_target' must ",
      "be a scalar value."
    )
  )
  # ---- Case: Defining length, width, and height -> overrides `isometric_body`
  expect_error(
    reforge(sardine,
      swimbladder_target = c(length = 0.2, width = 0.3, height = 0.1)
    ),
    paste0(
      "'swimbladder_target' contains more than 1 dimension while ",
      "'isometric_swimbladder' is TRUE. When TRUE 'swimbladder_target' must ",
      "be a scalar value."
    )
  )
  # ---- Case: Repeat previous but drop 'isometric_swimbladder'
  bladder_expand <- reforge(sardine,
    swimbladder_target = c(
      length = 0.2,
      width = 0.3,
      height = 0.1
    ),
    isometric_swimbladder = FALSE
  )

  new_length <- extract(bladder_expand, "shape_parameters")$body$length
  new_width <- max(extract(bladder_expand, "body")$rpos[2, ])
  z_pos <- extract(bladder_expand, "body")$rpos[3:4, ]
  new_height <- max(head(z_pos, 1) - tail(z_pos, 1))
  expect_equal(new_length, orig_length * scaling_length_ratio)
  expect_equal(new_width, orig_width * scaling_width_ratio)
  expect_equal(new_height, orig_height * scaling_height_ratio)
  expect_equal(
    extract(bladder_expand, "shape_parameters")$body$length,
    new_length
  )

  new_sblength <- extract(bladder_expand, "shape_parameters")$bladder$length
  new_sbwidth <- max(extract(bladder_expand, "bladder")$rpos[2, ])
  z_sbpos <- extract(bladder_expand, "bladder")$rpos[3:4, ]
  new_sbheight <- max(head(z_sbpos, 1) - tail(z_sbpos, 1))

  expect_equal(new_sblength, orig_sblength * scaling_length_ratio)
  expect_equal(new_sbwidth, orig_sbwidth * scaling_width_ratio)
  expect_equal(new_sbheight, orig_sbheight * scaling_height_ratio)
  expect_equal(
    extract(bladder_expand, "shape_parameters")$bladder$length,
    new_sblength
  )
  # ---- Case: Repeat previous but drop 'maintain_ratio'
  suppressWarnings({
    bladder_expand <- reforge(sardine,
      swimbladder_target = c(
        length = 0.2,
        width = 0.3,
        height = 0.1
      ),
      isometric_swimbladder = FALSE,
      maintain_ratio = FALSE
    )
  })
  expect_warning(
    reforge(sardine,
      swimbladder_target = c(length = 0.2, width = 0.3, height = 0.1),
      isometric_swimbladder = FALSE,
      maintain_ratio = FALSE
    ),
    "Swimbladder exceeds body bounds at some positions."
  )

  new_length <- extract(bladder_expand, "shape_parameters")$body$length
  new_width <- max(extract(bladder_expand, "body")$rpos[2, ])
  z_pos <- extract(bladder_expand, "body")$rpos[3:4, ]
  new_height <- max(head(z_pos, 1) - tail(z_pos, 1))
  expect_equal(new_length, orig_length)
  expect_equal(new_width, orig_width)
  expect_equal(new_height, orig_height)
  expect_equal(
    extract(bladder_expand, "shape_parameters")$body$length,
    new_length
  )

  new_sblength <- extract(bladder_expand, "shape_parameters")$bladder$length
  new_sbwidth <- max(extract(bladder_expand, "bladder")$rpos[2, ])
  z_sbpos <- extract(bladder_expand, "bladder")$rpos[3:4, ]
  new_sbheight <- max(head(z_sbpos, 1) - tail(z_sbpos, 1))

  expect_equal(new_sblength, orig_sblength * scaling_length_ratio)
  expect_equal(new_sbwidth, orig_sbwidth * scaling_width_ratio)
  expect_equal(new_sbheight, orig_sbheight * scaling_height_ratio)
  expect_equal(
    extract(bladder_expand, "shape_parameters")$bladder$length,
    orig_sblength * scaling_length_ratio
  )

  ##############################################################################
  # Mixed entries
  # ---- Case: body and swimbladder scaling
  expect_message(
    reforge(sardine, body_scale = 2, swimbladder_scale = 1.5),
    paste0(
      "Multiple axes specified for the body and swimbladder: 'maintain_ratio' ",
      "will be ignored for those axes."
    )
  )
  suppressMessages({
    body_bladder_expand <- reforge(
      sardine,
      body_scale = 2,
      swimbladder_scale = 1.5
    )
  })

  new_length <- extract(body_bladder_expand, "shape_parameters")$body$length
  new_width <- max(extract(body_bladder_expand, "body")$rpos[2, ])
  z_pos <- extract(body_bladder_expand, "body")$rpos[3:4, ]
  new_height <- max(head(z_pos, 1) - tail(z_pos, 1))
  expect_equal(new_length, orig_length * 2)
  expect_equal(new_width, orig_width * 2)
  expect_equal(new_height, orig_height * 2)
  expect_equal(
    extract(body_bladder_expand, "shape_parameters")$body$length,
    orig_length * 2
  )

  new_sblength <- extract(
    body_bladder_expand,
    "shape_parameters"
  )$bladder$length
  new_sbwidth <- max(extract(body_bladder_expand, "bladder")$rpos[2, ])
  z_sbpos <- extract(body_bladder_expand, "bladder")$rpos[3:4, ]
  new_sbheight <- max(head(z_sbpos, 1) - tail(z_sbpos, 1))

  expect_equal(new_sblength, orig_sblength * 1.5)
  expect_equal(new_sbwidth, orig_sbwidth * 1.5)
  expect_equal(new_sbheight, orig_sbheight * 1.5)
  expect_equal(
    extract(body_bladder_expand, "shape_parameters")$bladder$length,
    orig_sblength * 1.5
  )
  # ---- Case: body scaling and swimbladder target
  expect_message(
    reforge(
      sardine,
      body_scale = 20,
      swimbladder_target = c(length = 0.2, width = 0.3, height = 0.1),
      isometric_swimbladder = FALSE
    ),
    paste0(
      "Multiple axes specified for the body and swimbladder: 'maintain_ratio' ",
      "will be ignored for those axes."
    )
  )
  body_bladder_expand <- reforge(
    sardine,
    body_scale = 20,
    swimbladder_target = c(length = 0.2, width = 0.3, height = 0.1),
    isometric_swimbladder = FALSE,
    maintain_ratio = FALSE
  )

  new_length <- extract(body_bladder_expand, "shape_parameters")$body$length
  new_width <- max(extract(body_bladder_expand, "body")$rpos[2, ])
  z_pos <- extract(body_bladder_expand, "body")$rpos[3:4, ]
  new_height <- max(head(z_pos, 1) - tail(z_pos, 1))
  expect_equal(new_length, orig_length * 20)
  expect_equal(new_width, orig_width * 20)
  expect_equal(new_height, orig_height * 20)
  expect_equal(
    extract(body_bladder_expand, "shape_parameters")$body$length,
    orig_length * 20
  )

  new_sblength <- extract(
    body_bladder_expand,
    "shape_parameters"
  )$bladder$length
  new_sbwidth <- max(extract(body_bladder_expand, "bladder")$rpos[2, ])
  z_sbpos <- extract(body_bladder_expand, "bladder")$rpos[3:4, ]
  new_sbheight <- max(head(z_sbpos, 1) - tail(z_sbpos, 1))

  expect_equal(new_sblength, orig_sblength * scaling_length_ratio)
  expect_equal(new_sbwidth, orig_sbwidth * scaling_width_ratio)
  expect_equal(new_sbheight, orig_sbheight * scaling_height_ratio)
  expect_equal(
    extract(body_bladder_expand, "shape_parameters")$bladder$length,
    orig_sblength * scaling_length_ratio
  )

  ##############################################################################
  # Inflation factor
  # ---- Case: adjusting swimbladder inflaction factor
  bloated_bladder <- reforge(sardine, swimbladder_inflation_factor = 1.2)

  new_sblength <- extract(bloated_bladder, "shape_parameters")$bladder$length
  new_sbwidth <- max(extract(bloated_bladder, "bladder")$rpos[2, ])
  z_sbpos <- extract(bloated_bladder, "bladder")$rpos[3:4, ]
  new_sbheight <- max(head(z_sbpos, 1) - tail(z_sbpos, 1))

  expect_equal(new_sblength, orig_sblength)
  expect_equal(new_sbwidth, orig_sbwidth * 1.2)
  expect_equal(new_sbheight, orig_sbheight * 1.2)
  expect_equal(
    extract(bloated_bladder, "shape_parameters")$bladder$length,
    orig_sblength
  )
  # ---- Case: pinch the inflation factor
  pinched_bladder <- reforge(sardine, swimbladder_inflation_factor = 0.1)

  new_sblength <- extract(bloated_bladder, "shape_parameters")$bladder$length
  new_sbwidth <- max(extract(bloated_bladder, "bladder")$rpos[2, ])
  z_sbpos <- extract(bloated_bladder, "bladder")$rpos[3:4, ]
  new_sbheight <- max(head(z_sbpos, 1) - tail(z_sbpos, 1))

  expect_equal(new_sblength, orig_sblength)
  expect_equal(new_sbwidth, orig_sbwidth * 1.2)
  expect_equal(new_sbheight, orig_sbheight * 1.2)
  expect_equal(new_sblength, orig_sblength)
  # ---- Case: EXPLODE the swimbladder
  expect_warning(
    reforge(sardine, swimbladder_inflation_factor = 100),
    "Swimbladder exceeds body bounds at some positions."
  )

  ##############################################################################
  # Additional error states
  # ---- Case: Empty
  expect_error(
    reforge(sardine),
    paste0(
      "Must specify at least one scaling, target, inflation factor, or ",
      "segment count parameter."
    )
  )
  # ---- Case: Single invalid dimension
  expect_error(
    reforge(sardine, body_scale = c(invalid_dim = 2)),
    "'body_scale' has one or more invalid dimensions: 'invalid_dim'."
  )
  expect_error(
    reforge(sardine, body_target = c(invalid_dim = 2)),
    "'body_target' has one or more invalid dimensions: 'invalid_dim'."
  )
  # ---- Case: Multiple invalid dimensions
  expect_error(
    reforge(
      sardine,
      body_scale = c(length = 2, invalid_dim1 = 2, invalid_dim2 = 2),
      isometric_body = FALSE
    ),
    paste0(
      "'body_scale' has one or more invalid dimensions: ",
      "'invalid_dim1', 'invalid_dim2'."
    )
  )
  # ---- Case: Body scale & target
  expect_error(
    reforge(sardine, body_scale = 2, body_target = c(length = 1)),
    "Specify only one of body_scale or body_target, not both."
  )
  # ---- Case: Swimbladder scale & target
  expect_error(
    reforge(sardine, swimbladder_scale = 2, swimbladder_target = c(length = 1)),
    "Specify only one of swimbladder_scale or swimbladder_target, not both."
  )
  # ---- Case: Unspecified body dimensions
  expect_error(
    reforge(sardine, body_scale = c(1, 2, 3), isometric_body = FALSE),
    paste0(
      "'body_scale' must either be a scalar or a named vector with ",
      "dimensions: 'length', 'width', 'height'."
    )
  )
  # ---- Case: Unspecified swimbladder dimensions
  expect_error(
    reforge(
      sardine,
      swimbladder_scale = c(1, 2, 3),
      isometric_swimbladder = FALSE
    ),
    paste0(
      "'swimbladder_scale' must either be a scalar or a named vector with ",
      "dimensions: 'length', 'width', 'height'."
    )
  )
  # ---- Case: Unspecified body scale without isometry
  expect_error(
    reforge(sardine, body_scale = 1, isometric_body = FALSE),
    paste0(
      "'body_scale' must either be a scalar or a named vector with ",
      "dimensions: 'length', 'width', 'height'."
    )
  )
  # ---- Case: Unspecified swimbladder scale without isometry
  expect_error(
    reforge(sardine, swimbladder_scale = 1, isometric_swimbladder = FALSE),
    paste0(
      "'swimbladder_scale' must either be a scalar or a named vector ",
      "with dimensions: 'length', 'width', 'height'."
    )
  )
})

test_that("`reforge()` covers GAS, CAL, and ESS resize methods", {
  gas_obj <- gas_generate(
    shape = sphere(radius_body = 0.01, n_segments = 40),
    density_fluid = 1.24,
    sound_speed_fluid = 345
  )
  gas_reforged <- reforge(gas_obj, radius_target = 0.02, n_segments = 20)
  expect_equal(max(extract(gas_reforged, "shape_parameters")$radius), 0.02, tolerance = 1e-12)
  expect_equal(max(extract(gas_reforged, "body")$radius), 0.02, tolerance = 1e-12)
  expect_equal(extract(gas_reforged, "shape_parameters")$n_segments, 20)
  expect_equal(nrow(extract(gas_reforged, "body")$rpos), 21)
  expect_error(reforge(gas_obj), "Must specify at least one")
  expect_error(reforge(gas_obj, scale = 2, radius_target = 0.02), "Specify only one")
  expect_error(reforge(gas_obj, scale = -1), "'scale' must be a single positive number")
  expect_error(
    reforge(gas_obj, radius_target = -0.02),
    "'radius_target' must be a single positive number"
  )
  expect_error(
    reforge(gas_obj, n_segments = 0),
    "'n_segments' must be a single positive integer"
  )

  cal_obj <- cal_generate(diameter = 38.1e-3, n_segments = 60)
  cal_reforged <- reforge(cal_obj, diameter_target = 0.05, n_segments = 30)
  expect_equal(extract(cal_reforged, "shape_parameters")$diameter, 0.05, tolerance = 1e-12)
  expect_equal(extract(cal_reforged, "body")$diameter, 0.05, tolerance = 1e-12)
  expect_equal(extract(cal_reforged, "shape_parameters")$n_segments, 30)
  expect_equal(nrow(extract(cal_reforged, "body")$rpos), 31)
  expect_error(reforge(cal_obj), "Must specify at least one")
  expect_error(reforge(cal_obj, scale = 2, diameter_target = 0.05), "Specify only one")
  expect_error(reforge(cal_obj, scale = -1), "'scale' must be a single positive number")
  expect_error(
    reforge(cal_obj, diameter_target = -0.05),
    "'diameter_target' must be a single positive number"
  )
  expect_error(
    reforge(cal_obj, n_segments = 0),
    "'n_segments' must be a single positive integer"
  )

  ess_obj <- fixture_sphere("shelled_liquid")
  ess_reforged <- reforge(
    ess_obj,
    radius_target = 0.015,
    shell_thickness = 0.002,
    n_segments = 20
  )
  ess_shape <- extract(ess_reforged, "shape_parameters")
  expect_equal(max(ess_shape$shell$radius), 0.015, tolerance = 1e-12)
  expect_equal(max(ess_shape$fluid$radius), 0.013, tolerance = 1e-12)
  expect_equal(extract(ess_reforged, "shell")$shell_thickness, 0.002, tolerance = 1e-12)
  expect_equal(ess_shape$n_segments, 20)
  expect_equal(nrow(extract(ess_reforged, "shell")$rpos), 21)
  expect_equal(nrow(extract(ess_reforged, "fluid")$rpos), 21)
  expect_error(reforge(ess_obj), "Must specify at least one")
  expect_error(reforge(ess_obj, scale = 2, radius_target = 0.015), "Specify only one")
  expect_error(reforge(ess_obj, scale = -1), "'scale' must be a single positive number")
  expect_error(
    reforge(ess_obj, radius_target = -0.015),
    "'radius_target' must be a single positive number"
  )
  expect_error(
    reforge(ess_obj, shell_thickness = -0.002),
    "'shell_thickness' must be a single positive number"
  )
  expect_error(
    reforge(ess_obj, n_segments = 0),
    "'n_segments' must be a single positive integer"
  )
  expect_error(
    reforge(ess_obj, shell_thickness = 1),
    "shell_thickness exceeds shell radius"
  )

  ess_thickness_only <- reforge(ess_obj, shell_thickness = 0.0015)
  expect_equal(
    max(extract(ess_thickness_only, "shape_parameters")$fluid$radius),
    max(extract(ess_thickness_only, "shape_parameters")$shell$radius) - 0.0015,
    tolerance = 1e-12
  )
})
