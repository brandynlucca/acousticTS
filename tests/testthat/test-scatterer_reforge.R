library(acousticTS)

test_that("`.discover_reforge_params` returns expected keys.", {
  # Test FLS parameters
  expect_equal(
    .discover_reforge_params("FLS"),
    c("length", "radius", "length_radius_ratio_constant", "n_segments")
  )

  # Test SBF parameters
  expect_equal(
    .discover_reforge_params("SBF"),
    c(
      "body_scale", "swimbladder_scale", "body_target",
      "swimbladder_target", "swimbladder_inflation_factor",
      "isometric_body", "isometric_swimbladder", "maintain_ratio",
      "n_segments_body", "n_segments_swimbladder"
    )
  )

  # Test erroneous input
  expect_equal(
    .discover_reforge_params("OTHER"),
    character(0)
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
  suppressMessages({
    body_expand <- reforge(sardine, body_scale = c(length = 4, width = 1.5))
  })
  expect_message(
    reforge(sardine, body_scale = c(length = 4, width = 1.5)),
    paste0(
      "Multiple axes specified in body_scale/body_target: 'isometric_body' ",
      "will be ignored for those axes.."
    )
  )

  new_length <- extract(body_expand, "shape_parameters")$body$length
  new_width <- max(extract(body_expand, "body")$rpos[2, ])
  z_pos <- extract(body_expand, "body")$rpos[3:4, ]
  new_height <- max(head(z_pos, 1) - tail(z_pos, 1))
  expect_equal(new_length, orig_length * 4)
  expect_equal(new_width, orig_width * 1.5)
  expect_equal(new_height, orig_height)
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
  expect_equal(new_sbheight, orig_sbheight)
  expect_equal(
    extract(body_expand, "shape_parameters")$bladder$length,
    new_sblength
  )
  # ---- Case: Defining length, width, and height -> overrides `isometric_body`
  suppressMessages({
    body_expand <- reforge(sardine, body_scale = c(
      length = 4,
      width = 1.5, height = 2.0
    ))
  })
  expect_message(
    reforge(sardine, body_scale = c(length = 4, width = 1.5)),
    paste0(
      "Multiple axes specified in body_scale/body_target: 'isometric_body' ",
      "will be ignored for those axes."
    )
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
  # ---- Case: Repeat previous but drop 'isometric_body'
  suppressMessages({
    body_expand <- reforge(sardine,
      body_scale = c(length = 4, width = 1.5, height = 2.0),
      isometric_body = FALSE
    )
  })
  expect_no_message(
    reforge(sardine,
      body_scale = c(length = 4, width = 1.5, height = 2.0),
      isometric_body = FALSE
    )
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
  suppressMessages({
    body_expand <- reforge(sardine,
      body_scale = c(length = 4, width = 1.5, height = 2.0),
      isometric_body = FALSE,
      maintain_ratio = FALSE
    )
  })
  expect_no_message(
    reforge(sardine,
      body_scale = c(length = 4, width = 1.5, height = 2.0),
      isometric_body = FALSE,
      maintain_ratio = FALSE
    )
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
  suppressMessages({
    bladder_expand <- reforge(sardine,
      swimbladder_scale = c(length = 4, width = 1.5)
    )
  })
  expect_message(
    reforge(sardine,
      swimbladder_scale = c(length = 4, width = 1.5)
    ),
    paste0(
      "Multiple axes specified in swimbladder_scale/swimbladder_target: ",
      "'isometric_swimbladder' will be ignored for those axes."
    )
  )

  new_length <- extract(bladder_expand, "shape_parameters")$body$length
  new_width <- max(extract(bladder_expand, "body")$rpos[2, ])
  z_pos <- extract(bladder_expand, "body")$rpos[3:4, ]
  new_height <- max(head(z_pos, 1) - tail(z_pos, 1))
  expect_equal(new_length, orig_length * 4)
  expect_equal(new_width, orig_width * 1.5)
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
  expect_equal(new_sbheight, orig_sbheight)
  expect_equal(
    extract(bladder_expand, "shape_parameters")$bladder$length,
    new_sblength
  )
  # ---- Case: Defining length, width, and height -> overrides `isometric_body`
  suppressMessages({
    bladder_expand <- reforge(sardine,
      swimbladder_scale = c(
        length = 4,
        width = 1.5,
        height = 2.0
      )
    )
  })
  expect_message(
    reforge(sardine,
      swimbladder_scale = c(length = 4, width = 1.5, height = 2.0)
    ),
    paste0(
      "Multiple axes specified in swimbladder_scale/swimbladder_target: ",
      "'isometric_swimbladder' will be ignored for those axes."
    )
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
  # ---- Case: Repeat previous but drop 'isometric_swimbladder'
  suppressMessages({
    bladder_expand <- reforge(sardine,
      swimbladder_scale = c(
        length = 4,
        width = 1.5,
        height = 2.0
      ),
      isometric_swimbladder = FALSE
    )
  })
  expect_no_message(
    reforge(sardine,
      swimbladder_scale = c(length = 4, width = 1.5, height = 2.0),
      isometric_swimbladder = FALSE
    )
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
    bladder_expand <- reforge(sardine,
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
    "body_target must be a named numeric vector."
  )
  expect_error(
    reforge(sardine, swimbladder_target = 1),
    "swimbladder_target must be a named numeric vector."
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
  suppressMessages({
    body_expand <- reforge(sardine, body_target = c(length = 0.8, width = 0.4))
  })
  expect_message(
    reforge(sardine, body_target = c(length = 0.8, width = 0.4)),
    paste0(
      "Multiple axes specified in body_scale/body_target: 'isometric_body' ",
      "will be ignored for those axes.."
    )
  )

  new_length <- extract(body_expand, "shape_parameters")$body$length
  new_width <- max(extract(body_expand, "body")$rpos[2, ])
  scaling_width_ratio <- new_width / orig_width
  z_pos <- extract(body_expand, "body")$rpos[3:4, ]
  new_height <- max(head(z_pos, 1) - tail(z_pos, 1))
  expect_equal(new_length, orig_length * scaling_length_ratio)
  expect_equal(new_width, orig_width * scaling_width_ratio)
  expect_equal(new_height, orig_height)
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
  expect_equal(new_sbheight, orig_sbheight)
  expect_equal(
    extract(body_expand, "shape_parameters")$bladder$length,
    new_sblength
  )
  # # ---- Case: Defining length/width/height -> overrides `isometric_body`
  suppressMessages({
    body_expand <- reforge(sardine,
      body_target = c(
        length = 0.8,
        width = 0.4,
        height = 0.1
      )
    )
  })
  expect_message(
    reforge(sardine, body_target = c(length = 0.8, width = 0.4, height = 0.1)),
    paste0(
      "Multiple axes specified in body_scale/body_target: 'isometric_body' ",
      "will be ignored for those axes."
    )
  )

  new_length <- extract(body_expand, "shape_parameters")$body$length
  new_width <- max(extract(body_expand, "body")$rpos[2, ])
  z_pos <- extract(body_expand, "body")$rpos[3:4, ]
  new_height <- max(head(z_pos, 1) - tail(z_pos, 1))
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
  # ---- Case: Repeat previous but drop 'isometric_body'
  suppressMessages({
    body_expand <- reforge(sardine,
      body_target = c(length = 0.8, width = 0.4, height = 0.1),
      isometric_body = FALSE
    )
  })
  expect_no_message(
    reforge(sardine,
      body_target = c(length = 0.8, width = 0.4, height = 0.1),
      isometric_body = FALSE
    )
  )

  new_length <- extract(body_expand, "shape_parameters")$body$length
  new_width <- max(extract(body_expand, "body")$rpos[2, ])
  z_pos <- extract(body_expand, "body")$rpos[3:4, ]
  new_height <- max(head(z_pos, 1) - tail(z_pos, 1))
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
  suppressMessages({
    body_expand <- reforge(sardine,
      body_target = c(length = 0.8, width = 0.4, height = 0.1),
      isometric_body = FALSE,
      maintain_ratio = FALSE
    )
  })
  expect_no_message(
    reforge(sardine,
      body_target = c(length = 0.8, width = 0.4, height = 0.1),
      isometric_body = FALSE,
      maintain_ratio = FALSE
    )
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
  suppressMessages({
    bladder_expand <- reforge(sardine,
      swimbladder_target = c(length = 0.2, width = 0.3)
    )
  })
  expect_message(
    reforge(sardine,
      swimbladder_target = c(length = 0.2, width = 0.3)
    ),
    paste0(
      "Multiple axes specified in swimbladder_scale/swimbladder_target: ",
      "'isometric_swimbladder' will be ignored for those axes."
    )
  )

  new_length <- extract(bladder_expand, "shape_parameters")$body$length
  new_width <- max(extract(bladder_expand, "body")$rpos[2, ])
  z_pos <- extract(bladder_expand, "body")$rpos[3:4, ]
  new_height <- max(head(z_pos, 1) - tail(z_pos, 1))
  expect_equal(new_length, orig_length * scaling_length_ratio)
  expect_equal(new_width, orig_width * scaling_width_ratio)
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
  expect_equal(new_sbheight, orig_sbheight)
  expect_equal(
    extract(bladder_expand, "shape_parameters")$bladder$length,
    new_sblength
  )
  # ---- Case: Defining length, width, and height -> overrides `isometric_body`
  suppressMessages({
    bladder_expand <- reforge(sardine,
      swimbladder_target = c(
        length = 0.2,
        width = 0.3,
        height = 0.1
      )
    )
  })
  expect_message(
    reforge(sardine,
      swimbladder_target = c(length = 0.2, width = 0.3, height = 0.1)
    ),
    paste0(
      "Multiple axes specified in swimbladder_scale/swimbladder_target: ",
      "'isometric_swimbladder' will be ignored for those axes."
    )
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
  # ---- Case: Repeat previous but drop 'isometric_swimbladder'
  suppressMessages({
    bladder_expand <- reforge(sardine,
      swimbladder_target = c(
        length = 0.2,
        width = 0.3,
        height = 0.1
      ),
      isometric_swimbladder = FALSE
    )
  })
  expect_no_message(
    reforge(sardine,
      swimbladder_target = c(length = 0.2, width = 0.3, height = 0.1),
      isometric_swimbladder = FALSE
    )
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
    "Invalid body_scale names. Use: length, width, height."
  )
  expect_error(
    reforge(sardine, body_target = c(invalid_dim = 2)),
    "Invalid body_target names. Use: length, width, height."
  )
  # ---- Case: Multiple invalid dimensions
  expect_error(
    reforge(
      sardine,
      body_scale = c(length = 2, invalid_dim1 = 2, invalid_dim2 = 2),
      isometric_body = FALSE
    ),
    "Invalid body_scale names. Use: length, width, height."
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
      "body_scale must be a single number or named vector with ",
      "length/width/height"
    )
  )
  # ---- Case: Unspecified swimbladder dimensions
  expect_error(
    reforge(sardine,
      swimbladder_scale = c(1, 2, 3),
      isometric_swimbladder = FALSE
    ),
    paste0(
      "swimbladder_scale must be a single number or named vector with ",
      "length/width/height"
    )
  )
  # ---- Case: Unspecified body scale without isometry
  expect_error(
    reforge(sardine, body_scale = 1, isometric_body = FALSE),
    paste0(
      "body_scale must be a named vector with length/width/height when ",
      "isometric_body is FALSE."
    )
  )
  # ---- Case: Unspecified swimbladder scale without isometry
  expect_error(
    reforge(sardine, swimbladder_scale = 1, isometric_swimbladder = FALSE),
    paste0(
      "swimbladder_scale must be a named vector with length/width/height ",
      "when isometric_swimbladder is FALSE."
    )
  )
})
