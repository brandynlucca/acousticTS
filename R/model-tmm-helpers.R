################################################################################
# Transition matrix method (TMM) common helpers
################################################################################

# Resolve the default boundary from the scatterer class when the user omits it.
#' @noRd
.tmm_boundary_default <- function(object, boundary) {
  if (!is.null(boundary)) {
    return(boundary)
  }

  switch(
    class(object)[1],
    GAS = "gas_filled",
    FLS = "liquid_filled",
    stop(
      "Specify 'boundary' explicitly for TMM when the scatterer is not a ",
      "'FLS' or 'GAS' object.",
      call. = FALSE
    )
  )
}

# Guard the current public scope of TMM to the exact sphere and prolate
# spheroid geometries already supported by the package.
#' @noRd
.tmm_validate_shape <- function(shape_parameters) {
  supported <- c("Sphere", "ProlateSpheroid", "OblateSpheroid", "Cylinder")
  if (!(shape_parameters$shape %in% supported)) {
    stop(
      "The current TMM implementation supports the following shape-types: ",
      paste0("'", supported, "'", collapse = ", "),
      ". Input scatterer is shape-type '", shape_parameters$shape, "'.",
      call. = FALSE
    )
  }
}

# Decide whether the object should be routed through the spheroidal-coordinate
# backend rather than the spherical-coordinate backend.
#' @noRd
.tmm_is_spheroidal_branch <- function(shape_parameters) {
  identical(shape_parameters$shape, "ProlateSpheroid")
}

# Cylinders now have a dedicated geometry-matched cylindrical backend rather
# than the exploratory spherical stored-block branch.
#' @noRd
.tmm_is_cylindrical_branch <- function(shape_parameters) {
  identical(shape_parameters$shape, "Cylinder")
}

# The present penetrable implementation assumes a single homogeneous interior,
# so vector-valued contrasts/properties must be rejected here.
#' @noRd
.tmm_require_homogeneous_body <- function(body, boundary) {
  if (!(boundary %in% c("liquid_filled", "gas_filled"))) {
    return(invisible(TRUE))
  }

  for (nm in c("g", "h", "density", "sound_speed")) {
    value <- body[[nm]]
    if (!is.null(value) && length(value) > 1) {
      stop(
        "TMM currently assumes a homogeneous penetrable interior. The body ",
        "property '", nm, "' must therefore be scalar for boundary '",
        boundary, "'.",
        call. = FALSE
      )
    }
  }

  invisible(TRUE)
}

# Return the outer size scale used in the default spherical truncation rule.
#' @noRd
.tmm_bounding_radius <- function(shape_parameters) {
  switch(
    shape_parameters$shape,
    Sphere = as.numeric(shape_parameters$radius)[1],
    ProlateSpheroid = as.numeric(shape_parameters$semimajor_length)[1]
    ,
    OblateSpheroid = as.numeric(shape_parameters$semimajor_length)[1],
    Cylinder = sqrt(
      (as.numeric(shape_parameters$length)[1] / 2)^2 +
        max(as.numeric(shape_parameters$radius), na.rm = TRUE)^2
    )
  )
}

# Apply conservative hard overrides for slender prolates when a spherical-wave
# truncation is still needed.
#' @noRd
.tmm_prolate_nmax_override <- function(shape_parameters, boundary) {
  if (shape_parameters$shape != "ProlateSpheroid") {
    return(NA_integer_)
  }

  a <- as.numeric(shape_parameters$semimajor_length)[1]
  b <- as.numeric(shape_parameters$semiminor_length)[1]
  aspect_ratio <- a / b

  if (!is.finite(aspect_ratio) || aspect_ratio < 3) {
    return(NA_integer_)
  }

  switch(
    boundary,
    fixed_rigid = 30L,
    pressure_release = 24L,
    liquid_filled = 36L,
    gas_filled = 36L,
    NA_integer_
  )
}

# Cylinders are the least smooth geometry currently handled by the
# spherical-coordinate branch, so they benefit from a more conservative modal
# floor than the generic bounding-sphere rule at low and moderate ka.
#' @noRd
.tmm_cylinder_nmax_floor <- function(shape_parameters, boundary) {
  if (shape_parameters$shape != "Cylinder") {
    return(NA_integer_)
  }

  switch(
    boundary,
    fixed_rigid = 40L,
    pressure_release = 24L,
    liquid_filled = 36L,
    gas_filled = 24L,
    NA_integer_
  )
}

# Default truncation based on the exterior size parameter, with geometry-
# specific overrides when they are known to behave better.
#' @noRd
.tmm_default_n_max <- function(k_sw, shape_parameters, boundary) {
  n_override <- .tmm_prolate_nmax_override(shape_parameters, boundary)
  if (!is.na(n_override)) {
    return(rep.int(as.integer(n_override), length(k_sw)))
  }

  x <- abs(k_sw) * .tmm_bounding_radius(shape_parameters)
  n_default <- pmax(4L, as.integer(ceiling(x + 4 * x^(1 / 3) + 2)))
  n_floor <- .tmm_cylinder_nmax_floor(shape_parameters, boundary)
  if (!is.na(n_floor)) {
    n_default <- pmax(n_default, as.integer(n_floor))
  }
  n_default
}

# Normalize `n_max` into a per-frequency integer vector.
#' @noRd
.tmm_prepare_n_max <- function(n_max, frequency, k_sw, shape_parameters, boundary) {
  if (is.null(n_max)) {
    return(.tmm_default_n_max(k_sw, shape_parameters, boundary))
  }

  if (!is.numeric(n_max) || any(!is.finite(n_max)) || any(n_max < 1)) {
    stop("'n_max' must be NULL or a positive finite integer vector.", call. = FALSE)
  }

  if (length(n_max) == 1) {
    return(rep(as.integer(n_max), length(frequency)))
  }

  if (length(n_max) != length(frequency)) {
    stop(
      "'n_max' must be length 1 or match the length of 'frequency'.",
      call. = FALSE
    )
  }

  as.integer(n_max)
}

# Cylinder monostatic TMM uses the same modal truncation rule as FCMS so that
# the cylindrical-coordinate backend remains benchmark-compatible by default.
#' @noRd
.tmm_prepare_cylinder_n_max <- function(n_max, frequency, k_sw, shape_parameters) {
  if (is.null(n_max)) {
    return(as.integer(ceiling(k_sw * max(as.numeric(shape_parameters$radius), na.rm = TRUE)) + 10L))
  }

  if (!is.numeric(n_max) || any(!is.finite(n_max)) || any(n_max < 1)) {
    stop("'n_max' must be NULL or a positive finite integer vector.", call. = FALSE)
  }

  if (length(n_max) == 1) {
    return(rep(as.integer(n_max), length(frequency)))
  }

  if (length(n_max) != length(frequency)) {
    stop(
      "'n_max' must be length 1 or match the length of 'frequency'.",
      call. = FALSE
    )
  }

  as.integer(n_max)
}

# Shared collocation-node rule for the spherical-coordinate branch.
#' @noRd
.tmm_collocation_nodes <- function(shape_parameters, boundary, n_terms) {
  if (identical(shape_parameters$shape, "Cylinder") &&
      identical(boundary, "fixed_rigid")) {
    return(max(128L, 10L * n_terms))
  }

  max(64L, 4L * n_terms)
}

# Flatten the supported shape parameters into the compact numeric vector passed
# into the compiled spherical backend.
#' @noRd
.tmm_shape_values <- function(shape_parameters) {
  switch(
    shape_parameters$shape,
    Sphere = c(as.numeric(shape_parameters$radius)[1]),
    ProlateSpheroid = c(
      as.numeric(shape_parameters$semimajor_length)[1],
      as.numeric(shape_parameters$semiminor_length)[1]
    ),
    OblateSpheroid = c(
      as.numeric(shape_parameters$semiminor_length)[1],
      as.numeric(shape_parameters$semimajor_length)[1]
    ),
    Cylinder = c(
      as.numeric(shape_parameters$length)[1] / 2,
      max(as.numeric(shape_parameters$radius), na.rm = TRUE)
    ),
    stop("Unsupported TMM shape geometry.", call. = FALSE)
  )
}

# Evaluate the exact finite-cylinder modal coefficients inside the TMM cylinder
# branch so that the default monostatic cylinder response shares the same
# geometry-matched backend as FCMS.
#' @noRd
.tmm_run_cylindrical_branch <- function(shape_parameters,
                                        acoustics,
                                        body,
                                        boundary) {
  m_max <- max(acoustics$n_max)
  nu <- neumann(0:m_max)
  k1L <- shape_parameters$length * acoustics$k_sw
  k1a <- acoustics$k_sw * sin(body$theta_body) * max(as.numeric(shape_parameters$radius), na.rm = TRUE)
  k2a <- if (boundary %in% c("liquid_filled", "gas_filled")) {
    acoustics$k_sw * sin(body$theta_body) / body$h_body * max(as.numeric(shape_parameters$radius), na.rm = TRUE)
  } else {
    NA_real_
  }
  gh <- body$g_body * body$h_body

  Bm <- switch(
    boundary,
    liquid_filled = .fcms_bm_fluid(k1a, k2a, gh, nu, acoustics$n_max),
    gas_filled = .fcms_bm_fluid(k1a, k2a, gh, nu, acoustics$n_max),
    fixed_rigid = .fcms_bm_fixed_rigid(k1a, nu, acoustics$n_max),
    pressure_release = .fcms_bm_pressure_release(k1a, nu, acoustics$n_max),
    stop("Unsupported boundary for cylindrical TMM branch.", call. = FALSE)
  )

  if (!is.matrix(Bm)) {
    Bm <- t(as.matrix(Bm))
  }

  if (boundary %in% c("liquid_filled", "gas_filled")) {
    f_bs <- -shape_parameters$length / pi *
      sin(k1L * cos(body$theta_body)) / (k1L * cos(body$theta_body)) *
      colSums(Bm, na.rm = TRUE)
  } else {
    f_bs <- 1i * shape_parameters$length / pi *
      sin(k1L * cos(body$theta_body)) / (k1L * cos(body$theta_body)) *
      colSums(Bm, na.rm = TRUE)
  }

  list(
    model = data.frame(
      frequency = acoustics$frequency,
      f_bs = f_bs,
      sigma_bs = .sigma_bs(f_bs),
      TS = db(.sigma_bs(f_bs)),
      n_max = acoustics$n_max
    )
  )
}

# Build one lightweight retained-state placeholder for the cylindrical branch.
# The current cylindrical TMM post-processing reuses the exact-family modal
# formulas directly rather than a dense stored matrix, so the retained state is
# just enough to mark that this frequency was solved in the cylindrical family.
#' @noRd
.tmm_store_cylindrical_branch <- function(acoustics) {
  lapply(
    seq_len(nrow(acoustics)),
    function(i) {
      list(
        family = "cylindrical_mono",
        n_max = acoustics$n_max[i]
      )
    }
  )
}
