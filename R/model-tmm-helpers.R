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

# Enforce the current public TMM scope to homogeneous fluid/gas scatterers.
#' @noRd
.tmm_validate_object_scope <- function(object) {
  # Restrict the current TMM initializer to the supported scatterer classes ====
  if (!methods::is(object, "FLS") && !methods::is(object, "GAS")) {
    stop(
      "The current TMM implementation requires the scatterer to be either ",
      "'FLS' or 'GAS'. Input scatterer is type '", class(object)[1], "'.",
      call. = FALSE
    )
  }

  invisible(TRUE)
}

# Resolve the active geometric branch for one TMM initialization call.
#' @noRd
.tmm_branch_flags <- function(shape_parameters) {
  # Validate the supported geometry and identify the active backend ============
  .tmm_validate_shape(shape_parameters)
  use_spheroidal_branch <- .tmm_is_spheroidal_branch(shape_parameters)
  use_cylindrical_branch <- .tmm_is_cylindrical_branch(shape_parameters)
  if (use_cylindrical_branch &&
      isTRUE(getOption("acousticTS.warn_tmm_cylinder", TRUE))) {
    warning(
      "Cylinder 'TMM' support remains experimental. The default monostatic ",
      "branch is benchmark-matched to 'FCMS', but external BEM agreement is ",
      "not yet closed away from broadside, and stored cylinder TMM currently ",
      "supports only exact monostatic reuse plus orientation-averaged ",
      "monostatic products. Use 'FCMS' for production cylinder target ",
      "strength calculations when possible.",
      call. = FALSE
    )
  }

  list(
    use_spheroidal_branch = use_spheroidal_branch,
    use_cylindrical_branch = use_cylindrical_branch
  )
}

# Validate the public TMM boundary labels after applying the class default.
#' @noRd
.tmm_resolve_boundary <- function(object, boundary) {
  # Apply the class-specific default boundary and validate the final label =====
  boundary <- .tmm_boundary_default(object, boundary)
  if (!(boundary %in% c(
    "fixed_rigid",
    "pressure_release",
    "liquid_filled",
    "gas_filled"
  ))) {
    stop(
      "Only the following values for 'boundary' are available in TMM: ",
      "'fixed_rigid', 'pressure_release', 'liquid_filled', 'gas_filled'.",
      call. = FALSE
    )
  }

  boundary
}

# Validate the retained-block storage flag for one TMM initialization call.
#' @noRd
.tmm_validate_store_t_matrix <- function(store_t_matrix) {
  # Require an explicit scalar logical storage flag ============================
  if (!is.logical(store_t_matrix) || length(store_t_matrix) != 1 ||
      is.na(store_t_matrix)) {
    stop("'store_t_matrix' must be either TRUE or FALSE.", call. = FALSE)
  }

  invisible(TRUE)
}

# Disable unsupported spherical truncation overrides on the spheroidal branch.
#' @noRd
.tmm_branch_n_max <- function(n_max, use_spheroidal_branch) {
  # Ignore `n_max` on the prolate spheroidal branch ============================
  if (use_spheroidal_branch && !is.null(n_max)) {
    warning(
      "'n_max' is ignored for the current prolate-spheroidal TMM branch, ",
      "which uses a spheroidal-coordinate modal/T-matrix equivalent backend."
    )
    return(NULL)
  }

  n_max
}

# Validate the interior material properties required for penetrable TMM runs.
#' @noRd
.tmm_validate_penetrable_body <- function(body, boundary) {
  # Require scalar density and sound speed for penetrable targets ==============
  if (boundary %in% c("liquid_filled", "gas_filled") &&
      (is.null(body$density) || is.null(body$sound_speed))) {
    stop(
      "Penetrable TMM boundaries require scalar body density and sound speed ",
      "(or the corresponding scalar contrasts).",
      call. = FALSE
    )
  }

  invisible(TRUE)
}

# Resolve the homogeneous interior properties used by the TMM initializer.
#' @noRd
.tmm_prepare_body <- function(object, sound_speed_sw, density_sw, boundary) {
  # Complete the body properties against the surrounding medium ===============
  body <- .complete_material_props(
    acousticTS::extract(object, "body"),
    medium_sound_speed = sound_speed_sw,
    medium_density = density_sw
  )
  .tmm_require_homogeneous_body(body, boundary)
  .tmm_validate_penetrable_body(body, boundary)

  body
}

# Build the prolate-spheroidal geometry scalars used by the TMM backend.
#' @noRd
.tmm_spheroidal_geometry <- function(shape_parameters) {
  # Convert the stored prolate geometry to spheroidal coordinates =============
  xi <- 1 / sqrt(
    1 - (shape_parameters$radius / (shape_parameters$length / 2))^2
  )
  q <- shape_parameters$length / 2 / xi

  list(xi = xi, q = q)
}

# Build the common acoustics table shared by all TMM branches.
#' @noRd
.tmm_base_acoustics <- function(frequency, sound_speed_sw, body, boundary) {
  # Initialize the exterior and optional interior wavenumber columns ===========
  acoustics <- .init_acoustics_df(frequency, k_sw = sound_speed_sw)
  acoustics$k_body <- if (boundary %in% c("liquid_filled", "gas_filled")) {
    wavenumber(frequency, body$sound_speed)
  } else {
    NA_real_
  }

  acoustics
}

# Attach the spheroidal reduced-frequency bookkeeping to the acoustics table.
#' @noRd
.tmm_add_spheroidal_acoustics <- function(acoustics, boundary, geometry) {
  # Add reduced frequencies and spheroidal truncation heuristics ==============
  acoustics$chi_sw <- acoustics$k_sw * geometry$q
  acoustics$chi_body <- if (boundary %in% c("liquid_filled", "gas_filled")) {
    acoustics$k_body * geometry$q
  } else {
    acoustics$k_sw * geometry$q
  }
  acoustics$m_max <- ceiling(2 * acoustics$k_sw * geometry$radius)
  acoustics$n_max <- acoustics$m_max + ceiling(0.5 * acoustics$chi_sw)

  acoustics
}

# Apply the branch-specific truncation rules to the TMM acoustics table.
#' @noRd
.tmm_assign_branch_n_max <- function(acoustics,
                                     n_max,
                                     frequency,
                                     shape_parameters,
                                     boundary,
                                     use_spheroidal_branch,
                                     use_cylindrical_branch) {
  # Preserve the spheroidal truncation already attached to the acoustics table =
  if (use_spheroidal_branch) {
    return(acoustics)
  }

  # Apply the cylinder or spherical truncation rules as appropriate ============
  acoustics$n_max <- if (use_cylindrical_branch) {
    .tmm_prepare_cylinder_n_max(
      n_max = n_max,
      frequency = frequency,
      k_sw = acoustics$k_sw,
      shape_parameters = shape_parameters
    )
  } else {
    .tmm_prepare_n_max(
      n_max = n_max,
      frequency = frequency,
      k_sw = acoustics$k_sw,
      shape_parameters = shape_parameters,
      boundary = boundary
    )
  }

  acoustics
}

# Build the branch-specific acoustics table and optional spheroidal geometry.
#' @noRd
.tmm_prepare_acoustics <- function(frequency,
                                   sound_speed_sw,
                                   body,
                                   boundary,
                                   shape_parameters,
                                   use_spheroidal_branch,
                                   use_cylindrical_branch,
                                   n_max) {
  # Start from the shared wavenumber table ====================================
  acoustics <- .tmm_base_acoustics(frequency, sound_speed_sw, body, boundary)
  geometry <- NULL
  if (use_spheroidal_branch) {
    geometry <- .tmm_spheroidal_geometry(shape_parameters)
    geometry$radius <- shape_parameters$radius
    acoustics <- .tmm_add_spheroidal_acoustics(acoustics, boundary, geometry)
  }
  acoustics <- .tmm_assign_branch_n_max(
    acoustics = acoustics,
    n_max = n_max,
    frequency = frequency,
    shape_parameters = shape_parameters,
    boundary = boundary,
    use_spheroidal_branch = use_spheroidal_branch,
    use_cylindrical_branch = use_cylindrical_branch
  )

  list(acoustics = acoustics, geometry = geometry)
}

# Build the stored body metadata used by the TMM runtime helpers.
#' @noRd
.tmm_body_parameters <- function(body, geometry) {
  # Start from the common body-fixed incident and material properties ==========
  body_params <- list(
    theta_body = body$theta,
    theta_scatter = pi - body$theta,
    phi_body = pi,
    phi_scatter = 2 * pi,
    density = body$density,
    sound_speed = body$sound_speed,
    g_body = body$g,
    h_body = body$h
  )
  if (!is.null(geometry)) {
    body_params$xi <- geometry$xi
    body_params$q <- geometry$q
  }

  body_params
}

# Resolve the coordinate-system label stored with the initialized TMM object.
#' @noRd
.tmm_coordinate_system <- function(use_spheroidal_branch, use_cylindrical_branch) {
  # Label the active TMM backend coordinate system ============================
  if (use_spheroidal_branch) {
    return("spheroidal")
  }
  if (use_cylindrical_branch) {
    return("cylindrical")
  }

  "spherical"
}

# Resolve the stored precision label for one initialized TMM object.
#' @noRd
.tmm_precision_label <- function(use_spheroidal_branch, boundary) {
  # Penetrable spheroidal runs use the quad-precision backend =================
  if (use_spheroidal_branch && boundary %in% c("liquid_filled", "gas_filled")) {
    return("quad")
  }

  "double"
}

# Resolve the stored quadrature-node count for the spheroidal TMM branch.
#' @noRd
.tmm_n_integration_label <- function(use_spheroidal_branch, boundary) {
  # Only penetrable spheroidal runs currently keep an explicit quadrature count
  if (use_spheroidal_branch && boundary %in% c("liquid_filled", "gas_filled")) {
    return(96L)
  }

  NA_integer_
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
.tmm_prepare_n_max <- function(
    n_max, frequency, k_sw, shape_parameters, boundary
  ) {
  if (is.null(n_max)) {
    return(.tmm_default_n_max(k_sw, shape_parameters, boundary))
  }

  if (!is.numeric(n_max) || any(!is.finite(n_max)) || any(n_max < 1)) {
    stop("'n_max' must be NULL or a positive finite integer vector.",
         call. = FALSE)
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
    return(as.integer(ceiling(k_sw * max(as.numeric(shape_parameters$radius),
                                         na.rm = TRUE)) + 10L))
  }

  if (!is.numeric(n_max) || any(!is.finite(n_max)) || any(n_max < 1)) {
    stop("'n_max' must be NULL or a positive finite integer vector.",
         call. = FALSE)
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
  k1a <- acoustics$k_sw * sin(body$theta_body) *
    max(as.numeric(shape_parameters$radius), na.rm = TRUE)
  k2a <- if (boundary %in% c("liquid_filled", "gas_filled")) {
    acoustics$k_sw * sin(body$theta_body) / body$h_body *
      max(as.numeric(shape_parameters$radius), na.rm = TRUE)
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
