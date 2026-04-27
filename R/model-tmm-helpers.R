################################################################################
# Transition matrix method (TMM) common helpers
################################################################################

# Resolve the default boundary from the scatterer class when the user omits it.
#' @noRd
.tmm_boundary_default <- function(object, boundary) {
  if (!is.null(boundary)) {
    return(boundary)
  }

  switch(class(object)[1],
    ELA = "elastic",
    CAL = "elastic",
    GAS = "gas_filled",
    FLS = "liquid_filled",
    ESS = {
      shape_parameters <- acousticTS::extract(object, "shape_parameters")
      shell <- acousticTS::extract(object, "shell")
      shape_name <- as.character(shape_parameters[["shape"]])[1]
      if (identical(shape_name, "Sphere")) {
        if (.tmm_has_elastic_shell(shell)) {
          "elastic_shelled"
        } else {
          "shelled_liquid"
        }
      } else {
        stop(
          "Specify 'boundary' explicitly for TMM when the scatterer is not a ",
          "supported homogeneous fluid/gas sphere, spheroid, or supported ",
          "spherical shell ESS object.",
          call. = FALSE
        )
      }
    },
    stop(
      "Specify 'boundary' explicitly for TMM when the scatterer is not a ",
      "'FLS', 'GAS', 'ELA', or supported shell-sphere 'ESS' object.",
      call. = FALSE
    )
  )
}

# Guard the current public scope of TMM to the exact sphere and prolate
# spheroid geometries already supported by the package.
#' @noRd
.tmm_shape_name <- function(shape_parameters) {
  as.character(shape_parameters$shape)[1]
}

# Guard the current public scope of TMM to the exact sphere and prolate
# spheroid geometries already supported by the package.
#' @noRd
.tmm_validate_shape <- function(shape_parameters) {
  supported <- c("Sphere", "ProlateSpheroid", "OblateSpheroid", "Cylinder")
  shape_name <- .tmm_shape_name(shape_parameters)
  if (!(shape_name %in% supported)) {
    stop(
      "The current TMM implementation supports the following shape-types: ",
      paste0("'", supported, "'", collapse = ", "),
      ". Input scatterer is shape-type '", shape_name, "'.",
      call. = FALSE
    )
  }
}

# Decide whether the object should be routed through the spheroidal-coordinate
# backend rather than the spherical-coordinate backend.
#' @noRd
.tmm_is_spheroidal_branch <- function(shape_parameters, boundary = NULL) {
  identical(.tmm_shape_name(shape_parameters), "ProlateSpheroid") &&
    !identical(boundary, "elastic_shelled") &&
    !identical(boundary, "elastic")
}

# Identify spherical fluid-shell runs that should stay on the exact sphere
# modal path rather than the projected spherical T-matrix solve.
#' @noRd
.tmm_is_shell_sphere_branch <- function(object, shape_parameters, boundary) {
  methods::is(object, "ESS") &&
    identical(as.character(shape_parameters[["shape"]])[1], "Sphere") &&
    boundary %in% c(
      "shelled_pressure_release",
      "shelled_liquid",
      "shelled_gas"
    )
}

# Identify spherical elastic-shell runs that should stay on the exact elastic
# spherical modal path rather than the projected spherical T-matrix solve.
#' @noRd
.tmm_is_elastic_shell_sphere_branch <- function(object, shape_parameters, boundary) {
  methods::is(object, "ESS") &&
    identical(as.character(shape_parameters[["shape"]])[1], "Sphere") &&
    identical(boundary, "elastic_shelled")
}

# Identify elastic-solid TMM runs carried by plain ELA or CAL objects.
#' @noRd
.tmm_is_elastic_solid_branch <- function(object, boundary) {
  identical(boundary, "elastic") &&
    methods::is(object, "ELA") &&
    !methods::is(object, "ESS")
}

# Identify solid-elastic sphere runs that should stay on the exact calibration
# sphere modal path while still using the public TMM interface.
#' @noRd
.tmm_is_elastic_solid_sphere_branch <- function(object, shape_parameters, boundary) {
  .tmm_is_elastic_solid_branch(object, boundary) &&
    identical(as.character(shape_parameters[["shape"]])[1], "Sphere")
}

# Identify solid-elastic prolate runs carried by the public TMM interface.
#' @noRd
.tmm_is_elastic_solid_prolate_branch <- function(object, shape_parameters, boundary) {
  .tmm_is_elastic_solid_branch(object, boundary) &&
    identical(as.character(shape_parameters[["shape"]])[1], "ProlateSpheroid")
}

# Detect whether an ESS shell carries elastic properties.
#' @noRd
.tmm_has_elastic_shell <- function(shell) {
  any(vapply(
    c("E", "G", "K", "nu", "lambda"),
    function(nm) {
      value <- shell[[nm]]
      !is.null(value) && length(value) == 1L && is.finite(value)
    },
    logical(1)
  ))
}

# Identify all sphere-modal TMM branches that retain explicit spherical modal
# coefficients for downstream bistatic reconstruction.
#' @noRd
.tmm_is_sphere_modal_branch <- function(object, shape_parameters, boundary) {
  .tmm_is_shell_sphere_branch(object, shape_parameters, boundary) ||
    .tmm_is_elastic_shell_sphere_branch(object, shape_parameters, boundary)
}

# Resolve the active cylinder backend. The legacy branch keeps the
# FCMS-matched monostatic family, while the retained branch reuses the
# cylinder-native retained operator for stored general-angle work.
#' @noRd
.tmm_cylinder_backend <- function(shape_parameters = NULL,
                                  cylinder_backend = NULL,
                                  cylinder_endcap_fraction = NULL,
                                  boundary = NULL,
                                  store_t_matrix = FALSE) {
  if (!is.null(cylinder_backend)) {
    backend <- match.arg(
      cylinder_backend,
      c("auto", "legacy", "retained", "spherical")
    )
    if (identical(backend, "spherical")) {
      backend <- "retained"
    }
    if (!identical(backend, "auto")) {
      return(backend)
    }
  } else {
    backend <- match.arg(
      getOption("acousticTS.tmm_cylinder_backend", "legacy"),
      c("auto", "legacy", "retained", "spherical")
    )
    if (identical(backend, "spherical")) {
      backend <- "retained"
    }

    if (!identical(backend, "auto")) {
      return(backend)
    }
  }

  if (is.null(shape_parameters) || !identical(.tmm_shape_name(shape_parameters), "Cylinder")) {
    return("legacy")
  }

  if (isTRUE(store_t_matrix)) {
    return("retained")
  }

  taper_order <- if ("taper_order" %in% names(shape_parameters)) {
    as.numeric(shape_parameters$taper_order)[1]
  } else {
    NA_real_
  }

  if (is.finite(taper_order)) {
    return("spherical")
  }

  if (!is.null(cylinder_endcap_fraction) &&
    .tmm_resolve_cylinder_endcap_fraction(cylinder_endcap_fraction) > 0) {
    return("spherical")
  }

  "legacy"
}

# Cylinders keep a dedicated exact-family monostatic branch for plain
# `target_strength()` runs. Stored-cylinder TMM workflows keep that same
# cylindrical-family monostatic backend so exact monostatic reuse and
# orientation-averaged monostatic products stay tied to the FCMS-equivalent
# reference path.
#' @noRd
.tmm_is_cylindrical_branch <- function(shape_parameters,
                                       boundary = NULL,
                                       cylinder_backend = NULL,
                                       cylinder_endcap_fraction = NULL,
                                       store_t_matrix = FALSE) {
  identical(.tmm_shape_name(shape_parameters), "Cylinder") &&
    identical(
      .tmm_cylinder_backend(
        shape_parameters = shape_parameters,
        cylinder_backend = cylinder_backend,
        cylinder_endcap_fraction = cylinder_endcap_fraction,
        boundary = boundary,
        store_t_matrix = store_t_matrix
      ),
      "legacy"
    )
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
  if (!methods::is(object, "FLS") &&
    !methods::is(object, "GAS") &&
    !methods::is(object, "ELA")) {
    stop(
      "The current TMM implementation requires the scatterer to be either ",
      "'FLS', 'GAS', 'ELA', or a supported 'ESS'. Input scatterer is ",
      "type '", class(object)[1], "'.",
      call. = FALSE
    )
  }

  invisible(TRUE)
}

# Resolve the active geometric branch for one TMM initialization call.
#' @noRd
.tmm_branch_flags <- function(shape_parameters,
                              boundary,
                              store_t_matrix = FALSE,
                              cylinder_backend = NULL,
                              cylinder_endcap_fraction = NULL) {
  # Validate the supported geometry and identify the active backend ============
  .tmm_validate_shape(shape_parameters)
  use_spheroidal_branch <- .tmm_is_spheroidal_branch(shape_parameters, boundary)
  use_cylindrical_branch <- .tmm_is_cylindrical_branch(
    shape_parameters = shape_parameters,
    boundary = boundary,
    cylinder_backend = cylinder_backend,
    cylinder_endcap_fraction = cylinder_endcap_fraction,
    store_t_matrix = store_t_matrix
  )
  cylinder_backend <- .tmm_cylinder_backend(
    shape_parameters = shape_parameters,
    cylinder_backend = cylinder_backend,
    cylinder_endcap_fraction = cylinder_endcap_fraction,
    boundary = boundary,
    store_t_matrix = store_t_matrix
  )
  if (identical(.tmm_shape_name(shape_parameters), "Cylinder") &&
    isTRUE(getOption("acousticTS.warn_tmm_cylinder", TRUE))) {
    if (identical(cylinder_backend, "legacy")) {
      warning(
        "Cylinder 'TMM' support remains experimental. Plain monostatic ",
        "'target_strength()' calls and stored-cylinder monostatic reuse keep the ",
        "exact FCMS-matched cylinder branch. Public general-angle cylinder ",
        "bistatic and grid post-processing remain outside the current scope.",
        call. = FALSE
      )
    } else {
      warning(
        "Cylinder 'TMM' support remains experimental. The retained cylinder ",
        "backend is enabled, so stored or tapered cylinder TMM uses the ",
        "cylinder-native retained operator instead of the legacy FCMS-",
        "matched monostatic-only branch. Public general-angle cylinder ",
        "bistatic and grid post-processing remain outside the current scope.",
        call. = FALSE
      )
    }
  }

  list(
    use_spheroidal_branch = use_spheroidal_branch,
    use_cylindrical_branch = use_cylindrical_branch,
    cylinder_backend = cylinder_backend
  )
}

# Validate the public TMM boundary labels after applying the class default.
#' @noRd
.tmm_resolve_boundary <- function(object, boundary) {
  # Apply the class-specific default boundary and validate the final label =====
  boundary <- .tmm_boundary_default(object, boundary)
  if (methods::is(object, "ELA") &&
    !methods::is(object, "ESS") &&
    identical(boundary, "elastic")) {
    return(boundary)
  }
  if (methods::is(object, "ESS")) {
    shape_name <- as.character(acousticTS::extract(object, "shape_parameters")[["shape"]])[1]
    if ((identical(shape_name, "Sphere") &&
      boundary %in% c(
        "shelled_pressure_release",
        "shelled_liquid",
        "shelled_gas",
        "elastic_shelled"
      ))) {
      return(boundary)
    }

    if (identical(boundary, "elastic_shelled")) {
      return(boundary)
    }
  }

  if (identical(boundary, "elastic")) {
    stop(
      "TMM boundary 'elastic' is currently available only for solid elastic ",
      "ELA/CAL targets.",
      call. = FALSE
    )
  }

  if (!(boundary %in% c(
    "fixed_rigid",
    "pressure_release",
    "liquid_filled",
    "gas_filled"
  ))) {
    stop(
      "Only the following values for 'boundary' are available in TMM: ",
      "'fixed_rigid', 'pressure_release', 'liquid_filled', 'gas_filled', ",
      "'elastic', 'shelled_pressure_release', 'shelled_liquid', 'shelled_gas', ",
      "'elastic_shelled'.",
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

# Validate the optional retained-cylinder endcap smoothing fraction.
#' @noRd
.tmm_resolve_cylinder_endcap_fraction <- function(cylinder_endcap_fraction) {
  if (is.null(cylinder_endcap_fraction)) {
    return(NULL)
  }

  if (length(cylinder_endcap_fraction) != 1 ||
    is.na(cylinder_endcap_fraction) ||
    !is.finite(cylinder_endcap_fraction)) {
    stop(
      "'cylinder_endcap_fraction' must be NULL or one finite numeric value.",
      call. = FALSE
    )
  }

  cylinder_endcap_fraction <- as.numeric(cylinder_endcap_fraction)[1]
  if (cylinder_endcap_fraction < 0) {
    stop(
      "'cylinder_endcap_fraction' must be non-negative when supplied.",
      call. = FALSE
    )
  }

  if (cylinder_endcap_fraction == 0) {
    return(0)
  }

  min(cylinder_endcap_fraction, 0.45)
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
.tmm_prepare_body <- function(object,
                              sound_speed_sw,
                              density_sw,
                              boundary,
                              density_body = NULL,
                              sound_speed_longitudinal_body = NULL,
                              sound_speed_transversal_body = NULL) {
  shape_parameters <- acousticTS::extract(object, "shape_parameters")

  if (.tmm_is_elastic_solid_branch(object, boundary)) {
    body <- acousticTS::extract(object, "body")
    rho_body <- density_body %||% body$density
    cL_body <- sound_speed_longitudinal_body %||% body$sound_speed_longitudinal
    cT_body <- sound_speed_transversal_body %||% body$sound_speed_transversal
    if (is.null(rho_body) || is.null(cL_body) || is.null(cT_body)) {
      stop(
        "Elastic-solid TMM requires body density plus longitudinal and ",
        "transversal wave speeds, either stored on the scatterer or passed ",
        "directly to target_strength().",
        call. = FALSE
      )
    }

    shear_modulus <- rho_body * cT_body^2
    lame_lambda <- rho_body * cL_body^2 -
      2 * shear_modulus

    return(list(
      theta = body$theta,
      density = rho_body,
      radius = body$radius %||% shape_parameters$radius %||%
        shape_parameters$radius_body %||% NA_real_,
      diameter = body$diameter %||% shape_parameters$diameter %||%
        (2 * (shape_parameters$radius %||% shape_parameters$radius_body %||% NA_real_)),
      sound_speed_longitudinal = cL_body,
      sound_speed_transversal = cT_body,
      lambda = lame_lambda,
      mu = shear_modulus,
      medium_density = density_sw,
      medium_sound_speed = sound_speed_sw,
      semimajor_length = as.numeric(
        shape_parameters$semimajor_length %||%
          shape_parameters$length / 2
      )[1],
      semiminor_length = as.numeric(
        shape_parameters$semiminor_length %||%
          shape_parameters$radius
      )[1]
    ))
  }

  if (.tmm_is_shell_sphere_branch(
    object = object,
    shape_parameters = shape_parameters,
    boundary = boundary
  )) {
    shell <- .complete_material_props(
      acousticTS::extract(object, "shell"),
      medium_sound_speed = sound_speed_sw,
      medium_density = density_sw
    )
    fluid <- acousticTS::extract(object, "fluid")
    fluid <- .complete_material_props(
      fluid,
      medium_sound_speed = shell$sound_speed %||% (shell$h * sound_speed_sw),
      medium_density = shell$density %||% (shell$g * density_sw)
    )
    sph_body <- .sphms_body_parameters(
      object = object,
      exterior = shell,
      sound_speed_sw = sound_speed_sw,
      density_sw = density_sw
    )

    return(c(
      list(
        theta = shell$theta,
        density = shell$density,
        sound_speed = shell$sound_speed,
        g = shell$g,
        h = shell$h,
        shell_density = shell$density,
        shell_sound_speed = shell$sound_speed,
        shell_g = shell$g,
        shell_h = shell$h,
        fluid_density = fluid$density %||% NA_real_,
        fluid_sound_speed = fluid$sound_speed %||% NA_real_,
        fluid_g = fluid$g %||% NA_real_,
        fluid_h = fluid$h %||% NA_real_
      ),
      sph_body
    ))
  }

  if (.tmm_is_elastic_shell_sphere_branch(
    object = object,
    shape_parameters = shape_parameters,
    boundary = boundary
  )) {
    shell <- .extract_material_props(
      acousticTS::extract(object, "shell"),
      sound_speed_sw,
      density_sw
    )
    shell[names(.complete_elastic_moduli(
      E = shell$E,
      G = shell$G,
      K = shell$K,
      nu = shell$nu
    ))] <- .complete_elastic_moduli(
      E = shell$E,
      G = shell$G,
      K = shell$K,
      nu = shell$nu
    )
    fluid <- .extract_material_props(
      acousticTS::extract(object, "fluid"),
      sound_speed_sw,
      density_sw
    )
    shell_raw <- acousticTS::extract(object, "shell")
    fluid_raw <- acousticTS::extract(object, "fluid")

    return(list(
      theta = shell_raw$theta,
      density = shell$density,
      sound_speed = shell$sound_speed,
      radius_shell = shell_raw$radius,
      shell_thickness = shell_raw$shell_thickness,
      fluid_density = fluid$density,
      fluid_sound_speed = fluid$sound_speed,
      radius_fluid = fluid_raw$radius,
      shell_density = shell$density,
      shell_E = shell$E,
      shell_G = shell$G,
      shell_K = shell$K,
      shell_nu = shell$nu,
      shell_lambda = shell$lambda,
      medium_density = density_sw,
      medium_sound_speed = sound_speed_sw
    ))
  }

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
  shape_source <- shape_parameters$shell %||% shape_parameters
  a <- shape_source$semimajor_length %||% (shape_source$length / 2)
  b <- shape_source$semiminor_length %||% shape_source$radius
  xi <- 1 / sqrt(
    1 - (b / a)^2
  )
  q <- a / xi

  list(
    xi = xi,
    q = q,
    semimajor_length = a,
    semiminor_length = b
  )
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
  # Slender rigid prolates need a slightly higher retained ceiling than the
  # classic hard rule to stabilize broadside monostatic nulls in the stored
  # spheroidal TMM path.
  if (identical(boundary, "fixed_rigid")) {
    n_floor <- ceiling(abs(acoustics$k_sw) * geometry$semimajor_length + 4)
    m_floor <- ceiling(n_floor / 2)
    acoustics$m_max <- pmax(acoustics$m_max, m_floor)
    acoustics$n_max <- pmax(acoustics$n_max, n_floor)
  }

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
                                     use_cylindrical_branch,
                                     cylinder_backend = NULL) {
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
      boundary = boundary,
      cylinder_backend = cylinder_backend
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
                                   n_max,
                                   cylinder_backend = NULL,
                                   cylinder_endcap_fraction = NULL,
                                   store_t_matrix = FALSE) {
  if (identical(boundary, "elastic")) {
    acoustics <- .init_acoustics_df(
      frequency,
      k_sw = sound_speed_sw,
      k_l = body$sound_speed_longitudinal,
      k_t = body$sound_speed_transversal
    )
    if (identical(.tmm_shape_name(shape_parameters), "Sphere")) {
      acoustics$n_max <- NA_integer_
      return(list(acoustics = acoustics, geometry = NULL))
    }

    k_eff <- pmax(acoustics$k_sw, acoustics$k_l, acoustics$k_t)
    acoustics$n_max <- .tmm_prepare_n_max(
      n_max = n_max,
      frequency = frequency,
      k_sw = k_eff,
      shape_parameters = shape_parameters,
      boundary = boundary,
      cylinder_backend = cylinder_backend
    )

    return(list(acoustics = acoustics, geometry = NULL))
  }

  if (boundary %in% c(
    "shelled_pressure_release",
    "shelled_liquid",
    "shelled_gas"
  )) {
    acoustics <- .init_acoustics_df(
      frequency,
      k_sw = sound_speed_sw,
      k_shell = body$shell_sound_speed,
      k_fluid = if (is.finite(body$fluid_sound_speed)) body$fluid_sound_speed else NA_real_
    )
    acoustics$m_limit <- .sphms_m_limit(
      m_limit = n_max,
      acoustics = acoustics,
      body_params = body
    )
    acoustics$n_max <- acoustics$m_limit

    return(list(acoustics = acoustics, geometry = NULL))
  }

  if (identical(boundary, "elastic_shelled")) {
    if (identical(.tmm_shape_name(shape_parameters), "Sphere")) {
      acoustics <- .init_acoustics_df(
        frequency,
        k_sw = sound_speed_sw,
        k_fluid = body$fluid_sound_speed
      )
      acoustics$m_limit <- if (is.null(n_max)) {
        round(acoustics$k_sw * body$radius_shell) + 10L
      } else {
        as.integer(n_max)
      }
      acoustics$n_max <- acoustics$m_limit

      return(list(acoustics = acoustics, geometry = NULL))
    }

    acoustics <- .init_acoustics_df(
      frequency,
      k_sw = sound_speed_sw,
      k_body = body$fluid_sound_speed
    )
    acoustics$n_max <- .tmm_assign_branch_n_max(
      acoustics = acoustics,
      n_max = n_max,
      frequency = frequency,
      shape_parameters = shape_parameters,
      boundary = boundary,
      use_spheroidal_branch = FALSE,
      use_cylindrical_branch = FALSE,
      cylinder_backend = cylinder_backend
    )$n_max

    return(list(acoustics = acoustics, geometry = NULL))
  }

  # Start from the shared wavenumber table ====================================
  acoustics <- .tmm_base_acoustics(frequency, sound_speed_sw, body, boundary)
  geometry <- NULL
  if (use_spheroidal_branch) {
    geometry <- .tmm_spheroidal_geometry(shape_parameters)
    geometry$radius <- {
      shape_source <- shape_parameters$shell %||% shape_parameters
      shape_source$semiminor_length %||% shape_source$radius
    }
    acoustics <- .tmm_add_spheroidal_acoustics(acoustics, boundary, geometry)
  }
  acoustics <- .tmm_assign_branch_n_max(
    acoustics = acoustics,
    n_max = n_max,
    frequency = frequency,
    shape_parameters = shape_parameters,
    boundary = boundary,
    use_spheroidal_branch = use_spheroidal_branch,
    use_cylindrical_branch = use_cylindrical_branch,
    cylinder_backend = cylinder_backend
  )
  if (identical(.tmm_shape_name(shape_parameters), "Cylinder") &&
    identical(
      .tmm_cylinder_backend(
        shape_parameters = shape_parameters,
        cylinder_backend = cylinder_backend,
        cylinder_endcap_fraction = cylinder_endcap_fraction,
        boundary = boundary,
        store_t_matrix = store_t_matrix
      ),
      "retained"
    )) {
    acoustics$n_max <- .tmm_prepare_cylinder_n_max(
      n_max = n_max,
      frequency = frequency,
      k_sw = acoustics$k_sw,
      shape_parameters = shape_parameters
    )
  }

  list(acoustics = acoustics, geometry = geometry)
}

# Build the stored body metadata used by the TMM runtime helpers.
#' @noRd
.tmm_body_parameters <- function(body, geometry) {
  default_phi_body <- if (!is.null(geometry)) {
    pi / 2
  } else {
    pi
  }
  # Start from the common body-fixed incident and material properties ==========
  body_params <- list(
    theta_body = body$theta,
    theta_scatter = pi - body$theta,
    phi_body = body$phi_body %||% default_phi_body,
    phi_scatter = (body$phi_body %||% default_phi_body) + pi,
    density = body$density,
    sound_speed = body$sound_speed,
    g_body = body$g,
    h_body = body$h
  )
  if (!is.null(geometry)) {
    body_params$xi <- geometry$xi
    body_params$q <- geometry$q
  }
  for (nm in c(
    "radius_shell", "radius_fluid",
    "radius", "diameter",
    "semimajor_length", "semiminor_length", "shape_n_segments",
    "shell_density", "shell_sound_speed", "shell_g", "shell_h",
    "fluid_density", "fluid_sound_speed", "fluid_g", "fluid_h",
    "g21", "g31", "g32", "h21", "h31", "h32",
    "shell_thickness",
    "shell_E", "shell_G", "shell_K", "shell_nu", "shell_lambda",
    "sound_speed_longitudinal", "sound_speed_transversal",
    "lambda", "mu",
    "medium_density", "medium_sound_speed"
  )) {
    if (!is.null(body[[nm]])) {
      body_params[[nm]] <- body[[nm]]
    }
  }

  body_params
}

# Resolve the coordinate-system label stored with the initialized TMM object.
#' @noRd
.tmm_coordinate_system <- function(use_spheroidal_branch,
                                   use_cylindrical_branch,
                                   use_shell_sphere_branch = FALSE,
                                   use_elastic_solid_sphere_branch = FALSE,
                                   shape_parameters = NULL,
                                   cylinder_backend = NULL) {
  # Label the active TMM backend coordinate system ============================
  if (use_elastic_solid_sphere_branch) {
    return("sphere_elastic_exact")
  }
  if (use_shell_sphere_branch) {
    return("sphere_modal")
  }
  if (use_spheroidal_branch) {
    return("spheroidal")
  }
  if (use_cylindrical_branch) {
    return("cylindrical")
  }
  if (!is.null(shape_parameters) &&
    identical(.tmm_shape_name(shape_parameters), "Cylinder") &&
    identical(cylinder_backend, "retained")) {
    return("cylinder_native")
  }

  "spherical"
}

# Resolve the stored precision label for one initialized TMM object.
#' @noRd
.tmm_precision_label <- function(use_spheroidal_branch, boundary) {
  # Prolate spheroidal runs require the quad-precision backend for the
  # published benchmark regime; double precision drifts materially once the
  # reduced frequency is large enough.
  if (use_spheroidal_branch) {
    return("quad")
  }

  "double"
}

# Resolve the stored quadrature-node count for the spheroidal TMM branch.
#' @noRd
.tmm_n_integration_label <- function(use_spheroidal_branch, boundary) {
  # Only penetrable spheroidal runs currently keep an explicit quadrature count
  if (use_spheroidal_branch &&
    boundary %in% c("liquid_filled", "gas_filled", "elastic_shelled")) {
    return(96L)
  }

  NA_integer_
}

# Return the outer size scale used in the default spherical truncation rule.
#' @noRd
.tmm_bounding_radius <- function(shape_parameters) {
  shape_source <- shape_parameters$shell %||% shape_parameters
  switch(.tmm_shape_name(shape_parameters),
    Sphere = as.numeric(shape_source$radius)[1],
    ProlateSpheroid = as.numeric(shape_source$semimajor_length)[1],
    OblateSpheroid = as.numeric(shape_source$semimajor_length)[1],
    Cylinder = sqrt(
      (as.numeric(shape_source$length)[1] / 2)^2 +
        max(as.numeric(shape_source$radius), na.rm = TRUE)^2
    )
  )
}

# Apply conservative hard overrides for slender prolates when a spherical-wave
# truncation is still needed.
#' @noRd
.tmm_prolate_nmax_override <- function(shape_parameters, boundary) {
  if (.tmm_shape_name(shape_parameters) != "ProlateSpheroid") {
    return(NA_integer_)
  }

  shape_source <- shape_parameters$shell %||% shape_parameters
  a <- as.numeric(shape_source$semimajor_length)[1]
  b <- as.numeric(shape_source$semiminor_length)[1]
  aspect_ratio <- a / b

  if (!is.finite(aspect_ratio) || aspect_ratio < 3) {
    return(NA_integer_)
  }

  switch(boundary,
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
.tmm_cylinder_nmax_floor <- function(shape_parameters,
                                     boundary,
                                     cylinder_backend = NULL) {
  if (.tmm_shape_name(shape_parameters) != "Cylinder" ||
    identical(cylinder_backend, "retained")) {
    return(NA_integer_)
  }

  switch(boundary,
    fixed_rigid = 56L,
    pressure_release = 52L,
    liquid_filled = 64L,
    gas_filled = 56L,
    NA_integer_
  )
}

# Cylinder TMM in the spherical-coordinate branch behaves best when the
# retained degree is sized from the cylindrical radius rather than the bounding
# sphere. The offsets below were tuned against the finite-cylinder benchmark
# columns so the sharp-edged cylinder stays stable without the severe ringing
# caused by the old bounding-sphere floor. Pressure-release cylinders need a
# materially higher retained ceiling than the other boundaries, while the
# penetrable gas branch behaves best with a slightly leaner truncation.
#' @noRd
.tmm_prepare_cylinder_spherical_n_max <- function(n_max,
                                                  frequency,
                                                  k_sw,
                                                  shape_parameters,
                                                  boundary) {
  if (is.null(n_max)) {
    radius_body <- max(as.numeric(shape_parameters$radius), na.rm = TRUE)
    ka <- abs(k_sw) * radius_body
    offset <- switch(boundary,
      pressure_release = 36L,
      fixed_rigid = 46L,
      gas_filled = 15L,
      liquid_filled = 4L,
      10L
    )
    return(as.integer(pmax(4L, ceiling(ka) + offset)))
  }

  .tmm_prepare_n_max(
    n_max = n_max,
    frequency = frequency,
    k_sw = k_sw,
    shape_parameters = shape_parameters,
    boundary = boundary
  )
}

# Oblate spheroids need a more conservative retained spherical-wave floor than
# the generic bounding-sphere heuristic, especially for rigid and liquid cases
# where the default low-frequency truncation materially under-resolves broadside
# bistatic null structure.
#' @noRd
.tmm_oblate_nmax_floor <- function(shape_parameters, boundary) {
  if (.tmm_shape_name(shape_parameters) != "OblateSpheroid") {
    return(NA_integer_)
  }

  switch(boundary,
    fixed_rigid = 24L,
    liquid_filled = 24L,
    NA_integer_
  )
}

# Default truncation based on the exterior size parameter, with geometry-
# specific overrides when they are known to behave better.
#' @noRd
.tmm_default_n_max <- function(k_sw,
                               shape_parameters,
                               boundary,
                               cylinder_backend = NULL) {
  n_override <- .tmm_prolate_nmax_override(shape_parameters, boundary)
  if (!is.na(n_override)) {
    return(rep.int(as.integer(n_override), length(k_sw)))
  }

  x <- abs(k_sw) * .tmm_bounding_radius(shape_parameters)
  n_default <- pmax(4L, as.integer(ceiling(x + 4 * x^(1 / 3) + 2)))
  n_floor <- .tmm_cylinder_nmax_floor(
    shape_parameters,
    boundary,
    cylinder_backend = cylinder_backend
  )
  oblate_floor <- .tmm_oblate_nmax_floor(shape_parameters, boundary)
  if (!is.na(n_floor)) {
    n_default <- pmax(n_default, as.integer(n_floor))
  }
  if (!is.na(oblate_floor)) {
    n_default <- pmax(n_default, as.integer(oblate_floor))
  }
  n_default
}

# Normalize `n_max` into a per-frequency integer vector.
#' @noRd
.tmm_prepare_n_max <- function(
  n_max, frequency, k_sw, shape_parameters, boundary, cylinder_backend = NULL
) {
  if (is.null(n_max)) {
    return(.tmm_default_n_max(
      k_sw,
      shape_parameters,
      boundary,
      cylinder_backend = cylinder_backend
    ))
  }

  if (!is.numeric(n_max) || any(!is.finite(n_max)) || any(n_max < 1)) {
    stop("'n_max' must be NULL or a positive finite integer vector.",
      call. = FALSE
    )
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

# Build the candidate retained-degree set scanned by the elastic prolate TMM
# selector. The current elastic solid branch becomes unstable when the raw
# bounding-sphere rule drives `n_max` too aggressively, so the package scans a
# local neighborhood around that seed and selects a smooth modal path from the
# actual retained solves rather than trusting one closed-form cutoff.
#' @noRd
.tmm_elastic_prolate_candidate_n_max <- function(n_seed) {
  # Center the scan around the default retained-degree seed ====================
  n_seed <- as.integer(round(n_seed))[1]
  lower <- max(5L, n_seed - 12L)
  upper <- max(lower, n_seed + 6L)

  # Scan both parities because the elastic prolate branch can stabilize on
  # either side of the raw retained-degree seed.
  seq.int(lower, upper, by = 1L)
}

# Summarize the retained elastic-solid blocks with the diagnostics needed by
# the automatic `n_max` selector.
#' @noRd
.tmm_elastic_prolate_block_diagnostics <- function(blocks) {
  # Flatten the solved block bookkeeping across all retained azimuthal orders ===
  if (!length(blocks)) {
    return(list(
      max_coeff = Inf,
      low_max = Inf,
      tail_ratio = Inf
    ))
  }

  coeff_values <- unlist(
    lapply(blocks, function(block) Mod(block$coefficients)),
    use.names = FALSE
  )
  n_values <- unlist(
    lapply(blocks, function(block) block$n_seq),
    use.names = FALSE
  )
  m_values <- unlist(
    lapply(blocks, function(block) {
      rep.int(block$m, length(block$n_seq))
    }),
    use.names = FALSE
  )
  max_n <- max(n_values, na.rm = TRUE)
  low_mask <- m_values <= 2L & n_values <= (m_values + 2L)
  tail_mask <- n_values >= (max_n - 2L)

  max_coeff <- max(coeff_values, na.rm = TRUE)
  low_max <- if (any(low_mask)) {
    max(coeff_values[low_mask], na.rm = TRUE)
  } else {
    max_coeff
  }
  tail_max <- if (any(tail_mask)) {
    max(coeff_values[tail_mask], na.rm = TRUE)
  } else {
    max_coeff
  }

  list(
    max_coeff = max_coeff,
    low_max = low_max,
    tail_ratio = tail_max / max(max_coeff, sqrt(.Machine$double.eps))
  )
}

# Scan a local retained-degree neighborhood for the solid-elastic prolate TMM
# branch and record the modal diagnostics used by the automatic selector.
#' @noRd
.tmm_elastic_prolate_scan_n_path <- function(acoustics,
                                             body,
                                             medium,
                                             shape_parameters,
                                             boundary,
                                             svd_rel_tol = 1e-10) {
  # Keep the selector scoped to the current public elastic-prolate branch ======
  if (!identical(boundary, "elastic") ||
    !identical(.tmm_shape_name(shape_parameters), "ProlateSpheroid")) {
    return(NULL)
  }

  rows <- vector("list", 0L)
  theta_incident <- body$theta_body %||% body$theta
  for (i in seq_len(nrow(acoustics))) {
    candidates <- .tmm_elastic_prolate_candidate_n_max(acoustics$n_max[i])
    for (n_i in candidates) {
      fit_i <- .tmm_single_frequency_spherical(
        k_sw = acoustics$k_sw[i],
        k_body = acoustics$k_body[i],
        k_l = acoustics$k_l[i],
        k_t = acoustics$k_t[i],
        theta_body = theta_incident,
        boundary = boundary,
        shape_parameters = shape_parameters,
        rho_sw = medium$density[1],
        rho_body = body$density,
        lambda = body$lambda,
        mu = body$mu,
        n_max = n_i,
        store_t_matrix = FALSE,
        frequency_hz = acoustics$frequency[i],
        svd_rel_tol = svd_rel_tol
      )
      diag_i <- .tmm_elastic_prolate_block_diagnostics(fit_i$blocks)
      rows[[length(rows) + 1L]] <- data.frame(
        freq_index = i,
        frequency = acoustics$frequency[i],
        n_max = n_i,
        TS = db(.sigma_bs(fit_i$f_bs)),
        max_coeff = diag_i$max_coeff,
        low_max = diag_i$low_max,
        tail_ratio = diag_i$tail_ratio
      )
    }
  }

  do.call(rbind, rows)
}

# Select a smooth retained-degree path for the elastic prolate TMM branch from
# the scanned modal diagnostics.
#' @noRd
.tmm_elastic_prolate_select_n_path <- function(scan_df) {
  # Return the empty integer template when no scan was supplied ================
  if (is.null(scan_df) || !nrow(scan_df)) {
    return(integer(0))
  }

  scan_df <- scan_df[order(scan_df$freq_index, scan_df$n_max), , drop = FALSE]
  local_min_coeff <- tapply(scan_df$max_coeff, scan_df$freq_index, min)
  local_cost <- with(
    scan_df,
    log1p(low_max) +
      0.25 * log1p(max_coeff) +
      10 * sqrt(pmax(tail_ratio, 0))
  )
  unstable_mask <- scan_df$max_coeff >
    4 * local_min_coeff[as.character(scan_df$freq_index)]
  local_cost <- local_cost + ifelse(unstable_mask, 5, 0)
  scan_df$local_cost <- local_cost

  split_scan <- split(scan_df, scan_df$freq_index)
  dp_cost <- vector("list", length(split_scan))
  dp_prev <- vector("list", length(split_scan))
  dp_cost[[1L]] <- split_scan[[1L]]$local_cost
  dp_prev[[1L]] <- rep(NA_integer_, nrow(split_scan[[1L]]))

  if (length(split_scan) > 1L) {
    for (i in 2:length(split_scan)) {
      prev_df <- split_scan[[i - 1L]]
      curr_df <- split_scan[[i]]
      cost_i <- rep(Inf, nrow(curr_df))
      prev_i <- rep(NA_integer_, nrow(curr_df))

      for (j in seq_len(nrow(curr_df))) {
        transition_cost <- 0.35 * abs(curr_df$TS[j] - prev_df$TS) +
          0.05 * abs(curr_df$n_max[j] - prev_df$n_max)
        candidate_cost <- dp_cost[[i - 1L]] +
          transition_cost +
          curr_df$local_cost[j]
        best_idx <- which.min(candidate_cost)
        cost_i[j] <- candidate_cost[best_idx]
        prev_i[j] <- best_idx
      }

      dp_cost[[i]] <- cost_i
      dp_prev[[i]] <- prev_i
    }
  }

  path_idx <- integer(length(split_scan))
  path_idx[length(split_scan)] <- which.min(dp_cost[[length(split_scan)]])
  if (length(split_scan) > 1L) {
    for (i in seq(length(split_scan), 2L)) {
      path_idx[i - 1L] <- dp_prev[[i]][path_idx[i]]
    }
  }

  as.integer(vapply(
    seq_along(split_scan),
    function(i) as.integer(split_scan[[i]]$n_max[path_idx[i]]),
    integer(1)
  ))
}

# Cylinder monostatic TMM uses the same modal truncation rule as FCMS so that
# the cylindrical-coordinate backend remains benchmark-compatible by default.
#' @noRd
.tmm_prepare_cylinder_n_max <- function(n_max, frequency, k_sw, shape_parameters) {
  if (is.null(n_max)) {
    return(as.integer(ceiling(k_sw * max(as.numeric(shape_parameters$radius),
      na.rm = TRUE
    )) + 10L))
  }

  if (!is.numeric(n_max) || any(!is.finite(n_max)) || any(n_max < 1)) {
    stop("'n_max' must be NULL or a positive finite integer vector.",
      call. = FALSE
    )
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
  if (identical(.tmm_shape_name(shape_parameters), "Cylinder")) {
    return(max(128L, 10L * n_terms))
  }

  if (boundary %in% c("elastic_shelled", "elastic")) {
    return(max(128L, 8L * n_terms))
  }

  max(64L, 4L * n_terms)
}

# Flatten the supported shape parameters into the compact numeric vector passed
# into the compiled spherical backend.
#' @noRd
.tmm_shape_values <- function(shape_parameters) {
  shape_source <- shape_parameters$shell %||% shape_parameters
  switch(.tmm_shape_name(shape_parameters),
    Sphere = c(as.numeric(shape_source$radius)[1]),
    ProlateSpheroid = c(
      as.numeric(shape_source$semimajor_length)[1],
      as.numeric(shape_source$semiminor_length)[1]
    ),
    OblateSpheroid = c(
      as.numeric(shape_source$semiminor_length)[1],
      as.numeric(shape_source$semimajor_length)[1]
    ),
    Cylinder = c(
      as.numeric(shape_source$length)[1] / 2,
      max(as.numeric(shape_source$radius), na.rm = TRUE)
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

  Bm <- switch(boundary,
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
.tmm_store_cylindrical_branch <- function(acoustics,
                                         family = "cylindrical_mono") {
  lapply(
    seq_len(nrow(acoustics)),
    function(i) {
      list(
        family = family,
        n_max = acoustics$n_max[i]
      )
    }
  )
}
