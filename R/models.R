################################################################################
# MODEL REGISTRY
################################################################################

.model_registry_state <- local({
  env <- new.env(parent = emptyenv())
  env$user <- list()
  env$loaded <- FALSE
  env
})

#' Normalize model identifiers used by the model registry
#' @param x Character vector of model identifiers.
#' @param arg_name Argument name used in error messages.
#' @keywords internal
#' @noRd
.normalize_model_identifier <- function(x, arg_name = "model") {
  if (!is.character(x)) {
    stop("`", arg_name, "` must be character.", call. = FALSE)
  }

  x <- trimws(tolower(x))
  if (any(is.na(x) | !nzchar(x))) {
    stop(
      "`", arg_name, "` must contain non-empty model names only.",
      call. = FALSE
    )
  }

  x
}

#' Resolve the default result-slot name for one model identifier
#' @param name Canonical model identifier.
#' @keywords internal
#' @noRd
.default_model_slot <- function(name) {
  slot <- gsub("(_.*)", "\\L\\1", toupper(name), perl = TRUE)
  if (slot %in% c("CALIBRATION", "HIGH_pass_stanton")) {
    slot <- tolower(slot)
  }

  slot
}

#' Build one normalized model-registry entry
#' @keywords internal
#' @noRd
.new_model_registry_entry <- function(name,
                                      initialize,
                                      solver,
                                      slot = .default_model_slot(name),
                                      aliases = name,
                                      source = c("builtin", "user"),
                                      persistent = FALSE,
                                      initialize_ref = NULL,
                                      solver_ref = NULL) {
  source <- match.arg(source)
  canonical <- .normalize_model_identifier(name, "name")[[1]]

  if (!is.character(slot) || length(slot) != 1L || is.na(slot) || !nzchar(slot)) {
    stop("`slot` must be a single non-empty string.", call. = FALSE)
  }

  aliases <- unique(.normalize_model_identifier(c(aliases, canonical), "aliases"))

  list(
    canonical = canonical,
    slot = slot,
    initialize = initialize,
    solver = solver,
    aliases = aliases,
    source = source,
    persistent = isTRUE(persistent),
    initialize_ref = initialize_ref,
    solver_ref = solver_ref
  )
}

#' Built-in acousticTS model registry
#' @keywords internal
#' @noRd
.builtin_model_registry <- function() {
  list(
    dwba = .new_model_registry_entry(
      name = "dwba",
      initialize = "dwba_initialize",
      solver = "DWBA",
      slot = "DWBA"
    ),
    bbfm = .new_model_registry_entry(
      name = "bbfm",
      initialize = "bbfm_initialize",
      solver = "BBFM",
      slot = "BBFM"
    ),
    pcdwba = .new_model_registry_entry(
      name = "pcdwba",
      initialize = "pcdwba_initialize",
      solver = "PCDWBA",
      slot = "PCDWBA"
    ),
    dwba_curved = .new_model_registry_entry(
      name = "dwba_curved",
      initialize = "dwba_curved_initialize",
      solver = "DWBA_curved",
      slot = "DWBA_curved"
    ),
    sdwba = .new_model_registry_entry(
      name = "sdwba",
      initialize = "sdwba_initialize",
      solver = "SDWBA",
      slot = "SDWBA"
    ),
    sdwba_curved = .new_model_registry_entry(
      name = "sdwba_curved",
      initialize = "sdwba_curved_initialize",
      solver = "SDWBA_curved",
      slot = "SDWBA_curved"
    ),
    calibration = .new_model_registry_entry(
      name = "calibration",
      initialize = "calibration_initialize",
      solver = "calibration",
      slot = "calibration",
      aliases = c("calibration", "soems")
    ),
    essms = .new_model_registry_entry(
      name = "essms",
      initialize = "essms_initialize",
      solver = "ESSMS",
      slot = "ESSMS"
    ),
    vesms = .new_model_registry_entry(
      name = "vesms",
      initialize = "vesms_initialize",
      solver = "VESMS",
      slot = "VESMS"
    ),
    sphms = .new_model_registry_entry(
      name = "sphms",
      initialize = "sphms_initialize",
      solver = "SPHMS",
      slot = "SPHMS"
    ),
    tmm = .new_model_registry_entry(
      name = "tmm",
      initialize = "tmm_initialize",
      solver = "TMM",
      slot = "TMM"
    ),
    krm = .new_model_registry_entry(
      name = "krm",
      initialize = "krm_initialize",
      solver = "KRM",
      slot = "KRM"
    ),
    hpa = .new_model_registry_entry(
      name = "hpa",
      initialize = "hpa_initialize",
      solver = "HPA",
      slot = "HPA",
      aliases = c("hpa", "high_pass")
    ),
    psms = .new_model_registry_entry(
      name = "psms",
      initialize = "psms_initialize",
      solver = "PSMS",
      slot = "PSMS"
    ),
    fcms = .new_model_registry_entry(
      name = "fcms",
      initialize = "fcms_initialize",
      solver = "FCMS",
      slot = "FCMS"
    ),
    bcms = .new_model_registry_entry(
      name = "bcms",
      initialize = "bcms_initialize",
      solver = "BCMS",
      slot = "BCMS"
    ),
    ecms = .new_model_registry_entry(
      name = "ecms",
      initialize = "ecms_initialize",
      solver = "ECMS",
      slot = "ECMS"
    ),
    trcm = .new_model_registry_entry(
      name = "trcm",
      initialize = "trcm_initialize",
      solver = "TRCM",
      slot = "TRCM"
    )
  )
}

#' User-level registry path for persisted model registrations
#' @keywords internal
#' @noRd
.model_registry_user_path <- function() {
  file.path(
    tools::R_user_dir("acousticTS", which = "config"),
    "user-model-registry.rds"
  )
}

#' Resolve a function reference used in the model registry
#' @keywords internal
#' @noRd
.resolve_model_function_reference <- function(ref,
                                              envir = parent.frame(),
                                              arg_name = "reference") {
  if (!is.character(ref) || length(ref) != 1L || is.na(ref) || !nzchar(ref)) {
    stop(
      "`", arg_name, "` must be a single non-empty function reference.",
      call. = FALSE
    )
  }

  tryCatch(
    {
      if (grepl(":::", ref, fixed = TRUE)) {
        parts <- strsplit(ref, ":::", fixed = TRUE)[[1]]
        return(getFromNamespace(parts[[2]], parts[[1]]))
      }

      if (grepl("::", ref, fixed = TRUE)) {
        parts <- strsplit(ref, "::", fixed = TRUE)[[1]]
        return(getExportedValue(parts[[1]], parts[[2]]))
      }

      get(ref, mode = "function", envir = envir, inherits = TRUE)
    },
    error = function(e) {
      stop(
        "Could not resolve ", arg_name, " function reference '", ref, "'.",
        call. = FALSE
      )
    }
  )
}

#' Check whether a function reference is package-qualified
#' @keywords internal
#' @noRd
.is_qualified_model_reference <- function(ref) {
  is.character(ref) &&
    length(ref) == 1L &&
    !is.na(ref) &&
    nzchar(ref) &&
    (grepl(":::", ref, fixed = TRUE) || grepl("::", ref, fixed = TRUE))
}

#' Coerce one registry callable input into a function plus optional reference
#' @keywords internal
#' @noRd
.coerce_model_callable <- function(x,
                                   arg_name,
                                   envir = parent.frame()) {
  if (is.function(x)) {
    return(list(fn = x, ref = NULL))
  }

  if (is.character(x) && length(x) == 1L && !is.na(x) && nzchar(x)) {
    return(list(
      fn = .resolve_model_function_reference(x, envir = envir, arg_name = arg_name),
      ref = x
    ))
  }

  stop(
    "`", arg_name, "` must be either a function or a single function reference.",
    call. = FALSE
  )
}

#' Convert one user registry entry into its persisted representation
#' @keywords internal
#' @noRd
.serialize_model_registry_entry <- function(entry) {
  list(
    canonical = entry$canonical,
    slot = entry$slot,
    aliases = entry$aliases,
    initialize_ref = entry$initialize_ref,
    solver_ref = entry$solver_ref
  )
}

#' Rebuild one persisted user registry entry
#' @keywords internal
#' @noRd
.deserialize_model_registry_entry <- function(entry) {
  .new_model_registry_entry(
    name = entry$canonical,
    initialize = .resolve_model_function_reference(
      entry$initialize_ref,
      arg_name = paste0("initialize reference for model '", entry$canonical, "'")
    ),
    solver = .resolve_model_function_reference(
      entry$solver_ref,
      arg_name = paste0("solver reference for model '", entry$canonical, "'")
    ),
    slot = entry$slot,
    aliases = entry$aliases,
    source = "user",
    persistent = TRUE,
    initialize_ref = entry$initialize_ref,
    solver_ref = entry$solver_ref
  )
}

#' Persist the current user model registry
#' @keywords internal
#' @noRd
.write_persistent_model_registry <- function() {
  .ensure_model_registry_loaded()

  persistent_entries <- Filter(
    function(entry) isTRUE(entry$persistent),
    .model_registry_state$user
  )
  path <- .model_registry_user_path()

  if (length(persistent_entries) == 0L) {
    if (file.exists(path)) {
      unlink(path)
    }
    return(invisible(path))
  }

  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  saveRDS(
    lapply(persistent_entries, .serialize_model_registry_entry),
    file = path
  )

  invisible(path)
}

#' Load persisted user model registrations once per session
#' @keywords internal
#' @noRd
.ensure_model_registry_loaded <- function() {
  if (isTRUE(.model_registry_state$loaded)) {
    return(invisible(TRUE))
  }

  .model_registry_state$user <- list()
  path <- .model_registry_user_path()

  if (file.exists(path)) {
    stored <- tryCatch(
      readRDS(path),
      error = function(e) {
        warning(
          "Could not read persisted acousticTS user models from '",
          path, "'."
        )
        list()
      }
    )

    if (is.list(stored)) {
      for (entry in stored) {
        restored <- tryCatch(
          .deserialize_model_registry_entry(entry),
          error = function(e) {
            warning(conditionMessage(e))
            NULL
          }
        )

        if (!is.null(restored)) {
          .model_registry_state$user[[restored$canonical]] <- restored
        }
      }
    }
  }

  .model_registry_state$loaded <- TRUE
  invisible(TRUE)
}

#' Return the combined built-in and user model registry
#' @keywords internal
#' @noRd
.model_registry_entries <- function() {
  .ensure_model_registry_loaded()
  c(.builtin_model_registry(), .model_registry_state$user)
}

#' Resolve one model registry entry by canonical name or alias
#' @keywords internal
#' @noRd
.resolve_model_registry_entry <- function(name) {
  requested <- .normalize_model_identifier(name, "model")[[1]]
  entries <- .model_registry_entries()
  matches <- Filter(function(entry) requested %in% entry$aliases, entries)

  if (length(matches) == 1L) {
    return(matches[[1]])
  }

  available <- sort(unique(unlist(lapply(entries, function(entry) entry$aliases))))
  stop(
    "Unknown target strength model '", name, "'. Available models: ",
    paste(available, collapse = ", "),
    call. = FALSE
  )
}

#' Resolve a vector of model registry entries
#' @keywords internal
#' @noRd
.resolve_model_registry_entries <- function(model) {
  lapply(model, .resolve_model_registry_entry)
}

#' Resolve the initializer function for one registry entry
#' @keywords internal
#' @noRd
.model_registry_initializer <- function(entry) {
  if (is.function(entry$initialize)) {
    return(entry$initialize)
  }

  .resolve_model_function_reference(
    entry$initialize,
    envir = asNamespace("acousticTS"),
    arg_name = paste0("initialize function for model '", entry$canonical, "'")
  )
}

#' Resolve the solver function for one registry entry
#' @keywords internal
#' @noRd
.model_registry_solver <- function(entry) {
  if (is.function(entry$solver)) {
    return(entry$solver)
  }

  .resolve_model_function_reference(
    entry$solver,
    envir = asNamespace("acousticTS"),
    arg_name = paste0("solver function for model '", entry$canonical, "'")
  )
}

#' Model registry comprising target strength models currently supported by the
#' acousticTS R-package.
#' @keywords internal
#' @noRd
.get_models <- function() {
  entries <- .model_registry_entries()
  out <- list()

  for (entry in entries) {
    solver <- .model_registry_solver(entry)
    labels <- unique(c(entry$slot, entry$aliases))
    for (label in labels) {
      out[[label]] <- solver
    }
  }

  out
}

#' List available target-strength models
#'
#' @return
#' A data frame describing currently available built-in and user-registered
#' target-strength models.
#' @export
available_models <- function() {
  entries <- .model_registry_entries()
  out <- do.call(rbind, lapply(entries, function(entry) {
    data.frame(
      model = entry$canonical,
      slot = entry$slot,
      source = entry$source,
      persistent = isTRUE(entry$persistent),
      aliases = paste(setdiff(entry$aliases, entry$canonical), collapse = ", "),
      stringsAsFactors = FALSE
    )
  }))

  out[order(out$source != "builtin", out$model), , drop = FALSE]
}

#' Register a user-defined target-strength model
#'
#' @param name Canonical model name passed to [target_strength()].
#' @param initialize Initializer function or function reference.
#' @param solver Solver function or function reference.
#' @param slot Result-slot name stored on the scatterer. Defaults to the usual
#'   upper-case model label derived from `name`.
#' @param aliases Optional alternative model names.
#' @param persist Logical scalar. When `TRUE`, the model registration is written
#'   to the user's `R_user_dir()` config path and reloaded in later sessions.
#'   Persistent registrations require package-qualified function references such
#'   as `"mypkg::tsl_initialize"` and `"mypkg::TSL"`.
#' @param overwrite Logical scalar. When `TRUE`, an existing user registration
#'   with the same canonical name can be replaced. Built-in models cannot be
#'   overwritten.
#' @return Invisibly returns the normalized registry entry.
#' @export
register_model <- function(name,
                           initialize,
                           solver,
                           slot = NULL,
                           aliases = character(),
                           persist = FALSE,
                           overwrite = FALSE) {
  .ensure_model_registry_loaded()

  canonical <- .normalize_model_identifier(name, "name")[[1]]
  initialize_info <- .coerce_model_callable(
    initialize,
    arg_name = "initialize",
    envir = parent.frame()
  )
  solver_info <- .coerce_model_callable(
    solver,
    arg_name = "solver",
    envir = parent.frame()
  )
  slot <- if (is.null(slot)) .default_model_slot(canonical) else slot
  aliases <- unique(.normalize_model_identifier(c(aliases, canonical), "aliases"))

  builtin_entries <- .builtin_model_registry()
  if (canonical %in% names(builtin_entries)) {
    stop(
      "Built-in model '", canonical, "' cannot be overwritten.",
      call. = FALSE
    )
  }

  existing_user <- .model_registry_state$user[[canonical]]
  if (!is.null(existing_user) && !isTRUE(overwrite)) {
    stop(
      "A user model named '", canonical,
      "' is already registered. Use `overwrite = TRUE` to replace it.",
      call. = FALSE
    )
  }

  other_entries <- .model_registry_entries()
  if (!is.null(existing_user)) {
    other_entries <- other_entries[setdiff(names(other_entries), canonical)]
  }
  conflicts <- unique(unlist(lapply(other_entries, function(entry) {
    intersect(aliases, entry$aliases)
  })))
  if (length(conflicts) > 0L) {
    stop(
      "Cannot register model '", canonical,
      "' because these names are already in use: ",
      paste(conflicts, collapse = ", "),
      call. = FALSE
    )
  }

  if (isTRUE(persist) && (
    is.null(initialize_info$ref) ||
      is.null(solver_info$ref) ||
      !.is_qualified_model_reference(initialize_info$ref) ||
      !.is_qualified_model_reference(solver_info$ref)
  )) {
    stop(
      "Persistent model registrations require package-qualified function ",
      "references such as 'mypkg::tsl_initialize' and 'mypkg::TSL'.",
      call. = FALSE
    )
  }

  entry <- .new_model_registry_entry(
    name = canonical,
    initialize = initialize_info$fn,
    solver = solver_info$fn,
    slot = slot,
    aliases = aliases,
    source = "user",
    persistent = isTRUE(persist),
    initialize_ref = initialize_info$ref,
    solver_ref = solver_info$ref
  )

  .model_registry_state$user[[canonical]] <- entry
  .write_persistent_model_registry()

  invisible(entry)
}

#' Remove a user-defined target-strength model registration
#'
#' @param name Canonical model name or alias.
#' @param remove_persisted Logical scalar. When `TRUE`, any persisted registry
#'   entry is also removed from the user's config file.
#' @return Invisibly returns the removed canonical model name.
#' @export
unregister_model <- function(name, remove_persisted = TRUE) {
  .ensure_model_registry_loaded()

  entry <- .resolve_model_registry_entry(name)
  if (!identical(entry$source, "user")) {
    stop(
      "Only user-registered models can be removed with `unregister_model()`.",
      call. = FALSE
    )
  }

  .model_registry_state$user[[entry$canonical]] <- NULL

  if (isTRUE(remove_persisted) && isTRUE(entry$persistent)) {
    .write_persistent_model_registry()
  }

  invisible(entry$canonical)
}

#' Clear user-defined model registrations
#'
#' @param remove_persisted Logical scalar. When `TRUE`, persisted user model
#'   registrations are also deleted from the user's config file.
#' @return Invisibly returns `NULL`.
#' @export
reset_model_registry <- function(remove_persisted = FALSE) {
  .ensure_model_registry_loaded()
  .model_registry_state$user <- list()

  if (isTRUE(remove_persisted)) {
    path <- .model_registry_user_path()
    if (file.exists(path)) {
      unlink(path)
    }
  }

  invisible(NULL)
}
