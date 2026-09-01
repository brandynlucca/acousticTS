# Select implementation families, run each in an isolated sequential R process,
# validate manifest outputs, and write a provenance report.

# Locate the repository from the working directory, CI workspace, or script path.
find_repo_root <- function() {
  # A valid root contains both package metadata and the figure manifest.
  has_impl_manifest <- function(path) {
    file.exists(file.path(path, "DESCRIPTION")) &&
      file.exists(file.path(path, "tools", "implementation-figures", "manifest.csv"))
  }

  # Walk upward from one candidate until a valid root or filesystem root.
  search_upward <- function(start_path) {
    current <- normalizePath(start_path, winslash = "/", mustWork = FALSE)
    repeat {
      if (has_impl_manifest(current)) return(current)
      next_path <- dirname(current)
      if (identical(next_path, current)) return(NULL)
      current <- next_path
    }
  }

  search_roots <- c(Sys.getenv("GITHUB_WORKSPACE", unset = ""), getwd())
  script_arg <- commandArgs(trailingOnly = FALSE)
  file_arg <- script_arg[grep("^--file=", script_arg)]
  if (length(file_arg) == 1L) {
    search_roots <- c(search_roots, dirname(sub("^--file=", "", file_arg)))
  }

  for (root in unique(search_roots[nzchar(search_roots)])) {
    found_root <- search_upward(root)
    if (!is.null(found_root)) return(found_root)
  }

  stop(
    "Could not locate the repository root containing ",
    "'tools/implementation-figures/manifest.csv'.",
    call. = FALSE
  )
}

# Split semicolon-delimited manifest cells into normalized values.
split_manifest_values <- function(values) {
  values <- values[!is.na(values) & nzchar(values)]
  if (!length(values)) return(character())
  values <- trimws(unlist(strsplit(values, ";", fixed = TRUE)))
  unique(values[nzchar(values)])
}

# Fail when a completed family did not create every output it declares.
validate_outputs <- function(rows, family) {
  expected <- split_manifest_values(rows$outputs)
  missing <- expected[!file.exists(expected)]
  if (length(missing)) {
    stop(
      "Implementation figure family '", family,
      "' did not produce its declared outputs:\n- ",
      paste(missing, collapse = "\n- "),
      call. = FALSE
    )
  }
  expected
}

# Read the current commit for provenance without making Git a hard dependency.
git_commit <- function() {
  commit <- tryCatch(
    system2("git", c("rev-parse", "HEAD"), stdout = TRUE, stderr = FALSE),
    error = function(...) character()
  )
  if (length(commit) == 1L && nzchar(commit)) commit else NA_character_
}

# Write hashes, sizes, timing, and commit information for generated outputs.
write_provenance <- function(records, repo_root) {
  if (!length(records)) return(invisible(NULL))

  provenance_path <- Sys.getenv(
    "ACOUSTICTS_IMPL_PROVENANCE",
    unset = file.path(repo_root, ".tmp", "implementation-figures", "provenance.csv")
  )
  dir.create(dirname(provenance_path), recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(do.call(rbind, records), provenance_path, row.names = FALSE, na = "")
  message("Wrote implementation provenance to ", provenance_path)
  invisible(provenance_path)
}

# Resolve the requested scope and orchestrate one subprocess per family.
run_all <- function() {
  repo_root <- find_repo_root()
  setwd(repo_root)

  profile <- tolower(Sys.getenv("ACOUSTICTS_IMPL_PROFILE", unset = "light"))
  if (!identical(profile, "light")) {
    stop(
      "The package-side implementation runner supports only ",
      "'ACOUSTICTS_IMPL_PROFILE=light'. Heavy validation workflows live ",
      "outside this repository.",
      call. = FALSE
    )
  }

  manifest <- utils::read.csv(
    file.path("tools", "implementation-figures", "manifest.csv"),
    stringsAsFactors = FALSE
  )
  manifest <- manifest[tolower(manifest$profile) == profile, , drop = FALSE]

  families_csv <- Sys.getenv("ACOUSTICTS_IMPL_FAMILIES", unset = "")
  requested_families <- unique(tolower(trimws(
    strsplit(families_csv, ",", fixed = TRUE)[[1]]
  )))
  requested_families <- requested_families[nzchar(requested_families)]
  if (identical(requested_families, "all")) requested_families <- character()

  known_families <- sort(unique(tolower(manifest$family)))
  unknown_families <- setdiff(requested_families, known_families)
  if (length(unknown_families)) {
    stop(
      "Unknown implementation figure families: ",
      paste(unknown_families, collapse = ", "),
      call. = FALSE
    )
  }
  if (length(requested_families)) {
    manifest <- manifest[
      tolower(manifest$family) %in% requested_families,
      ,
      drop = FALSE
    ]
  }

  families <- unique(tolower(manifest$family))
  if (!length(families)) {
    message("No implementation figure families selected.")
    return(invisible(NULL))
  }

  worker <- normalizePath(
    file.path("tools", "implementation-figures", "run_family.R"),
    winslash = "/",
    mustWork = TRUE
  )
  rscript <- file.path(R.home("bin"), "Rscript")
  if (.Platform$OS.type == "windows") rscript <- paste0(rscript, ".exe")

  commit <- git_commit()
  generated_at <- format(Sys.time(), tz = "UTC", usetz = TRUE)
  records <- list()

  for (family in families) {
    rows <- manifest[tolower(manifest$family) == family, , drop = FALSE]
    scripts <- unique(file.path("tools", "implementation-figures", rows$script))
    missing_scripts <- scripts[!file.exists(scripts)]
    if (length(missing_scripts)) {
      stop(
        "Missing implementation builder scripts:\n- ",
        paste(missing_scripts, collapse = "\n- "),
        call. = FALSE
      )
    }

    message("Running implementation figure family: ", family)
    started <- proc.time()[["elapsed"]]
    status <- system2(
      rscript,
      c(
        "--vanilla",
        shQuote(worker),
        shQuote(repo_root),
        shQuote(family),
        shQuote(normalizePath(scripts, winslash = "/", mustWork = TRUE))
      )
    )
    elapsed <- proc.time()[["elapsed"]] - started

    if (!identical(status, 0L)) {
      stop(
        "Implementation figure family '", family,
        "' failed with exit status ", status, ".",
        call. = FALSE
      )
    }

    outputs <- validate_outputs(rows, family)
    output_info <- file.info(outputs)
    records[[family]] <- data.frame(
      family = family,
      output = outputs,
      bytes = output_info$size,
      md5 = unname(tools::md5sum(outputs)),
      elapsed_family_seconds = round(elapsed, 3),
      package_commit = commit,
      generated_at_utc = generated_at,
      stringsAsFactors = FALSE
    )
    message(
      "Completed ", family, " in ", round(elapsed, 2),
      " seconds (", length(outputs), " outputs)."
    )
  }

  write_provenance(records, repo_root)
  message("Implementation rebuild pass complete for profile: ", profile)
  invisible(TRUE)
}

run_all()
