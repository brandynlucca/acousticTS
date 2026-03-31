find_repo_root <- function(start_dir) {
  search_dir <- normalizePath(start_dir, winslash = "/", mustWork = FALSE)

  while (nzchar(search_dir) && dir.exists(search_dir)) {
    desc_path <- file.path(search_dir, "DESCRIPTION")
    manifest_path <- file.path(
      search_dir,
      "tools",
      "implementation-figures",
      "manifest.csv"
    )
    if (file.exists(desc_path) && file.exists(manifest_path)) {
      return(normalizePath(search_dir, winslash = "/", mustWork = TRUE))
    }

    parent_dir <- dirname(search_dir)
    if (identical(parent_dir, search_dir)) {
      break
    }
    search_dir <- parent_dir
  }

  NULL
}

repo_root <- find_repo_root(getwd())

if (is.null(repo_root)) {
  workspace <- Sys.getenv("GITHUB_WORKSPACE", unset = "")
  if (nzchar(workspace)) {
    repo_root <- find_repo_root(workspace)
  }
}

if (is.null(repo_root)) {
  script_arg <- commandArgs(trailingOnly = FALSE)
  script_file <- sub(
    "^--file=",
    "",
    script_arg[grep("^--file=", script_arg)][1]
  )
  if (!is.na(script_file) && nzchar(script_file)) {
    repo_root <- find_repo_root(dirname(script_file))
  }
}

if (is.null(repo_root)) {
  stop(
    "Could not locate the acousticTS repository root from the current ",
    "working directory or workflow environment.",
    call. = FALSE
  )
}

setwd(repo_root)

manifest <- utils::read.csv(
  file.path("tools", "implementation-figures", "manifest.csv"),
  stringsAsFactors = FALSE
)

light_families <- sort(unique(tolower(manifest$family[manifest$profile == "light"])))
requested_profile <- tolower(
  Sys.getenv("ACOUSTICTS_IMPL_PROFILE", unset = "light")
)
event_name <- tolower(Sys.getenv("GITHUB_EVENT_NAME", unset = ""))
base_ref <- Sys.getenv("GITHUB_BASE_REF", unset = "")
manual_families <- tolower(
  trimws(Sys.getenv("ACOUSTICTS_IMPL_REQUESTED_FAMILIES", unset = ""))
)
output_path <- Sys.getenv("GITHUB_OUTPUT", unset = "")

shared_patterns <- c(
  "^tools/implementation-figures/(helpers/|manifest\\.csv$|run_all\\.R$)",
  "^tools/implementation-figures/build_impl_benchmark",
  "^\\.github/workflows/implementation-figures\\.ya?ml$"
)

escape_regex <- function(x) {
  gsub("([][{}()+*^$|\\\\?.])", "\\\\\\1", x)
}

family_patterns <- lapply(
  split(manifest, tolower(manifest$family)),
  function(rows) {
    vignette_paths <- unique(file.path(
      "vignettes",
      tolower(rows$family),
      rows$vignette
    ))
    script_paths <- unique(file.path(
      "tools",
      "implementation-figures",
      rows$script
    ))

    c(
      paste0("^", escape_regex(vignette_paths), "$"),
      paste0("^", escape_regex(script_paths), "$")
    )
  }
)

write_output <- function(name, value, path) {
  if (!nzchar(path)) {
    return(invisible(NULL))
  }
  cat(
    paste0(name, "=", value, "\n"),
    file = path,
    append = TRUE
  )
}

if (event_name == "workflow_dispatch") {
  if (!nzchar(manual_families) || identical(manual_families, "all")) {
    write_output("scope", "all", output_path)
    write_output("families", "", output_path)
  } else {
    requested <- unique(strsplit(manual_families, ",", fixed = TRUE)[[1]])
    requested <- trimws(requested)
    requested <- requested[nzchar(requested)]

    unknown <- setdiff(requested, light_families)
    if (length(unknown) > 0L) {
      stop(
        "Unknown implementation figure families requested: ",
        paste(unknown, collapse = ", "),
        call. = FALSE
      )
    }

    write_output("scope", "selected", output_path)
    write_output("families", paste(requested, collapse = ","), output_path)
  }

  write_output("profile", requested_profile, output_path)
  quit(save = "no", status = 0)
}

if (event_name != "pull_request" || !nzchar(base_ref)) {
  write_output("scope", "all", output_path)
  write_output("families", "", output_path)
  write_output("profile", requested_profile, output_path)
  quit(save = "no", status = 0)
}

changed_files <- system2(
  "git",
  c("diff", "--name-only", paste0("origin/", base_ref, "...HEAD")),
  stdout = TRUE
)
changed_files <- changed_files[nzchar(changed_files)]

if (length(changed_files) == 0L) {
  write_output("scope", "skip", output_path)
  write_output("families", "", output_path)
  write_output("profile", requested_profile, output_path)
  quit(save = "no", status = 0)
}

if (any(vapply(
  shared_patterns,
  function(pattern) any(grepl(pattern, changed_files)),
  logical(1)
))) {
  write_output("scope", "all", output_path)
  write_output("families", "", output_path)
  write_output("profile", requested_profile, output_path)
  quit(save = "no", status = 0)
}

selected_families <- names(family_patterns)[vapply(
  family_patterns,
  function(patterns) {
    any(vapply(
      patterns,
      function(pattern) any(grepl(pattern, changed_files)),
      logical(1)
    ))
  },
  logical(1)
)]
selected_families <- sort(unique(selected_families))

if (length(selected_families) == 0L) {
  write_output("scope", "skip", output_path)
  write_output("families", "", output_path)
  write_output("profile", requested_profile, output_path)
  quit(save = "no", status = 0)
}

write_output("scope", "selected", output_path)
write_output("families", paste(selected_families, collapse = ","), output_path)
write_output("profile", requested_profile, output_path)
