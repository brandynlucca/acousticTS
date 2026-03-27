find_repo_root <- function() {
  has_impl_manifest <- function(path) {
    file.exists(file.path(path, "DESCRIPTION")) &&
      file.exists(file.path(
        path,
        "tools",
        "implementation-figures",
        "manifest.csv"
      ))
  }

  parent_dir <- function(path) {
    normalizePath(file.path(path, ".."), winslash = "/", mustWork = FALSE)
  }

  search_upward <- function(start_path) {
    current <- normalizePath(start_path, winslash = "/", mustWork = FALSE)
    repeat {
      if (has_impl_manifest(current)) {
        return(current)
      }

      next_path <- parent_dir(current)
      if (identical(next_path, current)) {
        return(NULL)
      }
      current <- next_path
    }
  }

  search_roots <- character()

  if (nzchar(Sys.getenv("GITHUB_WORKSPACE", unset = ""))) {
    search_roots <- c(search_roots, Sys.getenv("GITHUB_WORKSPACE"))
  }

  search_roots <- c(search_roots, getwd())

  script_arg <- commandArgs(trailingOnly = FALSE)
  file_arg <- script_arg[grep("^--file=", script_arg)]
  if (length(file_arg) == 1L) {
    search_roots <- c(search_roots, dirname(sub("^--file=", "", file_arg)))
  }

  search_roots <- unique(search_roots[nzchar(search_roots)])
  for (root in search_roots) {
    found_root <- search_upward(root)
    if (!is.null(found_root)) {
      return(found_root)
    }
  }

  stop(
    "Could not locate the repository root containing ",
    "'tools/implementation-figures/manifest.csv'.",
    call. = FALSE
  )
}

repo_root <- find_repo_root()
setwd(repo_root)

profile <- tolower(
  Sys.getenv("ACOUSTICTS_IMPL_PROFILE", unset = "light")
)
families_csv <- Sys.getenv("ACOUSTICTS_IMPL_FAMILIES", unset = "")
manifest <- utils::read.csv(
  file.path("tools", "implementation-figures", "manifest.csv"),
  stringsAsFactors = FALSE
)

allowed_profiles <- switch(
  profile,
  light = "light",
  all = c("light", "heavy"),
  stop("Unknown implementation figure profile: ", profile, call. = FALSE)
)

if (nzchar(families_csv)) {
  requested_families <- unique(
    tolower(trimws(strsplit(families_csv, ",", fixed = TRUE)[[1]]))
  )
  requested_families <- requested_families[nzchar(requested_families)]
  known_families <- sort(unique(tolower(manifest$family)))

  unknown_families <- setdiff(requested_families, known_families)
  if (length(unknown_families) > 0L) {
    stop(
      "Unknown implementation figure families: ",
      paste(unknown_families, collapse = ", "),
      call. = FALSE
    )
  }

  manifest <- manifest[
    tolower(manifest$family) %in% requested_families,
    ,
    drop = FALSE
  ]
}

scripts <- unique(
  file.path(
    "tools",
    "implementation-figures",
    manifest$script[manifest$profile %in% allowed_profiles]
  )
)

if (length(scripts) == 0L) {
  message(
    "No implementation figures selected for profile: ",
    profile,
    if (nzchar(families_csv)) {
      paste0(" (families: ", families_csv, ")")
    } else {
      ""
    }
  )
  quit(save = "no", status = 0)
}

for (script in scripts) {
  message("Running ", script)
  sys.source(script, envir = new.env(parent = globalenv()))
}

message("Implementation rebuild pass complete for profile: ", profile)
