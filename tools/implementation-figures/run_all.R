script_arg <- commandArgs(trailingOnly = FALSE)
script_file <- sub(
  "^--file=",
  "",
  script_arg[grep("^--file=", script_arg)][1]
)
if (is.na(script_file) || !nzchar(script_file)) {
  script_file <- file.path("tools", "implementation-figures", "run_all.R")
}
repo_root <- normalizePath(
  file.path(dirname(script_file), "..", ".."),
  winslash = "/",
  mustWork = FALSE
)
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
