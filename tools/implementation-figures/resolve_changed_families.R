if (file.exists("DESCRIPTION")) {
  repo_root <- normalizePath(".", winslash = "/", mustWork = TRUE)
} else {
  script_arg <- commandArgs(trailingOnly = FALSE)
  script_file <- sub(
    "^--file=",
    "",
    script_arg[grep("^--file=", script_arg)][1]
  )
  if (is.na(script_file) || !nzchar(script_file)) {
    script_file <- file.path(
      "tools",
      "implementation-figures",
      "resolve_changed_families.R"
    )
  }
  repo_root <- normalizePath(
    file.path(dirname(script_file), "..", ".."),
    winslash = "/",
    mustWork = FALSE
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
  "^src/",
  "^R/(acoustics|create_scatterer|create_shape|scatterer.*|utilities.*)\\.R$",
  "^tools/implementation-figures/(helpers/|manifest\\.csv$|run_all\\.R$)",
  "^tools/implementation-figures/build_impl_benchmark"
)

family_patterns <- list(
  bbfm = c(
    "^R/model-bbfm\\.R$",
    "^vignettes/bbfm/",
    "^tools/implementation-figures/bbfm/"
  ),
  bcms = c(
    "^R/model-bcms\\.R$",
    "^vignettes/bcms/",
    "^tools/implementation-figures/bcms/",
    "^tools/implementation-figures/bcms_compare_reference\\.R$"
  ),
  calibration = c(
    "^R/model-soems\\.R$",
    "^vignettes/calibration/",
    "^tools/implementation-figures/calibration/"
  ),
  dwba = c(
    "^R/model-dwba\\.R$",
    "^vignettes/dwba/",
    "^tools/implementation-figures/dwba/"
  ),
  ecms = c(
    "^R/model-ecms\\.R$",
    "^vignettes/ecms/",
    "^tools/implementation-figures/ecms/",
    "^tools/implementation-figures/ecms_compare_reference\\.R$"
  ),
  essms = c(
    "^R/model-essms\\.R$",
    "^vignettes/essms/",
    "^tools/implementation-figures/essms/"
  ),
  fcms = c(
    "^R/model-fcms\\.R$",
    "^vignettes/fcms/",
    "^tools/implementation-figures/fcms/"
  ),
  hpa = c(
    "^R/model-hpa\\.R$",
    "^vignettes/hpa/",
    "^tools/implementation-figures/hpa/"
  ),
  krm = c(
    "^R/model-krm\\.R$",
    "^vignettes/krm/",
    "^tools/implementation-figures/krm/"
  ),
  pcdwba = c(
    "^R/model-pcdwba\\.R$",
    "^vignettes/pcdwba/",
    "^tools/implementation-figures/pcdwba/",
    "^tools/implementation-figures/pcdwba_compare_summary\\.R$"
  ),
  psms = c(
    "^R/model-psms\\.R$",
    "^vignettes/psms/",
    "^tools/implementation-figures/psms/"
  ),
  sdwba = c(
    "^R/model-sdwba\\.R$",
    "^vignettes/sdwba/",
    "^tools/implementation-figures/sdwba/"
  ),
  sphms = c(
    "^R/model-sphms\\.R$",
    "^vignettes/sphms/",
    "^tools/implementation-figures/sphms/"
  ),
  tmm = c(
    "^R/model-tmm",
    "^R/tmm-",
    "^vignettes/tmm/",
    "^tools/implementation-figures/tmm/"
  ),
  trcm = c(
    "^R/model-trcm\\.R$",
    "^vignettes/trcm/",
    "^tools/implementation-figures/trcm/"
  ),
  vesm = c(
    "^R/model-vesms\\.R$",
    "^vignettes/vesm/",
    "^tools/implementation-figures/vesm/",
    "^tools/implementation-figures/vesm_compare_summary\\.R$"
  )
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
