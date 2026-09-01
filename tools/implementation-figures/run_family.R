# Execute all builders for one model family inside a disposable R process.
# This worker is invoked by run_all.R rather than called directly by users.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3L) {
  stop(
    "Usage: Rscript run_family.R <repo-root> <family> <script> [<script> ...]",
    call. = FALSE
  )
}

repo_root <- normalizePath(args[[1]], winslash = "/", mustWork = TRUE)
family <- args[[2]]
scripts <- args[-c(1L, 2L)]
setwd(repo_root)

for (script in scripts) {
  script <- normalizePath(script, winslash = "/", mustWork = TRUE)
  message("[", family, "] Running ", script)
  environment <- new.env(parent = globalenv())
  sys.source(script, envir = environment)
  rm(environment)
  invisible(gc())
}
