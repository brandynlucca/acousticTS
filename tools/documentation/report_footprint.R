# Report documentation source, generated-site, and figure-tooling footprints.
# Run before and after hosting changes to keep package and website costs distinct.

# Locate the package repository from the current working directory.
find_repo_root <- function(start = getwd()) {
  current <- normalizePath(start, winslash = "/", mustWork = FALSE)

  repeat {
    if (file.exists(file.path(current, "DESCRIPTION"))) {
      return(current)
    }

    parent <- dirname(current)
    if (identical(parent, current)) {
      stop("Could not locate the acousticTS repository root.", call. = FALSE)
    }
    current <- parent
  }
}

# Count files and total bytes beneath one documentation-related directory.
directory_summary <- function(path) {
  if (!dir.exists(path)) {
    return(data.frame(path = path, files = 0L, mebibytes = 0))
  }

  files <- list.files(path, recursive = TRUE, full.names = TRUE, all.files = TRUE)
  files <- files[file.info(files)$isdir %in% FALSE]
  bytes <- if (length(files)) sum(file.info(files)$size, na.rm = TRUE) else 0

  data.frame(
    path = path,
    files = length(files),
    mebibytes = round(bytes / 1024^2, 3)
  )
}

repo_root <- find_repo_root()
setwd(repo_root)

vignette_root_files <- list.files(
  "vignettes",
  full.names = TRUE,
  recursive = FALSE,
  all.files = TRUE
)
vignette_root_files <- vignette_root_files[
  file.info(vignette_root_files)$isdir %in% FALSE
]
root_bytes <- if (length(vignette_root_files)) {
  sum(file.info(vignette_root_files)$size, na.rm = TRUE)
} else {
  0
}

summary <- do.call(
  rbind,
  lapply(
    c("vignettes", "tools/implementation-figures", "docs"),
    directory_summary
  )
)
summary <- rbind(
  summary,
  data.frame(
    path = "vignettes (top-level files only)",
    files = length(vignette_root_files),
    mebibytes = round(root_bytes / 1024^2, 3)
  )
)

cat("Documentation footprint\n")
cat("Repository:", repo_root, "\n\n")
print(summary, row.names = FALSE)

vignette_files <- list.files("vignettes", recursive = TRUE, full.names = TRUE)
vignette_files <- vignette_files[file.info(vignette_files)$isdir %in% FALSE]
if (length(vignette_files)) {
  extensions <- tolower(tools::file_ext(vignette_files))
  extensions[!nzchar(extensions)] <- "(none)"
  sizes <- file.info(vignette_files)$size
  by_extension <- aggregate(
    sizes,
    by = list(extension = extensions),
    FUN = function(x) round(sum(x, na.rm = TRUE) / 1024^2, 3)
  )
  counts <- as.data.frame(table(extensions), stringsAsFactors = FALSE)
  names(counts) <- c("extension", "files")
  by_extension <- merge(counts, by_extension, by = "extension", all = TRUE)
  names(by_extension)[names(by_extension) == "x"] <- "mebibytes"
  by_extension <- by_extension[order(by_extension$mebibytes, decreasing = TRUE), ]

  cat("\nVignette files by extension\n")
  print(by_extension, row.names = FALSE)
}

cat(
  "\nUse this report before and after documentation changes. ",
  "Time R CMD build, R CMD INSTALL, pkgdown, and figure generation ",
  "separately so their costs are not conflated.\n",
  sep = ""
)
