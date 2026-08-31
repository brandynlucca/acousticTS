# Build a reduced pkgdown smoke site for pull requests.
# Site-wide inputs trigger a full build; otherwise only affected articles plus
# the home/reference surfaces are rendered.

base_ref <- Sys.getenv("GITHUB_BASE_REF", unset = "")
if (!nzchar(base_ref)) {
  stop("GITHUB_BASE_REF is required for a pull-request pkgdown build.", call. = FALSE)
}

changed_files <- system2(
  "git",
  c("diff", "--name-only", paste0("origin/", base_ref, "...HEAD")),
  stdout = TRUE
)
changed_files <- gsub("\\\\", "/", changed_files[nzchar(changed_files)])

full_site_patterns <- c(
  "^_pkgdown\\.ya?ml$",
  "^pkgdown/",
  "^vignettes/REFERENCES\\.bib$",
  "^R/pkgdown-helpers\\.R$",
  "^README\\.(Rmd|md)$",
  "^NEWS\\.md$"
)
requires_full_site <- any(vapply(
  full_site_patterns,
  function(pattern) any(grepl(pattern, changed_files)),
  logical(1)
))

if (requires_full_site) {
  message("Site-wide documentation input changed; building the complete site.")
  pkgdown::build_site_github_pages(
    new_process = FALSE,
    install = FALSE,
    clean = TRUE
  )
  quit(save = "no", status = 0)
}

article_sources <- list.files(
  "vignettes",
  pattern = "\\.Rmd$",
  recursive = TRUE,
  full.names = TRUE
)
article_sources <- gsub("\\\\", "/", article_sources)

# Convert a vignette source path to pkgdown's article identifier.
article_name <- function(path) {
  sub("\\.Rmd$", "", sub("^vignettes/", "", path))
}

articles <- article_name(changed_files[grepl("^vignettes/.*\\.Rmd$", changed_files)])

changed_vignette_assets <- changed_files[
  grepl("^vignettes/", changed_files) & !grepl("\\.Rmd$", changed_files)
]
if (length(changed_vignette_assets)) {
  for (asset in changed_vignette_assets) {
    asset_name <- basename(asset)
    referenced_by <- article_sources[vapply(
      article_sources,
      function(source) {
        any(grepl(asset_name, readLines(source, warn = FALSE), fixed = TRUE))
      },
      logical(1)
    )]
    articles <- c(articles, article_name(referenced_by))
  }
}

manifest <- utils::read.csv(
  file.path("tools", "implementation-figures", "manifest.csv"),
  stringsAsFactors = FALSE
)
# Split semicolon-delimited manifest inputs into normalized repository paths.
split_values <- function(values) {
  values <- values[!is.na(values) & nzchar(values)]
  if (!length(values)) return(character())
  values <- trimws(unlist(strsplit(values, ";", fixed = TRUE)))
  unique(gsub("\\\\", "/", values[nzchar(values)]))
}

for (i in seq_len(nrow(manifest))) {
  row <- manifest[i, , drop = FALSE]
  builder <- paste0("tools/implementation-figures/", row$script)
  inputs <- split_values(row$inputs)
  if (any(changed_files %in% c(builder, inputs))) {
    articles <- c(
      articles,
      paste(tolower(row$family), sub("\\.Rmd$", "", row$vignette), sep = "/")
    )
  }
}

if (any(grepl("^(R|src)/", changed_files))) {
  articles <- c(
    articles,
    "getting-started/getting-started",
    "model-library/model-library"
  )
}

articles <- sort(unique(articles[nzchar(articles)]))
articles <- intersect(articles, article_name(article_sources))

message("Building pkgdown home smoke test.")
pkgdown::build_home(quiet = FALSE)

reference_changed <- any(grepl(
  "^(R/|src/|man/|DESCRIPTION$|NAMESPACE$)",
  changed_files
))
if (reference_changed) {
  message("Package/reference input changed; building the function reference.")
  pkgdown::build_reference(lazy = FALSE, devel = FALSE)
}

if (length(articles)) {
  message("Building affected articles: ", paste(articles, collapse = ", "))
  for (article in articles) {
    pkgdown::build_article(
      article,
      lazy = FALSE,
      new_process = FALSE,
      quiet = FALSE
    )
  }
  pkgdown::build_articles_index()
} else {
  message("No affected articles detected; pull-request smoke build is complete.")
}
