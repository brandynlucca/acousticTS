impl_set_repo_root <- function() {
  has_repo_root <- function(path) {
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
      if (has_repo_root(current)) {
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

  args <- commandArgs(FALSE)
  file_arg <- args[grep("^--file=", args)]
  if (length(file_arg) == 1L) {
    search_roots <- c(search_roots, dirname(sub("^--file=", "", file_arg)))
  }

  search_roots <- unique(search_roots[nzchar(search_roots)])
  for (root in search_roots) {
    repo_root <- search_upward(root)
    if (!is.null(repo_root)) {
      setwd(repo_root)
      return(invisible(repo_root))
    }
  }

  stop(
    "Could not locate the repository root containing ",
    "'tools/implementation-figures/manifest.csv'.",
    call. = FALSE
  )
}

impl_load_all <- function() {
  impl_set_repo_root()
  devtools::load_all(".", quiet = TRUE)
  invisible(TRUE)
}

impl_data_path <- function(name) {
  file.path("tools", "implementation-figures", "data", name)
}

impl_output_path <- function(family, name) {
  file.path("vignettes", family, name)
}

impl_with_png <- function(path,
                          code,
                          width = 1600,
                          height = 1200,
                          res = 200) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  grDevices::png(path, width = width, height = height, res = res)
  on.exit(grDevices::dev.off(), add = TRUE)
  force(code)
  invisible(path)
}

impl_ts_label <- function() {
  expression(Target ~ strength ~ (dB ~ re. ~ 1 ~ m^2))
}

impl_delta_label <- function() {
  expression(Delta ~ TS ~ (dB))
}

impl_sort_compare_df <- function(df) {
  sort_cols <- intersect(
    c(
      "case",
      "frequency_hz",
      "frequency",
      "frequency_khz",
      "angle_deg",
      "theta_deg",
      "phi_deg",
      "ka"
    ),
    names(df)
  )

  if (!length(sort_cols)) {
    return(df)
  }

  order_args <- lapply(sort_cols, function(col) df[[col]])
  row_idx <- do.call(order, c(order_args, list(na.last = TRUE)))
  df[row_idx, , drop = FALSE]
}

impl_should_refresh_timings <- function() {
  refresh_override <- tolower(
    Sys.getenv("ACOUSTICTS_IMPL_REFRESH_TIMINGS", unset = "")
  )

  identical(refresh_override, "true") ||
    !nzchar(Sys.getenv("CI", unset = ""))
}

impl_round_timing_columns <- function(df, digits = 6) {
  timing_cols <- grep("elapsed", names(df), value = TRUE)

  for (col in timing_cols) {
    if (is.numeric(df[[col]])) {
      df[[col]] <- round(df[[col]], digits = digits)
    }
  }

  df
}
