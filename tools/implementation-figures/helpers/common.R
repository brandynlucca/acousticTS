impl_set_repo_root <- function() {
  if (file.exists("DESCRIPTION")) {
    return(invisible(normalizePath(".", winslash = "/", mustWork = FALSE)))
  }

  args <- commandArgs(FALSE)
  file_arg <- args[grep("^--file=", args)]

  if (length(file_arg) == 1L) {
    script_path <- normalizePath(sub("^--file=", "", file_arg))
    repo_root <- normalizePath(
      file.path(dirname(script_path), "..", "..", ".."),
      winslash = "/",
      mustWork = FALSE
    )
    setwd(repo_root)
  }

  invisible(normalizePath(".", winslash = "/", mustWork = FALSE))
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
