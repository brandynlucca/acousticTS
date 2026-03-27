source("tools/implementation-figures/helpers/common.R")
impl_set_repo_root()
data_dir <- impl_data_path("")
acousticts <- read.csv(file.path(data_dir, "vesm_reference_acousticts.csv"))
original <- read.csv(file.path(data_dir, "vesm_reference_original.csv"))

timing_acousticts <- read.csv(
  file.path(data_dir, "vesm_reference_acousticts_timing.csv")
)
timing_original <- read.csv(
  file.path(data_dir, "vesm_reference_original_timing.csv")
)

compare <- merge(
  acousticts[, c("frequency", "TS")],
  original,
  by = "frequency",
  suffixes = c("_acousticts", "_original")
)
compare <- impl_sort_compare_df(compare)

compare$delta_TS <- compare$TS_acousticts - compare$TS_original
compare$abs_delta_TS <- abs(compare$delta_TS)

write.csv(
  compare,
  file = file.path(data_dir, "vesm_reference_compare.csv"),
  row.names = FALSE
)

max_idx <- which.max(compare$abs_delta_TS)

summary_df <- data.frame(
  comparison = "acousticTS vs original VESM",
  n_frequency = nrow(compare),
  frequency_min_kHz = min(compare$frequency) * 1e-3,
  frequency_max_kHz = max(compare$frequency) * 1e-3,
  max_abs_delta_TS_dB = compare$abs_delta_TS[max_idx],
  mean_abs_delta_TS_dB = mean(compare$abs_delta_TS),
  frequency_at_max_delta_kHz = compare$frequency[max_idx] * 1e-3,
  elapsed_acousticts_s = timing_acousticts$elapsed_s[1],
  elapsed_original_s = timing_original$elapsed_s[1]
)

write.csv(
  summary_df,
  file = file.path(data_dir, "vesm_reference_compare_summary.csv"),
  row.names = FALSE
)
