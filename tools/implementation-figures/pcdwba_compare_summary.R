source("tools/implementation-figures/helpers/common.R")
impl_set_repo_root()
data_dir <- impl_data_path("")
pcdwba_acoustic <- read.csv(file.path(data_dir, "pcdwba_reference_acousticts.csv"))
pcdwba_echo <- read.csv(file.path(data_dir, "pcdwba_reference_echopop.csv"))
pcdwba_zoo <- read.csv(file.path(data_dir, "pcdwba_reference_zooscatr.csv"))

timing_acoustic <- read.csv(file.path(data_dir, "pcdwba_reference_acousticts_timing.csv"))
timing_echo <- read.csv(file.path(data_dir, "pcdwba_reference_echopop_timing.csv"))
timing_zoo <- read.csv(file.path(data_dir, "pcdwba_reference_zooscatr_timing.csv"))

merged <- Reduce(
  function(x, y) merge(x, y, by = c("frequency_hz", "frequency_khz", "ka")),
  list(
    setNames(
      pcdwba_acoustic,
      c("frequency_hz", "frequency_khz", "ka", "f_bs_real_acousticts",
        "f_bs_imag_acousticts", "sigma_bs_acousticts", "TS_acousticts")
    ),
    setNames(
      pcdwba_echo,
      c("frequency_hz", "frequency_khz", "ka", "f_bs_real_echopop",
        "f_bs_imag_echopop", "sigma_bs_echopop", "TS_echopop")
    ),
    setNames(
      pcdwba_zoo,
      c("frequency_hz", "frequency_khz", "ka", "f_bs_real_zooscatr",
        "f_bs_imag_zooscatr", "sigma_bs_zooscatr", "TS_zooscatr")
    )
  )
)
merged <- impl_sort_compare_df(merged)

write.csv(
  merged,
  file.path(data_dir, "pcdwba_reference_compare.csv"),
  row.names = FALSE
)

pairwise <- data.frame(
  comparison = c(
    "acousticTS vs echopop",
    "acousticTS vs ZooScatR-source",
    "echopop vs ZooScatR-source"
  ),
  max_abs_delta_TS_dB = c(
    max(abs(merged$TS_acousticts - merged$TS_echopop)),
    max(abs(merged$TS_acousticts - merged$TS_zooscatr)),
    max(abs(merged$TS_echopop - merged$TS_zooscatr))
  ),
  mean_abs_delta_TS_dB = c(
    mean(abs(merged$TS_acousticts - merged$TS_echopop)),
    mean(abs(merged$TS_acousticts - merged$TS_zooscatr)),
    mean(abs(merged$TS_echopop - merged$TS_zooscatr))
  )
)

timing_long <- rbind(timing_acoustic, timing_echo, timing_zoo)

write.csv(
  pairwise,
  file.path(data_dir, "pcdwba_reference_compare_summary.csv"),
  row.names = FALSE
)

write.csv(
  timing_long,
  file.path(data_dir, "pcdwba_reference_compare_timing.csv"),
  row.names = FALSE
)
