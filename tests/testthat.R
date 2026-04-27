library(testthat)
library(acousticTS)

local({
  close_test_devices <- function() {
    while (grDevices::dev.cur() > 1L) {
      try(grDevices::dev.off(), silent = TRUE)
    }
  }

  on.exit(
    {
      close_test_devices()
      invisible(gc())
    },
    add = TRUE
  )

  test_check("acousticTS", reporter = testthat::SummaryReporter$new())
  close_test_devices()
  invisible(gc())
})
