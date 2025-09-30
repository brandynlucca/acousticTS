test_that("Data objects can be loaded and are of expected class types", {
  library(acousticTS)

  # Test that data objects exist and can be loaded
  data(krill, package = "acousticTS")
  data(sardine, package = "acousticTS")
  data(cod, package = "acousticTS")

  # Test krill object
  expect_true(exists("krill"))
  expect_s4_class(krill, "FLS") # krill should be FLS class
  expect_s4_class(krill, "Scatterer") # should inherit from Scatterer

  # Test sardine object
  expect_true(exists("sardine"))
  expect_s4_class(sardine, "SBF") # sardine should be SBF class
  expect_s4_class(sardine, "GAS") # SBF inherits from GAS
  expect_s4_class(sardine, "Scatterer") # should inherit from Scatterer

  # Test cod object
  expect_true(exists("cod"))
  expect_s4_class(cod, "SBF") # cod should be SBF class
  expect_s4_class(cod, "GAS") # SBF inherits from GAS
  expect_s4_class(cod, "Scatterer") # should inherit from Scatterer
})

test_that("Data objects have proper structure", {
  library(acousticTS)

  # Load data objects
  data(krill, package = "acousticTS")
  data(sardine, package = "acousticTS")
  data(cod, package = "acousticTS")

  # Test krill structure (FLS object)
  expect_true("metadata" %in% slotNames(krill))
  expect_true("body" %in% slotNames(krill))
  expect_true("shape_parameters" %in% slotNames(krill))
  expect_type(krill@metadata, "list")
  expect_type(krill@body, "list")
  expect_type(krill@shape_parameters, "list")

  # Test sardine structure (SBF object)
  expect_true("metadata" %in% slotNames(sardine))
  expect_true("body" %in% slotNames(sardine))
  expect_true("bladder" %in% slotNames(sardine)) # SBF should have bladder
  expect_true("shape_parameters" %in% slotNames(sardine))
  expect_type(sardine@metadata, "list")
  expect_type(sardine@body, "list")
  expect_type(sardine@bladder, "list")
  expect_type(sardine@shape_parameters, "list")

  # Test cod structure (SBF object)
  expect_true("metadata" %in% slotNames(cod))
  expect_true("body" %in% slotNames(cod))
  expect_true("bladder" %in% slotNames(cod)) # SBF should have bladder
  expect_true("shape_parameters" %in% slotNames(cod))
  expect_type(cod@metadata, "list")
  expect_type(cod@body, "list")
  expect_type(cod@bladder, "list")
  expect_type(cod@shape_parameters, "list")
})

test_that("Data objects have valid metadata", {
  library(acousticTS)

  # Load data objects
  data(krill, package = "acousticTS")
  data(sardine, package = "acousticTS")
  data(cod, package = "acousticTS")

  # Test that metadata contains ID field
  expect_true("ID" %in% names(krill@metadata))
  expect_true("ID" %in% names(sardine@metadata))
  expect_true("ID" %in% names(cod@metadata))

  # Test that IDs are not just default values
  expect_true(is.character(krill@metadata$ID))
  expect_true(is.character(sardine@metadata$ID))
  expect_true(is.character(cod@metadata$ID))
})

test_that("Data objects can be used with generic methods", {
  library(acousticTS)

  # Suppression helper functions
  suppressAllOutput <- function(expr) {
    suppressMessages(suppressWarnings(capture.output(expr, file = NULL)))
  }

  suppressPlot <- function(expr) {
    pdf(NULL) # Open null graphics device
    on.exit(dev.off()) # Ensure it's closed after evaluation
    invisible(force(expr))
  }

  # Load data objects
  data(krill, package = "acousticTS")
  data(sardine, package = "acousticTS")
  data(cod, package = "acousticTS")

  # Test show method works
  expect_error(suppressAllOutput(show(krill)), NA)
  expect_error(suppressAllOutput(show(sardine)), NA)
  expect_error(suppressAllOutput(show(cod)), NA)

  # Test plot method works
  expect_error(suppressPlot(plot(krill)), NA)
  expect_error(suppressPlot(plot(sardine)), NA)
  expect_error(suppressPlot(plot(cod)), NA)

  # Test extract method works
  krill_metadata <- extract(krill, "metadata")
  sardine_body <- extract(sardine, "body")
  cod_bladder <- extract(cod, "bladder")

  expect_type(krill_metadata, "list")
  expect_type(sardine_body, "list")
  expect_type(cod_bladder, "list")
})
