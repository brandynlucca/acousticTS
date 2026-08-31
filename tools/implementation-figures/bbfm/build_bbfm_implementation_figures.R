# Build the BBFM shape and example-spectrum figures declared in the manifest.
# This wrapper is the supported entry point for the BBFM implementation article.
source("tools/implementation-figures/helpers/common.R")
impl_load_all()
source("tools/implementation-figures/bbfm/build_bbfm_shape_plot.R")
source("tools/implementation-figures/bbfm/build_bbfm_implementation_plot.R")
