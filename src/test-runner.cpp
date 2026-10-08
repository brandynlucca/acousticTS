#include <Rcpp.h>
#define TESTTHAT_TEST_RUNNER
#include <testthat.h>

// Keep the runner internal. The assertions live with the package tests and
// are included in the translation units containing the private kernels.
// [[Rcpp::export]]
bool native_kernel_tests_cpp() {
#ifdef TESTTHAT_DISABLED
    return true;
#else
    return testthat::run_tests(false);
#endif
}
