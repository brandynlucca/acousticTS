#include <Rcpp.h>
#include <cmath>
#include <complex>
#include <vector>
#include "bessel_helpers.h"

using namespace Rcpp;

// Evaluate complete cylindrical Bessel sequences once per frequency.
// [[Rcpp::export]]
ComplexMatrix fcms_bm_cpp(NumericVector ka, IntegerVector m_limit,
                           bool rigid) {
    if (m_limit.size() == 0 || ka.size() == 0)
        stop("Cylinder size parameters and modal limits must be non-empty.");
    int m_max = 0;
    for (int ml : m_limit) {
        if (ml == NA_INTEGER || ml < 0)
            stop("Cylinder modal limits must be non-negative integers.");
        m_max = std::max(m_max, ml);
    }
    int nf = std::max(ka.size(), m_limit.size());
    // Pad unused modes so callers retain their existing matrix dimensions.
    ComplexMatrix out(m_max + 1, nf);
    std::fill(out.begin(), out.end(), Rcomplex{NA_REAL, NA_REAL});
    std::vector<double> j(m_max + 2), y(m_max + 2);
    for (int i = 0; i < nf; ++i) {
        if (i % 64 == 0) Rcpp::checkUserInterrupt();
        int ml = m_limit[i % m_limit.size()];
        double x = ka[i % ka.size()];
        if (!(x > 0.0) || !std::isfinite(x)) {
            // Preserve the general helper's continuation at nonpositive inputs.
            for (int m = 0; m <= ml; ++m) {
                std::complex<double> jv = jc_single_impl({x, 0.0}, m);
                std::complex<double> yv = yc_single_impl({x, 0.0}, m);
                if (rigid) {
                    jv = 0.5 * (jc_single_impl({x, 0.0}, m - 1) -
                                jc_single_impl({x, 0.0}, m + 1));
                    yv = 0.5 * (yc_single_impl({x, 0.0}, m - 1) -
                                yc_single_impl({x, 0.0}, m + 1));
                }
                out(m,i) = to_Rcomplex(jv / (jv + std::complex<double>(0,1) * yv));
            }
            continue;
        }
        // Include one extra order for the centered derivative recurrence.
        R::bessel_j_ex(x, ml + 1.0, j.data());
        R::bessel_y_ex(x, ml + 1.0, y.data());
        for (int m = 0; m <= ml; ++m) {
            // Rigid boundaries use derivatives, soft boundaries use values.
            double jv = j[m], yv = y[m];
            if (rigid) {
                jv = m == 0 ? -j[1] : 0.5 * (j[m-1] - j[m+1]);
                yv = m == 0 ? -y[1] : 0.5 * (y[m-1] - y[m+1]);
            }
            out(m,i) = to_Rcomplex(jv / std::complex<double>(jv, yv));
        }
    }
    return out;
}

