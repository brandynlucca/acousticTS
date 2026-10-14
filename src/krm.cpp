#include <Rcpp.h>
#include <cmath>
#include <complex>
#include <vector>
#include "bessel_helpers.h"

using namespace Rcpp;

// Evaluate KRM segment sums without frequency-by-segment temporary matrices.
// [[Rcpp::export]]
ComplexVector krm_body_fbs_cpp(NumericMatrix rpos, NumericVector k_sw,
                               NumericVector k_body, double theta,
                               double reflection) {
    int nf = k_sw.size();
    int ns = rpos.ncol() - 1;
    if (rpos.nrow() != 4 || ns < 1 || k_body.size() != nf)
        stop("Invalid KRM body geometry or wavenumbers.");
    ComplexVector out(nf);
    double st = std::sin(theta), ct = std::cos(theta);
    double transmission = 1.0 - reflection * reflection;
    double factor = reflection / (2.0 * std::sqrt(M_PI));
    std::vector<double> a(ns), vu(ns), vl(ns), du(ns), dorsal(ns);
    // Rotate adjacent body nodes once for all frequencies.
    for (int j = 0; j < ns; ++j) {
        double x = rpos(0,j) + rpos(0,j+1);
        double u = rpos(2,j) + rpos(2,j+1);
        double l = rpos(3,j) + rpos(3,j+1);
        a[j] = (rpos(1,j) + rpos(1,j+1)) / 4.0;
        vu[j] = (x * ct + u * st) / 2.0;
        vl[j] = (x * ct + l * st) / 2.0;
        du[j] = (rpos(0,j+1) - rpos(0,j)) * st;
        dorsal[j] = u / 2.0;
    }
    for (int i = 0; i < nf; ++i) {
        if (i % 64 == 0) Rcpp::checkUserInterrupt();
        // Match R's extended-precision rowSums accumulation.
        long double re = 0.0L, im = 0.0L;
        for (int j = 0; j < ns; ++j) {
            // Combine dorsal reflection and the transmitted ventral return.
            double q = k_body[i] * dorsal[j];
            double psi = -M_PI * q / (2.0 * (q + 0.4));
            double p = -2.0 * k_sw[i] * vu[j];
            double p2 = p + 2.0 * k_body[i] * (vu[j] - vl[j]) + psi;
            double weight = factor * std::sqrt(k_sw[i] * a[j]) * du[j];
            re += weight * (std::sin(p) - transmission * std::sin(p2));
            im -= weight * (std::cos(p) - transmission * std::cos(p2));
        }
        out[i] = to_Rcomplex({static_cast<double>(re), static_cast<double>(im)});
    }
    return out;
}

// [[Rcpp::export]]
ComplexVector krm_bladder_fbs_cpp(NumericMatrix rpos, NumericVector k_sw,
                                  NumericVector k_medium, double theta,
                                  double reflection, double transmission) {
    int nf = k_sw.size();
    int ns = rpos.ncol() - 1;
    if (rpos.nrow() != 4 || ns < 1 || k_medium.size() != nf)
        stop("Invalid KRM bladder geometry or wavenumbers.");
    ComplexVector out(nf);
    double st = std::sin(theta), ct = std::cos(theta);
    double factor = reflection * transmission / (2.0 * std::sqrt(M_PI));
    std::vector<double> a(ns), v(ns), du(ns);
    // Rotate the swimbladder dorsal surface and projected segment lengths.
    for (int j = 0; j < ns; ++j) {
        a[j] = (rpos(1,j) + rpos(1,j+1)) / 4.0;
        v[j] = ((rpos(0,j) + rpos(0,j+1)) * ct +
                (rpos(2,j) + rpos(2,j+1)) * st) / 2.0;
        du[j] = (rpos(0,j+1) - rpos(0,j)) * st;
    }
    for (int i = 0; i < nf; ++i) {
        if (i % 64 == 0) Rcpp::checkUserInterrupt();
        long double re = 0.0L, im = 0.0L;
        for (int j = 0; j < ns; ++j) {
            // Apply the empirical amplitude and phase corrections.
            double ka = k_sw[i] * a[j];
            double amp = ka / (ka + 0.083);
            double psi = ka / (40.0 + ka) - 1.05;
            double p = -(2.0 * k_medium[i] * v[j] + psi);
            double weight = factor * amp *
                std::sqrt((k_medium[i] * a[j] + 1.0) * st) * du[j];
            re += weight * std::sin(p);
            im -= weight * std::cos(p);
        }
        out[i] = to_Rcomplex({static_cast<double>(re), static_cast<double>(im)});
    }
    return out;
}

