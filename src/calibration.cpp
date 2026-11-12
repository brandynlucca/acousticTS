#include <Rcpp.h>
#include <cmath>
#include <complex>
#include <vector>
#include "bessel_helpers.h"

using namespace Rcpp;

// [[Rcpp::export]]
ComplexVector calibration_modal_sum_cpp(NumericVector ka_sw,
                                         NumericVector ka_l,
                                         NumericVector ka_t,
                                         NumericVector theta,
                                         NumericVector density_ratio,
                                         bool adaptive) {
    int nf = ka_sw.size();
    if (ka_l.size() == 0 || ka_t.size() == 0 || theta.size() == 0 ||
        density_ratio.size() == 0) stop("Calibration inputs must be non-empty.");
    ComplexVector out(nf);
    for (int i = 0; i < nf; ++i) {
        if (i % 32 == 0) Rcpp::checkUserInterrupt();
        double x = ka_sw[i];
        double xl = ka_l[i % ka_l.size()];
        double xt = ka_t[i % ka_t.size()];
        double ct = std::cos(theta[i % theta.size()]);
        double g = density_ratio[i % density_ratio.size()];
        // Retain the existing initial modal cutoff and adaptive tail criterion.
        int ml = static_cast<int>(std::nearbyint(x)) + 10;
        int sequence_limit = ml + 1;
        // Regular waves use downward recurrence, irregular waves use upward.
        std::vector<double> j = js_sequence_miller_impl(sequence_limit, x);
        std::vector<double> jl = js_sequence_miller_impl(sequence_limit, xl);
        std::vector<double> jt = js_sequence_miller_impl(sequence_limit, xt);
        std::vector<double> y = ys_sequence_upward_impl(sequence_limit, x);
        double p_previous = 1.0, p = 1.0;
        long double re = 0.0L, im = 0.0L;
        for (int m = 0; ; ++m) {
            if (m > sequence_limit - 1) {
                // Extend radial sequences only when the adaptive sum needs them.
                sequence_limit = m + 16;
                j = js_sequence_miller_impl(sequence_limit, x);
                jl = js_sequence_miller_impl(sequence_limit, xl);
                jt = js_sequence_miller_impl(sequence_limit, xt);
                y = ys_sequence_upward_impl(sequence_limit, x);
            }
            // Advance the Legendre polynomial for the requested scattering angle.
            if (m == 1) p = ct;
            else if (m > 1) {
                double next = ((2.0 * m - 1.0) * ct * p -
                    (m - 1.0) * p_previous) / m;
                p_previous = p;
                p = next;
            }
            // Form logarithmic radial derivatives for the elastic boundary terms.
            double sw = x * j[m+1] / j[m] - m;
            double longitudinal = xl * jl[m+1] / jl[m] - m;
            double transversal = xt * jt[m+1] / jt[m] - m;
            double beta = x * y[m+1] / y[m] - m;
            double diff = -j[m] / y[m];
            double mm = m * (m + 1.0), t2 = xt * xt / 2.0;
            double la = longitudinal + 1.0;
            double td = mm - 1.0 - t2 + transversal;
            double numerator = longitudinal / la - mm / td;
            double denominator = (mm - t2 + 2.0 * longitudinal) / la -
                mm * (transversal + 1.0) / td;
            // Sum the MacLennan modal phase shifts with extended accumulation.
            double phi = t2 * numerator / denominator / g;
            double eta = diff * (phi + sw) / (phi + beta);
            double cosine = 1.0 / std::sqrt(1.0 + eta * eta);
            double sine = eta * cosine;
            double weight = (2.0 * m + 1.0) * p;
            std::complex<double> term = weight * sine *
                std::complex<double>(-sine, cosine);
            re += term.real();
            im += term.imag();
            if (m >= ml && (!adaptive || !(std::abs(term) > 1e-10) || m >= 500))
                break;
        }
        out[i] = to_Rcomplex({static_cast<double>(re), static_cast<double>(im)});
    }
    return out;
}

