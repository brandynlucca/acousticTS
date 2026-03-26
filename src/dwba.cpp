#include <Rcpp.h>
#include <R_ext/Applic.h>
#include <R_ext/Utils.h>
#include <cmath>
#include <complex>
#include <limits>
#include <vector>
#include "bessel_helpers.h"

using namespace Rcpp;

namespace {

// Carry the geometry-dependent quantities needed by the QUADPACK callback.
// Precomputing these terms once per frequency avoids rebuilding the same
// segment metadata every time the integrator samples the interval.
struct DwbaIntegrandData {
    std::vector<double> a0;
    std::vector<double> da;
    std::vector<double> pref;
    std::vector<double> phase0;
    std::vector<double> dphase;
    std::vector<double> cos_beta;
    bool imaginary = false;
    int segment = -1;
    double k_over_h = 0.0;
};

// Evaluate the cylindrical-Bessel factor that appears in the DWBA line
// integral. When cos(beta) is close to zero, use the removable limit instead
// of the raw expression to avoid catastrophic cancellation.
inline std::complex<double> dwba_bessel_over_cosbeta(double k_over_h,
                                                     double a_s,
                                                     double cos_beta) {
    if (std::abs(cos_beta) <= std::sqrt(std::numeric_limits<double>::epsilon())) {
        // Removable limit of J1(z) / cos(beta) as cos(beta) -> 0, with
        // z = 2 * (k / h) * a(s) * cos(beta).
        return std::complex<double>(k_over_h * a_s, 0.0);
    }

    double z_s = 2.0 * k_over_h * a_s * cos_beta;
    return jc_single_impl(std::complex<double>(z_s, 0.0), 1.0) / cos_beta;
}

inline double clamp_unit(double x) {
    if (x > 1.0) return 1.0;
    if (x < -1.0) return -1.0;
    return x;
}

// QUADPACK callback. The input vector x contains the current integration
// abscissae; on return it is overwritten with either the real or imaginary
// part of the complex DWBA integrand.
void dwba_integrand_eval(double *x, int n, void *ex) {
    auto *data = static_cast<DwbaIntegrandData*>(ex);
    int n_seg = static_cast<int>(data->a0.size());
    int j_start = (data->segment >= 0) ? data->segment : 0;
    int j_end = (data->segment >= 0) ? data->segment + 1 : n_seg;

    for (int i = 0; i < n; ++i) {
        double s = x[i];
        std::complex<double> total(0.0, 0.0);

        for (int j = j_start; j < j_end; ++j) {
            double a_s = data->a0[j] + s * data->da[j];
            double phase_s = data->phase0[j] + s * data->dphase[j];
            std::complex<double> bessel = dwba_bessel_over_cosbeta(
                data->k_over_h, a_s, data->cos_beta[j]
            );

            total += data->pref[j] * a_s *
                std::exp(std::complex<double>(0.0, 2.0 * phase_s)) *
                bessel;
        }

        x[i] = data->imaginary ? total.imag() : total.real();
    }
}

// Thin wrapper around R's QUADPACK entry point so the same adaptive integral
// can be reused for the real and imaginary parts of the DWBA amplitude.
double dwba_quadpack_integral(DwbaIntegrandData &data,
                              double rel_tol,
                              double abs_tol,
                              int subdivisions) {
    double lower = 0.0;
    double upper = 1.0;
    double result = NA_REAL;
    double abserr = NA_REAL;
    int neval = 0;
    int ier = 0;
    int limit = subdivisions;
    int lenw = 4 * limit;
    int last = 0;
    std::vector<int> iwork(limit);
    std::vector<double> work(lenw);

    // QUADPACK adaptive integrator that is used by R for `integrate()`
    Rdqags(
        dwba_integrand_eval,
        &data,
        &lower,
        &upper,
        &abs_tol,
        &rel_tol,
        &result,
        &abserr,
        &neval,
        &ier,
        &limit,
        &lenw,
        &last,
        iwork.data(),
        work.data()
    );

    if (ier != 0) {
        Rcpp::stop(
            "DWBA quadrature failed with QUADPACK error code %d.",
            ier
        );
    }

    return result;
}

}

// [[Rcpp::export]]
ComplexVector dwba_fbs_cpp(NumericMatrix rpos,
                           NumericVector k_sw,
                           double theta,
                           double h,
                           double R,
                           int subdivisions = 100,
                           double rel_tol = 0.0001220703125,
                           double abs_tol = 0.0001220703125) {
    if (rpos.nrow() != 4) {
        stop("'rpos' must be a 4 x n matrix with x, y, z, and radius rows.");
    }
    if (rpos.ncol() < 2) {
        stop("'rpos' must contain at least two points.");
    }
    if (subdivisions < 1) {
        stop("'subdivisions' must be a positive integer.");
    }
    if (h == 0.0) {
        stop("'h' must be non-zero.");
    }

    int n_freq = k_sw.size();
    int n_seg = rpos.ncol() - 1;
    ComplexVector out(n_freq);

    double cos_theta = std::cos(theta);
    double sin_theta = std::sin(theta);

    // Convert the sampled centerline/radius representation into segment-level
    // differences so the line integral can be evaluated over a unit interval.
    std::vector<double> dx(n_seg), dy(n_seg), dz(n_seg), da(n_seg), seglen(n_seg);
    std::vector<double> x0(n_seg), y0(n_seg), z0(n_seg), a0(n_seg);

    for (int j = 0; j < n_seg; ++j) {
        x0[j] = rpos(0, j);
        y0[j] = rpos(1, j);
        z0[j] = rpos(2, j);
        a0[j] = rpos(3, j);

        dx[j] = rpos(0, j + 1) - rpos(0, j);
        dy[j] = rpos(1, j + 1) - rpos(1, j);
        dz[j] = rpos(2, j + 1) - rpos(2, j);
        da[j] = rpos(3, j + 1) - rpos(3, j);

        seglen[j] = std::sqrt(dx[j] * dx[j] + dy[j] * dy[j] + dz[j] * dz[j]);
    }

    for (int i = 0; i < n_freq; ++i) {
        if (i % 32 == 0) R_CheckUserInterrupt();

        // Build the per-frequency phase and orientation terms that stay fixed
        // while QUADPACK samples the segment parameter.
        double k = k_sw[i];
        double kx = k * cos_theta;
        double ky = 0.0;
        double kz = k * sin_theta;

        DwbaIntegrandData data;
        data.a0 = a0;
        data.da = da;
        data.pref.resize(n_seg);
        data.phase0.resize(n_seg);
        data.dphase.resize(n_seg);
        data.cos_beta.resize(n_seg);
        data.k_over_h = k / h;

        for (int j = 0; j < n_seg; ++j) {
            double dot = dx[j] * kx + dy[j] * ky + dz[j] * kz;
            double denom = k * seglen[j];
            double cos_alpha = (denom == 0.0) ? 0.0 : clamp_unit(dot / denom);
            double cos_beta = std::sqrt(std::max(0.0, 1.0 - cos_alpha * cos_alpha));

            data.cos_beta[j] = cos_beta;
            data.pref[j] = (k / 4.0) * R * seglen[j];
            data.phase0[j] = (x0[j] * kx + y0[j] * ky + z0[j] * kz) / h;
            data.dphase[j] = (dx[j] * kx + dy[j] * ky + dz[j] * kz) / h;
        }

        // Integrate the real and imaginary parts separately because QUADPACK
        // works with real-valued callbacks.
        data.imaginary = false;
        double real_part = dwba_quadpack_integral(
            data, rel_tol, abs_tol, subdivisions
        );

        data.imaginary = true;
        double imag_part = dwba_quadpack_integral(
            data, rel_tol, abs_tol, subdivisions
        );

        out[i] = to_Rcomplex(std::complex<double>(real_part, imag_part));
    }

    return out;
}

// [[Rcpp::export]]
ComplexMatrix dwba_segment_integrals_cpp(NumericMatrix rpos,
                                         NumericVector k_sw,
                                         double theta,
                                         double h,
                                         double R,
                                         int subdivisions = 100,
                                         double rel_tol = 0.0001220703125,
                                         double abs_tol = 0.0001220703125) {
    if (rpos.nrow() != 4) {
        stop("'rpos' must be a 4 x n matrix with x, y, z, and radius rows.");
    }
    if (rpos.ncol() < 2) {
        stop("'rpos' must contain at least two points.");
    }
    if (subdivisions < 1) {
        stop("'subdivisions' must be a positive integer.");
    }
    if (h == 0.0) {
        stop("'h' must be non-zero.");
    }

    int n_freq = k_sw.size();
    int n_seg = rpos.ncol() - 1;
    ComplexMatrix out(n_freq, n_seg);

    double cos_theta = std::cos(theta);
    double sin_theta = std::sin(theta);

    std::vector<double> dx(n_seg), dy(n_seg), dz(n_seg), da(n_seg), seglen(n_seg);
    std::vector<double> x0(n_seg), y0(n_seg), z0(n_seg), a0(n_seg);

    for (int j = 0; j < n_seg; ++j) {
        x0[j] = rpos(0, j);
        y0[j] = rpos(1, j);
        z0[j] = rpos(2, j);
        a0[j] = rpos(3, j);

        dx[j] = rpos(0, j + 1) - rpos(0, j);
        dy[j] = rpos(1, j + 1) - rpos(1, j);
        dz[j] = rpos(2, j + 1) - rpos(2, j);
        da[j] = rpos(3, j + 1) - rpos(3, j);

        seglen[j] = std::sqrt(dx[j] * dx[j] + dy[j] * dy[j] + dz[j] * dz[j]);
    }

    for (int i = 0; i < n_freq; ++i) {
        if (i % 16 == 0) R_CheckUserInterrupt();

        // Rebuild the shared frequency-dependent prefactors once, then let the
        // callback isolate one segment at a time.
        double k = k_sw[i];
        double kx = k * cos_theta;
        double ky = 0.0;
        double kz = k * sin_theta;

        DwbaIntegrandData data;
        data.a0 = a0;
        data.da = da;
        data.pref.resize(n_seg);
        data.phase0.resize(n_seg);
        data.dphase.resize(n_seg);
        data.cos_beta.resize(n_seg);
        data.k_over_h = k / h;

        for (int j = 0; j < n_seg; ++j) {
            double dot = dx[j] * kx + dy[j] * ky + dz[j] * kz;
            double denom = k * seglen[j];
            double cos_alpha = (denom == 0.0) ? 0.0 : clamp_unit(dot / denom);
            double cos_beta = std::sqrt(std::max(0.0, 1.0 - cos_alpha * cos_alpha));

            data.cos_beta[j] = cos_beta;
            data.pref[j] = (k / 4.0) * R * seglen[j];
            data.phase0[j] = (x0[j] * kx + y0[j] * ky + z0[j] * kz) / h;
            data.dphase[j] = (dx[j] * kx + dy[j] * ky + dz[j] * kz) / h;
        }

        for (int j = 0; j < n_seg; ++j) {
            // Restrict the callback to a single segment so the returned matrix
            // exposes the contribution of each segment to the total amplitude.
            data.segment = j;
            data.imaginary = false;
            double real_part = dwba_quadpack_integral(
                data, rel_tol, abs_tol, subdivisions
            );

            data.imaginary = true;
            double imag_part = dwba_quadpack_integral(
                data, rel_tol, abs_tol, subdivisions
            );

            out(i, j) = to_Rcomplex(std::complex<double>(real_part, imag_part));
        }
    }

    return out;
}
