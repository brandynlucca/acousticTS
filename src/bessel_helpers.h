#ifndef BESSEL_HELPERS_H
#define BESSEL_HELPERS_H

#include <Rcpp.h>
#include <complex>

// Convert the internal std::complex representation used by the compiled
// solvers into the Rcomplex layout expected by Rcpp return values.
inline Rcomplex to_Rcomplex(const std::complex<double>& z) {
    Rcomplex rc;
    rc.r = z.real();
    rc.i = z.imag();
    return rc;
}

// Cylindrical Bessel helpers used by the weak-scattering and viscous-layer
// solvers. These evaluate one order/argument pair at a time.
std::complex<double> jc_single_impl(std::complex<double> zi, double nui);
std::complex<double> yc_single_impl(std::complex<double> zi, double nui);

// Real-argument spherical Bessel helpers used when the modal arguments stay on
// the real axis.
double js_single_impl(int li, double zi);
double ys_single_impl(int li, double zi);
std::complex<double> hs_single_impl(int li, double zi);

// Complex-argument spherical Bessel helpers used when attenuation or complex
// wavenumbers move the argument off the real axis.
std::complex<double> js_single_complex_impl(int li, std::complex<double> zi);
std::complex<double> ys_single_complex_impl(int li, std::complex<double> zi);
std::complex<double> hs_single_complex_impl(int li, std::complex<double> zi);

// Real-argument spherical Bessel derivatives. The derivative order k is
// explicit so callers can request higher derivatives directly.
double js_deriv_single_impl(int li, double zi, int k);
double ys_deriv_single_impl(int li, double zi, int k);
std::complex<double> hs_deriv_single_impl(int li, double zi, int k);

// Complex-argument spherical Bessel derivatives used by the elastic and
// viscous formulations.
std::complex<double> js_deriv_single_complex_impl(int li, std::complex<double> zi, int k);
std::complex<double> ys_deriv_single_complex_impl(int li, std::complex<double> zi, int k);
std::complex<double> hs_deriv_single_complex_impl(int li, std::complex<double> zi, int k);

#endif
