#ifndef BESSEL_HELPERS_H
#define BESSEL_HELPERS_H

#include <Rcpp.h>
#include <complex>

inline Rcomplex to_Rcomplex(const std::complex<double>& z) {
    Rcomplex rc;
    rc.r = z.real();
    rc.i = z.imag();
    return rc;
}

// Cylindrical Bessel helper functions
std::complex<double> jc_single_impl(std::complex<double> zi, double nui);
std::complex<double> yc_single_impl(std::complex<double> zi, double nui);

// Spherical Bessel helper functions
double js_single_impl(int li, double zi);
double ys_single_impl(int li, double zi);
std::complex<double> hs_single_impl(int li, double zi);

// Complex spherical Bessel helper functions
std::complex<double> js_single_complex_impl(int li, std::complex<double> zi);
std::complex<double> ys_single_complex_impl(int li, std::complex<double> zi);
std::complex<double> hs_single_complex_impl(int li, std::complex<double> zi);

// Spherical Bessel derivative helper functions
double js_deriv_single_impl(int li, double zi, int k);
double ys_deriv_single_impl(int li, double zi, int k);
std::complex<double> hs_deriv_single_impl(int li, double zi, int k);

// Complex spherical Bessel derivative helper functions
std::complex<double> js_deriv_single_complex_impl(int li, std::complex<double> zi, int k);
std::complex<double> ys_deriv_single_complex_impl(int li, std::complex<double> zi, int k);
std::complex<double> hs_deriv_single_complex_impl(int li, std::complex<double> zi, int k);

#endif
