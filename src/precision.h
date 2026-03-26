#pragma once
// ============================================================================
// INCLUDES
// ============================================================================
// [[Rcpp::depends(BH)]]
#include <cmath>
#include <complex>
#include <limits>

// Guard the optional quad-precision path so the same headers compile on
// toolchains that do not provide __float128 or libquadmath.
#ifdef __GNUC__
#include <boost/multiprecision/float128.hpp>
#include <quadmath.h>
#endif

// Keep the scalar conversions between double and __float128 in one place so
// the templated numerical code can switch precision cleanly.
#ifdef __GNUC__
inline __float128 double_to_quad(double x) {
    return static_cast<__float128>(x);
}

inline double quad_to_double(__float128 x) {
    return static_cast<double>(x);
}
#endif

// profcn returns mantissas plus base-10 exponents. These helpers reconstruct
// ordinary values without duplicating that logic in every caller.
template<typename T>
inline T pow10_typed(int exponent);

template<>
inline double pow10_typed<double>(int exponent) {
    return std::pow(10.0, static_cast<double>(exponent));
}

#ifdef __GNUC__
template<>
inline __float128 pow10_typed<__float128>(int exponent) {
    return powq(double_to_quad(10.0), double_to_quad(static_cast<double>(exponent)));
}
#endif

// Basic transcendental wrappers keep the templated code readable while still
// dispatching to quadmath when the quad backend is enabled.
template<typename T>
inline T preccos(T x) { return std::cos(x); }

#ifdef __GNUC__
template<>
inline __float128 preccos(__float128 x) { return cosq(x); }
#endif

// R uses sentinel values for missing integers. Mirror that test here so the
// compiled PSMS/TMM code can recognize missing profcn entries directly.
inline bool is_na_int(int x) {
    return x == NA_INTEGER;
}

// Floating-point missing values are represented as NaN on the C++ side.
template<typename T>
inline bool is_na_real(T x) {
    return std::isnan(x);
}

#ifdef __GNUC__
template<>
inline bool is_na_real<__float128>(__float128 x) {
    return isnanq(x);
}
#endif

// Shared conversion helper for returning complex values through Rcpp.
inline Rcomplex to_Rcomplex(const std::complex<double>& z) {
    Rcomplex out;
    out.r = z.real();
    out.i = z.imag();
    return out;
}

#ifdef __GNUC__
inline Rcomplex to_Rcomplex(const std::complex<__float128>& z) {
    Rcomplex out;
    out.r = static_cast<double>(z.real());
    out.i = static_cast<double>(z.imag());
    return out;
}
#endif
