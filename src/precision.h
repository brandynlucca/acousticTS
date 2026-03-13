#pragma once
// ============================================================================
// INCLUDES
// ============================================================================
// [[Rcpp::depends(BH)]]
#include <cmath>
#include <complex>
#include <limits>

// Guard for quad precision
#ifdef __GNUC__
#include <boost/multiprecision/float128.hpp>
#include <quadmath.h>
#endif

// Conversion from double-to-quad and quad-to-double
// Primarily for Rcpp interfacing; also required for Armadillo
#ifdef __GNUC__
inline __float128 double_to_quad(double x) {
    return static_cast<__float128>(x);
}

inline double quad_to_double(__float128 x) {
    return static_cast<double>(x);
}
#endif

// Specializations for power-10 exponentiation
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

// Specializations for cosine
template<typename T>
inline T preccos(T x) { return std::cos(x); }

#ifdef __GNUC__
template<>
inline __float128 preccos(__float128 x) { return cosq(x); }
#endif

// Check if integer is NA
inline bool is_na_int(int x) {
    return x == NA_INTEGER;
}

// Check if numeric is NA
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

// Convert C++ complex to required R-formatted complex
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
