#pragma once
// ============================================================================
// INCLUDES
// ============================================================================
// [[Rcpp::depends(BH)]]
#include <boost/multiprecision/cpp_bin_float.hpp>
#include <cmath>
#include <complex>
#include <limits>
#include <type_traits>

// Guard the optional quad-precision path so the same headers compile on
// toolchains that do not provide __float128 or libquadmath.
#ifndef ACOUSTICTS_HAVE_QUADMATH
#if defined(__has_include)
#if __has_include(<quadmath.h>) && defined(__SIZEOF_FLOAT128__)
#define ACOUSTICTS_HAVE_QUADMATH 1
#else
#define ACOUSTICTS_HAVE_QUADMATH 0
#endif
#else
#define ACOUSTICTS_HAVE_QUADMATH 0
#endif
#endif

#if ACOUSTICTS_HAVE_QUADMATH
#include <quadmath.h>
#endif

// Use the native GCC/libquadmath scalar when it exists, and otherwise fall
// back to Boost's IEEE-quad-equivalent binary float so clang/macOS can still
// compile the precision = "quad" path.
#if ACOUSTICTS_HAVE_QUADMATH
using acousticts_quad_t = __float128;
#else
using acousticts_quad_t = boost::multiprecision::cpp_bin_float_quad;
#endif

template<typename T>
struct is_quad_precision : std::false_type {};

template<>
struct is_quad_precision<acousticts_quad_t> : std::true_type {};

template<typename T>
inline constexpr bool is_quad_precision_v = is_quad_precision<T>::value;

// Keep the scalar conversions between double and the quad-precision backend in
// one place so the templated numerical code can switch precision cleanly.
inline acousticts_quad_t double_to_quad(double x) {
    return static_cast<acousticts_quad_t>(x);
}

inline double quad_to_double(const acousticts_quad_t& x) {
    return static_cast<double>(x);
}

// profcn returns mantissas plus base-10 exponents. These helpers reconstruct
// ordinary values without duplicating that logic in every caller.
template<typename T>
inline T pow10_typed(int exponent);

template<>
inline double pow10_typed<double>(int exponent) {
    return std::pow(10.0, static_cast<double>(exponent));
}

template<>
inline acousticts_quad_t pow10_typed<acousticts_quad_t>(int exponent) {
#if ACOUSTICTS_HAVE_QUADMATH
    return powq(double_to_quad(10.0), double_to_quad(static_cast<double>(exponent)));
#else
    using boost::multiprecision::pow;
    return pow(acousticts_quad_t(10), acousticts_quad_t(exponent));
#endif
}

// Basic transcendental wrappers keep the templated code readable while still
// dispatching to quadmath when the quad backend is enabled.
template<typename T>
inline T preccos(T x) { return std::cos(x); }

template<>
inline acousticts_quad_t preccos(acousticts_quad_t x) {
#if ACOUSTICTS_HAVE_QUADMATH
    return cosq(x);
#else
    using boost::multiprecision::cos;
    return cos(x);
#endif
}

template<typename T>
inline T precsqrt(T x) { return std::sqrt(x); }

template<>
inline acousticts_quad_t precsqrt(acousticts_quad_t x) {
#if ACOUSTICTS_HAVE_QUADMATH
    return sqrtq(x);
#else
    using boost::multiprecision::sqrt;
    return sqrt(x);
#endif
}

template<typename T>
inline T precabs(T x) { return std::abs(x); }

template<>
inline acousticts_quad_t precabs(acousticts_quad_t x) {
#if ACOUSTICTS_HAVE_QUADMATH
    return fabsq(x);
#else
    using boost::multiprecision::abs;
    return abs(x);
#endif
}

template<typename T>
inline T preceps() { return std::numeric_limits<T>::epsilon(); }

template<>
inline acousticts_quad_t preceps<acousticts_quad_t>() {
#if ACOUSTICTS_HAVE_QUADMATH
    static const acousticts_quad_t one = acousticts_quad_t(1);
    static const acousticts_quad_t two = acousticts_quad_t(2);
    static const acousticts_quad_t eps = nextafterq(one, two) - one;
    return eps;
#else
    return std::numeric_limits<acousticts_quad_t>::epsilon();
#endif
}

template<typename T>
inline T complex_abs_sq(const std::complex<T>& z) {
    return z.real() * z.real() + z.imag() * z.imag();
}

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

template<>
inline bool is_na_real<acousticts_quad_t>(acousticts_quad_t x) {
    return x != x;
}

// Shared conversion helper for returning complex values through Rcpp.
inline Rcomplex to_Rcomplex(const std::complex<double>& z) {
    Rcomplex out;
    out.r = z.real();
    out.i = z.imag();
    return out;
}

inline Rcomplex to_Rcomplex(const std::complex<acousticts_quad_t>& z) {
    Rcomplex out;
    out.r = static_cast<double>(z.real());
    out.i = static_cast<double>(z.imag());
    return out;
}
