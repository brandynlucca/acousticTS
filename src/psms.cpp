// ============================================================================
// INCLUDES
// ============================================================================
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <algorithm>
#include <cmath>
#include <complex>
#include <limits>
#include <map>
#include <stdexcept>
#include <type_traits>
#include <vector>
#ifdef __GNUC__
// #include <boost/multiprecision/float128.hpp>
#include <quadmath.h>
#endif

// ============================================================================
// FORTRAN INTERFACES - DOUBLE AND QUAD PRECISION
// ============================================================================
extern "C" {
    // Double precision
    void profcn_cpp_interface(
        double* c, int* m, int* lnum, int* ioprad, double* x1, int* iopang, int* iopnorm, 
        int* narg, double* arg, double* r1c, int* ir1e, double* r1dc, int* ir1de, double* r2c, 
        int* ir2e, double* r2dc, int* ir2de, int* naccr, double* s1c, int* is1e, double* s1dc, 
        int* is1de, int* naccs
    );
    void profcn_cpp_interface_batch(
        double* c, int* m_start, int* m_count, int* lnum, int* ioprad, double* x1,
        int* iopang, int* iopnorm, int* narg, double* arg, double* r1c, int* ir1e,
        double* r1dc, int* ir1de, double* r2c, int* ir2e, double* r2dc, int* ir2de,
        int* naccr, double* s1c, int* is1e, double* s1dc, int* is1de, int* naccs
    );

    // Quad precision
#ifdef __GNUC__
    // Quad precision interface
    void profcn_cpp_interface_quad(
        __float128* c, int* m, int* lnum, int* ioprad, __float128* x1, int* iopang, int* iopnorm, 
        int* narg, __float128* arg, __float128* r1c, int* ir1e, __float128* r1dc, int* ir1de, 
        __float128* r2c, int* ir2e, __float128* r2dc, int* ir2de, int* naccr,
        __float128* s1c, int* is1e, __float128* s1dc, int* is1de, int* naccs
    );
    void profcn_cpp_interface_batch_quad(
        __float128* c, int* m_start, int* m_count, int* lnum, int* ioprad, __float128* x1,
        int* iopang, int* iopnorm, int* narg, __float128* arg, __float128* r1c, int* ir1e,
        __float128* r1dc, int* ir1de, __float128* r2c, int* ir2e, __float128* r2dc, int* ir2de,
        int* naccr, __float128* s1c, int* is1e, __float128* s1dc, int* is1de, int* naccs
    );

    // Quad precision math functions
    __float128 powq(__float128, __float128);
#endif
}

// ============================================================================
// QUAD PRECISION HELPERS
// ============================================================================
// The spheroidal backend often stores values as mantissas plus base-10
// exponents to avoid overflow and underflow. These helpers provide the
// precision-aware scalar operations needed to reconstruct ordinary values in
// either double or quadruple precision.
template<typename T>
inline T pow10_typed(int exponent);

template<>
inline double pow10_typed<double>(int exponent) {
    if (exponent == 0) return 1.0;
    return std::pow(10.0, static_cast<double>(exponent));
}

#ifdef __GNUC__
// Keep the explicit double <-> quad conversions in one place so the rest of
// the PSMS implementation can remain templated and readable.
inline __float128 double_to_quad(double x) {
    return static_cast<__float128>(x);
}

inline double quad_to_double(__float128 x) {
    return static_cast<double>(x);
}

template<>
inline __float128 pow10_typed<__float128>(int exponent) {
    if (exponent == 0) return double_to_quad(1.0);
    return powq(double_to_quad(10.0), double_to_quad(static_cast<double>(exponent)));
}
#endif

template<typename T>
inline T preccos(T x) { return std::cos(x); }
#ifdef __GNUC__
template<>
inline __float128 preccos(__float128 x) { return cosq(x); }
#endif

template<typename T>
inline T precsqrt(T x) { return std::sqrt(x); }
#ifdef __GNUC__
template<>
inline __float128 precsqrt(__float128 x) { return sqrtq(x); }
#endif

template<typename T>
inline T precabs(T x) { return std::abs(x); }
#ifdef __GNUC__
template<>
inline __float128 precabs(__float128 x) { return fabsq(x); }
#endif

template<typename T>
inline T preceps() { return std::numeric_limits<T>::epsilon(); }
#ifdef __GNUC__
template<>
inline __float128 preceps<__float128>() {
    static const __float128 one = static_cast<__float128>(1.0);
    static const __float128 two = static_cast<__float128>(2.0);
    static const __float128 eps = nextafterq(one, two) - one;
    return eps;
}
#endif

template<typename T>
inline T complex_abs_sq(const std::complex<T>& z) {
    return z.real() * z.real() + z.imag() * z.imag();
}

template<typename T>
struct ProfcnResult;

template<typename T>
struct ProfcnBatchResult;

inline bool is_na_int(int x);

template<typename T>
inline bool is_na_real(T x);

template<typename T>
inline void scale_profcn_component(std::vector<T>& values, std::vector<int>& exponents);

// Internal PSMS support layer: special-function wrappers, profcn batching,
// exponent rescaling helpers, and precision-agnostic scaffolding.
#include "psms_support.h"

#include "psms_smn.h"

// C++-R INTERFACE FUNCTION [Smn]
// -----------------------------------------------------------
// [[Rcpp::export]]
Rcpp::List Smn_cpp(
    Rcpp::IntegerVector m,
    Rcpp::IntegerVector n, double c, Rcpp::NumericVector arg,
    bool normalize = false, std::string precision = "double"
) {
    std::vector<int> m_vec = Rcpp::as<std::vector<int>>(m);
    std::vector<int> n_vec = Rcpp::as<std::vector<int>>(n);
    int m_len = m_vec.size();
    int n_len = n_vec.size();
    int eta_len = arg.size();

    if (m_len == 0 || n_len == 0)
        Rcpp::stop("'m' and 'n' must have at least one element.");
    if (eta_len == 0)
        Rcpp::stop("'eta' must have at least one element.");

    if (precision == "quad") {
#ifndef __GNUC__
        Rcpp::stop("Quad precision requires GCC compilation.");
#else
        std::vector<__float128> eta_quad(eta_len);
        for (int i = 0; i < eta_len; ++i)
            eta_quad[i] = static_cast<__float128>(arg[i]);
        auto res = Smn_matrix<__float128>(
            m_vec, n_vec, static_cast<__float128>(c), eta_quad, normalize
        );
        return Smn_cpp_to_rcpp_list(res, m_len, n_len, eta_len);
#endif
    } else if (precision == "double") {
        std::vector<double> eta_vec = Rcpp::as<std::vector<double>>(arg);
        auto res = Smn_matrix<double>(m_vec, n_vec, c, eta_vec, normalize);
        return Smn_cpp_to_rcpp_list(res, m_len, n_len, eta_len);
    } else {
        Rcpp::stop("'precision' must be 'double' or 'quad'.");
    }
}

// ============================================================================
#include "psms_rmn.h"

// C++-R INTERFACE FUNCTION [Rmn]
// -----------------------------------------------------------
// [[Rcpp::export]]
Rcpp::List Rmn_cpp(
    Rcpp::IntegerVector m,
    Rcpp::IntegerVector n,
    double c,
    double x1,
    int kind = 1,
    std::string precision = "double"
) {
    std::vector<int> m_vec = Rcpp::as<std::vector<int>>(m);
    std::vector<int> n_vec = Rcpp::as<std::vector<int>>(n);
    int m_len = m_vec.size();
    int n_len = n_vec.size();

    if (precision == "quad") {
#ifndef __GNUC__
        Rcpp::stop("Quad precision requires GCC compilation.");
#else
        __float128 c_q = static_cast<__float128>(c);
        __float128 x1_q = static_cast<__float128>(x1);
        auto res = Rmn_matrix<__float128>(m_vec, n_vec, c_q, x1_q, kind);

        if (kind == 3 || kind == 4) {
            if (m_len == 1 || n_len == 1) {
                Rcpp::ComplexVector values(n_len), derivatives(n_len);
                for (int i = 0; i < n_len; ++i) {
                    values[i] = to_Rcomplex(res.value[0][i]);
                    derivatives[i] = to_Rcomplex(res.derivative[0][i]);
                }
                return Rcpp::List::create(
                    Rcpp::Named("value") = values,
                    Rcpp::Named("derivative") = derivatives
                );
            } else {
                Rcpp::ComplexMatrix values(m_len, n_len), derivatives(m_len, n_len);
                for (int i = 0; i < m_len; ++i)
                    for (int j = 0; j < n_len; ++j) {
                        values(i, j) = to_Rcomplex(res.value[i][j]);
                        derivatives(i, j) = to_Rcomplex(res.derivative[i][j]);
                    }
                return Rcpp::List::create(
                    Rcpp::Named("value") = values,
                    Rcpp::Named("derivative") = derivatives
                );
            }
        } else {
            if (m_len == 1 || n_len == 1) {
                Rcpp::NumericVector values(n_len), derivatives(n_len);
                for (int i = 0; i < n_len; ++i) {
                    values[i] = static_cast<double>(res.value[0][i].real());
                    derivatives[i] = static_cast<double>(res.derivative[0][i].real());
                }
                return Rcpp::List::create(
                    Rcpp::Named("value") = values,
                    Rcpp::Named("derivative") = derivatives
                );
            } else {
                Rcpp::NumericMatrix values(m_len, n_len), derivatives(m_len, n_len);
                for (int i = 0; i < m_len; ++i)
                    for (int j = 0; j < n_len; ++j) {
                        values(i, j) = static_cast<double>(res.value[i][j].real());
                        derivatives(i, j) = static_cast<double>(res.derivative[i][j].real());
                    }
                return Rcpp::List::create(
                    Rcpp::Named("value") = values,
                    Rcpp::Named("derivative") = derivatives
                );
            }
        }
#endif
    } else if (precision == "double") {
        auto res = Rmn_matrix<double>(m_vec, n_vec, c, x1, kind);
        if (kind == 3 || kind == 4) {
            if (m_len == 1 || n_len == 1) {
                Rcpp::ComplexVector values(n_len), derivatives(n_len);
                for (int i = 0; i < n_len; ++i) {
                    values[i] = to_Rcomplex(res.value[0][i]);
                    derivatives[i] = to_Rcomplex(res.derivative[0][i]);
                }
                return Rcpp::List::create(
                    Rcpp::Named("value") = values,
                    Rcpp::Named("derivative") = derivatives
                );
            } else {
                Rcpp::ComplexMatrix values(m_len, n_len), derivatives(m_len, n_len);
                for (int i = 0; i < m_len; ++i)
                    for (int j = 0; j < n_len; ++j) {
                        values(i, j) = to_Rcomplex(res.value[i][j]);
                        derivatives(i, j) = to_Rcomplex(res.derivative[i][j]);
                    }
                return Rcpp::List::create(
                    Rcpp::Named("value") = values,
                    Rcpp::Named("derivative") = derivatives
                );
            }
        } else {
            if (m_len == 1 || n_len == 1) {
                Rcpp::NumericVector values(n_len), derivatives(n_len);
                for (int i = 0; i < n_len; ++i) {
                    values[i] = res.value[0][i].real();
                    derivatives[i] = res.derivative[0][i].real();
                }
                return Rcpp::List::create(
                    Rcpp::Named("value") = values,
                    Rcpp::Named("derivative") = derivatives
                );
            } else {
                Rcpp::NumericMatrix values(m_len, n_len), derivatives(m_len, n_len);
                for (int i = 0; i < m_len; ++i)
                    for (int j = 0; j < n_len; ++j) {
                        values(i, j) = res.value[i][j].real();
                        derivatives(i, j) = res.derivative[i][j].real();
                    }
                return Rcpp::List::create(
                    Rcpp::Named("value") = values,
                    Rcpp::Named("derivative") = derivatives
                );
            }
        }
    } else {
        Rcpp::stop("'precision' must be 'double' or 'quad'.");
    }
    return R_NilValue;
}

// ============================================================================
#include "psms_backscatter.h"

// C++-R INTERFACE FUNCTION FOR PROLATE SPHEROID MODAL SERIES SOLUTION
// ----------------------------------------------------------------------------
// Based on the algorithm proposed by Furusawa (1988).
// ============================================================================
// [[Rcpp::export]]
Rcpp::ComplexVector prolate_spheroid_fbs(
    Rcpp::DataFrame acoustics,
    Rcpp::DataFrame body,
    Rcpp::DataFrame medium,
    Rcpp::List integration_pts,
    std::string precision = "double",
    std::string Amn_method = "Amn_fluid",
    bool adaptive = false,
    bool vectorized = false
) {
    // Extract acoustic parameters
    std::vector<double> chi_sw = Rcpp::as<std::vector<double>>(acoustics["chi_sw"]);
    std::vector<double> chi_body = Rcpp::as<std::vector<double>>(acoustics["chi_body"]);
    std::vector<int> m_max = Rcpp::as<std::vector<int>>(acoustics["m_max"]);
    std::vector<int> n_max = Rcpp::as<std::vector<int>>(acoustics["n_max"]);

    // Extract body parameters
    double xi = body["xi"];
    double theta_body = body["theta_body"];
    double phi_body = body["phi_body"];
    double phi_scatter = body["phi_scatter"];
    double density_body = body["density"];

    // Extract medium parameters
    double density_sw = medium["density"];
    
    // Extract the quadrature rule used for overlap integrals in the fluid/gas
    // pathways.
    Rcpp::NumericVector nodes = integration_pts["nodes"];
    Rcpp::NumericVector weights = integration_pts["weights"];

    int n_freq = acoustics.nrows();
    std::vector<Rcomplex> f_bs_out(n_freq);

    if (precision == "quad") {
#ifdef __GNUC__
        std::vector<__float128> nodes_q(nodes.size()), weights_q(weights.size());
        for (int i = 0; i < nodes.size(); ++i) nodes_q[i] = static_cast<__float128>(nodes[i]);
        for (int i = 0; i < weights.size(); ++i) weights_q[i] = static_cast<__float128>(weights[i]);
        for (int f = 0; f < n_freq; ++f) {
            // Each frequency is independent, so the wrapper simply dispatches
            // one retained PSMS solve per frequency and stores the resulting
            // complex backscatter amplitude.
            auto fbs = psms_fbs<__float128>(
                m_max[f], n_max[f],
                static_cast<__float128>(chi_sw[f]),
                static_cast<__float128>(chi_body[f]),
                static_cast<__float128>(xi),
                static_cast<__float128>(theta_body),
                static_cast<__float128>(phi_body),
                static_cast<__float128>(phi_scatter),
                static_cast<__float128>(density_body),
                static_cast<__float128>(density_sw),
                nodes_q, weights_q,
                Amn_method, adaptive, vectorized
            );
            f_bs_out[f] = to_Rcomplex(fbs);
        }
#else
        Rcpp::stop("Quad precision requires GCC and libquadmath.");
#endif
    } else {
        std::vector<double> nodes_vec = Rcpp::as<std::vector<double>>(nodes);
        std::vector<double> weights_vec = Rcpp::as<std::vector<double>>(weights);
        for (int f = 0; f < n_freq; ++f) {
            auto fbs = psms_fbs<double>(
                m_max[f], n_max[f],
                chi_sw[f], chi_body[f], xi,
                theta_body, phi_body, phi_scatter,
                density_body, density_sw,
                nodes_vec, weights_vec,
                Amn_method, adaptive, vectorized
            );
            f_bs_out[f] = to_Rcomplex(fbs);
        }
    }
    return Rcpp::wrap(f_bs_out);
}
