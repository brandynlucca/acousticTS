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
    // Normalize the R inputs once, then dispatch to the precision-specific
    // angular evaluator and convert the result back to the public R shape.
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
    // Mirror the public R interface while keeping the heavy radial evaluation
    // in the templated C++ backend.
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
    // Extract the frequency-dependent modal and acoustic state.
    std::vector<double> chi_sw = Rcpp::as<std::vector<double>>(acoustics["chi_sw"]);
    std::vector<double> chi_body = Rcpp::as<std::vector<double>>(acoustics["chi_body"]);
    std::vector<int> m_max = Rcpp::as<std::vector<int>>(acoustics["m_max"]);
    std::vector<int> n_max = Rcpp::as<std::vector<int>>(acoustics["n_max"]);

    // Extract the target geometry and scattering angles.
    double xi = body["xi"];
    double theta_body = body["theta_body"];
    double theta_scatter = body["theta_scatter"];
    double phi_body = body["phi_body"];
    double phi_scatter = body["phi_scatter"];
    double density_body = body["density"];

    // Extract the surrounding-medium density used in the boundary terms.
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
                static_cast<__float128>(theta_scatter),
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
                theta_body, theta_scatter, phi_body, phi_scatter,
                density_body, density_sw,
                nodes_vec, weights_vec,
                Amn_method, adaptive, vectorized
            );
            f_bs_out[f] = to_Rcomplex(fbs);
        }
    }
    return Rcpp::wrap(f_bs_out);
}

template<typename T>
Rcpp::List prolate_tmatrix_blocks_to_rcpp(
    int m_max,
    int n_max,
    const std::vector<std::vector<std::complex<T>>>& t_blocks
) {
    // Store each retained m-block as a square complex matrix together with the
    // corresponding degree sequence so the R-side post-processing helpers can
    // reconstruct the modal geometry later on.
    Rcpp::List blocks(m_max + 1);
    for (int m = 0; m <= m_max; ++m) {
        int size = n_max - m + 1;
        Rcpp::ComplexMatrix block_mat(size, size);
        for (int row = 0; row < size; ++row) {
            for (int col = 0; col < size; ++col) {
                block_mat(row, col) = to_Rcomplex(t_blocks[m][row * size + col]);
            }
        }
        blocks[m] = Rcpp::List::create(
            Rcpp::Named("m") = m,
            Rcpp::Named("n_seq") = Rcpp::seq(m, n_max),
            Rcpp::Named("T") = block_mat
        );
    }
    blocks.attr("coordinate_system") = "spheroidal";
    return blocks;
}

// Build and optionally retain the spheroidal T-matrix blocks used by the
// prolate branch, then evaluate the requested scattering geometry from those
// retained blocks.
// [[Rcpp::export]]
Rcpp::List prolate_spheroid_tmatrix_cpp(
    Rcpp::DataFrame acoustics,
    Rcpp::DataFrame body,
    Rcpp::DataFrame medium,
    Rcpp::List integration_pts,
    std::string precision = "double",
    std::string Amn_method = "Amn_fluid"
) {
    std::vector<double> chi_sw = Rcpp::as<std::vector<double>>(acoustics["chi_sw"]);
    std::vector<double> chi_body = Rcpp::as<std::vector<double>>(acoustics["chi_body"]);
    std::vector<int> m_max = Rcpp::as<std::vector<int>>(acoustics["m_max"]);
    std::vector<int> n_max = Rcpp::as<std::vector<int>>(acoustics["n_max"]);

    double xi = body["xi"];
    double theta_body = body["theta_body"];
    double theta_scatter = body["theta_scatter"];
    double phi_body = body["phi_body"];
    double phi_scatter = body["phi_scatter"];
    double density_body = body["density"];
    double density_sw = medium["density"];

    Rcpp::NumericVector nodes = integration_pts["nodes"];
    Rcpp::NumericVector weights = integration_pts["weights"];

    int n_freq = acoustics.nrows();
    Rcpp::ComplexVector f_scat(n_freq);
    Rcpp::List t_store(n_freq);

    if (precision == "quad") {
#ifdef __GNUC__
        std::vector<__float128> nodes_q(nodes.size()), weights_q(weights.size());
        for (int i = 0; i < nodes.size(); ++i) nodes_q[i] = static_cast<__float128>(nodes[i]);
        for (int i = 0; i < weights.size(); ++i) weights_q[i] = static_cast<__float128>(weights[i]);

        for (int f = 0; f < n_freq; ++f) {
            auto t_blocks = psms_tmatrix_blocks<__float128>(
                m_max[f], n_max[f],
                static_cast<__float128>(chi_sw[f]),
                static_cast<__float128>(chi_body[f]),
                static_cast<__float128>(xi),
                static_cast<__float128>(density_body),
                static_cast<__float128>(density_sw),
                nodes_q, weights_q,
                Amn_method
            );
            auto smn_inc = compute_smn_matrix<__float128>(
                m_max[f], n_max[f], static_cast<__float128>(chi_sw[f]),
                preccos(static_cast<__float128>(theta_body)), true
            );
            auto smn_scat = compute_smn_matrix<__float128>(
                m_max[f], n_max[f], static_cast<__float128>(chi_sw[f]),
                preccos(static_cast<__float128>(theta_scatter)), true
            );
            auto azimuth = compute_azimuth<__float128>(
                m_max[f],
                static_cast<__float128>(phi_body),
                static_cast<__float128>(phi_scatter)
            );

            auto f_val = compute_fbs_from_tmatrix_blocks<__float128>(
                m_max[f], n_max[f], azimuth, smn_inc, smn_scat, t_blocks
            );
            f_scat[f] = to_Rcomplex(f_val);
            t_store[f] = prolate_tmatrix_blocks_to_rcpp<__float128>(
                m_max[f], n_max[f], t_blocks
            );
        }
#else
        Rcpp::stop("Quad precision requires GCC and libquadmath.");
#endif
    } else {
        std::vector<double> nodes_vec = Rcpp::as<std::vector<double>>(nodes);
        std::vector<double> weights_vec = Rcpp::as<std::vector<double>>(weights);

        for (int f = 0; f < n_freq; ++f) {
            auto t_blocks = psms_tmatrix_blocks<double>(
                m_max[f], n_max[f], chi_sw[f], chi_body[f], xi,
                density_body, density_sw, nodes_vec, weights_vec, Amn_method
            );
            auto smn_inc = compute_smn_matrix<double>(
                m_max[f], n_max[f], chi_sw[f], preccos(theta_body), true
            );
            auto smn_scat = compute_smn_matrix<double>(
                m_max[f], n_max[f], chi_sw[f], preccos(theta_scatter), true
            );
            auto azimuth = compute_azimuth<double>(m_max[f], phi_body, phi_scatter);

            auto f_val = compute_fbs_from_tmatrix_blocks<double>(
                m_max[f], n_max[f], azimuth, smn_inc, smn_scat, t_blocks
            );
            f_scat[f] = to_Rcomplex(f_val);
            t_store[f] = prolate_tmatrix_blocks_to_rcpp<double>(
                m_max[f], n_max[f], t_blocks
            );
        }
    }

    return Rcpp::List::create(
        Rcpp::Named("f_scat") = f_scat,
        Rcpp::Named("t_matrix") = t_store
    );
}

namespace {

template<typename T>
std::vector<std::vector<std::complex<T>>> prolate_tmatrix_blocks_from_rcpp(
    const Rcpp::List& blocks
) {
    // Recover the retained T-matrix blocks from the R-side storage layout used
    // by prolate_tmatrix_blocks_to_rcpp().
    int m_max = blocks.size() - 1;
    std::vector<std::vector<std::complex<T>>> t_blocks(m_max + 1);

    for (int m = 0; m <= m_max; ++m) {
        Rcpp::List block = blocks[m];
        Rcpp::ComplexMatrix block_mat = block["T"];
        int n_row = block_mat.nrow();
        int n_col = block_mat.ncol();
        t_blocks[m].resize(n_row * n_col, std::complex<T>(0, 0));

        for (int row = 0; row < n_row; ++row) {
            for (int col = 0; col < n_col; ++col) {
                Rcomplex z = block_mat(row, col);
                t_blocks[m][row * n_col + col] = std::complex<T>(
                    static_cast<T>(z.r),
                    static_cast<T>(z.i)
                );
            }
        }
    }

    return t_blocks;
}

template<typename T>
std::complex<T> prolate_scattering_from_retained_blocks_single(
    const Rcpp::List& blocks,
    T chi_sw,
    T k_sw,
    T theta_body,
    T phi_body,
    T theta_scatter,
    T phi_scatter
) {
    // Rebuild the incident and scattered angular factors for the requested
    // geometry, then evaluate the stored spheroidal T-matrix blocks directly.
    int m_max = blocks.size() - 1;
    Rcpp::List last_block = blocks[m_max];
    Rcpp::IntegerVector n_seq_last = last_block["n_seq"];
    int n_max = n_seq_last[n_seq_last.size() - 1];

    auto t_blocks = prolate_tmatrix_blocks_from_rcpp<T>(blocks);
    auto smn_inc = compute_smn_matrix<T>(m_max, n_max, chi_sw, preccos(theta_body), true);
    auto smn_scat = compute_smn_matrix<T>(m_max, n_max, chi_sw, preccos(theta_scatter), true);
    auto azimuth = compute_azimuth<T>(m_max, phi_body, phi_scatter);
    auto raw_sum = compute_fbs_from_tmatrix_blocks<T>(
        m_max, n_max, azimuth, smn_inc, smn_scat, t_blocks
    );

    return std::complex<T>(0, -T(2) / k_sw) * raw_sum;
}

template<typename T>
Rcpp::ComplexMatrix prolate_scattering_grid_from_retained_blocks_single(
    const Rcpp::List& blocks,
    T chi_sw,
    T k_sw,
    T theta_body,
    T phi_body,
    const std::vector<T>& theta_scatter,
    const std::vector<T>& phi_scatter
) {
    // Precompute the scattered-angle angular basis for every theta slice, then
    // sweep over azimuth and evaluate the retained blocks pointwise.
    int m_max = blocks.size() - 1;
    Rcpp::List last_block = blocks[m_max];
    Rcpp::IntegerVector n_seq_last = last_block["n_seq"];
    int n_max = n_seq_last[n_seq_last.size() - 1];

    auto t_blocks = prolate_tmatrix_blocks_from_rcpp<T>(blocks);
    auto smn_inc = compute_smn_matrix<T>(m_max, n_max, chi_sw, preccos(theta_body), true);

    std::vector<std::vector<std::vector<T>>> smn_scat(theta_scatter.size());
    for (size_t i = 0; i < theta_scatter.size(); ++i) {
        smn_scat[i] = compute_smn_matrix<T>(
            m_max, n_max, chi_sw, preccos(theta_scatter[i]), true
        );
    }

    Rcpp::ComplexMatrix f_grid(theta_scatter.size(), phi_scatter.size());
    for (size_t i = 0; i < theta_scatter.size(); ++i) {
        for (size_t j = 0; j < phi_scatter.size(); ++j) {
            auto azimuth = compute_azimuth<T>(m_max, phi_body, phi_scatter[j]);
            auto raw_sum = compute_fbs_from_tmatrix_blocks<T>(
                m_max, n_max, azimuth, smn_inc, smn_scat[i], t_blocks
            );
            f_grid(i, j) = to_Rcomplex(std::complex<T>(0, -T(2) / k_sw) * raw_sum);
        }
    }

    return f_grid;
}

} // namespace

// Evaluate the stored prolate retained blocks in compiled code so the
// post-processing path uses the same spheroidal-function precision as the
// original boundary solve instead of rebuilding the angular factors in R.
// [[Rcpp::export]]
Rcpp::ComplexVector prolate_spheroid_scattering_from_tmatrix_cpp(
    Rcpp::DataFrame acoustics,
    Rcpp::List t_matrix,
    double theta_body,
    double phi_body,
    double theta_scatter,
    double phi_scatter,
    std::string precision = "double"
) {
    // Evaluate one retained block set per frequency at a single scattering
    // geometry.
    std::vector<double> chi_sw = Rcpp::as<std::vector<double>>(acoustics["chi_sw"]);
    std::vector<double> k_sw = Rcpp::as<std::vector<double>>(acoustics["k_sw"]);
    int n_freq = acoustics.nrows();
    Rcpp::ComplexVector f_scat(n_freq);

    if (precision == "quad") {
#ifdef __GNUC__
        for (int f = 0; f < n_freq; ++f) {
            f_scat[f] = to_Rcomplex(
                prolate_scattering_from_retained_blocks_single<__float128>(
                    t_matrix[f],
                    static_cast<__float128>(chi_sw[f]),
                    static_cast<__float128>(k_sw[f]),
                    static_cast<__float128>(theta_body),
                    static_cast<__float128>(phi_body),
                    static_cast<__float128>(theta_scatter),
                    static_cast<__float128>(phi_scatter)
                )
            );
        }
#else
        Rcpp::stop("Quad precision requires GCC and libquadmath.");
#endif
    } else {
        for (int f = 0; f < n_freq; ++f) {
            f_scat[f] = to_Rcomplex(
                prolate_scattering_from_retained_blocks_single<double>(
                    t_matrix[f],
                    chi_sw[f],
                    k_sw[f],
                    theta_body,
                    phi_body,
                    theta_scatter,
                    phi_scatter
                )
            );
        }
    }

    return f_scat;
}

// Evaluate one retained prolate T-matrix block set over a full
// theta-scatter/phi-scatter grid in compiled code.
// [[Rcpp::export]]
Rcpp::ComplexMatrix prolate_spheroid_scattering_grid_from_tmatrix_cpp(
    Rcpp::DataFrame acoustics,
    Rcpp::List t_matrix,
    double theta_body,
    double phi_body,
    Rcpp::NumericVector theta_scatter,
    Rcpp::NumericVector phi_scatter,
    std::string precision = "double"
) {
    // Evaluate the first retained-frequency block over a full
    // theta_scatter/phi_scatter grid.
    double chi_sw = Rcpp::as<std::vector<double>>(acoustics["chi_sw"])[0];
    double k_sw = Rcpp::as<std::vector<double>>(acoustics["k_sw"])[0];

    if (precision == "quad") {
#ifdef __GNUC__
        std::vector<__float128> theta_q(theta_scatter.size()), phi_q(phi_scatter.size());
        for (int i = 0; i < theta_scatter.size(); ++i) theta_q[i] = static_cast<__float128>(theta_scatter[i]);
        for (int j = 0; j < phi_scatter.size(); ++j) phi_q[j] = static_cast<__float128>(phi_scatter[j]);
        return prolate_scattering_grid_from_retained_blocks_single<__float128>(
            t_matrix,
            static_cast<__float128>(chi_sw),
            static_cast<__float128>(k_sw),
            static_cast<__float128>(theta_body),
            static_cast<__float128>(phi_body),
            theta_q,
            phi_q
        );
#else
        Rcpp::stop("Quad precision requires GCC and libquadmath.");
#endif
    }

    return prolate_scattering_grid_from_retained_blocks_single<double>(
        t_matrix,
        chi_sw,
        k_sw,
        theta_body,
        phi_body,
        Rcpp::as<std::vector<double>>(theta_scatter),
        Rcpp::as<std::vector<double>>(phi_scatter)
    );
}
