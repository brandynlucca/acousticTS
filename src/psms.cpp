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

    // Quad precision
#ifdef __GNUC__
    // Quad precision interface
    void profcn_cpp_interface_quad(
        __float128* c, int* m, int* lnum, int* ioprad, __float128* x1, int* iopang, int* iopnorm, 
        int* narg, __float128* arg, __float128* r1c, int* ir1e, __float128* r1dc, int* ir1de, 
        __float128* r2c, int* ir2e, __float128* r2dc, int* ir2de, int* naccr,
        __float128* s1c, int* is1e, __float128* s1dc, int* is1de, int* naccs
    );

    // Quad precision math functions
    __float128 powq(__float128, __float128);
#endif
}

// ============================================================================
// QUAD PRECISION HELPERS
// ============================================================================
template<typename T>
inline T pow10_typed(int exponent);

template<>
inline double pow10_typed<double>(int exponent) {
    return std::pow(10.0, static_cast<double>(exponent));
}

#ifdef __GNUC__
inline __float128 double_to_quad(double x) {
    return static_cast<__float128>(x);
}

inline double quad_to_double(__float128 x) {
    return static_cast<double>(x);
}

template<>
inline __float128 pow10_typed<__float128>(int exponent) {
    return powq(double_to_quad(10.0), double_to_quad(static_cast<double>(exponent)));
}
#endif

template<typename T>
inline T preccos(T x) { return std::cos(x); }
#ifdef __GNUC__
template<>
inline __float128 preccos(__float128 x) { return cosq(x); }
#endif

// ============================================================================
// HELPER FUNCTIONS
// ============================================================================
// Compute azimuth angles
template<typename T>
std::vector<T> compute_azimuth(int m_max, T phi_body, T phi_scatter) {
    std::vector<T> azimuth(m_max + 1);
    T dphi = phi_body - phi_scatter;
    for (int m = 0; m <= m_max; ++m)
        azimuth[m] = preccos(m * dphi);
    return azimuth;
}

inline bool is_na_int(int x) {
    return x == NA_INTEGER;
}

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

// ---- Helper: convert std::complex<double> to Rcomplex ----
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

// ============================================================================
// PROLATE SPHEROIDAL WAVE FUNCTION - PROFCN [DOUBLE]
// --------------------------------------------------
// This adds a safety net by padding the vector sizes to avoid segmentation 
// faults that can crash the program or the Rcpp interface
// ============================================================================
struct ProfcnResultDouble {
    std::vector<double> r1c, r1dc, r2c, r2dc;
    std::vector<int> ir1e, ir1de, ir2e, ir2de, naccr;
    std::vector<double> s1c, s1dc;
    std::vector<int> is1e, is1de, naccs;
};

ProfcnResultDouble cprofcn_double(
    double c,
    int m,
    int lnum,
    const std::vector<double>& arg,
    int ioprad,
    int iopnorm,
    int iopang,
    double x1
) {
    int narg = std::max<int>(arg.size(), 1);

    ProfcnResultDouble result;
    result.r1c.resize(lnum, 0.0);
    result.r1dc.resize(lnum, 0.0);
    result.r2c.resize(lnum, 0.0);
    result.r2dc.resize(lnum, 0.0);
    result.ir1e.resize(lnum, 0);
    result.ir1de.resize(lnum, 0);
    result.ir2e.resize(lnum, 0);
    result.ir2de.resize(lnum, 0);
    result.naccr.resize(lnum, 0);
    
    result.s1c.resize(lnum * narg, 0.0);
    result.s1dc.resize(lnum * narg, 0.0);
    result.is1e.resize(lnum * narg, 0);
    result.is1de.resize(lnum * narg, 0);
    result.naccs.resize(lnum * narg, 0);
    
    std::vector<double> arg_c = arg;

    for (int i = 0; i < narg; ++i) arg_c[i] = arg[i];

    // Interface
    profcn_cpp_interface(
        &c, &m, &lnum, &ioprad, &x1, &iopang, &iopnorm, &narg, arg_c.data(),
        result.r1c.data(), result.ir1e.data(), result.r1dc.data(), result.ir1de.data(),
        result.r2c.data(), result.ir2e.data(), result.r2dc.data(), result.ir2de.data(), 
        result.naccr.data(), result.s1c.data(), result.is1e.data(), result.s1dc.data(), 
        result.is1de.data(), result.naccs.data()
    );

    return result;
}

// ============================================================================
// PROLATE SPHEROIDAL WAVE FUNCTION - PROFCN [QUAD]
// --------------------------------------------------
// This adds a safety net by padding the vector sizes to avoid segmentation 
// faults that can crash the program or the Rcpp interface
// ============================================================================
#ifdef __GNUC__
struct ProfcnResultQuad {
    std::vector<__float128> r1c, r1dc, r2c, r2dc;
    std::vector<int> ir1e, ir1de, ir2e, ir2de, naccr;
    std::vector<__float128> s1c, s1dc;
    std::vector<int> is1e, is1de, naccs;
};

ProfcnResultQuad cprofcn_quad(
    double c,
    int m,
    int lnum,
    const std::vector<__float128>& arg,
    int ioprad,
    int iopnorm,
    int iopang,
    __float128 x1
) {
    int narg = std::max<int>(arg.size(), 1);
    
    ProfcnResultQuad result;
    result.r1c.resize(lnum, 0.0);
    result.r1dc.resize(lnum, 0.0);
    result.r2c.resize(lnum, 0.0);
    result.r2dc.resize(lnum, 0.0);
    result.ir1e.resize(lnum, 0);
    result.ir1de.resize(lnum, 0);
    result.ir2e.resize(lnum, 0);
    result.ir2de.resize(lnum, 0);
    result.naccr.resize(lnum, 0);
    
    result.s1c.resize(lnum * narg, 0.0);
    result.s1dc.resize(lnum * narg, 0.0);
    result.is1e.resize(lnum * narg, 0);
    result.is1de.resize(lnum * narg, 0);
    result.naccs.resize(lnum * narg, 0);
    
    std::vector<__float128> arg_c = arg;

    // Convert scalars to quad precision
    __float128 c_q = double_to_quad(c);
    __float128 x1_q = x1;
    int m_q = m;
    int lnum_q = lnum;
    int ioprad_q = ioprad;
    int iopnorm_q = iopnorm;
    int iopang_q = iopang;
    int narg_q = narg;

    // Interface
    profcn_cpp_interface_quad(
        &c_q, &m_q, &lnum_q, &ioprad_q, &x1_q, &iopang_q, &iopnorm_q, &narg_q, arg_c.data(),
        result.r1c.data(), result.ir1e.data(), result.r1dc.data(), result.ir1de.data(),
        result.r2c.data(), result.ir2e.data(), result.r2dc.data(), result.ir2de.data(), 
        result.naccr.data(), result.s1c.data(), result.is1e.data(), result.s1dc.data(), 
        result.is1de.data(), result.naccs.data()
    );

    return result;
}
#endif

// ============================================================================
// PROLATE SPHEROIDAL WAVE FUNCTION - PROFCN [PRIMARY WRAPPER]
// -----------------------------------------------------------
// This adds a safety net by padding the vector sizes to avoid segmentation 
// faults that can crash the program or the Rcpp interface
// ============================================================================
// Template version - works with double OR __float128
template<typename T>
struct ProfcnResult {
    std::vector<T> r1c, r1dc, r2c, r2dc;
    std::vector<int> ir1e, ir1de, ir2e, ir2de, naccr;
    std::vector<T> s1c, s1dc;
    std::vector<int> is1e, is1de, naccs;
};

// Forward declare the template interface
template<typename T>
ProfcnResult<T> cprofcn(
    T c, int m, int lnum,
    const std::vector<T>& arg,
    int ioprad, int iopnorm, int iopang, T x1
);

// -----------------------------------------------------------
// PRIMARY WRAPPER SPECIALIZATION - CPROFCN [DOUBLE]
// -----------------------------------------------------------
template<>
inline ProfcnResult<double> cprofcn<double>(
    double c, int m, int lnum,
    const std::vector<double>& arg,
    int ioprad, int iopnorm, int iopang, double x1
) {
    ProfcnResultDouble raw = cprofcn_double(c, m, lnum, arg, ioprad, iopnorm, iopang, x1);

    ProfcnResult<double> out;
    out.r1c   = std::move(raw.r1c);
    out.r1dc  = std::move(raw.r1dc);
    out.r2c   = std::move(raw.r2c);
    out.r2dc  = std::move(raw.r2dc);
    out.ir1e  = std::move(raw.ir1e);
    out.ir1de = std::move(raw.ir1de);
    out.ir2e  = std::move(raw.ir2e);
    out.ir2de = std::move(raw.ir2de);
    out.naccr = std::move(raw.naccr);
    out.s1c   = std::move(raw.s1c);
    out.s1dc  = std::move(raw.s1dc);
    out.is1e  = std::move(raw.is1e);
    out.is1de = std::move(raw.is1de);
    out.naccs = std::move(raw.naccs);

    return out;
}

// -----------------------------------------------------------
// PRIMARY WRAPPER SPECIALIZATION - CPROFCN [QUAD]
// -----------------------------------------------------------
#ifdef __GNUC__
template<>
inline ProfcnResult<__float128> cprofcn<__float128>(
    __float128 c, int m, int lnum,
    const std::vector<__float128>& arg,
    int ioprad, int iopnorm, int iopang, __float128 x1
) {
    ProfcnResultQuad raw = cprofcn_quad(c, m, lnum, arg, ioprad, iopnorm, iopang, x1);

    ProfcnResult<__float128> out;
    out.r1c   = std::move(raw.r1c);
    out.r1dc  = std::move(raw.r1dc);
    out.r2c   = std::move(raw.r2c);
    out.r2dc  = std::move(raw.r2dc);
    out.ir1e  = std::move(raw.ir1e);
    out.ir1de = std::move(raw.ir1de);
    out.ir2e  = std::move(raw.ir2e);
    out.ir2de = std::move(raw.ir2de);
    out.naccr = std::move(raw.naccr);
    out.s1c   = std::move(raw.s1c);
    out.s1dc  = std::move(raw.s1dc);
    out.is1e  = std::move(raw.is1e);
    out.is1de = std::move(raw.is1de);
    out.naccs = std::move(raw.naccs);

    return out;
}
#endif

// ============================================================================
// PROLATE SPHEROIDAL WAVE FUNCTION OF THE FIRST KIND [Smn - SCALAR]
// ============================================================================
template<typename T>
struct SmnResult {
    std::vector<T> value;
    std::vector<T> derivative;
};

template<typename T>
SmnResult<T> Smn_scalar(
    int m, int n, T c, const std::vector<T>& arg,
    bool normalize
) {
    int ioprad = 0;
    int iopnorm = normalize ? 1 : 0;
    int iopang = 2;
    T x1 = static_cast<T>(0);
    int lnum = n - m + 1;
    int narg = arg.size();
    int idx = n - m;
    
    ProfcnResult<T> result = cprofcn<T>(
        c, m, lnum, arg, ioprad, iopnorm, iopang, x1
    );
    
    SmnResult<T> out;
    out.value.resize(narg);
    out.derivative.resize(narg);
    
    for (int i = 0; i < narg; ++i) {
        int offset = idx + i * lnum;
        
        T pw = static_cast<T>(1);
        if (!is_na_int(result.is1e[offset])) {
            pw = pow10_typed<T>(result.is1e[offset]);
        }

        T pwde = static_cast<T>(1);
        if (!is_na_int(result.is1de[offset])) {
            pwde = pow10_typed<T>(result.is1de[offset]);
        }

        if (!is_na_real(result.s1c[offset])) {
            out.value[i] = result.s1c[offset] * pw;
        } else {
            out.value[i] = std::numeric_limits<T>::quiet_NaN();
        }

        if (!is_na_real(result.s1dc[offset])) {
            out.derivative[i] = result.s1dc[offset] * pwde;
        } else {
            out.derivative[i] = std::numeric_limits<T>::quiet_NaN();
        }
    }
    
    return out;
}

// ============================================================================
// PROLATE SPHEROIDAL WAVE FUNCTION OF THE FIRST KIND [Smn - VECTOR/MATRIX]
// ------------------------------------------------------------------------
// This outputs vectors of size n x m x eta
// ============================================================================
template<typename T>
struct SmnMatrixResult {
    std::vector<std::vector<T>> value;
    std::vector<std::vector<T>> derivative;
};

template<typename T>
SmnMatrixResult<T> Smn_matrix(
    const std::vector<int>& m,
    const std::vector<int>& n,
    T c,
    const std::vector<T>& eta,
    bool normalize = false
) {
    int m_len = m.size();
    int n_len = n.size();
    int eta_len = eta.size();

    // Validation
    if (m_len == 0 || n_len == 0)
        throw std::invalid_argument("'m' and 'n' must have at least one element.");
    if (eta_len == 0)
        throw std::invalid_argument("'eta' must have at least one element.");
    for (int mi : m)
        if (mi < 0) throw std::invalid_argument("All 'm' values must be >= 0.");
    for (int ni : n)
        if (ni < 0) throw std::invalid_argument("All 'n' values must be >= 0.");
    for (T e : eta)
        if (std::abs(e) > 1.0) {
            throw std::invalid_argument("|eta| must be <= 1 for the angular prolate domain.");
        }

    int ioprad = 0;
    int iopnorm = normalize ? 1 : 0;
    int iopang = 2;
    T x1 = static_cast<T>(0);

    SmnMatrixResult<T> out;

    // ------------------------------------------------------------------
    // CASE: Scalar (size[m] == size[n] == 1)
    if (m_len == 1 && n_len == 1) {
        if (n[0] < m[0])
            throw std::invalid_argument("'n' must be >= 'm'.");
        SmnResult<T> res = Smn_scalar<T>(m[0], n[0], c, eta, normalize);
        out.value.push_back(res.value);
        out.derivative.push_back(res.derivative);
        return out;
    }

    // ------------------------------------------------------------------
    // CASE: Scalar m, vector n (size[m] == 1 != size[n])
    if (m_len == 1 && n_len > 1) {
        int m_val = m[0];
        int n_min = *std::min_element(n.begin(), n.end());
        int n_max = *std::max_element(n.begin(), n.end());
        if (n_min < m_val)
            throw std::invalid_argument("All 'n' values must be >= 'm'.");
        int lnum = n_max - m_val + 1;

        ProfcnResult<T> result = cprofcn<T>(c, m_val, lnum, eta, ioprad, iopnorm, iopang, x1);

        out.value.resize(n_len, std::vector<T>(eta_len));
        out.derivative.resize(n_len, std::vector<T>(eta_len));
        for (int i = 0; i < n_len; ++i) {
            int idx = n[i] - m_val;
            for (int j = 0; j < eta_len; ++j) {
                int offset = idx + j * lnum;
                T pw = static_cast<T>(1), pwde = static_cast<T>(1);

                if (!is_na_int(result.is1e[offset])) pw = pow10_typed<T>(result.is1e[offset]);
                if (!is_na_int(result.is1de[offset])) pwde = pow10_typed<T>(result.is1de[offset]);

                bool valid_val = !is_na_real(result.s1c[offset]);
                if (valid_val) {
                    out.value[i][j] = result.s1c[offset] * pw;
                } else {
                    out.value[i][j] = std::numeric_limits<T>::quiet_NaN();
                }

                bool valid_der = !is_na_real(result.s1dc[offset]);
                if (valid_der) {
                    out.derivative[i][j] = result.s1dc[offset] * pwde;
                } else {
                    out.derivative[i][j] = std::numeric_limits<T>::quiet_NaN();
                }
            }
        }
        return out;
    }

    // ------------------------------------------------------------------
    // CASE: Pairwise vector m, vector n (size[m] > 1 == size[n])
    if (m_len == n_len) {
        for (int i = 0; i < m_len; ++i)
            if (n[i] < m[i])
                throw std::invalid_argument(
                    "For pairwise evaluation, each 'n[i]' must be >= 'm[i]'."
                );
        out.value.resize(m_len, std::vector<T>(eta_len));
        out.derivative.resize(m_len, std::vector<T>(eta_len));

        // Group by unique m values for efficiency
        std::map<int, std::vector<int>> m_to_indices;

        for (int i = 0; i < m_len; ++i) {
            m_to_indices[m[i]].push_back(i);
        }
        for (auto& kv : m_to_indices) {
            int m_val = kv.first;
            std::vector<int>& indices = kv.second;
            int n_max_local = n[indices[0]];
            for (size_t idx_i = 0; idx_i < indices.size(); ++idx_i) {
                if (n[indices[idx_i]] > n_max_local) n_max_local = n[indices[idx_i]];
            }
            int lnum = n_max_local - m_val + 1;
            ProfcnResult<T> result = cprofcn<T>(c, m_val, lnum, eta, ioprad, iopnorm, iopang, x1);
            for (size_t idx_i = 0; idx_i < indices.size(); ++idx_i) {
                int i = indices[idx_i];
                int idx = n[i] - m_val;
                for (int j = 0; j < eta_len; ++j) {
                    int offset = idx + j * lnum;

                    T pw = static_cast<T>(1), pwde = static_cast<T>(1);

                    if (!is_na_int(result.is1e[offset])) {
                        pw = pow10_typed<T>(result.is1e[offset]);
                    }
                    if (!is_na_int(result.is1de[offset])) {
                        pwde = pow10_typed<T>(result.is1de[offset]);
                    }

                    bool valid_val = !is_na_real(result.s1c[offset]);
                    if (valid_val) {
                        out.value[i][j] = result.s1c[offset] * pw;
                    } else {
                        out.value[i][j] = std::numeric_limits<T>::quiet_NaN();
                    }

                    bool valid_der = !is_na_real(result.s1dc[offset]);
                    if (valid_der) {
                        out.derivative[i][j] = result.s1dc[offset] * pwde;
                    } else {
                        out.derivative[i][j] = std::numeric_limits<T>::quiet_NaN();
                    }                    
                }
            }
        }
        return out;
    }

    // ------------------------------------------------------------------
    // CASE: Vector m, vector n (size[m] != size[n])
    // Outer product
    out.value.resize(m_len, std::vector<T>(n_len * eta_len));
    out.derivative.resize(m_len, std::vector<T>(n_len * eta_len));
    for (int i = 0; i < m_len; ++i) {
        int m_val = m[i];
        int n_max = *std::max_element(n.begin(), n.end());
        int lnum = n_max - m_val + 1;

        ProfcnResult<T> result = cprofcn<T>(c, m_val, lnum, eta, ioprad, iopnorm, iopang, x1);

        for (int j = 0; j < n_len; ++j) {
            if (n[j] < m_val) {
                for (int k = 0; k < eta_len; ++k) {
                    out.value[i][j * eta_len + k] = std::numeric_limits<T>::quiet_NaN();
                    out.derivative[i][j * eta_len + k] = std::numeric_limits<T>::quiet_NaN();
                }
            } else {
                int idx = n[j] - m_val;
                for (int k = 0; k < eta_len; ++k) {
                    int offset = idx + k * lnum;

                    T pw = static_cast<T>(1), pwde = static_cast<T>(1);

                    if (!is_na_int(result.is1e[offset])) {
                        pw = pow10_typed<T>(result.is1e[offset]);
                    }
                    if (!is_na_int(result.is1de[offset])) {
                        pwde = pow10_typed<T>(result.is1de[offset]);
                    }

                    bool valid_val = !is_na_real(result.s1c[offset]);
                    if (valid_val) {
                        out.value[i][j * eta_len + k] = result.s1c[offset] * pw;
                    } else {
                        out.value[i][j * eta_len + k] = std::numeric_limits<T>::quiet_NaN();
                    }

                    bool valid_der = !is_na_real(result.s1dc[offset]);
                    if (valid_der) {
                        out.derivative[i][j * eta_len + k] = result.s1dc[offset] * pwde;
                    } else {
                        out.derivative[i][j * eta_len + k] = std::numeric_limits<T>::quiet_NaN();
                    }    
                }
            }
        }
    }
    return out;
}

// -----------------------------------------------------------
// CONVERT C++ MATRIX TO R-COMPATIBLE LIST OF MATRICES
// -----------------------------------------------------------
template <typename T>
Rcpp::List Smn_cpp_to_rcpp_list(
    const SmnMatrixResult<T>& res,
    int m_len,
    int n_len,
    int eta_len
) {
    // Case: scalar m, vector n, scalar eta
    if (m_len == 1 && n_len > 1 && eta_len == 1) {
        Rcpp::NumericVector values(n_len), derivatives(n_len);
        for (int i = 0; i < n_len; ++i) {
            values[i] = static_cast<double>(res.value[i][0]);
            derivatives[i] = static_cast<double>(res.derivative[i][0]);
        }
        return Rcpp::List::create(
            Rcpp::Named("value") = values,
            Rcpp::Named("derivative") = derivatives
        );
    }
    // Case: scalar m, vector n, vector eta
    else if (m_len == 1 && n_len > 1) {
        Rcpp::NumericMatrix values(n_len, eta_len), derivatives(n_len, eta_len);
        for (int i = 0; i < n_len; ++i)
            for (int j = 0; j < eta_len; ++j) {
                values(i, j) = static_cast<double>(res.value[i][j]);
                derivatives(i, j) = static_cast<double>(res.derivative[i][j]);
            }
        return Rcpp::List::create(
            Rcpp::Named("value") = values,
            Rcpp::Named("derivative") = derivatives
        );
    }
    // Case: pairwise m/n, scalar eta
    else if (m_len == n_len && eta_len == 1) {
        Rcpp::NumericVector values(m_len), derivatives(m_len);
        for (int i = 0; i < m_len; ++i) {
            values[i] = static_cast<double>(res.value[i][0]);
            derivatives[i] = static_cast<double>(res.derivative[i][0]);
        }
        return Rcpp::List::create(
            Rcpp::Named("value") = values,
            Rcpp::Named("derivative") = derivatives
        );
    }
    // General case: if eta_len == 1, return m_len x n_len matrix as before
    else if (eta_len == 1) {
        Rcpp::NumericMatrix values(m_len, n_len), derivatives(m_len, n_len);
        for (int i = 0; i < m_len; ++i)
            for (int j = 0; j < n_len; ++j) {
                values(i, j) = static_cast<double>(res.value[i][j]);
                derivatives(i, j) = static_cast<double>(res.derivative[i][j]);
            }
        return Rcpp::List::create(
            Rcpp::Named("value") = values,
            Rcpp::Named("derivative") = derivatives
        );
    }
    // General case: eta_len > 1, return list of length eta_len, each with m x n matrix
    else {
        Rcpp::List out(eta_len);
        for (int k = 0; k < eta_len; ++k) {
            Rcpp::NumericMatrix values(m_len, n_len), derivatives(m_len, n_len);
            for (int i = 0; i < m_len; ++i)
                for (int j = 0; j < n_len; ++j) {
                    values(i, j) = static_cast<double>(res.value[i][j * eta_len + k]);
                    derivatives(i, j) = static_cast<double>(res.derivative[i][j * eta_len + k]);
                }
            out[k] = Rcpp::List::create(
                Rcpp::Named("value") = values,
                Rcpp::Named("derivative") = derivatives
            );
        }
        return out;
    }
}

// -----------------------------------------------------------
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
// RADIAL SPHEROIDAL WAVE FUNCTION OF THE FIRST KIND [Rmn - SCALAR]
// ============================================================================
template<typename T>
struct RmnResult {
    T val1;
    T der1;
    T val2;
    T der2;
    bool valid1;
    bool valid2;
};

template<typename T>
struct RadialValue {
    T val_real;
    T der_real;
    T val_imag;
    T der_imag;
};

// -----------------------------------------------------------
// HIGHER ORDER [3, 4] HELPER FUNCTION
// -----------------------------------------------------------
template<typename T>
RmnResult<T> Rmn_higher_order(int m, int n, int lnum, T c, T x1) {
    int iopnorm = 0, iopang = 2;
    std::vector<T> arg = {T(1)};
    int idx = n - m;

    auto result1 = cprofcn<T>(c, m, lnum, arg, 1, iopnorm, iopang, x1);
    auto result2 = cprofcn<T>(c, m, lnum, arg, 2, iopnorm, iopang, x1);

    T pw1 = 1, pwde1 = 1, pw2 = 1, pwde2 = 1;

    if (result1.ir1e.size() > static_cast<size_t>(idx) &&
        !is_na_int(result1.ir1e[idx]))
    {
        pw1 = pow10_typed<T>(result1.ir1e[idx]);
    }
    if (result1.ir1de.size() > static_cast<size_t>(idx) &&
        !is_na_int(result1.ir1de[idx]))
    {
        pwde1 = pow10_typed<T>(result1.ir1de[idx]);
    }
    if (result2.ir2e.size() > static_cast<size_t>(idx) &&
        !is_na_int(result2.ir2e[idx]))
    {
        pw2 = pow10_typed<T>(result2.ir2e[idx]);
    }
    if (result2.ir2de.size() > static_cast<size_t>(idx) &&
        !is_na_int(result2.ir2de[idx]))
    {
        pwde2 = pow10_typed<T>(result2.ir2de[idx]);
    }

    T R1 = std::numeric_limits<T>::quiet_NaN();
    if (result1.r1c.size() > static_cast<size_t>(idx) &&
        !is_na_real(result1.r1c[idx]))
    {
        R1 = result1.r1c[idx] * pw1;
    }

    T Rd1 = std::numeric_limits<T>::quiet_NaN();
    if (result1.r1dc.size() > static_cast<size_t>(idx) &&
        !is_na_real(result1.r1dc[idx]))
    {
        Rd1 = result1.r1dc[idx] * pwde1;
    }

    T R2 = std::numeric_limits<T>::quiet_NaN();
    if (result2.r2c.size() > static_cast<size_t>(idx) &&
        !is_na_real(result2.r2c[idx]))
    {
        R2 = result2.r2c[idx] * pw2;
    }

    T Rd2 = std::numeric_limits<T>::quiet_NaN();
    if (result2.r2dc.size() > static_cast<size_t>(idx) &&
        !is_na_real(result2.r2dc[idx]))
    {
        Rd2 = result2.r2dc[idx] * pwde2;
    }

    bool valid1 = !(result1.r1c.size() <= static_cast<size_t>(idx) || is_na_real(result1.r1c[idx]));
    bool valid2 = !(result2.r2c.size() <= static_cast<size_t>(idx) || is_na_real(result2.r2c[idx]));

    return {R1, Rd1, R2, Rd2, valid1, valid2};
}

template<typename T>
std::pair<std::complex<T>, std::complex<T>> Rmn_scalar(
    int m, int n, T c, T x1, int kind = 1
) {
    int lnum = n - m + 1;
    int lnum_safe = std::max(lnum, 2);
    T x1_r = x1 - T(1);
    int idx = n - m;

    if (kind == 1 || kind == 2) {
        int ioprad = kind;
        std::vector<T> arg = {T(1)};
        auto result = cprofcn<T>(c, m, lnum_safe, arg, ioprad, 0, 0, x1_r);
        const std::vector<T> *val, *der;
        const std::vector<int> *ie, *ide;
        if (kind == 1) {
            val = &result.r1c; der = &result.r1dc; ie = &result.ir1e; ide = &result.ir1de;
        } else {
            val = &result.r2c; der = &result.r2dc; ie = &result.ir2e; ide = &result.ir2de;
        }
        T pw = 1;
        if (ie->size() > static_cast<size_t>(idx) &&
            !is_na_int((*ie)[idx]))
        {
            pw = pow10_typed<T>((*ie)[idx]);
        }

        T pwde = 1;
        if (ide->size() > static_cast<size_t>(idx) &&
            !is_na_int((*ide)[idx]))
        {
            pwde = pow10_typed<T>((*ide)[idx]);
        }

        T R = std::numeric_limits<T>::quiet_NaN();
        if (val->size() > static_cast<size_t>(idx) &&
            !is_na_real((*val)[idx]))
        {
            R = (*val)[idx] * pw;
        }

        T Rd = std::numeric_limits<T>::quiet_NaN();
        if (der->size() > static_cast<size_t>(idx) &&
            !is_na_real((*der)[idx]))
        {
            Rd = (*der)[idx] * pwde;
        }

        // Specific numeric stability check for Rmn[2]
        if (kind == 2) {
            const T tiny = T(1e-300);
            bool R_bad = !std::isfinite(static_cast<double>(R)) ||
                        !std::isfinite(static_cast<double>(Rd)) ||
                        (std::abs(R) < tiny && std::abs(Rd) < tiny);
            if (R_bad) {
                R = std::numeric_limits<T>::quiet_NaN();
                Rd = std::numeric_limits<T>::quiet_NaN();
            }
        }
        return {std::complex<T>(R, 0), std::complex<T>(Rd, 0)};
    }

    // Higher order: Rmn[3] and Rmn[4]
    RmnResult<T> result = Rmn_higher_order<T>(m, n, lnum_safe, c, x1_r);
    std::complex<T> val, der;
    if (kind == 3) {
        val = std::complex<T>(result.val1, result.val2);
        der = std::complex<T>(result.der1, result.der2);
    } else {
        val = std::complex<T>(result.val1, -result.val2);
        der = std::complex<T>(result.der1, -result.der2);
    }
    return {val, der};
}

// ============================================================================
// RADIAL PROLATE SPHEROIDAL WAVE FUNCTIONS [Rmn - VECTOR/MATRIX]
// ------------------------------------------------------------------------
// This outputs vectors of size n x m with scalar xi
// ============================================================================
template<typename T>
struct RmnMatrixResult {
    std::vector<std::vector<std::complex<T>>> value;
    std::vector<std::vector<std::complex<T>>> derivative;
};

// -----------------------------------------------------------
// BATCH HELPER
// -----------------------------------------------------------
template<typename T>
RadialValue<T> extract_radial_from_batch(
    const ProfcnResult<T>& result1,
    const ProfcnResult<T>& result2, // ignored for real kinds
    int idx,
    int kind,
    int m,
    int n,
    T c,
    T x1
) {
    RadialValue<T> out;

    if (kind == 1 || kind == 2) {
        const std::vector<T> *rc, *rdc;
        const std::vector<int> *ie, *ide;
        if (kind == 1) {
            rc = &result1.r1c;
            rdc = &result1.r1dc;
            ie = &result1.ir1e;
            ide = &result1.ir1de;
        } else {
            rc = &result1.r2c;
            rdc = &result1.r2dc;
            ie = &result1.ir2e;
            ide = &result1.ir2de;
        }
        T pw = 1;
        if (ie->size() > static_cast<size_t>(idx) && !is_na_int((*ie)[idx])) {
            pw = pow10_typed<T>((*ie)[idx]);
        }
        T pwde = 1;
        if (ide->size() > static_cast<size_t>(idx) && !is_na_int((*ide)[idx])) {
            pwde = pow10_typed<T>((*ide)[idx]);
        }
        out.val_real = std::numeric_limits<T>::quiet_NaN();
        if (rc->size() > static_cast<size_t>(idx) && !is_na_real((*rc)[idx])) {
            out.val_real = (*rc)[idx] * pw;
        }
        out.der_real = std::numeric_limits<T>::quiet_NaN();
        if (rdc->size() > static_cast<size_t>(idx) && !is_na_real((*rdc)[idx])) {
            out.der_real = (*rdc)[idx] * pwde;
        }
        out.val_imag = 0;
        out.der_imag = 0;
        if (kind == 2) {
            const T tiny = T(1e-300);
            bool bad = !std::isfinite(static_cast<double>(out.val_real)) ||
                    !std::isfinite(static_cast<double>(out.der_real)) ||
                    (std::abs(out.val_real) < tiny && std::abs(out.der_real) < tiny);
            if (bad) {
                out.val_real = std::numeric_limits<T>::quiet_NaN();
                out.der_real = std::numeric_limits<T>::quiet_NaN();
            }
        }
        return out;
    }

    // Higher orders - Rmn[3] and Rmn[4]
    T pw1 = 1;
    if (result1.ir1e.size() > static_cast<size_t>(idx) &&
        !is_na_int(result1.ir1e[idx])) {
        pw1 = pow10_typed<T>(result1.ir1e[idx]);
    }
    T pwde1 = 1;
    if (result1.ir1de.size() > static_cast<size_t>(idx) &&
        !is_na_int(result1.ir1de[idx])) {
        pwde1 = pow10_typed<T>(result1.ir1de[idx]);
    }
    out.val_real = std::numeric_limits<T>::quiet_NaN();
    if (result1.r1c.size() > static_cast<size_t>(idx) &&
        !is_na_real(result1.r1c[idx])) {
        out.val_real = result1.r1c[idx] * pw1;
    }
    out.der_real = std::numeric_limits<T>::quiet_NaN();
    if (result1.r1dc.size() > static_cast<size_t>(idx) &&
        !is_na_real(result1.r1dc[idx])) {
        out.der_real = result1.r1dc[idx] * pwde1;
    }

    T pw2 = 1;
    if (result2.ir2e.size() > static_cast<size_t>(idx) &&
        !is_na_int(result2.ir2e[idx])) {
        pw2 = pow10_typed<T>(result2.ir2e[idx]);
    }
    T pwde2 = 1;
    if (result2.ir2de.size() > static_cast<size_t>(idx) &&
        !is_na_int(result2.ir2de[idx])) {
        pwde2 = pow10_typed<T>(result2.ir2de[idx]);
    }
    out.val_imag = std::numeric_limits<T>::quiet_NaN();
    if (result2.r2c.size() > static_cast<size_t>(idx) &&
        !is_na_real(result2.r2c[idx])) {
        out.val_imag = result2.r2c[idx] * pw2;
    }
    out.der_imag = std::numeric_limits<T>::quiet_NaN();
    if (result2.r2dc.size() > static_cast<size_t>(idx) &&
        !is_na_real(result2.r2dc[idx])) {
        out.der_imag = result2.r2dc[idx] * pwde2;
    }

    const T tiny = T(1e-300);
    bool bad = !std::isfinite(static_cast<double>(out.val_imag)) ||
            !std::isfinite(static_cast<double>(out.der_imag)) ||
            (std::abs(out.val_imag) < tiny && std::abs(out.der_imag) < tiny);
    if (bad) {
        out.val_imag = std::numeric_limits<T>::quiet_NaN();
        out.der_imag = std::numeric_limits<T>::quiet_NaN();
    }
    if (kind == 4) {
        out.val_imag = -out.val_imag;
        out.der_imag = -out.der_imag;
    }
    return out;
}

// -----------------------------------------------------------
// MAIN FUNCTION 
// -----------------------------------------------------------
template<typename T>
RmnMatrixResult<T> Rmn_matrix(
    const std::vector<int>& m,
    const std::vector<int>& n,
    T c,
    T x1,
    int kind = 1
) {
    int m_len = m.size();
    int n_len = n.size();
    bool is_complex = (kind == 3 || kind == 4);
    T x1_r = x1 - T(1);

    RmnMatrixResult<T> out;
    out.value.resize(m_len, std::vector<std::complex<T>>(n_len));
    out.derivative.resize(m_len, std::vector<std::complex<T>>(n_len));

    // ------------------------------------------------------------------
    // CASE: Scalar (size[m] == size[n] == 1)
    if (m_len == 1 && n_len == 1) {
        if (n[0] < m[0])
            throw std::invalid_argument("'n' must be >= 'm'.");
        auto res = Rmn_scalar<T>(m[0], n[0], c, x1, kind);
        out.value[0][0] = res.first;
        out.derivative[0][0] = res.second;
        return out;
    }

    // ------------------------------------------------------------------
    // CASE: Scalar m, vector n (size[m] == 1 != size[n])
    if (m_len == 1 && n_len > 1) {
        int m_val = m[0];
        int n_min = *std::min_element(n.begin(), n.end());
        int n_max = *std::max_element(n.begin(), n.end());
        if (n_min < m_val)
            throw std::invalid_argument("All 'n' values must be >= 'm'.");
        int lnum = n_max - m_val + 1;
        auto result1 = cprofcn<T>(c, m_val, lnum, {T(1)}, (is_complex ? 1 : kind), 0, 0, x1_r);
        ProfcnResult<T> result2;
        if (is_complex) result2 = cprofcn<T>(c, m_val, lnum, {T(1)}, 2, 0, 0, x1_r);
        for (int i = 0; i < n_len; ++i) {
            if (n[i] < m_val) {
                out.value[0][i] = std::complex<T>(std::numeric_limits<T>::quiet_NaN(), 0);
                out.derivative[0][i] = std::complex<T>(std::numeric_limits<T>::quiet_NaN(), 0);
                continue;
            }
            RadialValue<T> rv = is_complex
                ? extract_radial_from_batch<T>(
                    result1, result2, n[i] - m_val, kind, m_val, n[i], c, x1
                )
                : extract_radial_from_batch<T>(
                    result1, result1, n[i] - m_val, kind, m_val, n[i], c, x1
                );
            out.value[0][i] = std::complex<T>(rv.val_real, rv.val_imag);
            out.derivative[0][i] = std::complex<T>(rv.der_real, rv.der_imag);
        }
        return out;
    }

    // ------------------------------------------------------------------
    // CASE: Pairwise vector m, vector n (size[m] > 1 == size[n])
    if (m_len == n_len) {
        for (int i = 0; i < m_len; ++i)
            if (n[i] < m[i])
                throw std::invalid_argument(
                    "For pairwise evaluation, each 'n[i]' must be >= 'm[i]'."
                );
        std::map<int, std::vector<int>> m_to_indices;
        for (int i = 0; i < m_len; ++i) m_to_indices[m[i]].push_back(i);
        for (auto& kv : m_to_indices) {
            int m_val = kv.first;
            std::vector<int>& indices = kv.second;
            int n_max_local = n[indices[0]];
            for (size_t idx_i = 0; idx_i < indices.size(); ++idx_i)
                if (n[indices[idx_i]] > n_max_local) n_max_local = n[indices[idx_i]];
            int lnum = n_max_local - m_val + 1;
            auto result1 = cprofcn<T>(c, m_val, lnum, {T(1)}, (is_complex ? 1 : kind), 0, 0, x1_r);
            ProfcnResult<T> result2;
            if (is_complex) result2 = cprofcn<T>(c, m_val, lnum, {T(1)}, 2, 0, 0, x1_r);
            for (size_t idx_i = 0; idx_i < indices.size(); ++idx_i) {
                int i = indices[idx_i];
                RadialValue<T> rv = is_complex
                    ? extract_radial_from_batch<T>(
                        result1, result2, n[i] - m_val, kind, m_val, n[i], c, x1
                    )
                    : extract_radial_from_batch<T>(
                        result1, result1, n[i] - m_val, kind, m_val, n[i], c, x1
                    );
                out.value[i][i] = std::complex<T>(rv.val_real, rv.val_imag);
                out.derivative[i][i] = std::complex<T>(rv.der_real, rv.der_imag);
            }
        }
        return out;
    }

    // ------------------------------------------------------------------
    // CASE: Vector m, vector n (size[m] != size[n])
    // Outer product
    for (int i = 0; i < m_len; ++i) {
        int m_val = m[i];
        int n_max = *std::max_element(n.begin(), n.end());
        int lnum = n_max - m_val + 1;
        auto result1 = cprofcn<T>(c, m_val, lnum, {T(1)}, (is_complex ? 1 : kind), 0, 0, x1_r);
        ProfcnResult<T> result2;
        if (is_complex) result2 = cprofcn<T>(c, m_val, lnum, {T(1)}, 2, 0, 0, x1_r);
        for (int j = 0; j < n_len; ++j) {
            if (n[j] < m_val) {
                out.value[i][j] = std::complex<T>(std::numeric_limits<T>::quiet_NaN(), 0);
                out.derivative[i][j] = std::complex<T>(std::numeric_limits<T>::quiet_NaN(), 0);
                continue;
            }
            RadialValue<T> rv = is_complex
                ? extract_radial_from_batch<T>(
                    result1, result2, n[j] - m_val, kind, m_val, n[j], c, x1
                )
                : extract_radial_from_batch<T>(
                    result1, result1, n[j] - m_val, kind, m_val, n[j], c, x1
                );
            out.value[i][j] = std::complex<T>(rv.val_real, rv.val_imag);
            out.derivative[i][j] = std::complex<T>(rv.der_real, rv.der_imag);
        }
    }
    return out;
}

// -----------------------------------------------------------
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
// RADIAL PROLATE SPHEROIDAL WAVE FUNCTION MATRICES
// ----------------------------------------------------------------------------
template<typename T>
struct RadialMatrixResult {
    std::vector<std::vector<T>> value;
    std::vector<std::vector<T>> derivative;
};

template<typename T>
struct RadialComplexMatrixResult {
    std::vector<std::vector<std::complex<T>>> value;
    std::vector<std::vector<std::complex<T>>> derivative;
};

template<typename T>
RadialMatrixResult<T> radial_external_incident_matrix(
    int m_max,
    int n_max,
    T chi_sw,
    T xi,
    const std::vector<std::vector<T>>& smn_matrix
) {
    RadialMatrixResult<T> out;
    out.value.assign(m_max + 1, std::vector<T>(n_max + 1, std::numeric_limits<T>::quiet_NaN()));
    out.derivative.assign(
        m_max + 1, std::vector<T>(n_max + 1, std::numeric_limits<T>::quiet_NaN())
    );

    for (int m = 0; m <= m_max; ++m) {
        int size = n_max - m + 1;
        std::vector<int> n_vec;
        std::vector<int> valid_indices;
        for (int i = 0; i < size; ++i) {
            int n = m + i;
            T snm_val = smn_matrix[m][n];
            if (!std::isnan(static_cast<double>(snm_val)) && snm_val != T(0)) {
                n_vec.push_back(n);
                valid_indices.push_back(n);
            }
        }
        // Set zeros for smn=0 cases
        for (int n = m; n <= n_max; ++n) {
            T snm_val = smn_matrix[m][n];
            if (std::isnan(static_cast<double>(snm_val)) || snm_val == T(0)) {
                out.value[m][n] = T(0);
                out.derivative[m][n] = T(0);
            }
        }
        if (n_vec.empty()) continue;
        auto rmn_result = Rmn_matrix<T>({m}, n_vec, chi_sw, xi, 1);
        for (size_t idx = 0; idx < valid_indices.size(); ++idx) {
            int n = valid_indices[idx];
            out.value[m][n] = rmn_result.value[0][idx].real();
            out.derivative[m][n] = rmn_result.derivative[0][idx].real();
        }
    }
    return out;
}

template<typename T>
RadialComplexMatrixResult<T> radial_external_scattering_matrix(
    int m_max,
    int n_max,
    T chi_sw,
    T xi,
    const std::vector<std::vector<T>>& smn_matrix
) {
    RadialComplexMatrixResult<T> out;
    out.value.assign(
        m_max + 1, 
        std::vector<std::complex<T>>(
            n_max + 1, 
            std::complex<T>(
                std::numeric_limits<T>::quiet_NaN(), std::numeric_limits<T>::quiet_NaN()
            )
        )
    );
    out.derivative.assign(
        m_max + 1, 
        std::vector<std::complex<T>>(
            n_max + 1, 
            std::complex<T>(
                std::numeric_limits<T>::quiet_NaN(), std::numeric_limits<T>::quiet_NaN()
            )
        )
    );

    for (int m = 0; m <= m_max; ++m) {
        int size = n_max - m + 1;
        std::vector<int> n_vec;
        std::vector<int> valid_indices;
        for (int i = 0; i < size; ++i) {
            int n = m + i;
            T snm_val = smn_matrix[m][n];
            if (!std::isnan(static_cast<double>(snm_val)) && snm_val != T(0)) {
                n_vec.push_back(n);
                valid_indices.push_back(n);
            }
        }
        // Set zeros for smn=0 cases
        for (int n = m; n <= n_max; ++n) {
            T snm_val = smn_matrix[m][n];
            if (std::isnan(static_cast<double>(snm_val)) || snm_val == T(0)) {
                out.value[m][n] = std::complex<T>(0, 0);
                out.derivative[m][n] = std::complex<T>(0, 0);
            }
        }
        if (n_vec.empty()) continue;
        auto rmn_result = Rmn_matrix<T>({m}, n_vec, chi_sw, xi, 3);
        for (size_t idx = 0; idx < valid_indices.size(); ++idx) {
            int n = valid_indices[idx];
            out.value[m][n] = rmn_result.value[0][idx];
            out.derivative[m][n] = rmn_result.derivative[0][idx];
        }
    }
    return out;
}

template<typename T>
RadialMatrixResult<T> radial_internal_incident_matrix(
    int m_max,
    int n_max,
    T chi_body,
    T xi
) {
    RadialMatrixResult<T> out;
    out.value.assign(
        m_max + 1, std::vector<T>(n_max + 1, std::numeric_limits<T>::quiet_NaN())
    );
    out.derivative.assign(
        m_max + 1, std::vector<T>(n_max + 1, std::numeric_limits<T>::quiet_NaN())
    );

    for (int m = 0; m <= m_max; ++m) {
        int size = n_max - m + 1;
        std::vector<int> n_vec(size);
        for (int i = 0; i < size; ++i) n_vec[i] = m + i;
        auto rmn_result = Rmn_matrix<T>({m}, n_vec, chi_body, xi, 1);
        for (int i = 0; i < size; ++i) {
            int n = m + i;
            out.value[m][n] = rmn_result.value[0][i].real();
            out.derivative[m][n] = rmn_result.derivative[0][i].real();
        }
    }
    return out;
}
// ============================================================================
// BOUNDARY COUPLING MATRICES
// ----------------------------------------------------------------------------
// The incident and scattering boundary coupling matrices correspond to 
// E_{n,l}^{m(1)} and E_{n,l}^{m(3)} from Eq. (4) in Furusawa (1988)
// ============================================================================
template<typename T>
std::vector<std::vector<std::complex<T>>> boundary_coupling_incident_matrix(
    int m_max,
    int n_max,
    T density_body,
    T density_sw,
    const std::vector<std::vector<T>>& e1_val_mat,
    const std::vector<std::vector<T>>& e1_der_mat,
    const std::vector<std::vector<T>>& i1_val_mat,
    const std::vector<std::vector<T>>& i1_der_mat,
    const std::vector<std::vector<T>>& smn_matrix
) {
    std::vector<std::vector<std::complex<T>>> result(m_max + 1);
    for (int m = 0; m <= m_max; ++m) {
        int size = n_max - m + 1;
        std::vector<std::complex<T>> mat(size * size, std::complex<T>(0, 0));
        for (int li = 0; li < size; ++li) {
            int ell = m + li;
            T R_ml_1_val = i1_val_mat[m][ell];
            T R_ml_1_d   = i1_der_mat[m][ell];
            for (int ni = 0; ni < size; ++ni) {
                int n = m + ni;
                T snm_n = smn_matrix[m][n];
                if (std::isnan(static_cast<double>(snm_n)) || snm_n == T(0)) continue;
                T Rmn_i1_val = e1_val_mat[m][n];
                T Rmn_i1_d   = e1_der_mat[m][n];
                T factor = (R_ml_1_d != T(0))
                    ? (density_body * R_ml_1_val) / (density_sw * R_ml_1_d)
                    : T(0);
                T val = Rmn_i1_val - factor * Rmn_i1_d;
                mat[li * size + ni] = std::complex<T>(val, 0);
            }
        }
        result[m] = mat;
    }
    return result;
}

template<typename T>
std::vector<std::vector<std::complex<T>>> boundary_coupling_scattering_matrix(
    int m_max,
    int n_max,
    T density_body,
    T density_sw,
    const std::vector<std::vector<std::complex<T>>>& e3_val_mat,
    const std::vector<std::vector<std::complex<T>>>& e3_der_mat,
    const std::vector<std::vector<T>>& i1_val_mat,
    const std::vector<std::vector<T>>& i1_der_mat,
    const std::vector<std::vector<T>>& smn_matrix
) {
    std::vector<std::vector<std::complex<T>>> result(m_max + 1);
    for (int m = 0; m <= m_max; ++m) {
        int size = n_max - m + 1;
        std::vector<std::complex<T>> mat(size * size, std::complex<T>(0, 0));
        for (int li = 0; li < size; ++li) {
            int ell = m + li;
            T R_ml_1_val = i1_val_mat[m][ell];
            T R_ml_1_d   = i1_der_mat[m][ell];
            for (int ni = 0; ni < size; ++ni) {
                int n = m + ni;
                T snm_n = smn_matrix[m][n];
                if (std::isnan(static_cast<double>(snm_n)) || snm_n == T(0)) continue;
                std::complex<T> Rmn_i3_val = e3_val_mat[m][n];
                std::complex<T> Rmn_i3_d   = e3_der_mat[m][n];
                T factor = (R_ml_1_d != T(0))
                    ? (density_body * R_ml_1_val) / (density_sw * R_ml_1_d)
                    : T(0);
                std::complex<T> val = Rmn_i3_val - factor * Rmn_i3_d;
                mat[li * size + ni] = val;
            }
        }
        result[m] = mat;
    }
    return result;
}

// ============================================================================
// SIMPLIFIED BOUNDARY COUPLING MATRICES
// ----------------------------------------------------------------------------
// The incident and scattering boundary coupling matrices correspond to 
// E_{n,l}^{m(1)} and E_{n,l}^{m(3)} from Eq. (5) in Furusawa (1988)
// ============================================================================
template<typename T>
std::vector<std::vector<std::complex<T>>> simplified_boundary_coupling_incident_matrix(
    int m_max,
    int n_max,
    T density_body,
    T density_sw,
    const std::vector<std::vector<T>>& e1_val_mat,
    const std::vector<std::vector<T>>& e1_der_mat,
    const std::vector<std::vector<T>>& i1_val_mat,
    const std::vector<std::vector<T>>& i1_der_mat
) {
    // Compute lower trianglular result
    std::vector<std::vector<std::complex<T>>> lower(m_max + 1);
    for (int m = 0; m <= m_max; ++m) {
        int size = n_max - m + 1;
        lower[m].resize(size, std::complex<T>(std::numeric_limits<T>::quiet_NaN(), 0));
        for (int li = 0; li < size; ++li) {
            int ell = m + li;
            T R_ml_1_val = i1_val_mat[m][ell];
            T R_ml_1_d   = i1_der_mat[m][ell];
            int n = ell;
            T Rmn_i1_val = e1_val_mat[m][n];
            T Rmn_i1_d   = e1_der_mat[m][n];
            if (is_na_real(Rmn_i1_val)) {
                lower[m][li] = std::complex<T>(std::numeric_limits<T>::quiet_NaN(), 0);
                continue;
            }
            T factor = (R_ml_1_d != T(0))
                ? (density_body * R_ml_1_val) / (density_sw * R_ml_1_d)
                : T(0);
            lower[m][li] = std::complex<T>(Rmn_i1_val - factor * Rmn_i1_d, 0);
        }
    }

    // Convert to rectangular matrix (m_max+1) x (n_max+1), NA for n < m
    std::vector<std::vector<std::complex<T>>> result_mat(
        m_max + 1, 
        std::vector<std::complex<T>>(
            n_max + 1, std::complex<T>(std::numeric_limits<T>::quiet_NaN(), 0)
        )
    );
    for (int m = 0; m <= m_max; ++m) {
        int size = lower[m].size();
        for (int i = 0; i < size; ++i) {
            int n = m + i;
            result_mat[m][n] = lower[m][i];
        }
    }

    return result_mat;
}

template<typename T>
std::vector<std::vector<std::complex<T>>> simplified_boundary_coupling_scattering_matrix(
    int m_max,
    int n_max,
    T density_body,
    T density_sw,
    const std::vector<std::vector<std::complex<T>>>& e3_val_mat,
    const std::vector<std::vector<std::complex<T>>>& e3_der_mat,
    const std::vector<std::vector<T>>& i1_val_mat,
    const std::vector<std::vector<T>>& i1_der_mat
) {
    // Compute lower trianglular result
    std::vector<std::vector<std::complex<T>>> lower(m_max + 1);
    for (int m = 0; m <= m_max; ++m) {
        int size = n_max - m + 1;
        lower[m].resize(
            size, std::complex<T>(std::numeric_limits<T>::quiet_NaN(), 
            std::numeric_limits<T>::quiet_NaN())
        );
        for (int li = 0; li < size; ++li) {
            int ell = m + li;
            T R_ml_1_val = i1_val_mat[m][ell];
            T R_ml_1_d   = i1_der_mat[m][ell];
            int n = ell;
            std::complex<T> Rmn_i3_val = e3_val_mat[m][n];
            std::complex<T> Rmn_i3_d   = e3_der_mat[m][n];
            T factor = (R_ml_1_d != T(0))
                ? (density_body * R_ml_1_val) / (density_sw * R_ml_1_d)
                : T(0);
            lower[m][li] = Rmn_i3_val - factor * Rmn_i3_d;
        }
    }

    // Convert to rectangular matrix (m_max+1) x (n_max+1), NA for n < m
    std::vector<std::vector<std::complex<T>>> result_mat(
        m_max + 1, 
        std::vector<std::complex<T>>(
            n_max + 1, 
            std::complex<T>(std::numeric_limits<T>::quiet_NaN(), 
            std::numeric_limits<T>::quiet_NaN())
        )
    );
    for (int m = 0; m <= m_max; ++m) {
        int size = lower[m].size();
        for (int i = 0; i < size; ++i) {
            int n = m + i;
            result_mat[m][n] = lower[m][i];
        }
    }

    return result_mat;
}

// ============================================================================
// KERNEL MATRICES FOR BOUNDARY CONDITIONS
// ----------------------------------------------------------------------------
// The incident and scattering kernel matrices correspond to K_{n,l}^{m(1)} and
// K_{n,l}^{m(2)} from Eq. (4) in Furusawa (1988)
// ============================================================================
template<typename T>
std::vector<std::vector<std::complex<T>>> kernel_incident_matrix(
    int m_max,
    int n_max,
    const std::vector<std::vector<T>>& smn_matrix,
    const std::vector<std::vector<std::vector<T>>>& expansion_matrix,
    const std::vector<std::vector<std::complex<T>>>& E1_coupling
) {
    std::vector<std::vector<std::complex<T>>> result(m_max + 1);
    for (int m = 0; m <= m_max; ++m) {
        int size = n_max - m + 1;
        std::vector<std::complex<T>> mat(size * size, std::complex<T>(0, 0));
        for (int li = 0; li < size; ++li) {
            for (int ni = 0; ni < size; ++ni) {
                int n = m + ni;
                T snm_val = smn_matrix[m][n];
                if (std::isnan(static_cast<double>(snm_val)) || snm_val == T(0)) {
                    mat[li * size + ni] = std::complex<T>(0, 0);
                    continue;
                }
                T alpha = expansion_matrix[m][li][ni];
                std::complex<T> E1_nl = E1_coupling[m][li * size + ni];
                std::complex<T> jn = std::pow(std::complex<T>(0, 1), n);
                std::complex<T> K1_val = jn * snm_val * alpha * E1_nl;
                mat[li * size + ni] = K1_val;
            }
        }
        result[m] = mat;
    }
    return result;
}

template<typename T>
std::vector<std::vector<std::complex<T>>> kernel_scattering_matrix(
    int m_max,
    int n_max,
    const std::vector<std::vector<T>>& smn_matrix,
    const std::vector<std::vector<std::vector<T>>>& expansion_matrix,
    const std::vector<std::vector<std::complex<T>>>& E3_coupling
) {
    std::vector<std::vector<std::complex<T>>> result(m_max + 1);
    for (int m = 0; m <= m_max; ++m) {
        int size = n_max - m + 1;
        std::vector<std::complex<T>> mat(size * size, std::complex<T>(0, 0));
        for (int li = 0; li < size; ++li) {
            for (int ni = 0; ni < size; ++ni) {
                int n = m + ni;
                T snm_val = smn_matrix[m][n];
                if (std::isnan(static_cast<double>(snm_val)) || snm_val == T(0)) {
                    mat[li * size + ni] = std::complex<T>(0, 0);
                    continue;
                }
                T alpha = expansion_matrix[m][li][ni];
                std::complex<T> E3_nl = E3_coupling[m][li * size + ni];
                std::complex<T> jn = std::pow(std::complex<T>(0, 1), n);
                std::complex<T> K3_val = jn * snm_val * alpha * E3_nl;
                mat[li * size + ni] = K3_val;
            }
        }
        result[m] = mat;
    }
    return result;
}

// ============================================================================
// EXPANSION COEFFICIENT MATRIX FUNCTIONS FOR FLUID-FILLED PROLATE SPHEROIDS
// ----------------------------------------------------------------------------
// The expansion coefficient matrices for scattered waves of different 
// boundary conditions from Furusawa (1988) from Eq. (4)
// ============================================================================
// -----------------------------------------------------------
// PRE-COMPUTED SMN MATRIX USED FOR STREAMLINING COMPUTATIONAL EFFICIENCY
// -----------------------------------------------------------
template<typename T>
std::vector<std::vector<T>> compute_smn_matrix(
    int m_max,
    int n_max,
    T chi_sw,
    T eta_scalar,
    bool normalize = true
) {
    // Result is (m_max+1) x (n_max+1) matrix, NA (quiet_NaN) where n < m
    std::vector<std::vector<T>> result(
        m_max + 1, std::vector<T>(n_max + 1, std::numeric_limits<T>::quiet_NaN())
    );

    for (int m = 0; m <= m_max; ++m) {
        int size = n_max - m + 1;
        if (size <= 0) continue;
        std::vector<int> n_vec(size);
        for (int i = 0; i < size; ++i) n_vec[i] = m + i;
        std::vector<T> eta_vec = {eta_scalar};

        auto smn_result = Smn_matrix<T>(
            std::vector<int>{m}, n_vec, chi_sw, eta_vec, normalize
        );

        // smn_result.value is [n][eta], eta_len == 1
        for (int i = 0; i < size; ++i) {
            int n = m + i;
            result[m][n] = smn_result.value[i][0];
        }
    }
    return result;
}
// -----------------------------------------------------------
// ALPHA INTEGRATION COEFFICIENT
// -----------------------------------------------------------
template<typename T>
std::vector<std::vector<std::vector<T>>> alpha_integration_coefficient(
    int m_max,
    int n_max,
    T chi_sw,
    T chi_body,
    const std::vector<T>& nodes,
    const std::vector<T>& weights
) {
    std::vector<std::vector<std::vector<T>>> coef_list(m_max + 1);

    for (int m = 0; m <= m_max; ++m) {
        int size = n_max - m + 1;
        if (size <= 0) {
            coef_list[m] = std::vector<std::vector<T>>(); // empty
            continue;
        }
        // n_vec = [m, m+1, ..., n_max]
        std::vector<int> n_vec(size);
        for (int i = 0; i < size; ++i) n_vec[i] = m + i;

        // Batched calls for all n at all nodes
        auto smn_ext = Smn_matrix<T>(
            std::vector<int>{m}, n_vec, chi_sw, nodes, true
        ).value; // shape: [n][eta]
        auto smn_int = Smn_matrix<T>(
            std::vector<int>{m}, n_vec, chi_body, nodes, true
        ).value; // shape: [n][eta]

        // Build coef_mat (size x size)
        std::vector<std::vector<T>> coef_mat(size, std::vector<T>(size, T(0)));
        for (int li = 0; li < size; ++li) {
            for (int ni = 0; ni < size; ++ni) {
                T sum = T(0);
                for (size_t k = 0; k < nodes.size(); ++k) {
                    T v1 = smn_ext[ni][k];
                    T v2 = smn_int[li][k];
                    sum += v1 * v2 * weights[k];
                }
                coef_mat[li][ni] = sum;
            }
        }
        coef_list[m] = coef_mat;
    }
    return coef_list;
}

// -----------------------------------------------------------
// KERNEL BOUNDARY CONDITIONS FOR FLUID-FILLED SCATTERERS
// -----------------------------------------------------------
template<typename T>
struct KernelResult {
    std::vector<std::vector<std::complex<T>>> K1_kernel;
    std::vector<std::vector<std::complex<T>>> K3_kernel;
};

template<typename T>
KernelResult<T> fluid_boundary_kernels(
    int m_max,
    int n_max,
    T chi_sw,
    T chi_body,
    T xi,
    T density_body,
    T density_sw,
    const std::vector<T>& nodes,
    const std::vector<T>& weights,
    const std::vector<std::vector<T>>& smn_matrix
) {
    // Compute expansion matrix
    auto expansion_matrix = alpha_integration_coefficient<T>(
        m_max, n_max, chi_sw, chi_body, nodes, weights
    );

    // Radial function matrices
    auto radial_external_kind1 = radial_external_incident_matrix<T>(
        m_max, n_max, chi_sw, xi, smn_matrix
    );
    auto radial_external_kind3 = radial_external_scattering_matrix<T>(
        m_max, n_max, chi_sw, xi, smn_matrix
    );
    auto radial_internal_kind1 = radial_internal_incident_matrix<T>(
        m_max, n_max, chi_body, xi
    );

    // Coupling matrices
    auto E1_coupling = boundary_coupling_incident_matrix<T>(
        m_max, n_max, density_body, density_sw,
        radial_external_kind1.value, radial_external_kind1.derivative,
        radial_internal_kind1.value, radial_internal_kind1.derivative,
        smn_matrix
    );
    auto E3_coupling = boundary_coupling_scattering_matrix<T>(
        m_max, n_max, density_body, density_sw,
        radial_external_kind3.value, radial_external_kind3.derivative,
        radial_internal_kind1.value, radial_internal_kind1.derivative,
        smn_matrix
    );

    // Kernel matrices
    auto K1_kernel = kernel_incident_matrix<T>(
        m_max, n_max, smn_matrix, expansion_matrix, E1_coupling
    );
    auto K3_kernel = kernel_scattering_matrix<T>(
        m_max, n_max, smn_matrix, expansion_matrix, E3_coupling
    );
    KernelResult<T> out;
    out.K1_kernel = K1_kernel;
    out.K3_kernel = K3_kernel;
    return out;
}

// -----------------------------------------------------------
// EXPANSION COEFFICIENT MATRIX SOLVERS
// -----------------------------------------------------------
// Divide-and-conquer SVD
template<typename T>
std::vector<std::vector<std::complex<T>>> solve_fluid_Amn_divide_and_conquer(
    const std::vector<std::vector<std::complex<T>>>& K1_kernel,
    const std::vector<std::vector<std::complex<T>>>& K3_kernel
) {
    int m_max = K1_kernel.size() - 1;
    std::vector<std::vector<std::complex<T>>> result(m_max + 1);

    for (int m = 0; m <= m_max; ++m) {
        int size = std::sqrt(K1_kernel[m].size());
        // Convert to double precision for SVD
        arma::Col<std::complex<double>> b(size);
        arma::Mat<std::complex<double>> K3_arma(size, size);
        for (int i = 0; i < size; ++i) {
            std::complex<double> row_sum(0, 0);
            for (int j = 0; j < size; ++j) {
                auto val1 = K1_kernel[m][i * size + j];
                auto val2 = K3_kernel[m][i * size + j];
                row_sum += std::complex<double>(
                    static_cast<double>(val1.real()),
                    static_cast<double>(val1.imag())
                );
                K3_arma(i, j) = std::complex<double>(
                    static_cast<double>(val2.real()),
                    static_cast<double>(val2.imag())
                );
            }
            b(i) = -row_sum;
        }
        arma::Mat<std::complex<double>> U, V;
        arma::Col<double> s;
        bool svd_ok = arma::svd(U, s, V, K3_arma);
        if (!svd_ok) throw std::runtime_error("SVD failed");
        double tol = std::max(size, size) * s.max() * std::numeric_limits<double>::epsilon();
        arma::Col<double> d_inv(s.n_elem);
        for (arma::uword i = 0; i < s.n_elem; ++i)
            d_inv(i) = (s(i) > tol) ? (1.0 / s(i)) : 0.0;
        arma::Mat<std::complex<double>> diag_dinv = arma::diagmat(
            arma::conv_to<arma::Col<std::complex<double>>>::from(d_inv)
        );
        arma::Mat<std::complex<double>> K3_pinv = V * diag_dinv * U.t();
        arma::Col<std::complex<double>> A = K3_pinv * b;
        
        // Convert result to T (double or quad)
        std::vector<std::complex<T>> A_vec(A.n_elem);
        for (arma::uword i = 0; i < A.n_elem; ++i)
            A_vec[i] = std::complex<T>(
                static_cast<T>(A(i).real()),
                static_cast<T>(A(i).imag())
            );
        result[m] = A_vec;
    }
    return result;
}

// Two-sided Jacobi SVD
// template<typename T>
// std::vector<std::vector<std::complex<T>>> solve_fluid_Amn_jacobi(
//     const std::vector<std::vector<std::complex<T>>>& K1_kernel,
//     const std::vector<std::vector<std::complex<T>>>& K3_kernel
// ) {
//     int m_max = static_cast<int>(K1_kernel.size()) - 1;
//     std::vector<std::vector<std::complex<T>>> result(m_max + 1);

//     for (int m = 0; m <= m_max; ++m) {
//         int size = static_cast<int>(std::sqrt(K1_kernel[m].size()));

//         // Build b
//         std::vector<std::complex<T>> b(size, std::complex<T>(0));
//         for (int i = 0; i < size; ++i) {
//             std::complex<T> row_sum(0);
//             for (int j = 0; j < size; ++j)
//                 row_sum += K1_kernel[m][i * size + j];
//             b[i] = -row_sum;
//         }

//         // Solve using Jacobi SVD
//         result[m] = two_sided_jacobi(K3_kernel[m], b, size);
//     }

//     return result;
// }

// ============================================================================
// EXPANSION COEFFICIENT MATRICES FOR SCATTERED WAVES
// ----------------------------------------------------------------------------
// The expansion coefficient matrices for scattered waves of different 
// boundary conditions from Furusawa (1988):
// - Pressure release (soft): Eq. (3) 
// - Fixed rigid: Eq. (3)
// - Fluid-filled: Eq. (4)
// - Simplified fluid-filled: Eq. (5)
// ============================================================================
// -----------------------------------------------------------
// PRESSURE RELEASE
// -----------------------------------------------------------
template<typename T>
std::vector<std::vector<std::complex<T>>> pressure_release_Amn_expansion_matrix(
    int m_max,
    int n_max,
    T chi_sw,
    T xi
) {

    // Create a dummy Smn-based matrix that comprises all 1's
    std::vector<std::vector<T>> smn_matrix(m_max + 1, std::vector<T>(n_max + 1, 1.0));

    auto R_incident = radial_external_incident_matrix<T>(
        m_max, n_max, chi_sw, xi, smn_matrix
    );
    auto R_scattering = radial_external_scattering_matrix<T>(
        m_max, n_max, chi_sw, xi, smn_matrix
    );

    std::vector<std::vector<std::complex<T>>> amn(
        m_max + 1, std::vector<std::complex<T>>(n_max + 1, std::complex<T>(0, 0))
    );
    for (int m = 0; m <= m_max; ++m) {
        for (int n = m; n <= n_max; ++n) {
            auto v = R_incident.value[m][n];
            auto s = R_scattering.value[m][n];
            if (s != std::complex<T>(0, 0)) {
                amn[m][n] = -v / s;
            }
        }
    }
    
    return amn;
}

// -----------------------------------------------------------
// FIXED RIGID
// -----------------------------------------------------------
template<typename T>
std::vector<std::vector<std::complex<T>>> fixed_rigid_Amn_expansion_matrix(
    int m_max,
    int n_max,
    T chi_sw,
    T xi
) {

    // Create a dummy Smn-based matrix that comprises all 1's
    std::vector<std::vector<T>> smn_matrix(m_max + 1, std::vector<T>(n_max + 1, 1.0));

    auto R_incident = radial_external_incident_matrix<T>(
        m_max, n_max, chi_sw, xi, smn_matrix
    );
    auto R_scattering = radial_external_scattering_matrix<T>(
        m_max, n_max, chi_sw, xi, smn_matrix
    );

    std::vector<std::vector<std::complex<T>>> amn(
        m_max + 1, std::vector<std::complex<T>>(n_max + 1, std::complex<T>(0, 0))
    );
    for (int m = 0; m <= m_max; ++m) {
        for (int n = m; n <= n_max; ++n) {
            auto v = R_incident.derivative[m][n];
            auto s = R_scattering.derivative[m][n];
            if (s != std::complex<T>(0, 0)) {
                amn[m][n] = -v / s;
            }
        }
    }
    
    return amn;
}

// -----------------------------------------------------------
// FLUID-FILLED
// -----------------------------------------------------------
template<typename T>
std::vector<std::vector<std::complex<T>>> fluid_Amn_expansion_matrix(
    int m_max,
    int n_max,
    T chi_sw,
    T chi_body,
    T xi,
    T theta_body,
    T density_body,
    T density_sw,
    const std::vector<T>& nodes,
    const std::vector<T>& weights
) {
    // Compute S_mn matrix at incident angle
    auto smn_matrix = compute_smn_matrix<T>(
        m_max, n_max, chi_sw, preccos(theta_body), true
    );

    // Compute kernel matrices
    auto kernels = fluid_boundary_kernels<T>(
        m_max, n_max, chi_sw, chi_body, xi,
        density_body, density_sw, nodes, weights, smn_matrix
    );

    // Solve for expansion coefficients
    // std::vector<std::vector<std::complex<T>>> Amn;
    // if (solver_method == "divide_and_conquer") {
    auto Amn = solve_fluid_Amn_divide_and_conquer<T>(
        kernels.K1_kernel, kernels.K3_kernel
    );
    //} else if (solver_method == "two_sided_jacobi") {
    //     auto Amn = solve_fluid_Amn_jacobi<T>(
    //         kernels.K1_kernel, kernels.K3_kernel
    //     );
    // }

    // Convert to rectangular (m_max+1) x (n_max+1) matrix, NA for n < m
    std::vector<std::vector<std::complex<T>>> Amn_mat(
        m_max + 1, 
        std::vector<std::complex<T>>(
            n_max + 1, std::complex<T>(
                std::numeric_limits<T>::quiet_NaN(), std::numeric_limits<T>::quiet_NaN()
            )
        )
        );
    for (int m = 0; m <= m_max; ++m) {
        int len = Amn[m].size();
        for (int i = 0; i < len; ++i) {
            int n = m + i;
            if (n <= n_max) {
                auto val = Amn[m][i];
                // If val is nan, set to NA
                if (is_na_real(val.real()) || is_na_real(val.imag())) {
                    Amn_mat[m][n] = std::complex<T>(
                        std::numeric_limits<T>::quiet_NaN(), std::numeric_limits<T>::quiet_NaN()
                    );
                } else {
                    Amn_mat[m][n] = val;
                }
            }
        }
    }
    return Amn_mat;
}

// Simplified fluid-filled
template<typename T>
std::vector<std::vector<std::complex<T>>> simplified_fluid_Amn_expansion_matrix(
    int m_max,
    int n_max,
    T chi_sw,
    T chi_body,
    T xi,
    T density_body,
    T density_sw
) {
    // Create a dummy Smn-based matrix that comprises all 1's
    std::vector<std::vector<T>> smn_matrix(m_max + 1, std::vector<T>(n_max + 1, 1.0));

    // Radial function matrices
    auto radial_external_kind1 = radial_external_incident_matrix<T>(
        m_max, n_max, chi_sw, xi, smn_matrix
    );
    auto radial_external_kind3 = radial_external_scattering_matrix<T>(
        m_max, n_max, chi_sw, xi, smn_matrix
    );
    auto radial_internal_kind1 = radial_internal_incident_matrix<T>(
        m_max, n_max, chi_body, xi
    );

    // Coupling matrices
    auto E1_coupling = simplified_boundary_coupling_incident_matrix<T>(
        m_max, n_max, density_body, density_sw,
        radial_external_kind1.value, radial_external_kind1.derivative,
        radial_internal_kind1.value, radial_internal_kind1.derivative
    );
    auto E3_coupling = simplified_boundary_coupling_scattering_matrix<T>(
        m_max, n_max, density_body, density_sw,
        radial_external_kind3.value, radial_external_kind3.derivative,
        radial_internal_kind1.value, radial_internal_kind1.derivative
    );

    // Compute simplified Amn expansion matrix directly from diagonal of E1/E3.
    // The simplified method (Furusawa 1988, Eq. 5) is A_{mn} = -E1[m][n] / E3[m][n]
    // for each (m, n) independently. E1_coupling and E3_coupling are indexed as
    // [m][n] with n from 0 to n_max -- access them directly, no size*size needed.
    std::vector<std::vector<std::complex<T>>> Amn_mat(
        m_max + 1, 
        std::vector<std::complex<T>>(
            n_max + 1, std::complex<T>(
                std::numeric_limits<T>::quiet_NaN(), std::numeric_limits<T>::quiet_NaN()
            )
        )
    );
    for (int m = 0; m <= m_max; ++m) {
        for (int n = m; n <= n_max; ++n) {
            std::complex<T> e1 = E1_coupling[m][n];
            std::complex<T> e3 = E3_coupling[m][n];
            if (is_na_real(e1.real()) || is_na_real(e1.imag()) ||
                e1 == std::complex<T>(0, 0) ||
                is_na_real(e3.real()) || is_na_real(e3.imag()) ||
                e3 == std::complex<T>(0, 0)) {
                Amn_mat[m][n] = std::complex<T>(
                    std::numeric_limits<T>::quiet_NaN(), std::numeric_limits<T>::quiet_NaN()
                );
            } else {
                Amn_mat[m][n] = -e1 / e3;
            }
        }
    }
    return Amn_mat;
}

// ============================================================================
// LINEAR SCATTERING COEFFICIENT
// ----------------------------------------------------------------------------
// Linear scattering coefficient, f_inf(theta, phi|theta', phi'), based on 
// Eq. (2) from Furusawa (1988). The scattering coefficient is summed using 
// the Kaham summation algorithm.
// ============================================================================
template<typename T>
std::complex<T> compute_fbs(
    int m_max,
    int n_max,
    T chi_sw,
    T theta_body,
    T theta_scatter,
    const std::vector<T>& azimuth,
    const std::vector<std::vector<std::complex<T>>>& Amn
) {
    if (azimuth.size() != static_cast<size_t>(m_max + 1)){
        throw std::invalid_argument("azimuth must have length m_max + 1");
    }
    std::vector<T> nu(m_max + 1);
    for (int m = 0; m <= m_max; ++m)
        nu[m] = (m == 0) ? 1.0 : 2.0;

    auto smn_body_matrix = compute_smn_matrix<T>(
        m_max, n_max, chi_sw, preccos(theta_body), true
    );
    auto smn_scatter_matrix = compute_smn_matrix<T>(
        m_max, n_max, chi_sw, preccos(theta_scatter), true
    );

    std::complex<T> fbs_sum(0, 0), c(0, 0);

    // Loop over all (m, n) with n >= m, matching R logic
    for (int m = 0; m <= m_max; ++m) {
        for (int n = m; n <= n_max; ++n) {
            std::complex<T> Amn_val(0, 0);
            if (m < static_cast<int>(Amn.size()) && n < static_cast<int>(Amn[m].size()))
                Amn_val = Amn[m][n];

            // Compute Smn for incident and scattered field for (m, n)
            T smn_body = smn_body_matrix[m][n];
            T smn_scatter = smn_scatter_matrix[m][n];

            // Skip if any are NaN
            if (is_na_real(smn_body) || is_na_real(smn_scatter) ||
                is_na_real(Amn_val.real()) || is_na_real(Amn_val.imag()))
                continue;
            
            // Kahan summation
            std::complex<T> term = nu[m] * smn_body * smn_scatter * Amn_val * azimuth[m];
            std::complex<T> y = term - c;
            std::complex<T> t = fbs_sum + y;
            c = (t - fbs_sum) - y;
            fbs_sum = t;
        }
    }
    return fbs_sum;
}

// ============================================================================
// LINEAR SCATTERING COEFFICIENT FOR A PROLATE SPHEROID
// ----------------------------------------------------------------------------
// Linear scattering coefficient, f_inf(theta, phi|theta', phi'), based on 
// Eq. (2) from Furusawa (1988). 
// ============================================================================
template<typename T>
std::complex<T> psms_fbs(
    int m_max,
    int n_max,
    T chi_sw,
    T chi_body,
    T xi,
    T theta_body,
    T theta_scatter,
    T phi_body,
    T phi_scatter,
    T density_body,
    T density_sw,
    const std::vector<T>& nodes,
    const std::vector<T>& weights,
    const std::string& Amn_method 
) {
    // Determine appropriate Amn expansion matrix computation
    std::vector<std::vector<std::complex<T>>> Amn;

    // Fluid-filled
    if (Amn_method == "Amn_fluid") {
        Amn = fluid_Amn_expansion_matrix<T>(
            m_max, n_max, chi_sw, chi_body, xi,
            theta_body, density_body, density_sw, nodes, weights           
        );
    // Simplified fluid-filled
    } else if (Amn_method == "Amn_fluid_simplify") {
        Amn = simplified_fluid_Amn_expansion_matrix<T>(
            m_max, n_max, chi_sw, chi_body, xi, density_body, density_sw
        );
    // Fixed rigid
    } else if (Amn_method == "Amn_fixed_rigid") {
        Amn = fixed_rigid_Amn_expansion_matrix<T>(
            m_max, n_max, chi_sw, xi
        );
    // Pressure-release
    } else {
        Amn = pressure_release_Amn_expansion_matrix<T>(
            m_max, n_max, chi_sw, xi
        );
    }

    // Compute azimuth angle factors
    auto azimuth = compute_azimuth<T>(m_max, phi_body, phi_scatter);

    // Compute linear scattering coefficient
    return compute_fbs<T>(
        m_max, n_max, chi_sw, theta_body, theta_scatter, azimuth, Amn
    );
}

// ============================================================================
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
    std::string Amn_method = "Amn_fluid"
) {
    // Extract acoustic parameters
    Rcpp::NumericVector chi_sw = acoustics["chi_sw"];
    Rcpp::NumericVector chi_body = acoustics["chi_body"];
    Rcpp::IntegerVector m_max = acoustics["m_max"];
    Rcpp::IntegerVector n_max = acoustics["n_max"];

    // Extract body parameters
    double xi = body["xi"];
    double theta_body = body["theta_body"];
    double theta_scatter = body["theta_scatter"];
    double phi_body = body["phi_body"];
    double phi_scatter = body["phi_scatter"];
    double density_body = body["density"];

    // Extract medium parameters
    double density_sw = medium["density"];
    
    // Extract quadrature points
    Rcpp::NumericVector nodes = integration_pts["nodes"];
    Rcpp::NumericVector weights = integration_pts["weights"];

    int n_freq = acoustics.nrows();
    Rcpp::ComplexVector f_bs(n_freq);

    if (precision == "quad") {
#ifdef __GNUC__
        std::vector<__float128> nodes_q(nodes.size()), weights_q(weights.size());
        for (int i = 0; i < nodes.size(); ++i) nodes_q[i] = static_cast<__float128>(nodes[i]);
        for (int i = 0; i < weights.size(); ++i) weights_q[i] = static_cast<__float128>(weights[i]);

        for (int f = 0; f < n_freq; ++f) {
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
                Amn_method
            );
            f_bs[f].r = static_cast<double>(fbs.real());
            f_bs[f].i = static_cast<double>(fbs.imag());
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
                Amn_method
            );
            f_bs[f].r = fbs.real();
            f_bs[f].i = fbs.imag();
        }
    }
    return f_bs;
}