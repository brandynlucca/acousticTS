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

// ============================================================================
// HELPER FUNCTIONS
// ============================================================================
// Compute azimuthal factors cos(m(phi_body - phi_scatter)) for the general
// bistatic scattering sum.
template<typename T>
std::vector<T> compute_azimuth(int m_max, T phi_body, T phi_scatter) {
    std::vector<T> azimuth(m_max + 1);
    T dphi = phi_body - phi_scatter;
    for (int m = 0; m <= m_max; ++m)
        azimuth[m] = preccos(m * dphi);
    return azimuth;
}

template<typename T>
inline std::complex<T> imaginary_unit_power(int n) {
    switch (((n % 4) + 4) % 4) {
        case 0:
            return std::complex<T>(T(1), T(0));
        case 1:
            return std::complex<T>(T(0), T(1));
        case 2:
            return std::complex<T>(T(-1), T(0));
        default:
            return std::complex<T>(T(0), T(-1));
    }
}

// Extract one angular value from a single-m profcn call and apply the stored
// base-10 exponent immediately.
template<typename T>
inline T extract_angular_value_from_batch(
    const ProfcnResult<T>& result,
    int idx,
    int lnum = 0,
    int arg_idx = 0
) {
    int offset = idx;
    if (lnum > 0) {
        offset += arg_idx * lnum;
    }
    T scale = T(1);
    if (result.is1e.size() > static_cast<size_t>(offset) && !is_na_int(result.is1e[offset])) {
        scale = pow10_typed<T>(result.is1e[offset]);
    }
    if (result.s1c.size() > static_cast<size_t>(offset) && !is_na_real(result.s1c[offset])) {
        return result.s1c[offset] * scale;
    }
    return std::numeric_limits<T>::quiet_NaN();
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

template<typename T>
inline void scale_profcn_component(std::vector<T>& values, std::vector<int>& exponents) {
    size_t n = std::min(values.size(), exponents.size());
    for (size_t i = 0; i < n; ++i) {
        if (is_na_int(exponents[i])) {
            continue;
        }
        if (!is_na_real(values[i]) && exponents[i] != 0) {
            values[i] *= pow10_typed<T>(exponents[i]);
        }
        exponents[i] = 0;
    }
}

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
    // profcn writes into fixed-size arrays, so the wrapper allocates the full
    // expected Fortran layout before making the interface call.
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
    __float128 c,
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
    // The quad wrapper mirrors the double wrapper, but keeps the interface
    // scalars in __float128 so the Fortran call never falls back to double.
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
    __float128 c_q = c;
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

template<typename T>
struct ProfcnBatchResult {
    int lnum = 0;
    int narg = 0;
    int m_count = 0;
    std::vector<T> r1c, r1dc, r2c, r2dc;
    std::vector<int> ir1e, ir1de, ir2e, ir2de, naccr;
    std::vector<T> s1c, s1dc;
    std::vector<int> is1e, is1de, naccs;
};

// Batched profcn wrappers amortize the expensive spheroidal setup over a block
// of consecutive m values instead of repeating the full setup for each order.
inline ProfcnBatchResult<double> cprofcn_batch_double(
    double c,
    int m_start,
    int m_count,
    int lnum,
    const std::vector<double>& arg,
    int ioprad,
    int iopnorm,
    int iopang,
    double x1
) {
    int narg = std::max<int>(arg.size(), 1);
    ProfcnBatchResult<double> result;
    result.lnum = lnum;
    result.narg = narg;
    result.m_count = m_count;
    size_t radial_size = static_cast<size_t>(lnum) * static_cast<size_t>(m_count);
    size_t angular_size = radial_size * static_cast<size_t>(narg);

    result.r1c.resize(radial_size, 0.0);
    result.r1dc.resize(radial_size, 0.0);
    result.r2c.resize(radial_size, 0.0);
    result.r2dc.resize(radial_size, 0.0);
    result.ir1e.resize(radial_size, 0);
    result.ir1de.resize(radial_size, 0);
    result.ir2e.resize(radial_size, 0);
    result.ir2de.resize(radial_size, 0);
    result.naccr.resize(radial_size, 0);
    result.s1c.resize(angular_size, 0.0);
    result.s1dc.resize(angular_size, 0.0);
    result.is1e.resize(angular_size, 0);
    result.is1de.resize(angular_size, 0);
    result.naccs.resize(angular_size, 0);

    std::vector<double> arg_c = arg;
    profcn_cpp_interface_batch(
        &c, &m_start, &m_count, &lnum, &ioprad, &x1, &iopang, &iopnorm, &narg,
        arg_c.data(), result.r1c.data(), result.ir1e.data(), result.r1dc.data(),
        result.ir1de.data(), result.r2c.data(), result.ir2e.data(), result.r2dc.data(),
        result.ir2de.data(), result.naccr.data(), result.s1c.data(), result.is1e.data(),
        result.s1dc.data(), result.is1de.data(), result.naccs.data()
    );

    scale_profcn_component(result.r1c, result.ir1e);
    scale_profcn_component(result.r1dc, result.ir1de);
    scale_profcn_component(result.r2c, result.ir2e);
    scale_profcn_component(result.r2dc, result.ir2de);
    scale_profcn_component(result.s1c, result.is1e);
    scale_profcn_component(result.s1dc, result.is1de);
    return result;
}

#ifdef __GNUC__
inline ProfcnBatchResult<__float128> cprofcn_batch_quad(
    __float128 c,
    int m_start,
    int m_count,
    int lnum,
    const std::vector<__float128>& arg,
    int ioprad,
    int iopnorm,
    int iopang,
    __float128 x1
) {
    int narg = std::max<int>(arg.size(), 1);
    ProfcnBatchResult<__float128> result;
    result.lnum = lnum;
    result.narg = narg;
    result.m_count = m_count;
    size_t radial_size = static_cast<size_t>(lnum) * static_cast<size_t>(m_count);
    size_t angular_size = radial_size * static_cast<size_t>(narg);

    result.r1c.resize(radial_size, 0.0);
    result.r1dc.resize(radial_size, 0.0);
    result.r2c.resize(radial_size, 0.0);
    result.r2dc.resize(radial_size, 0.0);
    result.ir1e.resize(radial_size, 0);
    result.ir1de.resize(radial_size, 0);
    result.ir2e.resize(radial_size, 0);
    result.ir2de.resize(radial_size, 0);
    result.naccr.resize(radial_size, 0);
    result.s1c.resize(angular_size, 0.0);
    result.s1dc.resize(angular_size, 0.0);
    result.is1e.resize(angular_size, 0);
    result.is1de.resize(angular_size, 0);
    result.naccs.resize(angular_size, 0);

    std::vector<__float128> arg_c = arg;
    __float128 c_q = c;
    __float128 x1_q = x1;
    int m_start_q = m_start;
    int m_count_q = m_count;
    int lnum_q = lnum;
    int ioprad_q = ioprad;
    int iopnorm_q = iopnorm;
    int iopang_q = iopang;
    int narg_q = narg;

    profcn_cpp_interface_batch_quad(
        &c_q, &m_start_q, &m_count_q, &lnum_q, &ioprad_q, &x1_q, &iopang_q, &iopnorm_q,
        &narg_q, arg_c.data(), result.r1c.data(), result.ir1e.data(), result.r1dc.data(),
        result.ir1de.data(), result.r2c.data(), result.ir2e.data(), result.r2dc.data(),
        result.ir2de.data(), result.naccr.data(), result.s1c.data(), result.is1e.data(),
        result.s1dc.data(), result.is1de.data(), result.naccs.data()
    );

    scale_profcn_component(result.r1c, result.ir1e);
    scale_profcn_component(result.r1dc, result.ir1de);
    scale_profcn_component(result.r2c, result.ir2e);
    scale_profcn_component(result.r2dc, result.ir2de);
    scale_profcn_component(result.s1c, result.is1e);
    scale_profcn_component(result.s1dc, result.is1de);
    return result;
}
#endif

template<typename T>
inline size_t profcn_batch_radial_offset(
    int l_idx,
    int m_idx,
    int lnum
) {
    return static_cast<size_t>(l_idx) + static_cast<size_t>(lnum) * static_cast<size_t>(m_idx);
}

template<typename T>
inline size_t profcn_batch_angular_offset(
    int l_idx,
    int arg_idx,
    int m_idx,
    int lnum,
    int narg
) {
    return static_cast<size_t>(l_idx) +
        static_cast<size_t>(lnum) * (
            static_cast<size_t>(arg_idx) +
            static_cast<size_t>(narg) * static_cast<size_t>(m_idx)
        );
}

// Extract one angular value from a batched m-block result. The mantissa/exponent
// scaling has already been folded back into s1c/s1dc by this stage.
template<typename T>
inline T extract_angular_value_from_mblock(
    const ProfcnBatchResult<T>& result,
    int m_idx,
    int l_idx,
    int arg_idx = 0
) {
    size_t offset = profcn_batch_angular_offset<T>(
        l_idx, arg_idx, m_idx, result.lnum, result.narg
    );
    if (offset < result.s1c.size() && !is_na_real(result.s1c[offset])) {
        return result.s1c[offset];
    }
    return std::numeric_limits<T>::quiet_NaN();
}

template<typename T>
inline int psms_profcn_mblock_chunk_size() {
    return 8;
}

#ifdef __GNUC__
template<>
inline int psms_profcn_mblock_chunk_size<__float128>() {
    return 4;
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

template<typename T>
ProfcnBatchResult<T> cprofcn_mblock(
    T c,
    int m_start,
    int m_count,
    int lnum,
    const std::vector<T>& arg,
    int ioprad,
    int iopnorm,
    int iopang,
    T x1
);

template<typename T>
inline void normalize_profcn_result(ProfcnResult<T>& result) {
    scale_profcn_component(result.r1c, result.ir1e);
    scale_profcn_component(result.r1dc, result.ir1de);
    scale_profcn_component(result.r2c, result.ir2e);
    scale_profcn_component(result.r2dc, result.ir2de);
    scale_profcn_component(result.s1c, result.is1e);
    scale_profcn_component(result.s1dc, result.is1de);
}

// Dispatch to the double- or quad-precision batched profcn backend while
// preserving a common templated interface for the PSMS code.
template<>
inline ProfcnBatchResult<double> cprofcn_mblock<double>(
    double c,
    int m_start,
    int m_count,
    int lnum,
    const std::vector<double>& arg,
    int ioprad,
    int iopnorm,
    int iopang,
    double x1
) {
    return cprofcn_batch_double(c, m_start, m_count, lnum, arg, ioprad, iopnorm, iopang, x1);
}

#ifdef __GNUC__
template<>
inline ProfcnBatchResult<__float128> cprofcn_mblock<__float128>(
    __float128 c,
    int m_start,
    int m_count,
    int lnum,
    const std::vector<__float128>& arg,
    int ioprad,
    int iopnorm,
    int iopang,
    __float128 x1
) {
    return cprofcn_batch_quad(c, m_start, m_count, lnum, arg, ioprad, iopnorm, iopang, x1);
}
#endif

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
    normalize_profcn_result(out);

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
    normalize_profcn_result(out);

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
    // Angular-only request: no radial values are needed, so x1 can be fixed at
    // zero and ioprad = 0.
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

    // Validate the requested modal region before entering the expensive
    // spheroidal-function backend.
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

        // Group by unique m so each profcn call can be reused for all degrees
        // that share the same order.
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
// Rmn(3) and Rmn(4) are assembled from the first- and second-kind radial
// functions. This helper recovers those real-valued ingredients before they are
// combined into outgoing/incoming-wave conventions.
template<typename T>
RmnResult<T> Rmn_higher_order(int m, int n, int lnum, T c, T x1) {
    int iopnorm = 0, iopang = 0;
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
        // Kinds 1 and 2 are read directly from profcn without any extra
        // complex recombination.
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

    // Higher-order conventions use outgoing/incoming combinations built from
    // first- and second-kind radial functions.
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
// Translate packed profcn radial arrays into the real/imaginary pieces used by
// the PSMS algebra, applying the stored exponents and basic sanity checks.
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

template<typename T>
RadialValue<T> extract_radial_from_mblock(
    const ProfcnBatchResult<T>& result,
    int m_idx,
    int l_idx,
    int kind
) {
    RadialValue<T> out;
    size_t offset = profcn_batch_radial_offset<T>(l_idx, m_idx, result.lnum);

    // Kinds 1 and 2 are purely real. Kinds 3 and 4 use the first- and
    // second-kind outputs to form complex outgoing/incoming combinations.
    auto assign_kind12 = [&](const std::vector<T>& rc, const std::vector<T>& rdc) {
        out.val_real = (offset < rc.size() && !is_na_real(rc[offset]))
            ? rc[offset]
            : std::numeric_limits<T>::quiet_NaN();
        out.der_real = (offset < rdc.size() && !is_na_real(rdc[offset]))
            ? rdc[offset]
            : std::numeric_limits<T>::quiet_NaN();
        out.val_imag = T(0);
        out.der_imag = T(0);
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
    };

    if (kind == 1) {
        assign_kind12(result.r1c, result.r1dc);
        return out;
    }
    if (kind == 2) {
        assign_kind12(result.r2c, result.r2dc);
        return out;
    }

    out.val_real = (offset < result.r1c.size() && !is_na_real(result.r1c[offset]))
        ? result.r1c[offset]
        : std::numeric_limits<T>::quiet_NaN();
    out.der_real = (offset < result.r1dc.size() && !is_na_real(result.r1dc[offset]))
        ? result.r1dc[offset]
        : std::numeric_limits<T>::quiet_NaN();
    out.val_imag = (offset < result.r2c.size() && !is_na_real(result.r2c[offset]))
        ? result.r2c[offset]
        : std::numeric_limits<T>::quiet_NaN();
    out.der_imag = (offset < result.r2dc.size() && !is_na_real(result.r2dc[offset]))
        ? result.r2dc[offset]
        : std::numeric_limits<T>::quiet_NaN();

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
        bool batch_complex = is_complex && x1_r > T(0);
        int radial_kind = batch_complex ? 2 : (is_complex ? 1 : kind);
        auto result1 = cprofcn<T>(c, m_val, lnum, {T(1)}, radial_kind, 0, 0, x1_r);
        ProfcnResult<T> result2;
        if (batch_complex) {
            result2 = result1;
        } else if (is_complex) {
            result2 = cprofcn<T>(c, m_val, lnum, {T(1)}, 2, 0, 0, x1_r);
        }
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
            bool batch_complex = is_complex && x1_r > T(0);
            int radial_kind = batch_complex ? 2 : (is_complex ? 1 : kind);
            auto result1 = cprofcn<T>(c, m_val, lnum, {T(1)}, radial_kind, 0, 0, x1_r);
            ProfcnResult<T> result2;
            if (batch_complex) {
                result2 = result1;
            } else if (is_complex) {
                result2 = cprofcn<T>(c, m_val, lnum, {T(1)}, 2, 0, 0, x1_r);
            }
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
        bool batch_complex = is_complex && x1_r > T(0);
        int radial_kind = batch_complex ? 2 : (is_complex ? 1 : kind);
        auto result1 = cprofcn<T>(c, m_val, lnum, {T(1)}, radial_kind, 0, 0, x1_r);
        ProfcnResult<T> result2;
        if (batch_complex) {
            result2 = result1;
        } else if (is_complex) {
            result2 = cprofcn<T>(c, m_val, lnum, {T(1)}, 2, 0, 0, x1_r);
        }
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
struct ExternalRadialResult {
    RadialMatrixResult<T> incident;
    RadialComplexMatrixResult<T> scattering;
};

template<typename T>
struct ExternalBoundaryResult {
    std::vector<std::vector<T>> smn;
    ExternalRadialResult<T> radial;
};

template<typename T>
ExternalRadialResult<T> radial_external_matrices(
    int m_max,
    int n_max,
    T chi_sw,
    T xi,
    const std::vector<std::vector<T>>& smn_matrix
);

template<typename T>
ExternalBoundaryResult<T> external_boundary_matrices(
    int m_max,
    int n_max,
    T chi_sw,
    T xi,
    T eta_body,
    bool normalize = true
);

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
        T x1_r = xi - T(1);
        auto result1 = cprofcn<T>(chi_body, m, size, {T(1)}, 1, 0, 0, x1_r);
        for (int i = 0; i < size; ++i) {
            int n = m + i;
            auto rv = extract_radial_from_batch<T>(
                result1, result1, i, 1, m, n, chi_body, xi
            );
            out.value[m][n] = rv.val_real;
            out.derivative[m][n] = rv.der_real;
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
            T factor = (R_ml_1_d != T(0))
                ? (density_body * R_ml_1_val) / (density_sw * R_ml_1_d)
                : T(0);
            for (int ni = 0; ni < size; ++ni) {
                int n = m + ni;
                T snm_n = smn_matrix[m][n];
                if (std::isnan(static_cast<double>(snm_n)) || snm_n == T(0)) continue;
                T Rmn_i1_val = e1_val_mat[m][n];
                T Rmn_i1_d   = e1_der_mat[m][n];
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
            T factor = (R_ml_1_d != T(0))
                ? (density_body * R_ml_1_val) / (density_sw * R_ml_1_d)
                : T(0);
            for (int ni = 0; ni < size; ++ni) {
                int n = m + ni;
                T snm_n = smn_matrix[m][n];
                if (std::isnan(static_cast<double>(snm_n)) || snm_n == T(0)) continue;
                std::complex<T> Rmn_i3_val = e3_val_mat[m][n];
                std::complex<T> Rmn_i3_d   = e3_der_mat[m][n];
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
std::vector<std::vector<std::complex<T>>> kernel_incident_rhs(
    int m_max,
    int n_max,
    const std::vector<std::vector<T>>& smn_matrix,
    const std::vector<std::vector<std::vector<T>>>& expansion_matrix,
    const std::vector<std::vector<std::complex<T>>>& E1_coupling
) {
    std::vector<std::vector<std::complex<T>>> result(m_max + 1);
    for (int m = 0; m <= m_max; ++m) {
        int size = n_max - m + 1;
        std::vector<std::complex<T>> rhs(size, std::complex<T>(0, 0));
        std::vector<std::complex<T>> jn_smn(size, std::complex<T>(0, 0));
        for (int ni = 0; ni < size; ++ni) {
            int n = m + ni;
            T snm_val = smn_matrix[m][n];
            if (!std::isnan(static_cast<double>(snm_val)) && snm_val != T(0)) {
                jn_smn[ni] = imaginary_unit_power<T>(n) * snm_val;
            }
        }
        for (int li = 0; li < size; ++li) {
            std::complex<T> row_sum(0, 0);
            for (int ni = 0; ni < size; ++ni) {
                std::complex<T> phase_smn = jn_smn[ni];
                if (phase_smn == std::complex<T>(0, 0)) {
                    continue;
                }
                T alpha = expansion_matrix[m][li][ni];
                std::complex<T> E1_nl = E1_coupling[m][li * size + ni];
                row_sum += phase_smn * alpha * E1_nl;
            }
            rhs[li] = -row_sum;
        }
        result[m] = rhs;
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
        std::vector<std::complex<T>> jn_smn(size, std::complex<T>(0, 0));
        for (int ni = 0; ni < size; ++ni) {
            int n = m + ni;
            T snm_val = smn_matrix[m][n];
            if (!std::isnan(static_cast<double>(snm_val)) && snm_val != T(0)) {
                jn_smn[ni] = imaginary_unit_power<T>(n) * snm_val;
            }
        }
        for (int li = 0; li < size; ++li) {
            for (int ni = 0; ni < size; ++ni) {
                std::complex<T> phase_smn = jn_smn[ni];
                if (phase_smn == std::complex<T>(0, 0)) {
                    mat[li * size + ni] = std::complex<T>(0, 0);
                    continue;
                }
                T alpha = expansion_matrix[m][li][ni];
                std::complex<T> E3_nl = E3_coupling[m][li * size + ni];
                std::complex<T> K3_val = phase_smn * alpha * E3_nl;
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
        int ioprad = 0;
        int iopnorm = normalize ? 1 : 0;
        int iopang = 2;
        auto profcn_result = cprofcn<T>(
            chi_sw, m, size, {eta_scalar}, ioprad, iopnorm, iopang, T(0)
        );

        for (int i = 0; i < size; ++i) {
            int n = m + i;
            result[m][n] = extract_angular_value_from_batch<T>(profcn_result, i);
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
        // Batched calls for all n at all nodes
        auto ext_result = cprofcn<T>(chi_sw, m, size, nodes, 0, 1, 2, T(0));
        auto int_result = cprofcn<T>(chi_body, m, size, nodes, 0, 1, 2, T(0));

        std::vector<std::vector<T>> smn_ext(size, std::vector<T>(nodes.size(), T(0)));
        std::vector<std::vector<T>> smn_int(size, std::vector<T>(nodes.size(), T(0)));
        for (int ni = 0; ni < size; ++ni) {
            for (size_t k = 0; k < nodes.size(); ++k) {
                smn_ext[ni][k] = extract_angular_value_from_batch<T>(
                    ext_result, ni, size, static_cast<int>(k)
                );
                smn_int[ni][k] = extract_angular_value_from_batch<T>(
                    int_result, ni, size, static_cast<int>(k)
                );
            }
        }

        // Build coef_mat (size x size)
        std::vector<std::vector<T>> coef_mat(size, std::vector<T>(size, T(0)));
        for (int li = 0; li < size; ++li) {
            for (int ni = 0; ni < size; ++ni) {
                T sum = T(0);
                for (size_t k = 0; k < nodes.size(); ++k) {
                    sum += smn_ext[ni][k] * smn_int[li][k] * weights[k];
                }
                coef_mat[li][ni] = sum;
            }
        }
        coef_list[m] = coef_mat;
    }
    return coef_list;
}

template<typename T>
struct FullFluidPrepResult {
    std::vector<std::vector<T>> smn_body;
    ExternalRadialResult<T> radial_external;
    RadialMatrixResult<T> radial_internal;
    std::vector<std::vector<std::vector<T>>> expansion_matrix;
};

template<typename T>
FullFluidPrepResult<T> prepare_full_fluid_boundary_data(
    int m_max,
    int n_max,
    T chi_sw,
    T chi_body,
    T xi,
    T eta_body,
    const std::vector<T>& nodes,
    const std::vector<T>& weights,
    bool normalize = true
) {
    // Collect the full set of angular, radial, and overlap ingredients needed
    // by the batched fluid/gas PSMS pathway before kernel assembly.
    FullFluidPrepResult<T> out;
    out.smn_body.assign(
        m_max + 1, std::vector<T>(n_max + 1, std::numeric_limits<T>::quiet_NaN())
    );
    out.radial_external.incident.value.assign(
        m_max + 1, std::vector<T>(n_max + 1, std::numeric_limits<T>::quiet_NaN())
    );
    out.radial_external.incident.derivative.assign(
        m_max + 1, std::vector<T>(n_max + 1, std::numeric_limits<T>::quiet_NaN())
    );
    out.radial_external.scattering.value.assign(
        m_max + 1,
        std::vector<std::complex<T>>(
            n_max + 1,
            std::complex<T>(
                std::numeric_limits<T>::quiet_NaN(), std::numeric_limits<T>::quiet_NaN()
            )
        )
    );
    out.radial_external.scattering.derivative.assign(
        m_max + 1,
        std::vector<std::complex<T>>(
            n_max + 1,
            std::complex<T>(
                std::numeric_limits<T>::quiet_NaN(), std::numeric_limits<T>::quiet_NaN()
            )
        )
    );
    out.radial_internal.value.assign(
        m_max + 1, std::vector<T>(n_max + 1, std::numeric_limits<T>::quiet_NaN())
    );
    out.radial_internal.derivative.assign(
        m_max + 1, std::vector<T>(n_max + 1, std::numeric_limits<T>::quiet_NaN())
    );
    out.expansion_matrix.assign(m_max + 1, std::vector<std::vector<T>>());

    T x1_r = xi - T(1);
    int iopnorm = normalize ? 1 : 0;
    std::vector<T> ext_args(nodes.size() + 1);
    ext_args[0] = eta_body;
    for (size_t k = 0; k < nodes.size(); ++k) {
        ext_args[k + 1] = nodes[k];
    }

    // Process neighboring orders in small blocks so the expensive profcn setup
    // is shared across nearby m values.
    int chunk_size = psms_profcn_mblock_chunk_size<T>();
    for (int m_start = 0; m_start <= m_max; m_start += chunk_size) {
        int m_count = std::min(chunk_size, m_max - m_start + 1);
        int lnum = n_max - m_start + 1;
        auto ext_result_batch = cprofcn_mblock<T>(
            chi_sw, m_start, m_count, lnum, ext_args, 2, iopnorm, 2, x1_r
        );
        auto int_result_batch = cprofcn_mblock<T>(
            chi_body, m_start, m_count, lnum, nodes, 1, iopnorm, 2, x1_r
        );

        for (int local_m = 0; local_m < m_count; ++local_m) {
            int m = m_start + local_m;
            int size = n_max - m + 1;
            if (size <= 0) {
                continue;
            }

            std::vector<std::vector<T>> smn_ext(size, std::vector<T>(nodes.size(), T(0)));
            std::vector<std::vector<T>> smn_int(size, std::vector<T>(nodes.size(), T(0)));

            for (int i = 0; i < size; ++i) {
                int n = m + i;
                T smn_val = extract_angular_value_from_mblock<T>(ext_result_batch, local_m, i, 0);
                out.smn_body[m][n] = smn_val;

                if (is_na_real(smn_val) || smn_val == T(0)) {
                    out.radial_external.incident.value[m][n] = T(0);
                    out.radial_external.incident.derivative[m][n] = T(0);
                    out.radial_external.scattering.value[m][n] = std::complex<T>(0, 0);
                    out.radial_external.scattering.derivative[m][n] = std::complex<T>(0, 0);
                } else {
                    auto rv1 = extract_radial_from_mblock<T>(ext_result_batch, local_m, i, 1);
                    auto rv3 = extract_radial_from_mblock<T>(ext_result_batch, local_m, i, 3);
                    out.radial_external.incident.value[m][n] = rv1.val_real;
                    out.radial_external.incident.derivative[m][n] = rv1.der_real;
                    out.radial_external.scattering.value[m][n] =
                        std::complex<T>(rv3.val_real, rv3.val_imag);
                    out.radial_external.scattering.derivative[m][n] =
                        std::complex<T>(rv3.der_real, rv3.der_imag);
                }

                auto rv_int = extract_radial_from_mblock<T>(int_result_batch, local_m, i, 1);
                out.radial_internal.value[m][n] = rv_int.val_real;
                out.radial_internal.derivative[m][n] = rv_int.der_real;

                for (size_t k = 0; k < nodes.size(); ++k) {
                    smn_ext[i][k] = extract_angular_value_from_mblock<T>(
                        ext_result_batch, local_m, i, static_cast<int>(k + 1)
                    );
                    smn_int[i][k] = extract_angular_value_from_mblock<T>(
                        int_result_batch, local_m, i, static_cast<int>(k)
                    );
                }
            }

            std::vector<std::vector<T>> coef_mat(size, std::vector<T>(size, T(0)));
            for (int li = 0; li < size; ++li) {
                for (int ni = 0; ni < size; ++ni) {
                    T sum = T(0);
                    for (size_t k = 0; k < nodes.size(); ++k) {
                        sum += smn_ext[ni][k] * smn_int[li][k] * weights[k];
                    }
                    coef_mat[li][ni] = sum;
                }
            }
            out.expansion_matrix[m] = coef_mat;
        }
    }

    return out;
}

// -----------------------------------------------------------
// KERNEL BOUNDARY CONDITIONS FOR FLUID-FILLED SCATTERERS
// -----------------------------------------------------------
template<typename T>
struct KernelResult {
    std::vector<std::vector<std::complex<T>>> rhs;
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
    const std::vector<std::vector<T>>& smn_matrix,
    const ExternalRadialResult<T>& radial_external
) {
    // The full fluid/gas pathway factors naturally into:
    // 1. overlap coefficients between exterior and interior angular bases,
    // 2. boundary-coupling terms for the incident and scattered fields, and
    // 3. dense kernel systems whose solution gives the retained A_mn values.
    auto expansion_matrix = alpha_integration_coefficient<T>(
        m_max, n_max, chi_sw, chi_body, nodes, weights
    );

    auto radial_internal_kind1 = radial_internal_incident_matrix<T>(
        m_max, n_max, chi_body, xi
    );

    // E1 and E3 are the Furusawa boundary factors for the incident and
    // scattered parts of the field, respectively.
    auto E1_coupling = boundary_coupling_incident_matrix<T>(
        m_max, n_max, density_body, density_sw,
        radial_external.incident.value, radial_external.incident.derivative,
        radial_internal_kind1.value, radial_internal_kind1.derivative,
        smn_matrix
    );
    auto E3_coupling = boundary_coupling_scattering_matrix<T>(
        m_max, n_max, density_body, density_sw,
        radial_external.scattering.value, radial_external.scattering.derivative,
        radial_internal_kind1.value, radial_internal_kind1.derivative,
        smn_matrix
    );

    // Assemble the right-hand side and scattered kernel for the dense modal
    // solve carried out independently for each retained order m.
    auto rhs = kernel_incident_rhs<T>(
        m_max, n_max, smn_matrix, expansion_matrix, E1_coupling
    );
    auto K3_kernel = kernel_scattering_matrix<T>(
        m_max, n_max, smn_matrix, expansion_matrix, E3_coupling
    );
    KernelResult<T> out;
    out.rhs = rhs;
    out.K3_kernel = K3_kernel;
    return out;
}

// -----------------------------------------------------------
// EXPANSION COEFFICIENT MATRIX SOLVERS
// -----------------------------------------------------------
// Double precision keeps the mature Armadillo SVD path because it is robust for
// the moderately sized dense kernel systems that arise in the retained modal
// range.
template<typename T>
std::vector<std::vector<std::complex<T>>> solve_fluid_Amn_divide_and_conquer(
    const std::vector<std::vector<std::complex<T>>>& rhs,
    const std::vector<std::vector<std::complex<T>>>& K3_kernel
) {
    int m_max = K3_kernel.size() - 1;
    std::vector<std::vector<std::complex<T>>> result(m_max + 1);

    for (int m = 0; m <= m_max; ++m) {
        int size = static_cast<int>(rhs[m].size());
        // Convert each m-block into a dense matrix/vector pair for the least-
        // squares solve.
        arma::Col<std::complex<double>> b(size);
        arma::Mat<std::complex<double>> K3_arma(size, size);
        for (int i = 0; i < size; ++i) {
            auto rhs_val = rhs[m][i];
            for (int j = 0; j < size; ++j) {
                auto val2 = K3_kernel[m][i * size + j];
                K3_arma(i, j) = std::complex<double>(
                    static_cast<double>(val2.real()),
                    static_cast<double>(val2.imag())
                );
            }
            b(i) = std::complex<double>(
                static_cast<double>(rhs_val.real()),
                static_cast<double>(rhs_val.imag())
            );
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
        
        // Convert the solved coefficients back to the templated scalar type so
        // the calling code can stay precision-agnostic.
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

template<typename T>
std::vector<std::vector<T>> reflect_smn_matrix(
    const std::vector<std::vector<T>>& smn_matrix
) {
    std::vector<std::vector<T>> result = smn_matrix;
    // Exact backscatter replaces the scattered-angle angular function with its
    // parity-related counterpart at the incident angle.
    for (size_t m = 0; m < smn_matrix.size(); ++m) {
        for (size_t n = m; n < smn_matrix[m].size(); ++n) {
            T val = smn_matrix[m][n];
            if (is_na_real(val)) {
                result[m][n] = val;
            } else {
                result[m][n] = (((n - m) % 2) == 0) ? val : -val;
            }
        }
    }
    return result;
}

template<typename T>
std::vector<std::vector<std::complex<T>>> expand_Amn_triangular(
    int m_max,
    int n_max,
    const std::vector<std::vector<std::complex<T>>>& Amn_tri
) {
    std::vector<std::vector<std::complex<T>>> Amn_mat(
        m_max + 1,
        std::vector<std::complex<T>>(
            n_max + 1,
            std::complex<T>(
                std::numeric_limits<T>::quiet_NaN(),
                std::numeric_limits<T>::quiet_NaN()
            )
        )
    );
    for (int m = 0; m <= m_max; ++m) {
        int len = static_cast<int>(Amn_tri[m].size());
        for (int i = 0; i < len; ++i) {
            int n = m + i;
            if (n > n_max) continue;
            auto val = Amn_tri[m][i];
            if (is_na_real(val.real()) || is_na_real(val.imag())) {
                Amn_mat[m][n] = std::complex<T>(
                    std::numeric_limits<T>::quiet_NaN(),
                    std::numeric_limits<T>::quiet_NaN()
                );
            } else {
                Amn_mat[m][n] = val;
            }
        }
    }
    return Amn_mat;
}

template<typename T>
struct LUPFactorization {
    std::vector<std::complex<T>> lu;
    std::vector<int> row_piv;
    std::vector<int> col_piv;
    T tolerance;
    bool singular;
};

template<typename T>
LUPFactorization<T> lup_decompose(
    const std::vector<std::complex<T>>& A_in,
    int n
) {
    LUPFactorization<T> fac;
    fac.lu = A_in;
    fac.row_piv.resize(n);
    fac.col_piv.resize(n);
    for (int i = 0; i < n; ++i) {
        fac.row_piv[i] = i;
        fac.col_piv[i] = i;
    }
    fac.singular = false;

    T max_entry = T(0);
    // Scale the pivot tolerance to the largest matrix entry so the singularity
    // test remains meaningful across very different kernel magnitudes.
    for (const auto& val : fac.lu) {
        T val_abs_sq = complex_abs_sq(val);
        if (val_abs_sq > max_entry) {
            max_entry = val_abs_sq;
        }
    }
    max_entry = precsqrt(max_entry);
    fac.tolerance = std::max(T(1), max_entry) * T(64) * preceps<T>();

    for (int k = 0; k < n; ++k) {
        int pivot_row = k;
        int pivot_col = k;
        T pivot_abs_sq = complex_abs_sq(fac.lu[k * n + k]);
        for (int i = k + 1; i < n; ++i) {
            for (int j = k; j < n; ++j) {
                T cand_abs_sq = complex_abs_sq(fac.lu[i * n + j]);
                if (cand_abs_sq > pivot_abs_sq) {
                    pivot_row = i;
                    pivot_col = j;
                    pivot_abs_sq = cand_abs_sq;
                }
            }
        }
        for (int j = k; j < n; ++j) {
            T cand_abs_sq = complex_abs_sq(fac.lu[k * n + j]);
            if (cand_abs_sq > pivot_abs_sq) {
                pivot_row = k;
                pivot_col = j;
                pivot_abs_sq = cand_abs_sq;
            }
        }

        if (pivot_abs_sq <= fac.tolerance * fac.tolerance) {
            fac.singular = true;
            return fac;
        }

        if (pivot_row != k) {
            for (int j = 0; j < n; ++j) {
                std::swap(fac.lu[k * n + j], fac.lu[pivot_row * n + j]);
            }
            std::swap(fac.row_piv[k], fac.row_piv[pivot_row]);
        }

        if (pivot_col != k) {
            for (int i = 0; i < n; ++i) {
                std::swap(fac.lu[i * n + k], fac.lu[i * n + pivot_col]);
            }
            std::swap(fac.col_piv[k], fac.col_piv[pivot_col]);
        }

        std::complex<T> diag = fac.lu[k * n + k];
        for (int i = k + 1; i < n; ++i) {
            fac.lu[i * n + k] /= diag;
            std::complex<T> factor = fac.lu[i * n + k];
            for (int j = k + 1; j < n; ++j) {
                fac.lu[i * n + j] -= factor * fac.lu[k * n + j];
            }
        }
    }

    return fac;
}

template<typename T>
std::vector<std::complex<T>> lup_solve(
    const LUPFactorization<T>& fac,
    const std::vector<std::complex<T>>& b
) {
    int n = static_cast<int>(fac.row_piv.size());
    std::vector<std::complex<T>> x_perm(n), x(n), y(n);

    for (int i = 0; i < n; ++i) {
        y[i] = b[fac.row_piv[i]];
        for (int j = 0; j < i; ++j) {
            y[i] -= fac.lu[i * n + j] * y[j];
        }
    }

    for (int i = n - 1; i >= 0; --i) {
        std::complex<T> rhs = y[i];
        for (int j = i + 1; j < n; ++j) {
            rhs -= fac.lu[i * n + j] * x_perm[j];
        }
        x_perm[i] = rhs / fac.lu[i * n + i];
    }

    for (int i = 0; i < n; ++i) {
        x[fac.col_piv[i]] = x_perm[i];
    }

    return x;
}

template<typename T>
std::vector<std::complex<T>> matvec_product(
    const std::vector<std::complex<T>>& A,
    const std::vector<std::complex<T>>& x,
    int n
) {
    std::vector<std::complex<T>> out(n, std::complex<T>(0, 0));
    for (int i = 0; i < n; ++i) {
        std::complex<T> sum(0, 0);
        for (int j = 0; j < n; ++j) {
            sum += A[i * n + j] * x[j];
        }
        out[i] = sum;
    }
    return out;
}

template<typename T>
ExternalRadialResult<T> radial_external_matrices(
    int m_max,
    int n_max,
    T chi_sw,
    T xi,
    const std::vector<std::vector<T>>& smn_matrix
) {
    ExternalRadialResult<T> out;
    out.incident.value.assign(
        m_max + 1, std::vector<T>(n_max + 1, std::numeric_limits<T>::quiet_NaN())
    );
    out.incident.derivative.assign(
        m_max + 1, std::vector<T>(n_max + 1, std::numeric_limits<T>::quiet_NaN())
    );
    out.scattering.value.assign(
        m_max + 1,
        std::vector<std::complex<T>>(
            n_max + 1,
            std::complex<T>(
                std::numeric_limits<T>::quiet_NaN(), std::numeric_limits<T>::quiet_NaN()
            )
        )
    );
    out.scattering.derivative.assign(
        m_max + 1,
        std::vector<std::complex<T>>(
            n_max + 1,
            std::complex<T>(
                std::numeric_limits<T>::quiet_NaN(), std::numeric_limits<T>::quiet_NaN()
            )
        )
    );

    T x1_r = xi - T(1);

    // Only keep radial terms for modes whose angular factors are active on the
    // body boundary. Zeroing the inactive modes keeps later kernel assembly
    // smaller and avoids wasted special-function work.
    for (int m = 0; m <= m_max; ++m) {
        int size = n_max - m + 1;
        std::vector<int> valid_indices;
        for (int i = 0; i < size; ++i) {
            int n = m + i;
            T snm_val = smn_matrix[m][n];
            if (!std::isnan(static_cast<double>(snm_val)) && snm_val != T(0)) {
                valid_indices.push_back(n);
            } else {
                out.incident.value[m][n] = T(0);
                out.incident.derivative[m][n] = T(0);
                out.scattering.value[m][n] = std::complex<T>(0, 0);
                out.scattering.derivative[m][n] = std::complex<T>(0, 0);
            }
        }

        if (valid_indices.empty()) continue;

        int lnum = valid_indices.back() - m + 1;
        auto result_both = cprofcn<T>(chi_sw, m, lnum, {T(1)}, 2, 0, 0, x1_r);

        for (int n : valid_indices) {
            int idx = n - m;
            auto rv1 = extract_radial_from_batch<T>(
                result_both, result_both, idx, 1, m, n, chi_sw, xi
            );
            auto rv3 = extract_radial_from_batch<T>(
                result_both, result_both, idx, 3, m, n, chi_sw, xi
            );

            out.incident.value[m][n] = rv1.val_real;
            out.incident.derivative[m][n] = rv1.der_real;
            out.scattering.value[m][n] = std::complex<T>(rv3.val_real, rv3.val_imag);
            out.scattering.derivative[m][n] = std::complex<T>(rv3.der_real, rv3.der_imag);
        }
    }

    return out;
}

template<typename T>
ExternalBoundaryResult<T> external_boundary_matrices(
    int m_max,
    int n_max,
    T chi_sw,
    T xi,
    T eta_body,
    bool normalize
) {
    ExternalBoundaryResult<T> out;
    out.smn.assign(
        m_max + 1, std::vector<T>(n_max + 1, std::numeric_limits<T>::quiet_NaN())
    );
    out.radial.incident.value.assign(
        m_max + 1, std::vector<T>(n_max + 1, std::numeric_limits<T>::quiet_NaN())
    );
    out.radial.incident.derivative.assign(
        m_max + 1, std::vector<T>(n_max + 1, std::numeric_limits<T>::quiet_NaN())
    );
    out.radial.scattering.value.assign(
        m_max + 1,
        std::vector<std::complex<T>>(
            n_max + 1,
            std::complex<T>(
                std::numeric_limits<T>::quiet_NaN(), std::numeric_limits<T>::quiet_NaN()
            )
        )
    );
    out.radial.scattering.derivative.assign(
        m_max + 1,
        std::vector<std::complex<T>>(
            n_max + 1,
            std::complex<T>(
                std::numeric_limits<T>::quiet_NaN(), std::numeric_limits<T>::quiet_NaN()
            )
        )
    );

    T x1_r = xi - T(1);
    int iopnorm = normalize ? 1 : 0;

    // Evaluate the exterior angular basis on the body boundary together with
    // the matching incident and scattered radial terms.
    for (int m = 0; m <= m_max; ++m) {
        int size = n_max - m + 1;
        if (size <= 0) continue;

        auto result_both = cprofcn<T>(chi_sw, m, size, {eta_body}, 2, iopnorm, 2, x1_r);

        for (int i = 0; i < size; ++i) {
            int n = m + i;
            T smn_val = extract_angular_value_from_batch<T>(result_both, i);
            out.smn[m][n] = smn_val;

            if (is_na_real(smn_val) || smn_val == T(0)) {
                out.radial.incident.value[m][n] = T(0);
                out.radial.incident.derivative[m][n] = T(0);
                out.radial.scattering.value[m][n] = std::complex<T>(0, 0);
                out.radial.scattering.derivative[m][n] = std::complex<T>(0, 0);
                continue;
            }

            auto rv1 = extract_radial_from_batch<T>(
                result_both, result_both, i, 1, m, n, chi_sw, xi
            );
            auto rv3 = extract_radial_from_batch<T>(
                result_both, result_both, i, 3, m, n, chi_sw, xi
            );

            out.radial.incident.value[m][n] = rv1.val_real;
            out.radial.incident.derivative[m][n] = rv1.der_real;
            out.radial.scattering.value[m][n] = std::complex<T>(rv3.val_real, rv3.val_imag);
            out.radial.scattering.derivative[m][n] = std::complex<T>(rv3.der_real, rv3.der_imag);
        }
    }

    return out;
}

template<typename T>
std::vector<std::complex<T>> solve_linear_system_lup_refined(
    const std::vector<std::complex<T>>& A,
    const std::vector<std::complex<T>>& b,
    int n,
    int n_refine = 1
) {
    auto fac = lup_decompose<T>(A, n);
    if (fac.singular) {
        throw std::runtime_error("Fluid kernel system is singular or ill-conditioned");
    }

    auto x = lup_solve(fac, b);

    // Residual refinement is cheap relative to the decomposition and helps keep
    // the quad solve stable without introducing a second backend dependency.
    for (int iter = 0; iter < n_refine; ++iter) {
        auto Ax = matvec_product(A, x, n);
        std::vector<std::complex<T>> r(n, std::complex<T>(0, 0));
        T max_resid = T(0);
        for (int i = 0; i < n; ++i) {
            r[i] = b[i] - Ax[i];
            T resid_abs_sq = complex_abs_sq(r[i]);
            if (resid_abs_sq > max_resid) {
                max_resid = resid_abs_sq;
            }
        }
        if (precsqrt(max_resid) <= fac.tolerance) {
            break;
        }
        auto delta = lup_solve(fac, r);
        for (int i = 0; i < n; ++i) {
            x[i] += delta[i];
        }
    }

    return x;
}

template<typename T>
std::vector<std::vector<std::complex<T>>> solve_fluid_Amn_native(
    const std::vector<std::vector<std::complex<T>>>& rhs,
    const std::vector<std::vector<std::complex<T>>>& K3_kernel
) {
    int m_max = static_cast<int>(K3_kernel.size()) - 1;
    std::vector<std::vector<std::complex<T>>> result(m_max + 1);

    // Solve one dense truncated system per retained order m.
    for (int m = 0; m <= m_max; ++m) {
        int size = static_cast<int>(rhs[m].size());
        result[m] = solve_linear_system_lup_refined<T>(K3_kernel[m], rhs[m], size);
    }

    return result;
}

template<typename T>
std::vector<std::vector<std::complex<T>>> solve_fluid_Amn(
    const std::vector<std::vector<std::complex<T>>>& rhs,
    const std::vector<std::vector<std::complex<T>>>& K3_kernel
) {
    return solve_fluid_Amn_divide_and_conquer<T>(rhs, K3_kernel);
}

#ifdef __GNUC__
template<>
std::vector<std::vector<std::complex<__float128>>> solve_fluid_Amn<__float128>(
    const std::vector<std::vector<std::complex<__float128>>>& rhs,
    const std::vector<std::vector<std::complex<__float128>>>& K3_kernel
) {
    return solve_fluid_Amn_native<__float128>(rhs, K3_kernel);
}
#endif

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
    // For pressure release, each retained mode decouples after truncation, so
    // only the exterior radial functions are needed.
    // Create a dummy Smn-based matrix that comprises all 1's.
    std::vector<std::vector<T>> smn_matrix(m_max + 1, std::vector<T>(n_max + 1, 1.0));

    auto radial_external = radial_external_matrices<T>(
        m_max, n_max, chi_sw, xi, smn_matrix
    );

    std::vector<std::vector<std::complex<T>>> amn(
        m_max + 1, std::vector<std::complex<T>>(n_max + 1, std::complex<T>(0, 0))
    );
    for (int m = 0; m <= m_max; ++m) {
        for (int n = m; n <= n_max; ++n) {
            auto v = radial_external.incident.value[m][n];
            auto s = radial_external.scattering.value[m][n];
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
    // The rigid case is also diagonal after truncation, but the derivative
    // ratio replaces the pressure-release value ratio.
    // Create a dummy Smn-based matrix that comprises all 1's.
    std::vector<std::vector<T>> smn_matrix(m_max + 1, std::vector<T>(n_max + 1, 1.0));

    auto radial_external = radial_external_matrices<T>(
        m_max, n_max, chi_sw, xi, smn_matrix
    );

    std::vector<std::vector<std::complex<T>>> amn(
        m_max + 1, std::vector<std::complex<T>>(n_max + 1, std::complex<T>(0, 0))
    );
    for (int m = 0; m <= m_max; ++m) {
        for (int n = m; n <= n_max; ++n) {
            auto v = radial_external.incident.derivative[m][n];
            auto s = radial_external.scattering.derivative[m][n];
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
struct FluidAmnResult {
    std::vector<std::vector<T>> smn_body;
    std::vector<std::vector<std::complex<T>>> Amn_tri;
};

template<typename T>
FluidAmnResult<T> fluid_Amn_triangular_from_boundary(
    int m_max,
    int n_max,
    T chi_sw,
    T chi_body,
    T xi,
    T eta_body,
    T density_body,
    T density_sw,
    const std::vector<T>& nodes,
    const std::vector<T>& weights
) {
    // Build the full fluid/gas boundary machinery from explicit boundary
    // matrices, then solve the dense modal systems.
    auto external_boundary = external_boundary_matrices<T>(
        m_max, n_max, chi_sw, xi, eta_body, true
    );
    auto kernels = fluid_boundary_kernels<T>(
        m_max, n_max, chi_sw, chi_body, xi, density_body, density_sw,
        nodes, weights, external_boundary.smn, external_boundary.radial
    );

    FluidAmnResult<T> out;
    out.smn_body = std::move(external_boundary.smn);
    out.Amn_tri = solve_fluid_Amn<T>(kernels.rhs, kernels.K3_kernel);
    return out;
}

template<typename T>
FluidAmnResult<T> fluid_Amn_triangular_batched(
    int m_max,
    int n_max,
    T chi_sw,
    T chi_body,
    T xi,
    T eta_body,
    T density_body,
    T density_sw,
    const std::vector<T>& nodes,
    const std::vector<T>& weights
) {
    // The batched pathway precomputes the external, internal, and overlap data
    // together so the expensive spheroidal setup is shared across nearby m.
    auto prep = prepare_full_fluid_boundary_data<T>(
        m_max, n_max, chi_sw, chi_body, xi, eta_body, nodes, weights, true
    );

    auto E1_coupling = boundary_coupling_incident_matrix<T>(
        m_max, n_max, density_body, density_sw,
        prep.radial_external.incident.value, prep.radial_external.incident.derivative,
        prep.radial_internal.value, prep.radial_internal.derivative,
        prep.smn_body
    );
    auto E3_coupling = boundary_coupling_scattering_matrix<T>(
        m_max, n_max, density_body, density_sw,
        prep.radial_external.scattering.value, prep.radial_external.scattering.derivative,
        prep.radial_internal.value, prep.radial_internal.derivative,
        prep.smn_body
    );
    auto rhs = kernel_incident_rhs<T>(
        m_max, n_max, prep.smn_body, prep.expansion_matrix, E1_coupling
    );
    auto K3_kernel = kernel_scattering_matrix<T>(
        m_max, n_max, prep.smn_body, prep.expansion_matrix, E3_coupling
    );
    auto Amn = solve_fluid_Amn<T>(
        rhs, K3_kernel
    );

    FluidAmnResult<T> out;
    out.smn_body = std::move(prep.smn_body);
    out.Amn_tri = std::move(Amn);
    return out;
}

template<typename T>
FluidAmnResult<T> fluid_Amn_triangular(
    int m_max,
    int n_max,
    T chi_sw,
    T chi_body,
    T xi,
    T eta_body,
    T density_body,
    T density_sw,
    const std::vector<T>& nodes,
    const std::vector<T>& weights
) {
#ifdef __GNUC__
    if constexpr (std::is_same_v<T, __float128>) {
        // Quad precision benefits the most from the batched pathway because the
        // spheroidal backend is the dominant cost.
        return fluid_Amn_triangular_batched<T>(
            m_max, n_max, chi_sw, chi_body, xi, eta_body, density_body, density_sw,
            nodes, weights
        );
    }
#endif
    return fluid_Amn_triangular_from_boundary<T>(
        m_max, n_max, chi_sw, chi_body, xi, eta_body, density_body, density_sw,
        nodes, weights
    );
}

template<typename T>
std::vector<std::vector<std::complex<T>>> fluid_Amn_expansion_matrix(
    int m_max,
    int n_max,
    T chi_sw,
    T chi_body,
    T xi,
    T density_body,
    T density_sw,
    const std::vector<T>& nodes,
    const std::vector<T>& weights,
    const std::vector<std::vector<T>>& smn_matrix,
    const ExternalRadialResult<T>& radial_external
) {
    auto kernels = fluid_boundary_kernels<T>(
        m_max, n_max, chi_sw, chi_body, xi, density_body, density_sw,
        nodes, weights, smn_matrix, radial_external
    );
    auto Amn_tri = solve_fluid_Amn<T>(kernels.rhs, kernels.K3_kernel);
    return expand_Amn_triangular<T>(m_max, n_max, Amn_tri);
}

// Simplified fluid-filled
template<typename T>
std::vector<std::vector<std::complex<T>>> simplified_fluid_Amn_triangular(
    int m_max,
    int n_max,
    T chi_sw,
    T chi_body,
    T xi,
    T density_body,
    T density_sw,
    const ExternalRadialResult<T>& radial_external
) {
    auto radial_internal_kind1 = radial_internal_incident_matrix<T>(
        m_max, n_max, chi_body, xi
    );

    // Coupling matrices
    auto E1_coupling = simplified_boundary_coupling_incident_matrix<T>(
        m_max, n_max, density_body, density_sw,
        radial_external.incident.value, radial_external.incident.derivative,
        radial_internal_kind1.value, radial_internal_kind1.derivative
    );
    auto E3_coupling = simplified_boundary_coupling_scattering_matrix<T>(
        m_max, n_max, density_body, density_sw,
        radial_external.scattering.value, radial_external.scattering.derivative,
        radial_internal_kind1.value, radial_internal_kind1.derivative
    );

    // The simplified fluid approximation discards degree coupling and keeps
    // only the diagonal E1/E3 ratio for each retained (m, n).
    std::vector<std::vector<std::complex<T>>> Amn_tri(
        m_max + 1,
        std::vector<std::complex<T>>()
    );
    for (int m = 0; m <= m_max; ++m) {
        Amn_tri[m].assign(
            n_max - m + 1,
            std::complex<T>(
                std::numeric_limits<T>::quiet_NaN(),
                std::numeric_limits<T>::quiet_NaN()
            )
        );
        for (int n = m; n <= n_max; ++n) {
            std::complex<T> e1 = E1_coupling[m][n];
            std::complex<T> e3 = E3_coupling[m][n];
            if (is_na_real(e1.real()) || is_na_real(e1.imag()) ||
                e1 == std::complex<T>(0, 0) ||
                is_na_real(e3.real()) || is_na_real(e3.imag()) ||
                e3 == std::complex<T>(0, 0)) {
                Amn_tri[m][n - m] = std::complex<T>(
                    std::numeric_limits<T>::quiet_NaN(), std::numeric_limits<T>::quiet_NaN()
                );
            } else {
                Amn_tri[m][n - m] = -e1 / e3;
            }
        }
    }
    return Amn_tri;
}

template<typename T>
std::vector<std::vector<std::complex<T>>> simplified_fluid_Amn_expansion_matrix(
    int m_max,
    int n_max,
    T chi_sw,
    T chi_body,
    T xi,
    T density_body,
    T density_sw,
    const ExternalRadialResult<T>& radial_external
) {
    return expand_Amn_triangular<T>(
        m_max,
        n_max,
        simplified_fluid_Amn_triangular<T>(
            m_max, n_max, chi_sw, chi_body, xi, density_body, density_sw,
            radial_external
        )
    );
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
    const std::vector<T>& azimuth,
    const std::vector<std::vector<T>>& smn_body_matrix,
    const std::vector<std::vector<T>>& smn_scatter_matrix,
    const std::vector<std::vector<std::complex<T>>>& Amn
) {
    if (azimuth.size() != static_cast<size_t>(m_max + 1)){
        throw std::invalid_argument("azimuth must have length m_max + 1");
    }
    std::vector<T> nu(m_max + 1);
    for (int m = 0; m <= m_max; ++m)
        nu[m] = (m == 0) ? 1.0 : 2.0;

    std::complex<T> fbs_sum(0, 0), c(0, 0);

    // Sum the retained modal field using compensated summation because the
    // physically relevant answer can come from cancellation among sizable
    // complex terms.
    for (int m = 0; m <= m_max; ++m) {
        for (int n = m; n <= n_max; ++n) {
            std::complex<T> Amn_val(0, 0);
            if (m < static_cast<int>(Amn.size()) && n < static_cast<int>(Amn[m].size()))
                Amn_val = Amn[m][n];

            // Pull the incident and scattered angular factors for the same
            // retained mode.
            T smn_body = smn_body_matrix[m][n];
            T smn_scatter = smn_scatter_matrix[m][n];

            // Skip if any are NaN
            if (is_na_real(smn_body) || is_na_real(smn_scatter) ||
                is_na_real(Amn_val.real()) || is_na_real(Amn_val.imag()))
                continue;
            
            // Kahan-style compensated summation limits cancellation error in
            // the long modal tail.
            std::complex<T> term = nu[m] * smn_body * smn_scatter * Amn_val * azimuth[m];
            std::complex<T> y = term - c;
            std::complex<T> t = fbs_sum + y;
            c = (t - fbs_sum) - y;
            fbs_sum = t;
        }
    }
    return fbs_sum;
}

template<typename T>
inline bool psms_modal_band_is_negligible(
    int m,
    int m_max,
    T band_abs,
    T prev_band_abs,
    T max_band
);

template<typename T>
inline bool psms_modal_term_is_negligible(
    int idx,
    int max_idx,
    T term_abs,
    T prev_term_abs,
    T max_term
);

template<typename T>
std::complex<T> compute_fbs_backscatter(
    int m_max,
    int n_max,
    const std::vector<std::vector<T>>& smn_body_matrix,
    const std::vector<std::vector<std::complex<T>>>& Amn,
    bool adaptive = false
) {
    std::complex<T> fbs_sum(0, 0), c(0, 0);
    T max_band = T(0);
    T prev_band = T(-1);
    int n_small_bands = 0;

    for (int m = 0; m <= m_max; ++m) {
        T nu_m = (m == 0) ? T(1) : T(2);
        std::complex<T> band_sum(0, 0), band_c(0, 0);
        T max_term = T(0);
        T prev_term = T(-1);
        int n_small_terms = 0;

        for (int n = m; n <= n_max; ++n) {
            T smn_body = smn_body_matrix[m][n];
            std::complex<T> Amn_val = Amn[m][n];

            if (is_na_real(smn_body) ||
                is_na_real(Amn_val.real()) || is_na_real(Amn_val.imag())) {
                continue;
            }

            T parity_n = ((n % 2) == 0) ? T(1) : T(-1);
            std::complex<T> term = nu_m * parity_n * smn_body * smn_body * Amn_val;
            T term_abs = precsqrt(complex_abs_sq(term));
            max_term = std::max(max_term, term_abs);
            if (adaptive &&
                psms_modal_term_is_negligible<T>(
                    n - m, n_max - m, term_abs, prev_term, max_term
                )) {
                ++n_small_terms;
            } else {
                n_small_terms = 0;
            }
            prev_term = term_abs;

            std::complex<T> y_band = term - band_c;
            std::complex<T> t_band = band_sum + y_band;
            band_c = (t_band - band_sum) - y_band;
            band_sum = t_band;

            if (adaptive && n_small_terms >= 2) {
                break;
            }
        }

        std::complex<T> y = band_sum - c;
        std::complex<T> t = fbs_sum + y;
        c = (t - fbs_sum) - y;
        fbs_sum = t;

        T band_abs = precsqrt(complex_abs_sq(band_sum));
        max_band = std::max(max_band, band_abs);
        if (adaptive &&
            psms_modal_band_is_negligible<T>(m, m_max, band_abs, prev_band, max_band)) {
            ++n_small_bands;
        } else {
            n_small_bands = 0;
        }
        prev_band = band_abs;
        if (adaptive && n_small_bands >= 2) {
            return fbs_sum;
        }
    }

    return fbs_sum;
}

template<typename T>
std::complex<T> compute_fbs_backscatter_triangular(
    int m_max,
    const std::vector<std::vector<T>>& smn_body_matrix,
    const std::vector<std::vector<std::complex<T>>>& Amn_tri,
    bool adaptive = false
) {
    std::complex<T> fbs_sum(0, 0), c(0, 0);
    T max_band = T(0);
    T prev_band = T(-1);
    int n_small_bands = 0;

    for (int m = 0; m <= m_max; ++m) {
        T nu_m = (m == 0) ? T(1) : T(2);
        int len = static_cast<int>(Amn_tri[m].size());
        std::complex<T> band_sum(0, 0), band_c(0, 0);
        T max_term = T(0);
        T prev_term = T(-1);
        int n_small_terms = 0;

        for (int i = 0; i < len; ++i) {
            int n = m + i;
            T smn_body = smn_body_matrix[m][n];
            std::complex<T> Amn_val = Amn_tri[m][i];

            if (is_na_real(smn_body) ||
                is_na_real(Amn_val.real()) || is_na_real(Amn_val.imag())) {
                continue;
            }

            T parity_n = ((n % 2) == 0) ? T(1) : T(-1);
            std::complex<T> term = nu_m * parity_n * smn_body * smn_body * Amn_val;
            T term_abs = precsqrt(complex_abs_sq(term));
            max_term = std::max(max_term, term_abs);
            if (adaptive &&
                psms_modal_term_is_negligible<T>(i, len - 1, term_abs, prev_term, max_term)) {
                ++n_small_terms;
            } else {
                n_small_terms = 0;
            }
            prev_term = term_abs;

            std::complex<T> y_band = term - band_c;
            std::complex<T> t_band = band_sum + y_band;
            band_c = (t_band - band_sum) - y_band;
            band_sum = t_band;

            if (adaptive && n_small_terms >= 2) {
                break;
            }
        }

        std::complex<T> y = band_sum - c;
        std::complex<T> t = fbs_sum + y;
        c = (t - fbs_sum) - y;
        fbs_sum = t;

        T band_abs = precsqrt(complex_abs_sq(band_sum));
        max_band = std::max(max_band, band_abs);
        if (adaptive &&
            psms_modal_band_is_negligible<T>(m, m_max, band_abs, prev_band, max_band)) {
            ++n_small_bands;
        } else {
            n_small_bands = 0;
        }
        prev_band = band_abs;
        if (adaptive && n_small_bands >= 2) {
            return fbs_sum;
        }
    }

    return fbs_sum;
}

template<typename T>
inline T psms_modal_rel_tol() {
    return static_cast<T>(1e-8L);
}

#ifdef __GNUC__
template<>
inline __float128 psms_modal_rel_tol<__float128>() {
    return static_cast<__float128>(1e-14L);
}
#endif

template<typename T>
inline T psms_modal_abs_tol() {
    return static_cast<T>(1e-12L);
}

#ifdef __GNUC__
template<>
inline __float128 psms_modal_abs_tol<__float128>() {
    return static_cast<__float128>(1e-18L);
}
#endif

template<typename T>
inline bool psms_modal_band_is_negligible(
    int m,
    int m_max,
    T band_abs,
    T prev_band_abs,
    T max_band
) {
    if (m < 2 || m >= m_max || prev_band_abs <= T(0)) {
        return false;
    }
    T rel_cutoff = psms_modal_rel_tol<T>() * std::max(T(1), max_band);
    T cutoff = std::max(psms_modal_abs_tol<T>(), rel_cutoff);
    T gradient = precabs(band_abs - prev_band_abs) / std::max(prev_band_abs, T(1));
    return band_abs <= cutoff && gradient <= static_cast<T>(1.0L);
}

template<typename T>
inline bool psms_modal_term_is_negligible(
    int idx,
    int max_idx,
    T term_abs,
    T prev_term_abs,
    T max_term
) {
    if (idx < 4 || idx >= max_idx || prev_term_abs <= T(0)) {
        return false;
    }
    T rel_cutoff = psms_modal_rel_tol<T>() * std::max(T(1), max_term);
    T cutoff = std::max(psms_modal_abs_tol<T>(), rel_cutoff);
    T gradient = precabs(term_abs - prev_term_abs) / std::max(prev_term_abs, T(1));
    return term_abs <= cutoff && gradient <= static_cast<T>(1.0L);
}

// Backscatter for pressure-release PSMS can be evaluated directly from the
// exterior angular and radial functions, with no dense kernel solve.
template<typename T>
std::complex<T> psms_backscatter_pressure_release(
    int m_max,
    int n_max,
    T chi_sw,
    T xi,
    T eta_body,
    bool adaptive = false
) {
    std::complex<T> fbs_sum(0, 0), c(0, 0);
    T x1_r = xi - T(1);
    int chunk_size = psms_profcn_mblock_chunk_size<T>();
    T max_band = T(0);
    T prev_band = T(-1);
    int n_small_bands = 0;

    for (int m_start = 0; m_start <= m_max; m_start += chunk_size) {
        int m_count = std::min(chunk_size, m_max - m_start + 1);
        int lnum = n_max - m_start + 1;
        auto result_both = cprofcn_mblock<T>(
            chi_sw, m_start, m_count, lnum, {eta_body}, 2, 1, 2, x1_r
        );

        for (int local_m = 0; local_m < m_count; ++local_m) {
            int m = m_start + local_m;
            int size = n_max - m + 1;
            T nu_m = (m == 0) ? T(1) : T(2);
            std::complex<T> band_sum(0, 0), band_c(0, 0);
            T max_term = T(0);
            T prev_term = T(-1);
            int n_small_terms = 0;

            for (int i = 0; i < size; ++i) {
                int n = m + i;
                T smn = extract_angular_value_from_mblock<T>(result_both, local_m, i, 0);
                if (is_na_real(smn)) continue;

                auto rv1 = extract_radial_from_mblock<T>(result_both, local_m, i, 1);
                auto rv3 = extract_radial_from_mblock<T>(result_both, local_m, i, 3);
                std::complex<T> denom(rv3.val_real, rv3.val_imag);
                if (is_na_real(denom.real()) || is_na_real(denom.imag()) ||
                    denom == std::complex<T>(0, 0)) {
                    continue;
                }

                std::complex<T> amn = -std::complex<T>(rv1.val_real, T(0)) / denom;
                T parity_n = ((n % 2) == 0) ? T(1) : T(-1);
                std::complex<T> term = nu_m * parity_n * smn * smn * amn;
                T term_abs = precsqrt(complex_abs_sq(term));
                max_term = std::max(max_term, term_abs);
                if (adaptive &&
                    psms_modal_term_is_negligible<T>(i, size - 1, term_abs, prev_term, max_term)) {
                    ++n_small_terms;
                } else {
                    n_small_terms = 0;
                }
                prev_term = term_abs;
                std::complex<T> y_band = term - band_c;
                std::complex<T> t_band = band_sum + y_band;
                band_c = (t_band - band_sum) - y_band;
                band_sum = t_band;

                if (adaptive && n_small_terms >= 2) {
                    break;
                }
            }

            std::complex<T> y = band_sum - c;
            std::complex<T> t = fbs_sum + y;
            c = (t - fbs_sum) - y;
            fbs_sum = t;

            T band_abs = precsqrt(complex_abs_sq(band_sum));
            max_band = std::max(max_band, band_abs);
            if (adaptive &&
                psms_modal_band_is_negligible<T>(m, m_max, band_abs, prev_band, max_band)) {
                ++n_small_bands;
            } else {
                n_small_bands = 0;
            }
            prev_band = band_abs;
            // Stop once several successive m-bands are numerically silent.
            if (adaptive && n_small_bands >= 2) {
                return fbs_sum;
            }
        }
    }

    return fbs_sum;
}

// The rigid backscatter path is analogous to the pressure-release case, except
// the derivative ratio enforces the vanishing normal-velocity condition.
template<typename T>
std::complex<T> psms_backscatter_fixed_rigid(
    int m_max,
    int n_max,
    T chi_sw,
    T xi,
    T eta_body,
    bool adaptive = false
) {
    std::complex<T> fbs_sum(0, 0), c(0, 0);
    T x1_r = xi - T(1);
    int chunk_size = psms_profcn_mblock_chunk_size<T>();
    T max_band = T(0);
    T prev_band = T(-1);
    int n_small_bands = 0;

    for (int m_start = 0; m_start <= m_max; m_start += chunk_size) {
        int m_count = std::min(chunk_size, m_max - m_start + 1);
        int lnum = n_max - m_start + 1;
        auto result_both = cprofcn_mblock<T>(
            chi_sw, m_start, m_count, lnum, {eta_body}, 2, 1, 2, x1_r
        );

        for (int local_m = 0; local_m < m_count; ++local_m) {
            int m = m_start + local_m;
            int size = n_max - m + 1;
            T nu_m = (m == 0) ? T(1) : T(2);
            std::complex<T> band_sum(0, 0), band_c(0, 0);
            T max_term = T(0);
            T prev_term = T(-1);
            int n_small_terms = 0;

            for (int i = 0; i < size; ++i) {
                int n = m + i;
                T smn = extract_angular_value_from_mblock<T>(result_both, local_m, i, 0);
                if (is_na_real(smn)) continue;

                auto rv1 = extract_radial_from_mblock<T>(result_both, local_m, i, 1);
                auto rv3 = extract_radial_from_mblock<T>(result_both, local_m, i, 3);
                std::complex<T> denom(rv3.der_real, rv3.der_imag);
                if (is_na_real(denom.real()) || is_na_real(denom.imag()) ||
                    denom == std::complex<T>(0, 0)) {
                    continue;
                }

                std::complex<T> amn = -std::complex<T>(rv1.der_real, T(0)) / denom;
                T parity_n = ((n % 2) == 0) ? T(1) : T(-1);
                std::complex<T> term = nu_m * parity_n * smn * smn * amn;
                T term_abs = precsqrt(complex_abs_sq(term));
                max_term = std::max(max_term, term_abs);
                if (adaptive &&
                    psms_modal_term_is_negligible<T>(i, size - 1, term_abs, prev_term, max_term)) {
                    ++n_small_terms;
                } else {
                    n_small_terms = 0;
                }
                prev_term = term_abs;
                std::complex<T> y_band = term - band_c;
                std::complex<T> t_band = band_sum + y_band;
                band_c = (t_band - band_sum) - y_band;
                band_sum = t_band;

                if (adaptive && n_small_terms >= 2) {
                    break;
                }
            }

            std::complex<T> y = band_sum - c;
            std::complex<T> t = fbs_sum + y;
            c = (t - fbs_sum) - y;
            fbs_sum = t;

            T band_abs = precsqrt(complex_abs_sq(band_sum));
            max_band = std::max(max_band, band_abs);
            if (adaptive &&
                psms_modal_band_is_negligible<T>(m, m_max, band_abs, prev_band, max_band)) {
                ++n_small_bands;
            } else {
                n_small_bands = 0;
            }
            prev_band = band_abs;
            if (adaptive && n_small_bands >= 2) {
                return fbs_sum;
            }
        }
    }

    return fbs_sum;
}

// Solve a single dense fluid/gas kernel system for one retained order m.
// Double precision uses the Armadillo SVD path here; quad precision dispatches
// to the native LUP-refined solve below.
template<typename T>
std::vector<std::complex<T>> solve_fluid_mode_system(
    const std::vector<std::complex<T>>& K3_kernel,
    const std::vector<std::complex<T>>& rhs,
    int size
) {
    arma::Col<std::complex<double>> b(size);
    arma::Mat<std::complex<double>> K3_arma(size, size);
    for (int i = 0; i < size; ++i) {
        auto rhs_val = rhs[i];
        for (int j = 0; j < size; ++j) {
            auto val = K3_kernel[i * size + j];
            K3_arma(i, j) = std::complex<double>(
                static_cast<double>(val.real()),
                static_cast<double>(val.imag())
            );
        }
        b(i) = std::complex<double>(
            static_cast<double>(rhs_val.real()),
            static_cast<double>(rhs_val.imag())
        );
    }
    arma::Mat<std::complex<double>> U, V;
    arma::Col<double> s;
    bool svd_ok = arma::svd(U, s, V, K3_arma);
    if (!svd_ok) {
        throw std::runtime_error("SVD failed");
    }
    double tol =
        std::max(size, size) * s.max() * std::numeric_limits<double>::epsilon();
    arma::Col<double> d_inv(s.n_elem);
    for (arma::uword i = 0; i < s.n_elem; ++i) {
        d_inv(i) = (s(i) > tol) ? (1.0 / s(i)) : 0.0;
    }
    arma::Mat<std::complex<double>> diag_dinv =
        arma::diagmat(arma::conv_to<arma::Col<std::complex<double>>>::from(d_inv));
    arma::Col<std::complex<double>> A = V * diag_dinv * U.t() * b;

    std::vector<std::complex<T>> out(size);
    for (int i = 0; i < size; ++i) {
        out[i] = std::complex<T>(
            static_cast<T>(A(i).real()),
            static_cast<T>(A(i).imag())
        );
    }
    return out;
}

#ifdef __GNUC__
template<>
std::vector<std::complex<__float128>> solve_fluid_mode_system<__float128>(
    const std::vector<std::complex<__float128>>& K3_kernel,
    const std::vector<std::complex<__float128>>& rhs,
    int size
) {
    return solve_linear_system_lup_refined<__float128>(K3_kernel, rhs, size);
}
#endif

template<typename T>
std::complex<T> psms_backscatter_fluid_adaptive(
    int m_max,
    int n_max,
    T chi_sw,
    T chi_body,
    T xi,
    T eta_body,
    T density_body,
    T density_sw,
    const std::vector<T>& nodes,
    const std::vector<T>& weights,
    bool normalize = true,
    bool adaptive = false
) {
    std::complex<T> fbs_sum(0, 0), c(0, 0);
    T max_band = T(0);
    T prev_band = T(-1);
    int n_small_bands = 0;
    T x1_r = xi - T(1);
    int iopnorm = normalize ? 1 : 0;
    std::vector<T> ext_args(nodes.size() + 1);
    ext_args[0] = eta_body;
    for (size_t k = 0; k < nodes.size(); ++k) {
        ext_args[k + 1] = nodes[k];
    }

    // Batch the exterior and interior spheroidal evaluations over small blocks
    // of m so the high setup cost is shared across nearby orders.
    int chunk_size = psms_profcn_mblock_chunk_size<T>();
    for (int m_start = 0; m_start <= m_max; m_start += chunk_size) {
        int m_count = std::min(chunk_size, m_max - m_start + 1);
        int lnum = n_max - m_start + 1;
        auto ext_result_batch = cprofcn_mblock<T>(
            chi_sw, m_start, m_count, lnum, ext_args, 2, iopnorm, 2, x1_r
        );
        auto int_result_batch = cprofcn_mblock<T>(
            chi_body, m_start, m_count, lnum, nodes, 1, iopnorm, 2, x1_r
        );

        for (int local_m = 0; local_m < m_count; ++local_m) {
            int m = m_start + local_m;
            int size = n_max - m + 1;
            if (size <= 0) {
                continue;
            }

            if (adaptive) {
                T max_proxy = T(0);
                T prev_proxy = T(-1);
                int n_small_proxy = 0;
                int full_size = size;
                for (int i = 0; i < full_size; ++i) {
                    T smn_val = extract_angular_value_from_mblock<T>(
                        ext_result_batch, local_m, i, 0
                    );
                    if (is_na_real(smn_val)) {
                        continue;
                    }
                    auto rv1 = extract_radial_from_mblock<T>(ext_result_batch, local_m, i, 1);
                    auto rv3 = extract_radial_from_mblock<T>(ext_result_batch, local_m, i, 3);
                    auto rv_int = extract_radial_from_mblock<T>(int_result_batch, local_m, i, 1);
                    T ext_scale = std::max(
                        precabs(rv1.val_real),
                        std::max(
                            precsqrt(rv3.val_real * rv3.val_real + rv3.val_imag * rv3.val_imag),
                            precabs(rv_int.val_real)
                        )
                    );
                    T proxy = precabs(smn_val) * std::max(
                        T(1),
                        std::max(ext_scale, precabs(rv_int.der_real))
                    );
                    max_proxy = std::max(max_proxy, proxy);
                    if (psms_modal_term_is_negligible<T>(
                        i, full_size - 1, proxy, prev_proxy, max_proxy
                    )) {
                        ++n_small_proxy;
                    } else {
                        n_small_proxy = 0;
                    }
                    prev_proxy = proxy;
                    if (n_small_proxy >= 2) {
                        size = std::max(1, i - 1);
                        break;
                    }
                }
            }

            std::vector<T> smn_body(size, T(0));
            std::vector<T> ext_inc_val(size, T(0));
            std::vector<T> ext_inc_der(size, T(0));
            std::vector<std::complex<T>> ext_scat_val(size, std::complex<T>(0, 0));
            std::vector<std::complex<T>> ext_scat_der(size, std::complex<T>(0, 0));
            std::vector<T> int_val(size, T(0));
            std::vector<T> int_der(size, T(0));
            std::vector<std::vector<T>> smn_ext(size, std::vector<T>(nodes.size(), T(0)));
            std::vector<std::vector<T>> smn_int(size, std::vector<T>(nodes.size(), T(0)));

            for (int i = 0; i < size; ++i) {
                T smn_val = extract_angular_value_from_mblock<T>(
                    ext_result_batch, local_m, i, 0
                );
                smn_body[i] = smn_val;

                if (!(is_na_real(smn_val) || smn_val == T(0))) {
                    auto rv1 = extract_radial_from_mblock<T>(ext_result_batch, local_m, i, 1);
                    auto rv3 = extract_radial_from_mblock<T>(ext_result_batch, local_m, i, 3);
                    ext_inc_val[i] = rv1.val_real;
                    ext_inc_der[i] = rv1.der_real;
                    ext_scat_val[i] = std::complex<T>(rv3.val_real, rv3.val_imag);
                    ext_scat_der[i] = std::complex<T>(rv3.der_real, rv3.der_imag);
                }

                auto rv_int = extract_radial_from_mblock<T>(int_result_batch, local_m, i, 1);
                int_val[i] = rv_int.val_real;
                int_der[i] = rv_int.der_real;

                for (size_t k = 0; k < nodes.size(); ++k) {
                    smn_ext[i][k] = extract_angular_value_from_mblock<T>(
                        ext_result_batch, local_m, i, static_cast<int>(k + 1)
                    );
                    smn_int[i][k] = extract_angular_value_from_mblock<T>(
                        int_result_batch, local_m, i, static_cast<int>(k)
                    );
                }
            }

            // The overlap matrix couples interior and exterior degrees for the
            // current order m.
            std::vector<std::vector<T>> coef_mat(size, std::vector<T>(size, T(0)));
            // Form the dense right-hand side and scattered kernel for this m.
            for (int li = 0; li < size; ++li) {
                for (int ni = 0; ni < size; ++ni) {
                    T sum = T(0);
                    for (size_t k = 0; k < nodes.size(); ++k) {
                        sum += smn_ext[ni][k] * smn_int[li][k] * weights[k];
                    }
                    coef_mat[li][ni] = sum;
                }
            }

            std::vector<std::complex<T>> rhs(size, std::complex<T>(0, 0));
            std::vector<std::complex<T>> K3(size * size, std::complex<T>(0, 0));
            std::vector<std::complex<T>> jn_smn(size, std::complex<T>(0, 0));
            for (int ni = 0; ni < size; ++ni) {
                int n = m + ni;
                T smn_val = smn_body[ni];
                if (!is_na_real(smn_val) && smn_val != T(0)) {
                    jn_smn[ni] = imaginary_unit_power<T>(n) * smn_val;
                }
            }

            for (int li = 0; li < size; ++li) {
                T factor = (int_der[li] != T(0))
                    ? (density_body * int_val[li]) / (density_sw * int_der[li])
                    : T(0);
                std::complex<T> row_sum(0, 0);
                for (int ni = 0; ni < size; ++ni) {
                    if (jn_smn[ni] == std::complex<T>(0, 0)) {
                        continue;
                    }
                    std::complex<T> E1 =
                        std::complex<T>(ext_inc_val[ni] - factor * ext_inc_der[ni], T(0));
                    std::complex<T> E3 = ext_scat_val[ni] - factor * ext_scat_der[ni];
                    T alpha = coef_mat[li][ni];
                    row_sum += jn_smn[ni] * alpha * E1;
                    K3[li * size + ni] = jn_smn[ni] * alpha * E3;
                }
                rhs[li] = -row_sum;
            }

            auto Amn = solve_fluid_mode_system<T>(K3, rhs, size);
            T nu_m = (m == 0) ? T(1) : T(2);
            std::complex<T> band_sum(0, 0), band_c(0, 0);
            T max_term = T(0);
            T prev_term = T(-1);
            int n_small_terms = 0;
            for (int i = 0; i < size; ++i) {
                int n = m + i;
                T smn_val = smn_body[i];
                auto amn = Amn[i];
                if (is_na_real(smn_val) ||
                    is_na_real(amn.real()) || is_na_real(amn.imag())) {
                    continue;
                }
                T parity_n = ((n % 2) == 0) ? T(1) : T(-1);
                std::complex<T> term = nu_m * parity_n * smn_val * smn_val * amn;
                T term_abs = precsqrt(complex_abs_sq(term));
                max_term = std::max(max_term, term_abs);
                if (adaptive &&
                    psms_modal_term_is_negligible<T>(i, size - 1, term_abs, prev_term, max_term)) {
                    ++n_small_terms;
                } else {
                    n_small_terms = 0;
                }
                prev_term = term_abs;
                std::complex<T> y_band = term - band_c;
                std::complex<T> t_band = band_sum + y_band;
                band_c = (t_band - band_sum) - y_band;
                band_sum = t_band;

                if (adaptive && n_small_terms >= 2) {
                    break;
                }
            }

            std::complex<T> y = band_sum - c;
            std::complex<T> t = fbs_sum + y;
            c = (t - fbs_sum) - y;
            fbs_sum = t;

            T band_abs = precsqrt(complex_abs_sq(band_sum));
            max_band = std::max(max_band, band_abs);
            if (adaptive &&
                psms_modal_band_is_negligible<T>(m, m_max, band_abs, prev_band, max_band)) {
                ++n_small_bands;
            } else {
                n_small_bands = 0;
            }
            prev_band = band_abs;
            if (adaptive && n_small_bands >= 2) {
                return fbs_sum;
            }
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
    T phi_body,
    T phi_scatter,
    T density_body,
    T density_sw,
    const std::vector<T>& nodes,
    const std::vector<T>& weights,
    const std::string& Amn_method,
    bool adaptive = false
) {
    T eta_body = preccos(theta_body);
    bool is_backscatter =
        precabs(preccos(phi_body - phi_scatter) + T(1)) <= T(64) * preceps<T>();

    if (Amn_method == "Amn_fixed_rigid") {
        return psms_backscatter_fixed_rigid<T>(
            m_max, n_max, chi_sw, xi, eta_body, adaptive
        );
    }

    if (Amn_method == "Amn_pressure_release") {
        return psms_backscatter_pressure_release<T>(
            m_max, n_max, chi_sw, xi, eta_body, adaptive
        );
    }

#ifdef __GNUC__
    if constexpr (std::is_same_v<T, __float128>) {
        if (is_backscatter && Amn_method == "Amn_fluid") {
            return psms_backscatter_fluid_adaptive<T>(
                m_max, n_max, chi_sw, chi_body, xi, eta_body,
                density_body, density_sw, nodes, weights, true, adaptive
            );
        }
    }
#endif

    // Determine appropriate Amn expansion matrix computation
    std::vector<std::vector<std::complex<T>>> Amn;
    std::vector<std::vector<T>> smn_body_matrix;

    // Fluid-filled
    if (Amn_method == "Amn_fluid") {
        auto fluid_result = fluid_Amn_triangular<T>(
            m_max, n_max, chi_sw, chi_body, xi,
            eta_body, density_body, density_sw, nodes, weights
        );
        smn_body_matrix = std::move(fluid_result.smn_body);
        if (is_backscatter) {
            // Exact backscatter can stay in the triangular modal storage and
            // avoid expanding to a rectangular matrix.
            return compute_fbs_backscatter_triangular<T>(
                m_max, smn_body_matrix, fluid_result.Amn_tri, adaptive
            );
        }
        Amn = expand_Amn_triangular<T>(m_max, n_max, fluid_result.Amn_tri);
    // Simplified fluid-filled
    } else if (Amn_method == "Amn_fluid_simplify") {
        auto external_boundary = external_boundary_matrices<T>(
            m_max, n_max, chi_sw, xi, eta_body, true
        );
        smn_body_matrix = external_boundary.smn;
        auto Amn_tri = simplified_fluid_Amn_triangular<T>(
            m_max, n_max, chi_sw, chi_body, xi, density_body, density_sw,
            external_boundary.radial
        );
        if (is_backscatter) {
            return compute_fbs_backscatter_triangular<T>(
                m_max, smn_body_matrix, Amn_tri, adaptive
            );
        }
        Amn = expand_Amn_triangular<T>(m_max, n_max, Amn_tri);
    // Fixed rigid
    } else if (Amn_method == "Amn_fixed_rigid") {
        auto external_boundary = external_boundary_matrices<T>(
            m_max, n_max, chi_sw, xi, eta_body, true
        );
        smn_body_matrix = external_boundary.smn;
        Amn = fixed_rigid_Amn_expansion_matrix<T>(
            m_max, n_max, chi_sw, xi
        );
    // Pressure-release
    } else {
        auto external_boundary = external_boundary_matrices<T>(
            m_max, n_max, chi_sw, xi, eta_body, true
        );
        smn_body_matrix = external_boundary.smn;
        Amn = pressure_release_Amn_expansion_matrix<T>(
            m_max, n_max, chi_sw, xi
        );
    }

    // The exported PSMS pathway is primarily backscatter. In that case the
    // angular parity collapses the scattered-angle factor to the exact sign
    // (-1)^n, which is cheaper and numerically cleaner than rebuilding the
    // second angular matrix explicitly.
    if (is_backscatter) {
        return compute_fbs_backscatter<T>(
            m_max, n_max, smn_body_matrix, Amn, adaptive
        );
    }

    // Compute azimuth angle factors
    auto azimuth = compute_azimuth<T>(m_max, phi_body, phi_scatter);
    auto smn_scatter_matrix = reflect_smn_matrix<T>(smn_body_matrix);

    // Compute linear scattering coefficient
    return compute_fbs<T>(
        m_max, n_max, azimuth, smn_body_matrix, smn_scatter_matrix, Amn
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
    std::string Amn_method = "Amn_fluid",
    bool adaptive = false
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
                Amn_method, adaptive
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
                Amn_method, adaptive
            );
            f_bs_out[f] = to_Rcomplex(fbs);
        }
    }
    return Rcpp::wrap(f_bs_out);
}
