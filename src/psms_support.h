#pragma once

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
    // Exact powers of i appear throughout the spheroidal modal sums. Evaluate
    // them algebraically instead of with std::pow.
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
    // profcn stores angular values plus base-10 exponents. Apply the exponent
    // immediately so callers always see the physical-scale quantity.
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

template<typename T>
inline void scale_profcn_component(std::vector<T>& values, std::vector<int>& exponents) {
    // Normalize mantissa/exponent pairs in place so downstream code can treat
    // the arrays as ordinary floating-point values.
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

// ============================================================================
// PROLATE SPHEROIDAL WAVE FUNCTION - PROFCN [DOUBLE]
// --------------------------------------------------
// This adds a safety net by padding the vector sizes to avoid segmentation 
// faults that can crash the program or the Rcpp interface
// ============================================================================
struct ProfcnResultDouble {
    // Direct mirror of the double-precision profcn output arrays.
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
#if ACOUSTICTS_HAVE_QUADMATH
struct ProfcnResultQuad {
    // Quad-precision mirror of the profcn output arrays.
    std::vector<acousticts_quad_t> r1c, r1dc, r2c, r2dc;
    std::vector<int> ir1e, ir1de, ir2e, ir2de, naccr;
    std::vector<acousticts_quad_t> s1c, s1dc;
    std::vector<int> is1e, is1de, naccs;
};

ProfcnResultQuad cprofcn_quad(
    acousticts_quad_t c,
    int m,
    int lnum,
    const std::vector<acousticts_quad_t>& arg,
    int ioprad,
    int iopnorm,
    int iopang,
    acousticts_quad_t x1
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
    
    std::vector<acousticts_quad_t> arg_c = arg;

    // Convert scalars to quad precision
    acousticts_quad_t c_q = c;
    acousticts_quad_t x1_q = x1;
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
    // Batched profcn layout for a block of neighboring m values.
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
    double x1,
    bool scale = true
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

    if (scale) {
        scale_profcn_component(result.r1c, result.ir1e);
        scale_profcn_component(result.r1dc, result.ir1de);
        scale_profcn_component(result.r2c, result.ir2e);
        scale_profcn_component(result.r2dc, result.ir2de);
        scale_profcn_component(result.s1c, result.is1e);
        scale_profcn_component(result.s1dc, result.is1de);
    }
    return result;
}

#if ACOUSTICTS_HAVE_QUADMATH
inline ProfcnBatchResult<acousticts_quad_t> cprofcn_batch_quad(
    acousticts_quad_t c,
    int m_start,
    int m_count,
    int lnum,
    const std::vector<acousticts_quad_t>& arg,
    int ioprad,
    int iopnorm,
    int iopang,
    acousticts_quad_t x1
) {
    int narg = std::max<int>(arg.size(), 1);
    ProfcnBatchResult<acousticts_quad_t> result;
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

    std::vector<acousticts_quad_t> arg_c = arg;
    acousticts_quad_t c_q = c;
    acousticts_quad_t x1_q = x1;
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
    // Radial terms are packed as contiguous m-blocks.
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
    // Angular terms are packed as (degree, argument, order-within-block).
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
    // Default chunk size for double precision.
    return 8;
}

template<>
inline int psms_profcn_mblock_chunk_size<acousticts_quad_t>() {
    // Smaller quad chunks limit temporary-buffer growth.
    return 4;
}

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
    // Convert any mantissa/exponent storage returned by profcn into plain
    // values so the rest of the code can ignore exponent bookkeeping.
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

template<>
inline ProfcnBatchResult<acousticts_quad_t> cprofcn_mblock<acousticts_quad_t>(
    acousticts_quad_t c,
    int m_start,
    int m_count,
    int lnum,
    const std::vector<acousticts_quad_t>& arg,
    int ioprad,
    int iopnorm,
    int iopang,
    acousticts_quad_t x1
) {
#if ACOUSTICTS_HAVE_QUADMATH
    return cprofcn_batch_quad(c, m_start, m_count, lnum, arg, ioprad, iopnorm, iopang, x1);
#else
    std::vector<double> arg_double(arg.size());
    for (size_t i = 0; i < arg.size(); ++i) {
        arg_double[i] = quad_to_double(arg[i]);
    }

    auto raw = cprofcn_batch_double(
        quad_to_double(c),
        m_start,
        m_count,
        lnum,
        arg_double,
        ioprad,
        iopnorm,
        iopang,
        quad_to_double(x1),
        false
    );

    ProfcnBatchResult<acousticts_quad_t> out;
    out.lnum = raw.lnum;
    out.narg = raw.narg;
    out.m_count = raw.m_count;
    out.ir1e = raw.ir1e;
    out.ir1de = raw.ir1de;
    out.ir2e = raw.ir2e;
    out.ir2de = raw.ir2de;
    out.naccr = raw.naccr;
    out.is1e = raw.is1e;
    out.is1de = raw.is1de;
    out.naccs = raw.naccs;
    out.r1c.reserve(raw.r1c.size());
    out.r1dc.reserve(raw.r1dc.size());
    out.r2c.reserve(raw.r2c.size());
    out.r2dc.reserve(raw.r2dc.size());
    out.s1c.reserve(raw.s1c.size());
    out.s1dc.reserve(raw.s1dc.size());
    for (double value : raw.r1c) out.r1c.push_back(double_to_quad(value));
    for (double value : raw.r1dc) out.r1dc.push_back(double_to_quad(value));
    for (double value : raw.r2c) out.r2c.push_back(double_to_quad(value));
    for (double value : raw.r2dc) out.r2dc.push_back(double_to_quad(value));
    for (double value : raw.s1c) out.s1c.push_back(double_to_quad(value));
    for (double value : raw.s1dc) out.s1dc.push_back(double_to_quad(value));
    scale_profcn_component(out.r1c, out.ir1e);
    scale_profcn_component(out.r1dc, out.ir1de);
    scale_profcn_component(out.r2c, out.ir2e);
    scale_profcn_component(out.r2dc, out.ir2de);
    scale_profcn_component(out.s1c, out.is1e);
    scale_profcn_component(out.s1dc, out.is1de);
    return out;
#endif
}

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
    // Double-precision dispatcher for the scalar profcn wrapper.
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
template<>
inline ProfcnResult<acousticts_quad_t> cprofcn<acousticts_quad_t>(
    acousticts_quad_t c, int m, int lnum,
    const std::vector<acousticts_quad_t>& arg,
    int ioprad, int iopnorm, int iopang, acousticts_quad_t x1
) {
#if ACOUSTICTS_HAVE_QUADMATH
    // Quad-precision dispatcher for the scalar profcn wrapper.
    ProfcnResultQuad raw = cprofcn_quad(c, m, lnum, arg, ioprad, iopnorm, iopang, x1);

    ProfcnResult<acousticts_quad_t> out;
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
#else
    std::vector<double> arg_double(arg.size());
    for (size_t i = 0; i < arg.size(); ++i) {
        arg_double[i] = quad_to_double(arg[i]);
    }

    auto raw = cprofcn<double>(
        quad_to_double(c),
        m,
        lnum,
        arg_double,
        ioprad,
        iopnorm,
        iopang,
        quad_to_double(x1)
    );

    ProfcnResult<acousticts_quad_t> out;
    out.ir1e = raw.ir1e;
    out.ir1de = raw.ir1de;
    out.ir2e = raw.ir2e;
    out.ir2de = raw.ir2de;
    out.naccr = raw.naccr;
    out.is1e = raw.is1e;
    out.is1de = raw.is1de;
    out.naccs = raw.naccs;
    out.r1c.reserve(raw.r1c.size());
    out.r1dc.reserve(raw.r1dc.size());
    out.r2c.reserve(raw.r2c.size());
    out.r2dc.reserve(raw.r2dc.size());
    out.s1c.reserve(raw.s1c.size());
    out.s1dc.reserve(raw.s1dc.size());
    for (double value : raw.r1c) out.r1c.push_back(double_to_quad(value));
    for (double value : raw.r1dc) out.r1dc.push_back(double_to_quad(value));
    for (double value : raw.r2c) out.r2c.push_back(double_to_quad(value));
    for (double value : raw.r2dc) out.r2dc.push_back(double_to_quad(value));
    for (double value : raw.s1c) out.s1c.push_back(double_to_quad(value));
    for (double value : raw.s1dc) out.s1dc.push_back(double_to_quad(value));
    return out;
#endif
}

// ============================================================================
