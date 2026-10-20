#pragma once
#include <climits>
#include <cstddef>
#include <cstdint>
#include <array>
#include <cstring>
#include <new>

inline constexpr size_t psms_retained_cache_budget_bytes = 64 * 1024 * 1024;

struct ProfcnWorkspaceDimensions {
    int narg = 0, maxm = 0, maxint = 0, maxj = 0, maxp = 0, maxn = 0;
    int maxpdr = 0, maxq = 0, maxdr = 0, maxd = 0, maxlp = 0, maxmp = 0;
    int ngau = 200, maxt = 1;
    size_t radial_size = 0, angular_size = 0;
};

inline int checked_profcn_integer(long long value, const char* expression) {
    if (value < 0 || value > INT_MAX)
        throw std::length_error(std::string("Spheroidal native integer range exceeded by ") + expression + ".");
    return static_cast<int>(value);
}

template<typename T> int checked_profcn_real_integer(T value, const char* expression) {
    // INT truncates toward zero, so test its result range before any native conversion.
    if (!is_finite_real(value) || value <= T(INT_MIN) - T(1) || value >= T(INT_MAX) + T(1))
        throw std::length_error(std::string("Spheroidal native integer conversion exceeded by ") + expression + ".");
    return static_cast<int>(value);
}

inline size_t checked_profcn_product(size_t first, size_t second, size_t maximum, const char* expression) {
    if (second != 0 && first > maximum / second)
        throw std::length_error(std::string("Spheroidal workspace range exceeded by ") + expression + ".");
    return first * second;
}

template<typename T> T profcn_log10(T value) {
#if ACOUSTICTS_HAVE_QUADMATH
    if constexpr (std::is_same_v<T, acousticts_quad_t>) return log10q(value);
    else
#endif
    return std::log10(value);
}

template<typename T> ProfcnWorkspaceDimensions preflight_profcn_dimensions(
        T c, int m_start, int m_count, int lnum, size_t arguments,
        int ioprad, int iopang, T x1) {
    if (m_start < 0 || m_count < 1 || lnum < 1 || arguments > size_t(INT_MAX) ||
            ioprad < 0 || ioprad > 2 || iopang < 0 || iopang > 2 ||
            !is_finite_real(c) || c < T(0) || (ioprad != 0 &&
            (!is_finite_real(x1) || x1 < T(0))) || (ioprad == 2 && x1 == T(0)))
        throw std::invalid_argument("Spheroidal dimensions or radial inputs are outside the native domain.");
    ProfcnWorkspaceDimensions dims;
    dims.narg = static_cast<int>(std::max(size_t(1), arguments));
    dims.maxm = checked_profcn_integer(static_cast<long long>(m_start) + m_count - 1, "maximum order");
    const int digits = is_quad_precision_v<T> ? 33 : 15;
    const int integer_c = checked_profcn_real_integer(c, "int(c)");
    const long long maximum_order = dims.maxm;
    // These products occur in geteig and pleg before conversion to native real.
    checked_profcn_integer(maximum_order * maximum_order, "m*m in geteig");
    if (iopang != 0 || ioprad == 2)
        checked_profcn_integer(2 * maximum_order * maximum_order + 1, "2*m*m+1 in pleg");
    checked_profcn_integer(2 * (maximum_order + lnum - 1) + 1, "twice the maximum degree plus one");
    dims.maxint = checked_profcn_integer(static_cast<long long>(lnum) + 3 * digits + integer_c + 5, "maxint");
    dims.maxj = checked_profcn_integer(static_cast<long long>(dims.maxint) + maximum_order, "maxj");
    dims.maxp = dims.maxint;
    dims.maxn = dims.maxj;
    dims.maxpdr = 4 * digits + 5;
    if (ioprad == 2) {
        const int degrees = std::max(checked_profcn_integer(static_cast<long long>(lnum) + maximum_order,
            "lnum+maxm"), 64);
        const T logarithm = profcn_log10(x1);
        // Match every executed workspace conversion, including subsequently replaced estimates.
        if (x1 >= T(13) / T(20000))
            dims.maxn = checked_profcn_real_integer(T(2) * (T(degrees) * (-T(37) / T(2) - T(20) * logarithm) +
                T(5 * digits) + T(4 * maximum_order) + c + T(5000)) + T(maximum_order) + T(5), "maxn near xi=1");
        if (x1 > T(2) / T(25))
            dims.maxn = checked_profcn_real_integer(T(2) * (T(degrees) * (T(1) / T(2) - T(3) * logarithm) +
                T(5 * digits) + T(4 * maximum_order) + c + T(1000)) + T(maximum_order) + T(5), "maxn intermediate xi");
        if (x1 > T(1))
            dims.maxn = checked_profcn_real_integer(T(2) * (T(degrees) * T(0.5) + T(5 * digits) +
                T(4 * maximum_order) + c + T(500)) + T(maximum_order) + T(5), "maxn large xi");
        dims.maxn = checked_profcn_integer(dims.maxn, "maxn");
        dims.maxp = std::max(dims.maxn, dims.maxp);
        if (x1 < T(1) / T(1000))
            dims.ngau = checked_profcn_integer(200LL - 50LL * checked_profcn_real_integer(
                logarithm - pow10_typed<T>(-30), "log10(x1) for quadrature"), "ngau");
        if (x1 < pow10_typed<T>(-10))
            dims.ngau = checked_profcn_integer(250LL - 50LL * checked_profcn_real_integer(
                logarithm - pow10_typed<T>(-30), "log10(x1) for quadrature"), "ngau");
        if (x1 < pow10_typed<T>(-11)) dims.ngau = 1200;
        if (x1 < pow10_typed<T>(-12)) dims.ngau = 2400;
        if (x1 <= T(0.5))
            dims.maxpdr = checked_profcn_integer(static_cast<long long>(dims.maxpdr) +
                checked_profcn_real_integer(T(2) * c + T(100) * x1, "radial Legendre support") + 400, "maxpdr");
        // The integration path forms this window count before visiting any quadrature nodes.
        if (c > T(2000))
            checked_profcn_real_integer(T(10) * (c / T(2000)) * (c / T(2000)), "radial integration windows");
    }
    dims.maxq = checked_profcn_integer(static_cast<long long>(dims.maxint) + 2 * maximum_order, "maxq");
    dims.maxdr = dims.maxpdr / 2 + 1;
    dims.maxp = std::max(dims.maxp, dims.maxpdr);
    dims.maxd = dims.maxp / 2 + 1;
    checked_profcn_integer(2LL * dims.maxd, "normalization-factor support");
    dims.maxlp = checked_profcn_integer(static_cast<long long>(lnum) + maximum_order + 5, "maxlp");
    dims.maxmp = checked_profcn_integer(maximum_order + 5, "maxmp");
    if (iopang != 0) dims.maxt = dims.narg;
    if (ioprad != 0 && x1 != T(0)) {
        const T bessel_argument = c * precsqrt(x1 * (x1 + T(2)));
        checked_profcn_real_integer(bessel_argument + bessel_argument + T(20), "spherical Bessel start index");
    }
    // Output offsets include default-integer products in the legacy scalar accessors.
    dims.radial_size = checked_profcn_product(size_t(lnum), size_t(m_count), size_t(INT_MAX), "radial output");
    dims.angular_size = checked_profcn_product(dims.radial_size, size_t(dims.narg), size_t(INT_MAX), "angular output");
    const size_t addressable_reals = size_t(std::numeric_limits<std::ptrdiff_t>::max()) / sizeof(T);
    checked_profcn_product(size_t(dims.maxt), size_t(dims.maxp), addressable_reals, "Legendre workspace");
    checked_profcn_product(size_t(10), size_t(dims.maxn), addressable_reals, "radial recurrence workspace");
    checked_profcn_product(size_t(10), size_t(dims.maxp), addressable_reals, "Legendre recurrence workspace");
    return dims;
}

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
    return precnan<T>();
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

inline ProfcnResultDouble cprofcn_double(
    double c,
    int m,
    int lnum,
    const std::vector<double>& arg,
    int ioprad,
    int iopnorm,
    int iopang,
    double x1
) {
    const auto dimensions = preflight_profcn_dimensions(c, m, 1, lnum, arg.size(), ioprad, iopang, x1);
    int narg = dimensions.narg;

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
    
    result.s1c.resize(dimensions.angular_size, 0.0);
    result.s1dc.resize(dimensions.angular_size, 0.0);
    result.is1e.resize(dimensions.angular_size, 0);
    result.is1de.resize(dimensions.angular_size, 0);
    result.naccs.resize(dimensions.angular_size, 0);
    
    std::vector<double> arg_c = arg;
    arg_c.resize(narg, 0.0);

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

inline ProfcnResultQuad cprofcn_quad(
    acousticts_quad_t c,
    int m,
    int lnum,
    const std::vector<acousticts_quad_t>& arg,
    int ioprad,
    int iopnorm,
    int iopang,
    acousticts_quad_t x1
) {
    const auto dimensions = preflight_profcn_dimensions(c, m, 1, lnum, arg.size(), ioprad, iopang, x1);
    int narg = dimensions.narg;

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
    
    result.s1c.resize(dimensions.angular_size, 0.0);
    result.s1dc.resize(dimensions.angular_size, 0.0);
    result.is1e.resize(dimensions.angular_size, 0);
    result.is1de.resize(dimensions.angular_size, 0);
    result.naccs.resize(dimensions.angular_size, 0);
    
    std::vector<acousticts_quad_t> arg_c = arg;
    arg_c.resize(narg, acousticts_quad_t(0));

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
    const auto dimensions = preflight_profcn_dimensions(c, m_start, m_count, lnum,
        arg.size(), ioprad, iopang, x1);
    int narg = dimensions.narg;
    ProfcnBatchResult<double> result;
    result.lnum = lnum;
    result.narg = narg;
    result.m_count = m_count;
    size_t radial_size = dimensions.radial_size;
    size_t angular_size = dimensions.angular_size;

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
    arg_c.resize(narg, 0.0);
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
    const auto dimensions = preflight_profcn_dimensions(c, m_start, m_count, lnum,
        arg.size(), ioprad, iopang, x1);
    int narg = dimensions.narg;
    ProfcnBatchResult<acousticts_quad_t> result;
    result.lnum = lnum;
    result.narg = narg;
    result.m_count = m_count;
    size_t radial_size = dimensions.radial_size;
    size_t angular_size = dimensions.angular_size;

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
    arg_c.resize(narg, acousticts_quad_t(0));
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
struct ProfcnCoefficientBatchResult {
    std::uint64_t cache_identity = 0;
    ProfcnBatchResult<T> functions;
    int capacity = 0;
    std::vector<T> coefficients;
    std::vector<int> lengths;
};

template<typename T>
ProfcnCoefficientBatchResult<T> cprofcn_mblock_coefficients(
    T c, int m_start, int m_count, int lnum, const std::vector<T>& arguments,
    int ioprad, int iopang, T x1, bool packed_radial = false
) {
    const auto dimensions = preflight_profcn_dimensions(c, m_start, m_count, lnum,
        arguments.size(), ioprad, iopang, x1);
    // Check native indexing before allocating any of the output rectangles.
    const int initial_capacity = dimensions.maxint;
    const int argument_count = dimensions.narg;
    const size_t radial_size = dimensions.radial_size;
    if (radial_size > size_t(INT_MAX) / size_t(initial_capacity) ||
            radial_size > size_t(INT_MAX) / size_t(argument_count)) {
        throw std::length_error("Spheroidal coefficient workspace exceeds the native indexing range.");
    }
    ProfcnCoefficientBatchResult<T> output;
    auto& result = output.functions;
    result.lnum = lnum;
    result.narg = argument_count;
    result.m_count = m_count;
    const size_t angular_size = radial_size * result.narg;
    // Allocate the unchanged native radial and angular output layout.
    result.r1c.resize(radial_size);
    result.r1dc.resize(radial_size);
    result.r2c.resize(radial_size);
    result.r2dc.resize(radial_size);
    result.ir1e.resize(radial_size);
    result.ir1de.resize(radial_size);
    result.ir2e.resize(radial_size);
    result.ir2de.resize(radial_size);
    result.naccr.resize(radial_size);
    result.s1c.resize(angular_size);
    result.s1dc.resize(angular_size);
    result.is1e.resize(angular_size);
    result.is1de.resize(angular_size);
    result.naccs.resize(angular_size);
    std::vector<T> args = arguments;
    args.resize(result.narg, T(0));
    int iopnorm = 1, status = 0;
    // Match the native coefficient workspace, retrying any larger reported support.
    output.capacity = initial_capacity;
    output.lengths.resize(radial_size);
    for (;;) {
        if (radial_size > size_t(INT_MAX) / size_t(output.capacity) ||
                radial_size > output.coefficients.max_size() / size_t(output.capacity)) {
            throw std::length_error("Spheroidal coefficient workspace exceeds the native indexing range.");
        }
        output.coefficients.assign(radial_size * output.capacity, T(0));
        if constexpr (std::is_same_v<T, double>) {
            profcn_cpp_interface_batch_coeff(
                &c, &m_start, &m_count, &lnum, &ioprad, &x1, &iopang, &iopnorm,
                &result.narg, args.data(), result.r1c.data(), result.ir1e.data(),
                result.r1dc.data(), result.ir1de.data(), result.r2c.data(),
                result.ir2e.data(), result.r2dc.data(), result.ir2de.data(),
                result.naccr.data(), result.s1c.data(), result.is1e.data(),
                result.s1dc.data(), result.is1de.data(), result.naccs.data(),
                &output.capacity, output.coefficients.data(), output.lengths.data(), &status
            );
        } else {
#if ACOUSTICTS_HAVE_QUADMATH
            profcn_cpp_interface_batch_coeff_quad(
                &c, &m_start, &m_count, &lnum, &ioprad, &x1, &iopang, &iopnorm,
                &result.narg, args.data(), result.r1c.data(), result.ir1e.data(),
                result.r1dc.data(), result.ir1de.data(), result.r2c.data(),
                result.ir2e.data(), result.r2dc.data(), result.ir2de.data(),
                result.naccr.data(), result.s1c.data(), result.is1e.data(),
                result.s1dc.data(), result.is1de.data(), result.naccs.data(),
                &output.capacity, output.coefficients.data(), output.lengths.data(), &status
            );
#else
            throw std::runtime_error("Native quad spheroidal coefficients are unavailable.");
#endif
        }
        if (status == 0) {
            if (std::any_of(output.lengths.begin(), output.lengths.end(),
                    [&](int length) { return length < 1 || length > output.capacity; })) {
                throw std::runtime_error("Native spheroidal coefficient support is invalid.");
            }
            break;
        }
        const int required = *std::max_element(output.lengths.begin(), output.lengths.end());
        if (status != 1 || required <= output.capacity) {
            throw std::runtime_error("Native spheroidal expansion coefficients could not be resolved.");
        }
        output.capacity = required;
    }
    // Restore the existing mantissa/exponent representation in native precision.
    if (!packed_radial) {
        scale_profcn_component(result.r1c, result.ir1e);
        scale_profcn_component(result.r1dc, result.ir1de);
        scale_profcn_component(result.r2c, result.ir2e);
        scale_profcn_component(result.r2dc, result.ir2de);
    }
    scale_profcn_component(result.s1c, result.is1e);
    scale_profcn_component(result.s1dc, result.is1de);
    return output;
}

template<typename T>
T spheroidal_coefficient_overlap(
    const ProfcnCoefficientBatchResult<T>& exterior,
    const ProfcnCoefficientBatchResult<T>& interior,
    int local_order, int external_degree_index, int internal_degree_index
) {
    // Apply Ferrers expansion and orthogonality (DLMF 30.8.1 and 14.17.6).
    if ((external_degree_index - internal_degree_index) % 2 != 0) return T(0);
    const size_t external_index = external_degree_index +
        static_cast<size_t>(exterior.functions.lnum) * local_order;
    const size_t internal_index = internal_degree_index +
        static_cast<size_t>(interior.functions.lnum) * local_order;
    const int length = std::min(exterior.lengths[external_index],
        interior.lengths[internal_index]);
    const int parity = external_degree_index % 2;
#if ACOUSTICTS_HAVE_EXACT_QUAD_DOT
    if constexpr (std::is_same_v<T, acousticts_quad_t>) {
        // Round the complete coefficient dot product once in native binary128.
        static thread_local quad_dot_product::Accumulator accumulator;
        const size_t count = length > parity ? (size_t(length - parity) + 1) / 2 : 0;
        return quad_dot_product::dot(
            exterior.coefficients.data() + exterior.capacity * external_index + parity,
            interior.coefficients.data() + interior.capacity * internal_index + parity,
            count, accumulator, 2, 2);
    }
#endif
    T sum = 0, correction = 0;
    for (int r = parity; r < length; r += 2) {
        const T product = exterior.coefficients[r + exterior.capacity * external_index] *
            interior.coefficients[r + interior.capacity * internal_index];
        const T next = product - correction;
        const T updated = sum + next;
        correction = (updated - sum) - next;
        sum = updated;
    }
    return sum;
}

template<typename T> struct SpheroidalCoefficientOverlapCache {
    struct Columns {
        std::vector<size_t> offsets;
        std::vector<T> values;

        size_t bytes() const {
            return offsets.capacity() * sizeof(size_t) + values.capacity() * sizeof(T);
        }

        void capture(const ProfcnCoefficientBatchResult<T>& batch, int parity, int count) {
            offsets.resize(size_t(count) + 1);
            for (int col = 0; col < count; ++col)
                offsets[col + 1] = offsets[col] + size_t(batch.lengths[parity + 2 * col]);
            values.resize(offsets.back());
            for (int col = 0; col < count; ++col) {
                const int index = parity + 2 * col;
                std::memcpy(values.data() + offsets[col],
                    batch.coefficients.data() + size_t(batch.capacity) * index,
                    (offsets[col + 1] - offsets[col]) * sizeof(T));
            }
        }

        bool identical(const Columns& other, int col) const {
            if (size_t(col) + 1 >= offsets.size() || size_t(col) + 1 >= other.offsets.size()) return false;
            const size_t length = offsets[col + 1] - offsets[col];
            return length == other.offsets[col + 1] - other.offsets[col] &&
                std::memcmp(values.data() + offsets[col], other.values.data() + other.offsets[col],
                    length * sizeof(T)) == 0;
        }
    };

    struct Block {
        int columns = 0, rows = 0;
        Columns exterior, interior;
        std::vector<T> overlaps;

        size_t bytes() const {
            return exterior.bytes() + interior.bytes() + overlaps.capacity() * sizeof(T);
        }
    };

    std::array<Block, 2> blocks;
    int order = -1;
    size_t budget, peak_bytes = 0, reused_entries = 0, computed_entries = 0;

    explicit SpheroidalCoefficientOverlapCache(size_t maximum_bytes = psms_retained_cache_budget_bytes)
        : budget(maximum_bytes) {}

    size_t bytes() const { return sizeof(*this) + blocks[0].bytes() + blocks[1].bytes(); }

    void clear() {
        for (auto& block : blocks) block = Block();
        order = -1;
    }

    static bool reserve_bytes(size_t count, size_t width, size_t limit, size_t& used) {
        if (used > limit || count > (limit - used) / width) return false;
        used += count * width;
        return true;
    }

    static bool column_bytes(const ProfcnCoefficientBatchResult<T>& batch, int parity,
            int count, size_t limit, size_t& used) {
        if (!reserve_bytes(size_t(count) + 1, sizeof(size_t), limit, used)) return false;
        for (int col = 0; col < count; ++col)
            if (!reserve_bytes(size_t(batch.lengths[parity + 2 * col]), sizeof(T), limit, used)) return false;
        return true;
    }

    void fill(int m, int parity, const ProfcnCoefficientBatchResult<T>& exterior,
            const ProfcnCoefficientBatchResult<T>& interior, std::vector<T>& overlaps) {
        if (m < 0 || parity < 0 || parity > 1 || exterior.functions.lnum < 1 ||
                interior.functions.lnum < 1 || exterior.functions.m_count != 1 || interior.functions.m_count != 1)
            throw std::invalid_argument("Coefficient-overlap reuse requires one valid azimuthal order.");
        const int columns = int((static_cast<long long>(exterior.functions.lnum) + 1 - parity) / 2);
        const int rows = int((static_cast<long long>(interior.functions.lnum) + 1 - parity) / 2);
        const size_t entries = checked_profcn_product(size_t(rows), size_t(columns),
            overlaps.max_size(), "coefficient-overlap cache rectangle");
        if (overlaps.size() != entries)
            throw std::invalid_argument("Coefficient-overlap output has an inconsistent rectangle.");
        auto validate = [parity](const ProfcnCoefficientBatchResult<T>& batch, int count) {
            if (batch.capacity < 1 || batch.lengths.size() < size_t(batch.functions.lnum) ||
                    batch.coefficients.size() / size_t(batch.capacity) < size_t(batch.functions.lnum))
                throw std::invalid_argument("Coefficient-overlap input has inconsistent native storage.");
            for (int col = 0; col < count; ++col) {
                const int length = batch.lengths[parity + 2 * col];
                if (length < 1 || length > batch.capacity)
                    throw std::invalid_argument("Coefficient-overlap input has invalid native support.");
            }
        };
        validate(exterior, columns);
        validate(interior, rows);
        if (order != m) { clear(); order = m; }
        peak_bytes = bytes();
        reused_entries = computed_entries = 0;
        Block candidate;
        std::vector<unsigned char> same_exterior, same_interior;
        size_t required = sizeof(*this) + blocks[1 - parity].bytes();
        bool retain = column_bytes(exterior, parity, columns, budget, required) &&
            column_bytes(interior, parity, rows, budget, required) &&
            reserve_bytes(entries, sizeof(T), budget, required);
        // Retention is optional and never reduces the complete coefficient product.
        if (retain) {
            try {
                candidate.columns = columns;
                candidate.rows = rows;
                candidate.exterior.capture(exterior, parity, columns);
                candidate.interior.capture(interior, parity, rows);
                candidate.overlaps.resize(entries);
                same_exterior.resize(columns);
                same_interior.resize(rows);
                peak_bytes = bytes() + sizeof(candidate) + candidate.bytes() +
                    sizeof(same_exterior) + sizeof(same_interior) +
                    same_exterior.capacity() + same_interior.capacity();
                retain = sizeof(*this) + blocks[1 - parity].bytes() + candidate.bytes() <= budget;
            } catch (const std::bad_alloc&) {
                peak_bytes = bytes() + sizeof(candidate) + candidate.bytes() +
                    sizeof(same_exterior) + sizeof(same_interior) +
                    same_exterior.capacity() + same_interior.capacity();
                retain = false;
            }
        }
        if (!retain) {
            clear();
            candidate = Block();
            std::vector<unsigned char>().swap(same_exterior);
            std::vector<unsigned char>().swap(same_interior);
        } else {
            // Compare every exported coefficient and its length before reusing an intersection.
            for (int col = 0; col < columns; ++col)
                same_exterior[col] = candidate.exterior.identical(blocks[parity].exterior, col);
            for (int row = 0; row < rows; ++row)
                same_interior[row] = candidate.interior.identical(blocks[parity].interior, row);
        }
        for (int row = 0; row < rows; ++row) {
            for (int col = 0; col < columns; ++col) {
                const size_t index = size_t(row) * columns + col;
                if (retain && same_interior[row] && same_exterior[col]) {
                    overlaps[index] = blocks[parity].overlaps[size_t(row) * blocks[parity].columns + col];
                    ++reused_entries;
                } else {
                    overlaps[index] = spheroidal_coefficient_overlap(exterior, interior,
                        0, parity + 2 * col, parity + 2 * row);
                    ++computed_entries;
                }
                if (!is_finite_real(overlaps[index]))
                    throw std::runtime_error("Angular overlap is nonfinite.");
                if (retain) candidate.overlaps[index] = overlaps[index];
            }
        }
        if (retain) blocks[parity] = std::move(candidate);
    }
};

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
    return precnan<T>();
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
