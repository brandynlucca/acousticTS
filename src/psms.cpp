// ============================================================================
// INCLUDES
// ============================================================================
// [[Rcpp::depends(RcppArmadillo, BH)]]
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
#include "precision.h"

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
    void profcn_cpp_interface_batch_coeff(
        double* c, int* m_start, int* m_count, int* lnum, int* ioprad, double* x1,
        int* iopang, int* iopnorm, int* narg, double* arg, double* r1c, int* ir1e,
        double* r1dc, int* ir1de, double* r2c, int* ir2e, double* r2dc, int* ir2de,
        int* naccr, double* s1c, int* is1e, double* s1dc, int* is1de, int* naccs,
        int* capacity, double* coefficients, int* lengths, int* status
    );

    // Quad precision
#if ACOUSTICTS_HAVE_QUADMATH
    // Quad precision interface
    void profcn_cpp_interface_quad(
        acousticts_quad_t* c, int* m, int* lnum, int* ioprad, acousticts_quad_t* x1,
        int* iopang, int* iopnorm, int* narg, acousticts_quad_t* arg,
        acousticts_quad_t* r1c, int* ir1e, acousticts_quad_t* r1dc, int* ir1de,
        acousticts_quad_t* r2c, int* ir2e, acousticts_quad_t* r2dc, int* ir2de,
        int* naccr, acousticts_quad_t* s1c, int* is1e, acousticts_quad_t* s1dc,
        int* is1de, int* naccs
    );
    void profcn_cpp_interface_batch_quad(
        acousticts_quad_t* c, int* m_start, int* m_count, int* lnum, int* ioprad,
        acousticts_quad_t* x1, int* iopang, int* iopnorm, int* narg,
        acousticts_quad_t* arg, acousticts_quad_t* r1c, int* ir1e,
        acousticts_quad_t* r1dc, int* ir1de, acousticts_quad_t* r2c, int* ir2e,
        acousticts_quad_t* r2dc, int* ir2de, int* naccr, acousticts_quad_t* s1c,
        int* is1e, acousticts_quad_t* s1dc, int* is1de, int* naccs
    );
    void profcn_cpp_interface_batch_coeff_quad(
        acousticts_quad_t* c, int* m_start, int* m_count, int* lnum, int* ioprad,
        acousticts_quad_t* x1, int* iopang, int* iopnorm, int* narg,
        acousticts_quad_t* arg, acousticts_quad_t* r1c, int* ir1e,
        acousticts_quad_t* r1dc, int* ir1de, acousticts_quad_t* r2c, int* ir2e,
        acousticts_quad_t* r2dc, int* ir2de, int* naccr, acousticts_quad_t* s1c,
        int* is1e, acousticts_quad_t* s1dc, int* is1de, int* naccs,
        int* capacity, acousticts_quad_t* coefficients, int* lengths, int* status
    );

    // Quad precision math functions
    acousticts_quad_t powq(acousticts_quad_t, acousticts_quad_t);
#endif
}
template<typename T>
struct ProfcnResult;

template<typename T>
struct ProfcnBatchResult;

// Internal PSMS support layer: special-function wrappers, profcn batching,
// exponent rescaling helpers, and precision-agnostic scaffolding.
#include "psms_support.h"

#if ACOUSTICTS_HAVE_QUADMATH
// Keep unsupported arithmetic and rounding modes on the original integral path.
extern "C" int pint_exact_sums_available_quad() {
#if ACOUSTICTS_HAVE_EXACT_QUAD_DOT
    return std::fegetround() == FE_TONEAREST;
#else
    return 0;
#endif
}

static acousticts_quad_t pint_original_dot_quad(const acousticts_quad_t* basis,
        const acousticts_quad_t* kernel, int nodes, int basis_stride) {
    acousticts_quad_t sum = 0;
    for (int node = 0; node < nodes; ++node)
        sum += basis[size_t(node) * basis_stride] * kernel[size_t(node) * 3];
    return sum;
}

static acousticts_quad_t pint_exact_dot_quad(const acousticts_quad_t* basis,
        const acousticts_quad_t* kernel, int nodes, int basis_stride) {
#if ACOUSTICTS_HAVE_EXACT_QUAD_DOT
    if (std::fegetround() == FE_TONEAREST) {
        static thread_local quad_dot_product::Accumulator accumulator;
        try {
            accumulator.reset();
            bool finite = true;
            for (int node = 0; node < nodes; ++node) {
                if (!accumulator.add(basis[size_t(node) * basis_stride], kernel[size_t(node) * 3])) {
                    finite = false;
                    break;
                }
            }
            if (finite) return accumulator.rounded();
        } catch (...) {
            // Preserve native arithmetic and its existing failure checks on fallback.
        }
    }
#endif
    return pint_original_dot_quad(basis, kernel, nodes, basis_stride);
}

extern "C" void pint_exact_sums_quad(const acousticts_quad_t* basis,
        const acousticts_quad_t* kernels, int degree_capacity, int nodes,
        int first, int last, int extra, acousticts_quad_t* first_integrals,
        acousticts_quad_t* second_integrals, acousticts_quad_t* fourth_integrals) {
    // Accumulate the same stored binary128 products and round each moment once.
    for (int degree = first - 1; degree < last; degree += 2) {
        first_integrals[degree] = pint_exact_dot_quad(basis + degree, kernels, nodes, degree_capacity);
        second_integrals[degree] = pint_exact_dot_quad(basis + degree, kernels + 1, nodes, degree_capacity);
        fourth_integrals[degree] = pint_exact_dot_quad(basis + degree, kernels + 2, nodes, degree_capacity);
    }
    fourth_integrals[extra - 1] = pint_exact_dot_quad(
        basis + extra - 1, kernels + 2, nodes, degree_capacity);
}

extern "C" void profcn_set_integral_sums_quad(
    int (*available)(), decltype(&pint_exact_sums_quad) sum);

namespace {
struct RegisterIntegralSums {
    RegisterIntegralSums() {
        // Register optional acceleration without adding a standalone Fortran link dependency.
        profcn_set_integral_sums_quad(&pint_exact_sums_available_quad, &pint_exact_sums_quad);
    }
};
const RegisterIntegralSums registered_integral_sums;
}
#endif


#if ACOUSTICTS_HAVE_QUADMATH
namespace psms_checked_integrals {
static thread_local bool enabled = false;
}

extern "C" int pint_set_checked_moments_quad(int enabled) {
    const int previous = psms_checked_integrals::enabled;
    psms_checked_integrals::enabled = enabled != 0;
    return previous;
}

extern "C" int pint_checked_moments_available_quad() {
#if ACOUSTICTS_HAVE_EXACT_QUAD_DOT && !defined(__FAST_MATH__) && \
    !defined(__ASSOCIATIVE_MATH__) && !defined(__RECIPROCAL_MATH__) && \
    (!defined(__FINITE_MATH_ONLY__) || __FINITE_MATH_ONLY__ == 0) && \
    (!defined(__GCC_IEC_559) || __GCC_IEC_559 > 0)
    return psms_checked_integrals::enabled && std::numeric_limits<double>::is_iec559 &&
        std::numeric_limits<double>::digits == 53 && std::numeric_limits<double>::max_exponent == 1024 &&
        std::fegetround() == FE_TONEAREST;
#else
    return 0;
#endif
}

#if ACOUSTICTS_HAVE_EXACT_QUAD_DOT
namespace psms_checked_integrals {
// The Taylor identity follows DLMF 10.51.2 and Jacobi Rodrigues integration.
// DLMF 10.50.3 supplies the discrete Green factor for the scaled Neumann recurrence.
// Bounds concern the real expression at fixed stored native inputs. Values and signed
// products remain binary128. Only outward positive bounds use scaled binary64 mantissas.
using Quad = acousticts_quad_t;
using quad_dot_product::Accumulator;

// Binary128 unit roundoff, and the smallest double step used for outward rounding.
static const Quad unit_roundoff = ldexpq(Quad(1), -113);
static const double positive_upper_factor = 1.0 + 0x1p-52, positive_lower_factor = 1.0 - 0x1p-52;

// ---- Bit-level binary128 helpers ------------------------------------------
inline void words(Quad x, std::uint64_t w[2]) { std::memcpy(w, &x, 16); }
inline Quad fromw(std::uint64_t lo, std::uint64_t hi) { const std::uint64_t w[2]{lo, hi}; Quad x; std::memcpy(&x, w, 16); return x; }
// A rounded result is normal exactly when its exponent field is neither 0 nor 0x7fff.
inline bool normal_quad(Quad x) { std::uint64_t w[2]; words(x, w); return ((w[1] >> 48) & 0x7fff) - 1u < 0x7ffeu; }
inline bool finite_quad(Quad x) { std::uint64_t w[2]; words(x, w); return ((w[1] >> 48) & 0x7fff) != 0x7fff; }
// Next value toward +inf by integer increment of the sign-magnitude encoding.
inline Quad next_upper(Quad x) {
    std::uint64_t w[2]; words(x, w);
    if (((w[1] >> 48) & 0x7fff) == 0x7fff) return x;
    if ((w[1] << 1 | w[0]) == 0) return fromw(1, 0);
    if (w[1] >> 63) { if (w[0]-- == 0) --w[1]; }
    else { if (++w[0] == 0) ++w[1]; }
    return fromw(w[0], w[1]);
}
inline Quad next_lower(Quad x) { return -next_upper(-x); }
inline Quad quad_abs(Quad x) { std::uint64_t w[2]; words(x, w); w[1] &= ~(UINT64_C(1) << 63); return fromw(w[0], w[1]); }

// ---- Outward positive arithmetic with a free integer exponent --------------
// A positive mantissa times an integer power of two. Results are upper bounds unless named lower.
struct PositiveBound { double mantissa; int exponent; };
inline double pow2d(int d) { const std::uint64_t b = std::uint64_t(d + 1023) << 52; double v; std::memcpy(&v, &b, 8); return v; }
inline PositiveBound pnorm_(double m, int e) { if (m >= 2.0) { m *= 0.5; e += 1; } return {m, e}; }
inline PositiveBound pmul(PositiveBound a, PositiveBound b) { return pnorm_(a.mantissa * b.mantissa * positive_upper_factor, a.exponent + b.exponent); }
// Lower bound for a*b, used only to build denominators.
inline PositiveBound pmul_lo(PositiveBound a, PositiveBound b) { return pnorm_(a.mantissa * b.mantissa * positive_lower_factor, a.exponent + b.exponent); }
inline PositiveBound padd(PositiveBound a, PositiveBound b) {
    if (a.mantissa == 0) return b;
    if (b.mantissa == 0) return a;
    if (a.exponent < b.exponent) std::swap(a, b);
    const int d = a.exponent - b.exponent;
    const double s = d < 1000 ? a.mantissa + b.mantissa * pow2d(-d) : a.mantissa + 0x1p-998;
    return pnorm_(s * positive_upper_factor, a.exponent);
}
// Upper bound for a/b where b.mantissa>0 is a lower bound of the true denominator.
inline PositiveBound pdiv(PositiveBound a, PositiveBound b) {
    double m = a.mantissa / b.mantissa * positive_upper_factor; int e = a.exponent - b.exponent;
    if (m == 0) return {0, 0};
    if (m < 1.0) { m *= 2.0; e -= 1; }
    return {m, e};
}
// Exact outward conversion of a nonzero binary128 subnormal magnitude.
inline PositiveBound psubnormal(std::uint64_t lo, std::uint64_t hi, bool upper) {
    const __uint128_t v = (static_cast<__uint128_t>(hi & UINT64_C(0xffffffffffff)) << 64) | lo;
    const std::uint64_t top = std::uint64_t(v >> 64);
    const int bits = top ? 128 - __builtin_clzll(top) : 64 - __builtin_clzll(lo);
    const int shift = std::max(0, bits - 53);
    std::uint64_t kept = std::uint64_t(v >> shift);
    if (upper && shift > 0) ++kept;
    return pnorm_(double(kept) * pow2d(-(bits - shift - 1)), bits - 1 - 16494);
}
inline PositiveBound pfromq_up(Quad x, bool& ok) {
    std::uint64_t w[2]; words(x, w);
    const unsigned E = unsigned(w[1] >> 48) & 0x7fff;
    if (E == 0x7fff) { ok = false; return {0, 0}; }
    if (E == 0) return ((w[1] << 16) | w[0]) == 0 ? PositiveBound{0, 0} : psubnormal(w[0], w[1], true);
    const std::uint64_t f = (w[1] & UINT64_C(0xffffffffffff)) | (UINT64_C(1) << 48);
    return pnorm_(double(f + 1) * 0x1p-48, int(E) - 16383);
}
// A zero lower bound stays 0, which every division rejects.
inline PositiveBound pfromq_dn(Quad x, bool& ok) {
    std::uint64_t w[2]; words(x, w);
    const unsigned E = unsigned(w[1] >> 48) & 0x7fff;
    if (E == 0x7fff) { ok = false; return {0, 0}; }
    if (E == 0) return ((w[1] << 16) | w[0]) == 0 ? PositiveBound{0, 0} : psubnormal(w[0], w[1], false);
    const std::uint64_t f = (w[1] & UINT64_C(0xffffffffffff)) | (UINT64_C(1) << 48);
    return {double(f) * 0x1p-48, int(E) - 16383};
}
inline PositiveBound pfromd(double x) {
    if (x == 0) return {0, 0};
    int e; double m = std::frexp(x, &e) * 2.0;
    return {m, e - 1};
}
inline PositiveBound pint(unsigned long long n) { return pfromd(double(n) * (n >> 53 ? positive_upper_factor : 1.0)); }
// Upper double, saturating to +inf above and to 2^-1000 below.
inline double pdbl_up(PositiveBound a) {
    if (a.mantissa == 0) return 0;
    if (a.exponent > 1000) return HUGE_VAL;
    if (a.exponent < -1000) return 0x1p-999;
    return a.mantissa * pow2d(a.exponent);
}
inline bool ple(PositiveBound a, PositiveBound b) {
    if (a.mantissa == 0) return true;
    if (b.mantissa == 0) return false;
    const int d = a.exponent - b.exponent;
    if (d < -2) return true;
    if (d > 2) return false;
    return a.mantissa * pow2d(d) <= b.mantissa;
}
inline PositiveBound psqrt(PositiveBound a) {
    if (a.mantissa == 0) return a;
    double m = a.mantissa; int e = a.exponent;
    if (e & 1) { m *= 2.0; e -= 1; }
    return pnorm_(std::sqrt(m) * positive_upper_factor, e / 2);
}
// Exact conversion to binary128 storage. Reject values above the format.
inline Quad ptoq_up(PositiveBound a, bool& ok) {
    if (a.mantissa == 0) return 0;
    if (!(a.mantissa < 4.0)) { ok = false; return 0; }
    double m = a.mantissa; int e = a.exponent;
    if (m >= 2.0) { m *= 0.5; e += 1; }
    if (e > 16382) { ok = false; return 0; }
    if (e < -16381) return fromw(0, UINT64_C(2) << 48);
    std::uint64_t b; std::memcpy(&b, &m, 8);
    b &= UINT64_C(0xfffffffffffff);
    return fromw(b << 60, (std::uint64_t(e + 16383) << 48) | (b >> 4));
}

// ---- Outward binary128 interval helpers used only in setup and fallback -----
struct Interval { Quad lower, upper; };
inline Interval ivadd(Interval a, Interval b) { return {next_lower(a.lower + b.lower), next_upper(a.upper + b.upper)}; }
inline Interval ivsub(Interval a, Interval b) { return {next_lower(a.lower - b.upper), next_upper(a.upper - b.lower)}; }
inline Interval ivmul(Interval a, Interval b) {
    const Quad p[4]{a.lower * b.lower, a.lower * b.upper, a.upper * b.lower, a.upper * b.upper};
    Quad lo = p[0], hi = p[0];
    for (int i = 0; i < 4; ++i) { if (!finite_quad(p[i])) return {-HUGE_VALQ, HUGE_VALQ}; if (p[i] < lo) lo = p[i]; if (p[i] > hi) hi = p[i]; }
    return {next_lower(lo), next_upper(hi)};
}
inline Interval ivdivp(Interval a, Interval b) {
    if (!(b.lower > 0) || !finite_quad(b.upper)) return {-HUGE_VALQ, HUGE_VALQ};
    return ivmul(a, Interval{next_lower(1 / b.upper), next_upper(1 / b.lower)});
}

// ---- Rigorous relative bound for r rounded steps ---------------------------
// Q_S=S/(1-S). The lemma needs S<1, so S>=1 returns failure.
struct CoefficientBounds {
    std::vector<PositiveBound> q; std::vector<double> onep; PositiveBound sigma; bool ok = true;
    void init(PositiveBound sig, int cap) {
        sigma = sig; q.assign(cap + 1, PositiveBound{0, 0}); onep.assign(cap + 1, 1.0);
        for (int r = 1; r <= cap; ++r) {
            const PositiveBound S = pmul(pint(r), sig);
            const double sd = pdbl_up(S), den = (1.0 - sd) * positive_lower_factor;
            if (!(den > 0)) { ok = false; return; }
            q[r] = pdiv(S, pfromd(den));
            onep[r] = (1.0 / den) * positive_upper_factor;
        }
    }
};
inline PositiveBound relq(PositiveBound S, bool& ok) {
    const double den = (1.0 - pdbl_up(S)) * positive_lower_factor;
    if (!(den > 0)) { ok = false; return {0, 0}; }
    return pdiv(S, pfromd(den));
}

// Sign/cos on a binary128 point with outward radii. The series and rounding model replace
// the interval recurrence. Alternating remainders are DLMF 4.19 with |x|<=1.
// Fix the nearest binary128 value independently of the rounding mode at DLL load.
static const Quad half_pi = fromw(UINT64_C(0x8469898cc51701b8), UINT64_C(0x3fff921fb54442d1));
static bool validated_sincos(Quad z, Quad& sval, Quad& cval, PositiveBound& srad, PositiveBound& crad) {
    bool ok = true;
    const Quad quot = roundq(z / half_pi);
    if (fabsq(quot) > Quad(1e18)) return false;
    const long long j = (long long)quot;
    const Quad prod = quot * half_pi, x = z - prod;
    if (fabsq(x) > 1 || !normal_quad(x)) return false;
    const PositiveBound pu{1.0, -113};
    // |pi/2-stored|<=2^-112, the verified adjacent-endpoint enclosure.
    PositiveBound dr = padd(padd(pmul(pint(j < 0 ? -j : j), PositiveBound{1.0, -112}), pmul(pfromq_up(prod, ok), pu)), pmul(pfromq_up(x, ok), pu));
    const Quad x2 = x * x;
    if (!normal_quad(x2)) return false;
    const PositiveBound stepS = pmul(pint(3), pu);
    Quad ts = x, tc = 1, ss = x, sc = 1;
    PositiveBound es{0, 0}, ec{0, 0};
    PositiveBound lasts = pfromq_up(ts, ok), lastc = PositiveBound{1.0, 0};
    for (int j2 = 1; j2 <= 24; ++j2) {
        const Quad ds = Quad(2 * j2 * (2 * j2 + 1)), dc = Quad((2 * j2 - 1) * 2 * j2);
        ts = (ts * x2) / ds; tc = (tc * x2) / dc;
        if (!normal_quad(ts) || !normal_quad(tc)) return false;
        const PositiveBound S = pmul(pint(j2), stepS), qr = relq(S, ok);
        const PositiveBound as = pfromq_up(ts, ok), ac = pfromq_up(tc, ok);
        if (j2 & 1) { ss -= ts; sc -= tc; } else { ss += ts; sc += tc; }
        es = padd(es, padd(pmul(as, qr), pmul(pfromq_up(ss, ok), pu)));
        ec = padd(ec, padd(pmul(ac, qr), pmul(pfromq_up(sc, ok), pu)));
        lasts = as; lastc = ac;
    }
    // One further term dominates the alternating truncation remainder.
    srad = padd(padd(es, lasts), dr); crad = padd(padd(ec, lastc), dr);
    switch (((j % 4) + 4) % 4) {
        case 0: sval = ss; cval = sc; break;
        case 1: sval = sc; cval = -ss; std::swap(srad, crad); break;
        case 2: sval = -ss; cval = -sc; break;
        default: sval = -sc; cval = ss; std::swap(srad, crad); break;
    }
    return ok;
}

// Interval route with the same series, for arguments outside the normal contract.
static bool interval_sincos(Quad z, Quad& sval, Quad& cval, PositiveBound& srad, PositiveBound& crad) {
    bool ok = true;
    const Interval pib{next_lower(half_pi), next_upper(half_pi)};
    const Quad quot = roundq(z / half_pi);
    if (fabsq(quot) > Quad(1e18)) return false;
    const long long j = (long long)quot;
    const Interval x = ivsub(Interval{z, z}, ivmul(Interval{quot, quot}, pib));
    if (!finite_quad(x.lower) || !finite_quad(x.upper) || fmaxq(fabsq(x.lower), fabsq(x.upper)) > 1) return false;
    Interval sq = ivmul(x, x);
    if (x.lower <= 0 && x.upper >= 0) sq.lower = 0;
    const Interval negsq{-sq.upper, -sq.lower};
    Interval ts = x, tc{1, 1}, ss = x, sc{1, 1};
    for (int i = 1; i <= 24; ++i) {
        const Quad ds = Quad(2 * i * (2 * i + 1)), dc = Quad((2 * i - 1) * 2 * i);
        ts = ivdivp(ivmul(ts, negsq), Interval{ds, ds}); ss = ivadd(ss, ts);
        tc = ivdivp(ivmul(tc, negsq), Interval{dc, dc}); sc = ivadd(sc, tc);
    }
    const Quad rs = fmaxq(fabsq(ts.lower), fabsq(ts.upper)), rc = fmaxq(fabsq(tc.lower), fabsq(tc.upper));
    const Interval s2{next_lower(ss.lower - rs), next_upper(ss.upper + rs)}, c2{next_lower(sc.lower - rc), next_upper(sc.upper + rc)};
    const Quad smid = (s2.lower + s2.upper) / 2, cmid = (c2.lower + c2.upper) / 2;
    const Quad sr = next_upper(fmaxq(next_upper(smid - s2.lower), next_upper(s2.upper - smid))), cr = next_upper(fmaxq(next_upper(cmid - c2.lower), next_upper(c2.upper - cmid)));
    srad = pfromq_up(sr, ok); crad = pfromq_up(cr, ok);
    switch (((j % 4) + 4) % 4) {
        case 0: sval = smid; cval = cmid; break;
        case 1: sval = cmid; cval = -smid; std::swap(srad, crad); break;
        case 2: sval = -smid; cval = -cmid; break;
        default: sval = -cmid; cval = smid; std::swap(srad, crad); break;
    }
    return ok;
}

// ---- Green table in outward positive arithmetic ----------------------------
struct MomentRecurrence {
    int nmax = 0, last = -1, avail = -1;
    Quad c, q, s, origin, root, z, h, inverse;
    Interval zb, hb, rootb, originb;
    PositiveBound zL, zU, dz, green, asum, bsum, eps_h, eps_root, eps_origin, abss;
    PositiveBound pref_step;
    bool inverse_normal = true;
    std::vector<Quad> y;
    std::vector<PositiveBound> ay, E, Y, J, Bm;
    CoefficientBounds rel;
    Accumulator defect_acc;
};

static PositiveBound h_norm(int order, PositiveBound zLinv2) {
    PositiveBound term = zLinv2, total = term;
    for (int k = 0; k < order; ++k) {
        const double ratio = (double(2 * k + 1) / double(2 * k + 2)) * positive_upper_factor * double(order - k) * positive_upper_factor * double(order + k + 1) * positive_upper_factor;
        term = pmul(term, pmul(pfromd(ratio), zLinv2));
        total = padd(total, term);
    }
    return psqrt(total);
}

static bool initialize_recurrence(MomentRecurrence& T, Quad c, Quad q, Quad s, int center, int nmax, int cap) {
    bool ok = true;
    if (std::fegetround() != FE_TONEAREST) return false;
    T.c = c; T.q = q; T.s = s; T.nmax = nmax;
    if (!(c > 0) || !(q > 0) || s == 0 || !finite_quad(c) || !finite_quad(q) || !finite_quad(s)) return false;
    T.origin = q + Quad(center);
    T.originb = ivadd(Interval{q, q}, Interval{Quad(center), Quad(center)});
    if (!(T.originb.lower > 1) || !finite_quad(T.originb.upper)) return false;
    T.root = sqrtq(T.origin);
    if (!(T.root > 0) || !finite_quad(T.root)) return false;
    // A residual encloses the root without an accuracy assumption on sqrt.
    const Interval r2 = ivmul(Interval{T.root, T.root}, Interval{T.root, T.root});
    const Interval res = ivsub(T.originb, r2);
    const Quad rerr = next_upper(fmaxq(fabsq(res.lower), fabsq(res.upper)) / T.root);
    T.rootb = Interval{next_lower(T.root - rerr), next_upper(T.root + rerr)};
    if (!(T.rootb.lower > 0)) return false;
    T.z = c * T.root;
    T.h = c / (2 * T.root);
    T.zb = ivmul(Interval{c, c}, T.rootb);
    T.hb = ivdivp(Interval{c, c}, ivmul(Interval{2, 2}, T.rootb));
    if (!(T.zb.lower > 0) || !finite_quad(T.zb.upper) || !(T.hb.lower > 0) || !finite_quad(T.hb.upper)) return false;
    T.inverse = 1 / T.z;
    T.zL = pfromq_dn(T.zb.lower, ok); T.zU = pfromq_up(T.zb.upper, ok);
    const Quad dzq = next_upper(fmaxq(fabsq(T.z - T.zb.lower), fabsq(T.zb.upper - T.z)));
    T.dz = pfromq_up(dzq, ok);
    if (T.zL.mantissa == 0) return false;
    T.abss = pfromq_up(quad_abs(s), ok);
    const PositiveBound absL = pfromq_dn(quad_abs(s), ok);
    if (absL.mantissa == 0) return false;
    T.green = pdiv(pmul(T.zU, T.zU), absL);
    // |hhat/h-1|<=max(hhat-hL,hU-hhat)/hL.
    const Quad hrad = next_upper(fmaxq(next_upper(T.h - T.hb.lower), next_upper(T.hb.upper - T.h)));
    T.eps_h = pdiv(pfromq_up(hrad, ok), pfromq_dn(T.hb.lower, ok));
    const PositiveBound pu{1.0, -113};
    T.eps_root = pdiv(pfromq_up(rerr, ok), pfromq_dn(T.rootb.lower, ok));
    T.eps_origin = pdiv(pfromq_up(next_upper(fmaxq(fabsq(T.origin - T.originb.lower), fabsq(T.originb.upper - T.origin))), ok), pfromq_dn(T.originb.lower, ok));
    T.rel.init(padd(T.eps_h, pmul(pint(3), pu)), cap);
    if (!T.rel.ok) return false;
    T.pref_step = padd(T.eps_h, pmul(pint(3), pu));
    // Seed values use the same library sin/cos as the native y table, checked against validated values.
    Quad sl, cl, sv, cv; PositiveBound srad, crad;
    sincosq(T.z, &sl, &cl);
    if (!validated_sincos(T.z, sv, cv, srad, crad)) {
        if (!interval_sincos(T.z, sv, cv, srad, crad)) return false;
    }
    const PositiveBound Se = padd(padd(pfromq_up(quad_abs(sl - sv), ok), pmul(pfromq_up(quad_abs(sl - sv), ok), pu)), padd(srad, T.dz));
    const PositiveBound Ce = padd(padd(pfromq_up(quad_abs(cl - cv), ok), pmul(pfromq_up(quad_abs(cl - cv), ok), pu)), padd(crad, T.dz));
    T.y.assign(nmax + 1, Quad(0));
    T.ay.assign(nmax + 1, PositiveBound{0, 0}); T.E = T.ay; T.Y = T.ay; T.J = T.ay; T.Bm = T.ay;
    T.inverse_normal = normal_quad(T.inverse);
    const Quad a0 = cl * T.inverse;
    T.y[0] = -a0 * s;
    const Quad b1 = sl * T.inverse, a1 = a0 * T.inverse, d1 = -a1 - b1;
    T.y[1] = d1 * s;
    const PositiveBound zL2 = pmul_lo(T.zL, T.zL), acl = pfromq_up(quad_abs(cl), ok), asl = pfromq_up(quad_abs(sl), ok);
    const PositiveBound Q2 = relq(pmul(pint(2), pu), ok), Q3 = relq(pmul(pint(3), pu), ok), Q4 = relq(pmul(pint(4), pu), ok);
    T.ay[0] = pfromq_up(quad_abs(T.y[0]), ok);
    T.ay[1] = pfromq_up(quad_abs(T.y[1]), ok);
    PositiveBound e0 = padd(pmul(T.ay[0], Q3), pmul(T.abss, padd(pdiv(Ce, T.zL), pdiv(pmul(acl, T.dz), zL2))));
    const PositiveBound zL4 = pmul_lo(zL2, zL2), two{1.0, 1};
    const PositiveBound inner = padd(padd(pmul(pfromq_up(quad_abs(d1), ok), pu), pmul(pfromq_up(quad_abs(a1), ok), Q4)),
        padd(pmul(pfromq_up(quad_abs(b1), ok), Q2), padd(pdiv(Ce, zL2), padd(pdiv(pmul(pmul(acl, pmul(two, T.zU)), T.dz), zL4),
        padd(pdiv(Se, T.zL), pdiv(pmul(asl, T.dz), zL2))))));
    PositiveBound e1 = padd(pmul(T.ay[1], pu), pmul(T.abss, inner));
    const bool seed_normal = T.inverse_normal && normal_quad(cl) && normal_quad(sl) && normal_quad(a0) && normal_quad(T.y[0]) && normal_quad(a1) && normal_quad(b1) && normal_quad(T.y[1]);
    if (!seed_normal) {
        // Non-normal seeds take the absolute interval route.
        const Quad Ceq = ptoq_up(Ce, ok), Seq = ptoq_up(Se, ok);
        const Interval cosb{next_lower(cl - Ceq), next_upper(cl + Ceq)}, sinb{next_lower(sl - Seq), next_upper(sl + Seq)};
        const Interval un0 = ivdivp(Interval{-cosb.upper, -cosb.lower}, T.zb), un1 = ivdivp(ivsub(un0, sinb), T.zb);
        const Interval sb{s, s}, y0b = ivmul(un0, sb), y1b = ivmul(un1, sb);
        if (!finite_quad(y0b.lower) || !finite_quad(y0b.upper) || !finite_quad(y1b.lower) || !finite_quad(y1b.upper)) return false;
        e0 = pfromq_up(next_upper(fmaxq(fabsq(T.y[0] - y0b.lower), fabsq(y0b.upper - T.y[0]))), ok);
        e1 = pfromq_up(next_upper(fmaxq(fabsq(T.y[1] - y1b.lower), fabsq(y1b.upper - T.y[1]))), ok);
    }
    T.E[0] = e0; T.E[1] = e1;
    T.Y[0] = padd(T.ay[0], e0); T.Y[1] = padd(T.ay[1], e1);
    // Sparse order anchors from DLMF 10.49.17-18, then the backward ratio sweep above z.
    const PositiveBound zLinv2 = pdiv(PositiveBound{1.0, 0}, zL2);
    const double zUd = pdbl_up(T.zU);
    int n0 = nmax;
    if (zUd < double(nmax - 1)) n0 = int(std::ceil(zUd)) + 1;
    int previous = -1;
    for (int sample = 0; sample <= 8; ++sample) {
        const int anchor = (sample * n0 + 7) / 8;
        const PositiveBound b = h_norm(anchor, zLinv2);
        for (int i = previous + 1; i <= anchor; ++i) T.J[i] = b;
        previous = anchor;
    }
    if (n0 < nmax) {
        std::vector<double> ratios(nmax + 1, 1.0);
        for (int i = nmax; i > n0; --i) {
            const double den = (double(2 * i + 1) * positive_lower_factor - zUd * ratios[i] * positive_upper_factor) * positive_lower_factor;
            if (!(den > 0)) return false;
            ratios[i - 1] = std::min(1.0, (zUd / den) * positive_upper_factor);
        }
        for (int i = n0 + 1; i <= nmax; ++i) T.J[i] = pmul(T.J[i - 1], pfromd(ratios[i - 1]));
    }
    T.asum = padd(pmul(T.J[1], e0), pmul(T.J[0], e1));
    T.bsum = padd(pmul(T.Y[1], e0), pmul(T.Y[0], e1));
    T.Bm[0] = pdiv(T.abss, T.zL);
    T.Bm[1] = padd(T.Bm[0], pdiv(T.Bm[0], T.zL));
    T.last = 1; T.avail = 1;
    return ok;
}

// Enclose each Neumann recurrence defect with exact products before propagating it.
static inline bool extend_recurrence(MomentRecurrence& T, int target) {
    bool ok = true;
    while (T.last < target) {
        const int n = T.last + 1, i = n - 1;
        const Quad w = Quad(2 * i + 1), kappa = w * T.inverse, p = kappa * T.y[i];
        const Quad yn = p - T.y[i - 1];
        T.y[n] = yn;
        if (!finite_quad(yn)) return false;
        const PositiveBound an = pfromq_up(quad_abs(yn), ok), ap1 = T.ay[i - 1];
        T.ay[n] = an;
        const PositiveBound sum = padd(an, ap1);
        PositiveBound R;
            T.defect_acc.reset();
            if (!T.defect_acc.add(T.z, yn) || !T.defect_acc.add(T.z, T.y[i - 1]) || !T.defect_acc.add(-w, T.y[i])) return false;
            const Quad D = T.defect_acc.rounded();
            if (!finite_quad(D)) return false;
            R = pdiv(padd(pfromq_up(next_upper(quad_abs(D)), ok), pmul(T.dz, sum)), T.zL);
        T.asum = padd(T.asum, pmul(T.J[i], R));
        T.bsum = padd(T.bsum, pmul(T.Y[i], R));
        const PositiveBound prod = pmul(T.green, T.asum);
        const double den = (1.0 - pdbl_up(prod)) * positive_lower_factor;
        if (!(den > 0)) return false;
        const PositiveBound numer = padd(pmul(an, T.asum), pmul(T.J[n], T.bsum));
        T.E[n] = pdiv(pmul(T.green, numer), pfromd(den));
        T.Y[n] = padd(an, T.E[n]);
        T.last = n;
    }
    return ok;
}

static inline void majorant_extend(MomentRecurrence& T, int target) {
    while (T.avail < target) {
        const int j = T.avail;
        T.Bm[j + 1] = padd(pmul(pdiv(pint(2 * j + 1), T.zL), T.Bm[j]), T.Bm[j - 1]);
        T.avail = j + 1;
    }
}

// ---- Enclosed moment table -------------------------------------------------
struct MomentResult {
    int status = 0;
    std::vector<Quad> values, bounds;
};

// Exceeding this work budget rejects the candidate without truncating the original solver.
static const int moment_term_budget = 1024;

static void enclose_moments(MomentRecurrence& T, int m, int center, int kmax, MomentResult& out, int first_k=0) {
    bool ok = true;
    const PositiveBound pu{1.0, -113}, eta{1.0, -16495};
    if (std::fegetround() != FE_TONEAREST) { out.status = 9; return; }
    out.status = 0;
    out.values.assign(3 * (kmax + 1), Quad(0)); out.bounds = out.values;
    static thread_local Accumulator acc[3];
    std::vector<Quad> coef(moment_term_budget);
    const Quad origin = T.origin, root = T.root, h = T.h;
    const Quad step = center == 0 ? -h : h;
    const Interval stepb = center == 0 ? Interval{-T.hb.upper, -T.hb.lower} : T.hb;
    const Quad mq = Quad(m);
    // Leading factor and its rigorous starting relative bound.
    Quad leading = Quad(2) / Quad(2 * m + 1);
    Interval leadb = ivdivp(Interval{2, 2}, Interval{Quad(2 * m + 1), Quad(2 * m + 1)});
    if (center == 1) {
        const Quad ratio = T.q / origin;
        const Interval ratb = ivdivp(Interval{T.q, T.q}, T.originb);
        Quad p = 1; Interval pb{1, 1}, fb = ratb; Quad f = ratio;
        int ex = m / 2;
        while (ex > 0) { if (ex & 1) { p *= f; pb = ivmul(pb, fb); } ex >>= 1; if (ex) { f *= f; fb = ivmul(fb, fb); } }
        leading *= p; leadb = ivmul(leadb, pb);
        if (m & 1) {
            const Quad rt = sqrtq(ratio);
            const Interval r2 = ivmul(Interval{rt, rt}, Interval{rt, rt}), rs = ivsub(ratb, r2);
            const Quad e = next_upper(fmaxq(fabsq(rs.lower), fabsq(rs.upper)) / rt);
            leading *= rt; leadb = ivmul(leadb, Interval{next_lower(rt - e), next_upper(rt + e)});
        }
    }
    if (!(leadb.lower > 0) || !finite_quad(leadb.upper)) { out.status = 1; return; }
    PositiveBound lead_rel0 = pdiv(pfromq_up(next_upper(fmaxq(next_upper(leading - leadb.lower), next_upper(leadb.upper - leading))), ok), pfromq_dn(leadb.lower, ok));
    bool lead_fast = true;
    const PositiveBound invoL = pdiv(PositiveBound{1.0, 0}, pfromq_dn(T.originb.lower, ok));
    const PositiveBound c2_upper = pmul(pfromq_up(T.c, ok), pfromq_up(T.c, ok));
    for (int k = 0; k <= kmax; ++k) {
        const PositiveBound Sk = padd(lead_rel0, pmul(pint(k), T.pref_step));
        if(k>=first_k) {
        Quad a = 1;
        Interval ab{1, 1};
        bool fast = true;
        PositiveBound Cint{0, 0};
        PositiveBound Ae[3] = {{0, 0}, {0, 0}, {0, 0}}, Ce[3] = {{0, 0}, {0, 0}, {0, 0}}, tail[3] = {{0, 0}, {0, 0}, {0, 0}};
        int r = 0;
        const double fam = double(std::max(0, 2 * (m + k + 2) - 1));
        bool stopped = false;
        for (; r < moment_term_budget; ++r) {
            const int n = m + k + r;
            if (n + 2 > T.nmax) { out.status = 2; return; }
            if (n + 2 > T.last) { if (!extend_recurrence(T, n + 2)) { out.status = 2; return; } }
            if (n + 2 > T.avail) majorant_extend(T, n + 2);
            coef[r] = a;
            const PositiveBound aU = pfromq_up(quad_abs(a), ok);
            PositiveBound aC, aTail;
            double onep;
            if (fast) { aC = pmul(aU, T.rel.q[r]); onep = T.rel.onep[r]; aTail = pmul(aU, pfromd(onep)); }
            else { aC = Cint; aTail = pfromq_up(fmaxq(quad_abs(ab.lower), quad_abs(ab.upper)), ok); }
            for (int f = 0; f < 3; ++f) {
                Ae[f] = padd(Ae[f], pmul(aU, T.E[n + f]));
                Ce[f] = padd(Ce[f], pmul(aC, T.Y[n + f]));
            }
            if (n >= 1) {
                const PositiveBound rho_bound=padd(invoL,padd(
                    pdiv(pmul(invoL,pint(static_cast<unsigned long long>(fam))),pint(2*(r+1))),
                    pdiv(c2_upper,pint(2ULL*(2*n-1)*(r+1)))));
                const double rho=pdbl_up(rho_bound);
                if (rho < 1.0) {
                    const double den=(1.0-rho)*positive_lower_factor;
                    const PositiveBound frac=pdiv(rho_bound,pfromd(den));
                    bool done = true;
                    for (int f = 0; f < 3; ++f) {
                        tail[f] = pmul(pmul(aTail, T.Bm[n + f]), frac);
                        if (!ple(tail[f], padd(Ae[f], Ce[f]))) done = false;
                    }
                    if (done) { stopped = true; break; }
                }
            }
            const Quad num = center == 0 ? Quad(k + r) + Quad(0.5) : Quad(k + m + r + 1);
            const Quad den = Quad(r + 1) * (Quad(2 * k + m + r) + Quad(1.5));
            const Quad t1 = a * step, t2 = t1 * num, anext = t2 / den;
            if (fast && !(normal_quad(t1) && normal_quad(t2) && normal_quad(anext))) {
                // Switch to the absolute interval from the proved bound at the last normal step.
                const Quad cq = ptoq_up(pmul(aU, T.rel.q[r]), ok);
                ab = Interval{next_lower(a - cq), next_upper(a + cq)};
                fast = false;
            }
            if (!fast) {
                ab = ivdivp(ivmul(ivmul(ab, stepb), Interval{num, num}), Interval{den, den});
                if (!finite_quad(ab.lower) || !finite_quad(ab.upper)) { out.status = 6; return; }
                Cint = pfromq_up(next_upper(fmaxq(fabsq(anext - ab.lower), fabsq(ab.upper - anext))), ok);
            }
            a = anext;
        }
        if (!stopped) { out.status = 3; return; }
        const int count = r + 1;
        // Rounded exact signed sums use the same stored products as the native dot.
        for (int f = 0; f < 3; ++f) acc[f].reset();
        for (int i = 0; i < count; ++i) {
            const int n = m + k + i;
            for (int f = 0; f < 3; ++f) if (!acc[f].add(coef[i], T.y[n + f])) { out.status = 4; return; }
        }
        Quad sums[3];
        for (int f = 0; f < 3; ++f) { sums[f] = acc[f].rounded(); if (!finite_quad(sums[f])) { out.status = 4; return; } }
        PositiveBound arith[3];
        for (int f = 0; f < 3; ++f) arith[f] = padd(Ae[f], Ce[f]);
        // Family prefactors and final rounded moments.
        for (int f = 0; f < 3; ++f) {
            const PositiveBound Es = padd(padd(arith[f], tail[f]), padd(pmul(pfromq_up(quad_abs(sums[f]), ok), pu), eta));
            const Quad P = f == 0 ? leading : (f == 1 ? leading / root : leading / origin);
            const Quad value = P * sums[f];
            out.values[3 * k + f] = value;
            PositiveBound relP = Sk;
            if (f == 1) relP = padd(relP, padd(T.eps_root, pu));
            if (f == 2) relP = padd(relP, padd(T.eps_origin, pu));
            const PositiveBound qP = relq(relP, ok), aP = pfromq_up(quad_abs(P), ok);
            if (!lead_fast || !normal_quad(P) || !(normal_quad(value) || sums[f] == 0)) {
                // Absolute product of the prefactor interval and the sum interval.
                Interval Lb = leadb;
                if (lead_fast) { const Quad rq = ptoq_up(pmul(pfromq_up(quad_abs(leading), ok), relq(Sk, ok)), ok); Lb = Interval{next_lower(leading - rq), next_upper(leading + rq)}; }
                const Interval Pb = f == 0 ? Lb : (f == 1 ? ivdivp(Lb, T.rootb) : ivdivp(Lb, T.originb));
                const Quad Eq = ptoq_up(Es, ok);
                const Interval mb = ivmul(Pb, Interval{next_lower(sums[f] - Eq), next_upper(sums[f] + Eq)});
                if (!finite_quad(mb.lower) || !finite_quad(mb.upper)) { out.status = 5; return; }
                out.bounds[3 * k + f] = next_upper(fmaxq(fabsq(value - mb.lower), fabsq(mb.upper - value)));
                if (!finite_quad(out.bounds[3 * k + f])) { out.status = 5; return; }
                continue;
            }
            const PositiveBound b = padd(padd(pmul(pfromq_up(quad_abs(value), ok), pu), pmul(pmul(aP, pfromq_up(quad_abs(sums[f]), ok)), qP)),
                pmul(pmul(aP, padd(PositiveBound{1.0, 0}, qP)), Es));
            out.bounds[3 * k + f] = ptoq_up(b, ok);
        }
        }
        const Quad num = mq + Quad(k) + Quad(0.5), d1 = Quad(2 * k + m) + Quad(1.5), d2 = Quad(2 * k + m) + Quad(2.5);
        const Quad den = d1 * d2, t1 = leading * (-h), t2 = t1 * num, next = t2 / den;
        if (lead_fast && !(normal_quad(t1) && normal_quad(t2) && normal_quad(next))) {
            // Hand the proved relative radius to the interval recurrence at the last normal step.
            const Quad rq = ptoq_up(pmul(pfromq_up(quad_abs(leading), ok), relq(Sk, ok)), ok);
            leadb = Interval{next_lower(leading - rq), next_upper(leading + rq)};
            lead_fast = false;
        }
        if (!lead_fast) {
            leadb = ivdivp(ivmul(ivmul(leadb, Interval{-T.hb.upper, -T.hb.lower}), Interval{num, num}), Interval{den, den});
            if (k < kmax && (!finite_quad(leadb.lower) || !finite_quad(leadb.upper))) { out.status = 5; return; }
        }
        leading = next;
        if (!ok) { out.status = 8; return; }
    }
}


static int enclose_stream_rounding(Quad* out) {
    bool ok=true;
    const PositiveBound pu{1,-112},eta{1,-16494};
    for(int f=0;f<2;++f) {
        const PositiveBound radius=padd(pfromq_up(next_upper(out[2+f]),ok),
            padd(pmul(pfromq_up(quad_abs(out[f]),ok),pu),eta));
        const PositiveBound last=padd(pfromq_up(next_upper(out[6+f]),ok),
            padd(pmul(pfromq_up(quad_abs(out[4+f]),ok),pu),eta));
        out[2+f]=ptoq_up(radius,ok);out[6+f]=ptoq_up(last,ok);
    }
    return ok?1:6;
}


// Fixed-input propagation only. Successful returns preserve the supplied values.
// Inputs and radii are exact stored binary128 numbers. Every result encloses the
// corresponding real expression, not errors in creating its stored factors.
static bool input_interval(Quad value,Quad radius,Interval& interval) noexcept {
    if(!finite_quad(value) || !finite_quad(radius) || radius<0) return false;
    interval=radius==0 ? Interval{value,value} : Interval{next_lower(value-radius),next_upper(value+radius)};
    return finite_quad(interval.lower) && finite_quad(interval.upper) && interval.lower<=interval.upper;
}

static bool enclosing_radius(Quad actual,Interval interval,Quad& radius) noexcept {
    if(!finite_quad(actual) || !finite_quad(interval.lower) || !finite_quad(interval.upper) || interval.lower>interval.upper) return false;
    if(actual==interval.lower && actual==interval.upper) {radius=0;return true;}
    const Quad left=next_upper(quad_abs(actual-interval.lower)),right=next_upper(quad_abs(interval.upper-actual));
    radius=left>right ? left : right;
    return finite_quad(radius) && radius>=0;
}


} // namespace psms_checked_integrals
#endif
extern "C" int pint_checked_moments_quad(int m,const acousticts_quad_t* c,const acousticts_quad_t* q,const acousticts_quad_t* s,int capacity,
        int required,int first,acousticts_quad_t* p1,acousticts_quad_t* p2,acousticts_quad_t* p4,acousticts_quad_t* errors) {
#if ACOUSTICTS_HAVE_EXACT_QUAD_DOT
    using namespace psms_checked_integrals;
    if(!pint_checked_moments_available_quad()) return 0;
    if(!c || !q || !s || !p1 || !p2 || !p4 || !errors || capacity>INT_MAX/6) return 0;
    if(m<0 || required<2 || required>capacity || first<1 || first>required ||
       !(first&1) || m>1022 || std::fegetround()!=FE_TONEAREST) return 0;
    const int kmax=(required-2)/2,first_k=(first-1)/2;
    if(kmax>1022-m || first_k>kmax) return 0;
    try {
        MomentRecurrence table;
        if(!initialize_recurrence(table,*c,*q,*s,0,m+kmax+moment_term_budget+2,moment_term_budget)) return 0;
        MomentResult result;
        enclose_moments(table,m,0,kmax,result,first_k);
        if(result.status) return 0;
        for(int k=first_k;k<=kmax;++k) for(int f=0;f<3;++f)
            if(!finite_quad(result.values[3*k+f]) || !finite_quad(result.bounds[3*k+f]) ||
               result.bounds[3*k+f]<0) return 0;
        for(int k=first_k;k<=kmax;++k) {
            const int i=2*k;
            p1[i]=result.values[3*k];p2[i]=result.values[3*k+1];p4[i]=result.values[3*k+2];
            errors[2*capacity+i]=result.bounds[3*k];
            errors[3*capacity+i]=result.bounds[3*k+1];
            errors[5*capacity+i]=result.bounds[3*k+2];
        }
        return 1;
    } catch(...) { return 0; }

#else
    return 0;
#endif
}
extern "C" int pint_checked_parity_quad(const acousticts_quad_t* alpha,const acousticts_quad_t* beta,const acousticts_quad_t* minus,const acousticts_quad_t* plus,
        const acousticts_quad_t* eminus,const acousticts_quad_t* eplus,const acousticts_quad_t* returned,acousticts_quad_t* error) {
#if ACOUSTICTS_HAVE_EXACT_QUAD_DOT
    using namespace psms_checked_integrals;
    if(!pint_checked_moments_available_quad()) return 0;
    if(!alpha || !beta || !minus || !plus || !eminus || !eplus || !returned || !error) return 0;
    if(std::fegetround()!=FE_TONEAREST || !finite_quad(*alpha) || *alpha==0 || !finite_quad(*beta) ||
       !finite_quad(*returned)) return 0;
    Interval low,high;
    if(!input_interval(*minus,*eminus,low) || !input_interval(*plus,*eplus,high)) return 0;
    Interval target=ivsub(high,ivmul(Interval{*beta,*beta},low));
    if(!finite_quad(target.lower) || !finite_quad(target.upper)) return 0;
    if(*alpha<0) target=Interval{-target.upper,-target.lower};
    const Quad denominator=quad_abs(*alpha);
    target=ivdivp(target,Interval{denominator,denominator});
    return enclosing_radius(*returned,target,*error) ? 1 : 0;

#else
    return 0;
#endif
}
extern "C" int pint_checked_radial_quad(const acousticts_quad_t* F,const acousticts_quad_t* D,const acousticts_quad_t* C,const acousticts_quad_t* s1,const acousticts_quad_t* s2,
        const acousticts_quad_t* u1,const acousticts_quad_t* u2,const acousticts_quad_t* full,const acousticts_quad_t* first_deriv,const acousticts_quad_t* correction,
        const acousticts_quad_t* derivsum,acousticts_quad_t* error,acousticts_quad_t* error_d) {
#if ACOUSTICTS_HAVE_EXACT_QUAD_DOT
    using namespace psms_checked_integrals;
    if(!pint_checked_moments_available_quad()) return 0;
    if(!F || !D || !C || !s1 || !s2 || !u1 || !u2 || !full || !first_deriv || !correction ||
       !derivsum || !error || !error_d) return 0;
    if(std::fegetround()!=FE_TONEAREST) return 0;
    for(Quad value:{*F,*D,*C,*s1,*s2,*u1,*u2,*full,*first_deriv,*correction,*derivsum})
        if(!finite_quad(value)) return 0;
    if(*u1<0 || *u2<0) return 0;
    // s1/s2 and u1/u2 are rounded exact signed and positive dot sums.
    const Quad upper1=next_upper(*u1),upper2=next_upper(*u2);
    if(!finite_quad(upper1) || !finite_quad(upper2)) return 0;
    const Interval x1{next_lower(next_lower(*s1)-upper1),next_upper(next_upper(*s1)+upper1)};
    const Interval x2{next_lower(next_lower(*s2)-upper2),next_upper(next_upper(*s2)+upper2)};
    if(!finite_quad(x1.lower) || !finite_quad(x1.upper) || !finite_quad(x2.lower) || !finite_quad(x2.upper)) return 0;
    const Interval radial=ivmul(Interval{*F,*F},x1);
    const Interval derivative=ivadd(ivmul(Interval{*D,*D},x2),ivmul(Interval{*C,*C},x1));
    Quad radial_error=0,derivative_error=0;
    if(!enclosing_radius(*full,radial,radial_error) || !enclosing_radius(*derivsum,derivative,derivative_error)) return 0;
    *error=radial_error;*error_d=derivative_error;return 1;

#else
    return 0;
#endif
}
extern "C" int pint_checked_pack_quad(const acousticts_quad_t* value,const acousticts_quad_t* uncertainty,int decimal_shift,
        const acousticts_quad_t* returned_packed,acousticts_quad_t* error) {
#if ACOUSTICTS_HAVE_EXACT_QUAD_DOT
    using namespace psms_checked_integrals;
    if(!pint_checked_moments_available_quad()) return 0;
    if(!value || !uncertainty || !returned_packed || !error) return 0;
    if(std::fegetround()!=FE_TONEAREST || decimal_shift < -5000 || decimal_shift > 5000 ||
       !finite_quad(*returned_packed)) return 0;
    Interval target;
    if(!input_interval(*value,*uncertainty,target)) return 0;
    if(target.lower==0 && target.upper==0) return enclosing_radius(*returned_packed,target,*error) ? 1 : 0;
    // 5**48 has 112 bits, so every integer power 10**j through j=48 is exact.
    static const std::array<Quad,49> powers=[] {
        std::array<Quad,49> p{};p[0]=1;
        for(int i=1;i<=48;++i) p[i]=p[i-1]*Quad(10);
        return p;
    }();
    int remaining=decimal_shift<0 ? -decimal_shift : decimal_shift;
    while(remaining>0) {
        const int step=remaining<48 ? remaining : 48;
        const Quad factor=powers[step];
        if(decimal_shift>0) target=Interval{next_lower(target.lower*factor),next_upper(target.upper*factor)};
        else target=Interval{next_lower(target.lower/factor),next_upper(target.upper/factor)};
        if(!finite_quad(target.lower) || !finite_quad(target.upper)) return 0;
        remaining-=step;
    }
    return enclosing_radius(*returned_packed,target,*error) ? 1 : 0;

#else
    return 0;
#endif
}
#endif

#if ACOUSTICTS_HAVE_QUADMATH
// Preserve exact binary128 products in radial sums and panel discrepancies.

using PintQuad=acousticts_quad_t;
extern "C" int pint_radial_stream_quad(int m,int parity,int seed,int lim,int maxint,
    const PintQuad* enr,const PintQuad* p1,const PintQuad* p2,const PintQuad* p3,const PintQuad* p4,const PintQuad* errors,
    PintQuad dcon,PintQuad* out,int* jint) {
#if ACOUSTICTS_HAVE_EXACT_QUAD_DOT
  if(std::fegetround()!=FE_TONEAREST) return 6;
  if(m<0 || (parity!=0 && parity!=1) || seed<0 || lim<=seed ||
      2LL*lim+parity>=maxint) return 3;
  try {
    static thread_local std::array<quad_dot_product::Accumulator,4> accumulators;
    for(auto& accumulator:accumulators) accumulator.reset();
    std::fill(out,out+8,PintQuad(0));
    *jint=lim;
    PintQuad partial1=0,partial2=0,correction1=0,correction2=0;
    bool stopped=false;
    auto multiply=[](PintQuad a,PintQuad b,PintQuad& result) {
      result=a*b;
      if(!finiteq(a)||!finiteq(b)||!finiteq(result)) return 6;
      return a!=0 && b!=0 && result==0 ? 5 : 1;
    };
    auto add=[](PintQuad term,PintQuad& total,PintQuad& correction) {
      const PintQuad adjusted=term-correction;
      const PintQuad updated=total+adjusted;
      correction=(updated-total)-adjusted;
      total=updated;
      return finiteq(total)&&finiteq(correction);
    };
    for(int pass=0;pass<2;++pass) {
      PintQuad weight=1;
      const int first=pass==0?seed:seed-1;
      const int last=pass==0?lim:0;
      const int increment=pass==0?1:-1;
      for(int j=first;pass==0?j<=last:j>=last;j+=increment) {
        if((pass==0&&j>seed)||pass==1) {
          const int ratio_index=pass==0?j:j+1;
          const int i=2*ratio_index+parity+1;
          const PintQuad ratio=(PintQuad(2)*m+i-1)/(i-1);
          PintQuad factor=0,next=0;
          int status=multiply(enr[ratio_index-1],ratio,factor);
          if(status!=1) return status;
          if(pass==0) status=multiply(weight,factor,next);
          else {
            if(factor==0) return 6;
            next=weight/factor;
            if(!finiteq(next)) return 6;
            if(weight!=0&&next==0) return 5;
          }
          if(status!=1) return status;
          weight=next;
        }
        const int index=2*j+parity;
        const PintQuad value1=parity==0?p1[index]:p3[index];
        const PintQuad value2=parity==0?p2[index]:p4[index];
        const PintQuad error1=errors[index],error2=errors[maxint+index];
        if(error1<0||error2<0) return 6;
        PintQuad term1=0,term2=0,term_error1=0,term_error2=0;
        int status=multiply(weight,value1,term1);
        if(status!=1) return status;
        status=multiply(weight,value2,term2);
        if(status!=1) return status;
        status=multiply(fabsq(weight),error1,term_error1);
        if(status!=1) return status;
        status=multiply(fabsq(weight),error2,term_error2);
        if(status!=1) return status;
        if(!accumulators[0].add(weight,value1)||!accumulators[1].add(weight,value2)||
           !accumulators[2].add(fabsq(weight),error1)||!accumulators[3].add(fabsq(weight),error2)) return 6;
        if((pass==0&&!stopped)||j==lim) {
          if(pass==0&&!stopped) {
            if(!add(term1,partial1,correction1)||!add(term2,partial2,correction2)) return 6;
            if(j>seed&&partial1!=0&&partial2!=0) {
              const PintQuad endpoint1=fabsq(term1/partial1),endpoint2=fabsq(term2/partial2);
              if(finiteq(endpoint1)&&finiteq(endpoint2)&&endpoint1+endpoint2<dcon) {
                *jint=j;
                stopped=true;
              }
            }
          }
          if(j==lim) {
            out[4]=fabsq(term1);
            out[5]=fabsq(term2);
            out[6]=term_error1;
            out[7]=term_error2;
          }
        }
      }
    }
    for(int i=0;i<4;++i) {
      out[i]=accumulators[i].rounded();
      if(!finiteq(out[i])) return 6;
    }
    return pint_checked_moments_available_quad() ? psms_checked_integrals::enclose_stream_rounding(out) : 1;
  } catch(...) { return 6; }
#else
  return 6;
#endif
}
using PintRadialStream=decltype(&pint_radial_stream_quad);
extern "C" void profcn_set_radial_stream_quad(PintRadialStream);
namespace {
struct RegisterRadialStream {
  RegisterRadialStream() { profcn_set_radial_stream_quad(&pint_radial_stream_quad); }
};
const RegisterRadialStream registered_radial_stream;
}

extern "C" void profcn_set_checked_radial_callbacks_quad(
    decltype(&pint_checked_moments_available_quad), decltype(&pint_checked_moments_quad),
    decltype(&pint_checked_parity_quad), decltype(&pint_checked_radial_quad),
    decltype(&pint_checked_pack_quad));
namespace {
struct RegisterCheckedRadials {
  RegisterCheckedRadials() {
    // Keep checked evaluation optional in standalone native capability probes.
    profcn_set_checked_radial_callbacks_quad(&pint_checked_moments_available_quad,
        &pint_checked_moments_quad, &pint_checked_parity_quad,
        &pint_checked_radial_quad, &pint_checked_pack_quad);
  }
};
const RegisterCheckedRadials registered_checked_radials;
}


extern "C" int pint_panel_error_quad(const PintQuad* basis,const PintQuad* kernels,
    const PintQuad* coarse,int capacity,int nodes,int panel_nodes,int first,int last,int extra,
    PintQuad* errors1,PintQuad* errors2,PintQuad* errors4) {
#if ACOUSTICTS_HAVE_EXACT_QUAD_DOT
  if(std::fegetround()!=FE_TONEAREST||capacity<1||nodes<1||panel_nodes<1||
      nodes%panel_nodes!=0||first<1||last<first||last>capacity||extra<1||extra>capacity) return 6;
  try {
    static thread_local quad_dot_product::Accumulator panel,total;
    const PintQuad infinity=HUGE_VALQ;
    auto error=[&](int degree,int family,PintQuad& result) {
      total.reset();
      for(int start=0;start<nodes;start+=panel_nodes) {
        panel.reset();
        for(int node=start;node<start+panel_nodes;++node) {
          const PintQuad value=basis[size_t(node)*capacity+degree];
          if(!panel.add(value,kernels[size_t(node)*3+family])||
              !panel.add(value,-coarse[size_t(node)*3+family])) return false;
        }
        const PintQuad discrepancy=fabsq(panel.rounded());
        if(!finiteq(discrepancy)) return false;
        if(!total.add(nextafterq(discrepancy,infinity),PintQuad(1))) return false;
      }
      result=nextafterq(total.rounded(),infinity);
      return finiteq(result)&&result>=0;
    };
    for(int degree=first-1;degree<last;degree+=2) {
      if(!error(degree,0,errors1[degree])||!error(degree,1,errors2[degree])||
          !error(degree,2,errors4[degree])) return 6;
    }
    if(!error(extra-1,2,errors4[extra-1])) return 6;
    return 1;
  } catch(...) { return 6; }
#else
  return 6;
#endif
}
using PintPanelError=decltype(&pint_panel_error_quad);
extern "C" void profcn_set_panel_error_quad(PintPanelError);
namespace {
struct RegisterPanelError {
  RegisterPanelError() { profcn_set_panel_error_quad(&pint_panel_error_quad); }
};
const RegisterPanelError registered_panel_error;
}

#endif


#if ACOUSTICTS_HAVE_QUADMATH
namespace psms_degree_replay {
using Real = acousticts_quad_t;
struct PrimitiveSnapshot {
    std::vector<Real> bessel_ratios, bessel_derivative_ratios, bessel_values, bessel_log_derivatives, neumann_ratios,
        neumann_derivative_ratios, neumann_values, neumann_log_derivatives, legendre_ratios;
    std::vector<int> bessel_exponents, neumann_exponents;
    bool ready = false;
    size_t bytes() const {
        return sizeof(Real) * (bessel_ratios.capacity() + bessel_derivative_ratios.capacity() +
            bessel_values.capacity() + bessel_log_derivatives.capacity() +
            neumann_ratios.capacity() + neumann_derivative_ratios.capacity() + neumann_values.capacity() +
                neumann_log_derivatives.capacity() + legendre_ratios.capacity()) +
            sizeof(int) * (bessel_exponents.capacity() + neumann_exponents.capacity());
    }
};
struct DegreeCheckpoint {
    std::vector<Real> raw_ratios, normalized_ratios, coefficients;
    std::array<Real, 5> input_scalars{};
    std::array<int, 28> input_state{};
    std::array<Real, 18> real_outputs{};
    std::array<int, 31> integer_outputs{};
    int coefficient_length = 0, coefficient_prefix = 0;
    bool eligible = false;
    size_t bytes() const {
        return sizeof(*this) + sizeof(Real) * (raw_ratios.capacity() + normalized_ratios.capacity() +
            coefficients.capacity());
    }
};
struct BlockCheckpoint {
    int m = -1, count = 0;
    Real c = 0, x1 = 0;
    PrimitiveSnapshot primitive;
    std::map<int, DegreeCheckpoint> degrees;
    size_t bytes() const {
        size_t result = sizeof(*this) + primitive.bytes();
        for (const auto& degree : degrees) result += degree.second.bytes() + 64;
        return result;
    }
};
static BlockCheckpoint previous, current;
static bool active = false, nearest_history = true;
static constexpr size_t budget = 2 * 1024 * 1024;

static void clear() noexcept {
    previous = BlockCheckpoint{};
    current = BlockCheckpoint{};
}
static void invalidate() noexcept {
    active = false;
    clear();
}
static bool available() noexcept {
#if defined(__FAST_MATH__) || defined(__ASSOCIATIVE_MATH__) || defined(__RECIPROCAL_MATH__) || \
    (defined(__FINITE_MATH_ONLY__) && __FINITE_MATH_ONLY__ != 0) || \
    (defined(__GCC_IEC_559) && __GCC_IEC_559 <= 0)
    nearest_history = false;
    invalidate();
    return false;
#else
    if (std::fegetround() != FE_TONEAREST) {
        nearest_history = false;
        invalidate();
    }
    return active && nearest_history;
#endif
}
static bool admit(size_t additional) noexcept {
    if (additional > budget || previous.bytes() + current.bytes() > budget - additional) {
        invalidate();
        return false;
    }
    return true;
}
static bool equal(const Real* left, const Real* right, size_t count) noexcept {
    return count == 0 || std::memcmp(left, right, count * sizeof(Real)) == 0;
}
template<class T> static bool prefix(const std::vector<T>& left, const std::vector<T>& right, int count) noexcept {
    return count >= 0 && size_t(count) <= left.size() && size_t(count) <= right.size() &&
        (count == 0 || std::memcmp(left.data(), right.data(), size_t(count) * sizeof(T)) == 0);
}
Scope::Scope() noexcept {
    clear();
    active = true;
    available();
}
Scope::~Scope() noexcept { invalidate(); }
void Scope::disable() noexcept { invalidate(); }

extern "C" int acousticts_degree_replay_available_quad() noexcept { return available() ? 1 : 0; }
extern "C" void acousticts_degree_replay_clear_quad() noexcept { invalidate(); }
extern "C" void acousticts_degree_replay_begin_quad(int m, int count, const Real* c, const Real* x1) noexcept {
    if (!available()) return;
    if (current.m == m && equal(c, &current.c, 1) && equal(x1, &current.x1, 1)) previous = std::move(current);
    else previous = BlockCheckpoint{};
    current = BlockCheckpoint{};
    current.m = m;
    current.count = count;
    current.c = *c;
    current.x1 = *x1;
}
extern "C" int acousticts_degree_replay_probe_quad(int m, int count, int l, int ratio_count,
        int bessel_support, int neumann_support, int bessel_ratio_capacity, int radial_value_capacity,
            int neumann_ratio_capacity, int legendre_ratio_capacity, const int* input_state,
            const Real* input_scalars,
        Real* raw_ratios, const Real* bessel_ratios, const Real* bessel_derivative_ratios, const Real* bessel_values,
            const Real* bessel_log_derivatives, const int* bessel_exponents,
        const Real* neumann_ratios, const Real* neumann_derivative_ratios, const Real* neumann_values,
            const Real* neumann_log_derivatives, const int* neumann_exponents,
        const Real* legendre_ratios, int coefficient_capacity, Real* coefficients, int* coefficient_length,
            Real* real_outputs, int* integer_outputs) noexcept {
    if (!available()) return 0;
    const int parity = (l - m) % 2, degree_seed = (l - m) / 2, coefficient_prefix = 2 * ratio_count + parity + 1;
    if (m != current.m || count != current.count || ratio_count <= 0 || degree_seed < 0 ||
        degree_seed > ratio_count || coefficient_prefix > coefficient_capacity) return 0;
    const int bessel_last_term = std::max(degree_seed, bessel_support / 2 - parity),
        neumann_last_term = std::max(degree_seed, neumann_support / 2 - parity);
    const int bessel_ratio_prefix = m + 2 * bessel_last_term + parity, neumann_ratio_prefix = m +
        2 * neumann_last_term + parity, legendre_ratio_prefix = 2 * neumann_last_term + parity + 1;
    const int bessel_value_prefix = std::max(l + 1, 3), neumann_value_prefix = l + 1,
        bessel_value_count = std::min(radial_value_capacity, bessel_ratio_capacity),
        neumann_value_count = std::min(radial_value_capacity, neumann_ratio_capacity + 1);
    if (bessel_ratio_prefix > bessel_ratio_capacity || bessel_value_prefix > bessel_value_count ||
        neumann_ratio_prefix > neumann_ratio_capacity || neumann_value_prefix > neumann_value_count ||
        legendre_ratio_prefix > legendre_ratio_capacity) return 0;
    try {
        // These native tables are immutable within the single-order block.
        auto& primitive = current.primitive;
        if (!primitive.ready) {
            const size_t bytes = sizeof(Real) * (2ULL * bessel_ratio_capacity + 2ULL * bessel_value_count +
                2ULL * neumann_ratio_capacity + 2ULL * neumann_value_count + legendre_ratio_capacity) +
                sizeof(int) * (size_t(bessel_value_count) + neumann_value_count);
            if (!admit(bytes)) return 0;
            primitive.bessel_ratios.assign(bessel_ratios, bessel_ratios + bessel_ratio_capacity);
            primitive.bessel_derivative_ratios.assign(bessel_derivative_ratios, bessel_derivative_ratios +
                bessel_ratio_capacity);
            primitive.bessel_values.assign(bessel_values, bessel_values + bessel_value_count);
            primitive.bessel_log_derivatives.assign(bessel_log_derivatives, bessel_log_derivatives +
                bessel_value_count);
            primitive.bessel_exponents.assign(bessel_exponents, bessel_exponents + bessel_value_count);
            primitive.neumann_ratios.assign(neumann_ratios, neumann_ratios + neumann_ratio_capacity);
            primitive.neumann_derivative_ratios.assign(neumann_derivative_ratios, neumann_derivative_ratios +
                neumann_ratio_capacity);
            primitive.neumann_values.assign(neumann_values, neumann_values + neumann_value_count);
            primitive.neumann_log_derivatives.assign(neumann_log_derivatives, neumann_log_derivatives +
                neumann_value_count);
            primitive.neumann_exponents.assign(neumann_exponents, neumann_exponents + neumann_value_count);
            primitive.legendre_ratios.assign(legendre_ratios, legendre_ratios + legendre_ratio_capacity);
            primitive.ready = true;
        }
        if (!admit(sizeof(DegreeCheckpoint) + 64 + sizeof(Real) * (size_t(2) * ratio_count +
            coefficient_prefix))) return 0;
        auto& checkpoint = current.degrees[l];
        checkpoint.raw_ratios.assign(raw_ratios, raw_ratios + ratio_count);
        checkpoint.coefficient_prefix = coefficient_prefix;
        std::copy(input_state, input_state + 28, checkpoint.input_state.begin());
        std::copy(input_scalars, input_scalars + 5, checkpoint.input_scalars.begin());
        const auto found = previous.degrees.find(l);
        if (found == previous.degrees.end() || !found->second.eligible || !previous.primitive.ready) return 0;
        const auto& saved_checkpoint = found->second;
        const auto& old_primitive = previous.primitive;
        // Fresh support and eigenpairs must reproduce every downstream input bit.
        const bool matches = checkpoint.raw_ratios.size() == saved_checkpoint.raw_ratios.size() &&
            checkpoint.coefficient_prefix == saved_checkpoint.coefficient_prefix &&
            equal(checkpoint.raw_ratios.data(), saved_checkpoint.raw_ratios.data(), checkpoint.raw_ratios.size()) &&
                checkpoint.input_state == saved_checkpoint.input_state &&
            equal(checkpoint.input_scalars.data(), saved_checkpoint.input_scalars.data(), 5) &&
            prefix(primitive.bessel_ratios, old_primitive.bessel_ratios, bessel_ratio_prefix) &&
                prefix(primitive.bessel_derivative_ratios, old_primitive.bessel_derivative_ratios,
                bessel_ratio_prefix) &&
            prefix(primitive.bessel_values, old_primitive.bessel_values, bessel_value_prefix) &&
                prefix(primitive.bessel_log_derivatives, old_primitive.bessel_log_derivatives, bessel_value_prefix) &&
            prefix(primitive.bessel_exponents, old_primitive.bessel_exponents, bessel_value_prefix) &&
                prefix(primitive.neumann_ratios, old_primitive.neumann_ratios, neumann_ratio_prefix) &&
            prefix(primitive.neumann_derivative_ratios, old_primitive.neumann_derivative_ratios,
                neumann_ratio_prefix) && prefix(primitive.neumann_values, old_primitive.neumann_values,
                neumann_value_prefix) &&
            prefix(primitive.neumann_log_derivatives, old_primitive.neumann_log_derivatives, neumann_value_prefix) &&
                prefix(primitive.neumann_exponents, old_primitive.neumann_exponents, neumann_value_prefix) &&
            prefix(primitive.legendre_ratios, old_primitive.legendre_ratios, legendre_ratio_prefix);
        if (!matches) return 0;
        checkpoint = saved_checkpoint;
        std::copy(checkpoint.normalized_ratios.begin(), checkpoint.normalized_ratios.end(), raw_ratios);
        std::fill(coefficients, coefficients + coefficient_capacity, Real(0));
        std::copy(checkpoint.coefficients.begin(), checkpoint.coefficients.end(), coefficients);
        *coefficient_length = checkpoint.coefficient_length;
        std::copy(checkpoint.real_outputs.begin(), checkpoint.real_outputs.end(), real_outputs);
        std::copy(checkpoint.integer_outputs.begin(), checkpoint.integer_outputs.end(), integer_outputs);
        previous.degrees.erase(found);
        return 1;
    } catch (const std::bad_alloc&) {
        invalidate();
        return 0;
    }
}
extern "C" void acousticts_degree_replay_finish_quad(int l, int eligible, int normalization_terms, int ratio_count,
        const Real* normalized_ratios, int coefficient_capacity, const Real* coefficients, int coefficient_length,
            int status,
        const Real* real_outputs, const int* integer_outputs) noexcept {
    if (!available()) return;
    const auto found = current.degrees.find(l);
    if (found == current.degrees.end()) return;
    auto& checkpoint = found->second;
    if (!eligible || normalization_terms > ratio_count || status != 0 || coefficient_length <= 0 ||
        coefficient_length > checkpoint.coefficient_prefix || checkpoint.coefficient_prefix > coefficient_capacity) {
        previous.degrees.erase(l);
        return;
    }
    try {
        if (!admit(sizeof(Real) * (size_t(ratio_count) + checkpoint.coefficient_prefix))) return;
        checkpoint.normalized_ratios.assign(normalized_ratios, normalized_ratios + ratio_count);
        checkpoint.coefficients.assign(coefficients, coefficients + checkpoint.coefficient_prefix);
        checkpoint.coefficient_length = coefficient_length;
        std::copy(real_outputs, real_outputs + 18, checkpoint.real_outputs.begin());
        std::copy(integer_outputs, integer_outputs + 31, checkpoint.integer_outputs.begin());
        checkpoint.eligible = true;
        previous.degrees.erase(l);
    } catch (const std::bad_alloc&) {
        invalidate();
    }
}
extern "C" void profcn_set_degree_replay_callbacks_quad(
    decltype(&acousticts_degree_replay_available_quad), decltype(&acousticts_degree_replay_begin_quad),
    decltype(&acousticts_degree_replay_probe_quad), decltype(&acousticts_degree_replay_finish_quad),
    decltype(&acousticts_degree_replay_clear_quad));
struct RegisterCallbacks {
    RegisterCallbacks() {
        // Standalone Fortran retains the original path when these callbacks are absent.
        profcn_set_degree_replay_callbacks_quad(&acousticts_degree_replay_available_quad,
            &acousticts_degree_replay_begin_quad, &acousticts_degree_replay_probe_quad,
            &acousticts_degree_replay_finish_quad, &acousticts_degree_replay_clear_quad);
    }
};
static const RegisterCallbacks registered_callbacks;
}
#endif

#include "psms_smn.h"

// Report whether configure enabled the complete native C++/Fortran binary128
// backend.  R-side validation uses this before accepting precision = "quad".
// [[Rcpp::export]]
bool quad_precision_available_cpp() {
    return ACOUSTICTS_HAVE_QUADMATH != 0;
}

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
        std::vector<acousticts_quad_t> eta_quad(eta_len);
        for (int i = 0; i < eta_len; ++i)
            eta_quad[i] = static_cast<acousticts_quad_t>(arg[i]);
        auto res = Smn_matrix<acousticts_quad_t>(
            m_vec, n_vec, static_cast<acousticts_quad_t>(c), eta_quad, normalize
        );
        return Smn_cpp_to_rcpp_list(res, m_len, n_len, eta_len);
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
        acousticts_quad_t c_q = static_cast<acousticts_quad_t>(c);
        acousticts_quad_t x1_q = static_cast<acousticts_quad_t>(x1);
        auto res = Rmn_matrix<acousticts_quad_t>(m_vec, n_vec, c_q, x1_q, kind);

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
#include "psms_rigid_release.h"
#include "psms_gas_support.h"

// Dispatch modal support selection and full boundary solves.
template<typename T>
std::complex<T> psms_public_field(
    int m_seed, int n_seed, T chi_sw, T chi_body, T xi,
    T theta_body, T theta_scatter, T phi_body, T phi_scatter,
    T density_body, T density_sw, const std::vector<T>& nodes,
    const std::vector<T>& weights, const std::string& method,
    bool backscatter, bool adaptive, bool vectorized, bool retain_geometry_orders = true
) {
    // Full solutions share packed native arithmetic in fixed and adaptive evaluations.
    const bool fluid = method == "Amn_fluid_gas" ||
        method == "Amn_fluid_coefficients" || method == "Amn_fluid";
    const bool rigid_or_release = method == "Amn_pressure_release" ||
        method == "Amn_fixed_rigid";
    if (fluid || rigid_or_release) {
        const auto incidence = psms_direction_from_angle(theta_body);
        const auto reception = psms_direction_from_angle(theta_scatter);
#if ACOUSTICTS_HAVE_QUADMATH
        if constexpr (is_quad_precision_v<T>) {
            if (adaptive && backscatter && (rigid_or_release || method == "Amn_fluid_gas")) {
                if (m_seed < 0 || n_seed < m_seed)
                    throw std::invalid_argument("Invalid spheroidal modal support.");
                if (rigid_or_release) return psms_adaptive_rigid_release_field(chi_sw, xi, incidence,
                    phi_scatter - phi_body, method == "Amn_fixed_rigid", retain_geometry_orders);
                return psms_gas_support::fields(chi_sw, chi_body, xi,
                    std::vector<int>{n_seed}, density_body / density_sw,
                    std::vector<SpheroidalDirection<T>>{incidence},
                    phi_scatter - phi_body).front();
            }
            // Adapt full-liquid backscatter without changing fixed integration or other boundaries.
            if (adaptive && backscatter &&
                    (method == "Amn_fluid" || method == "Amn_fluid_coefficients"))
                return psms_adaptive_liquid_field(m_seed, n_seed, chi_sw, chi_body, xi,
                    incidence, reception, phi_scatter - phi_body, density_body / density_sw);
        }
#endif
        return psms_fixed_support_field<T>(m_seed, n_seed, chi_sw, chi_body, xi,
            incidence, reception, phi_scatter - phi_body, density_body / density_sw,
            backscatter, rigid_or_release, method == "Amn_fixed_rigid", adaptive,
            method == "Amn_fluid" ? &nodes : nullptr,
            method == "Amn_fluid" ? &weights : nullptr,
            method == "Amn_fluid_gas" && is_quad_precision_v<T>);
    }
    // Retain the selected diagonal approximation and its optional modal stopping.
    return psms_fbs<T>(m_seed, n_seed, chi_sw, chi_body, xi,
        theta_body, theta_scatter, phi_body, phi_scatter,
        density_body, density_sw, nodes, weights, method,
        backscatter, adaptive, vectorized);
}

// C++-R INTERFACE FUNCTION FOR PROLATE SPHEROID MODAL SERIES SOLUTION
// ----------------------------------------------------------------------------
// Based on the algorithm proposed by Furusawa (1988).
// ============================================================================
inline bool psms_public_backscatter_angles(double theta_body, double theta_scatter,
        double phi_body, double phi_scatter) {
    const double pi = std::acos(-1.0);
    // Match the public monostatic construction exactly, preserving neighboring receiver angles.
    return theta_scatter == pi - theta_body &&
        (phi_scatter == phi_body + pi || phi_scatter == phi_body - pi);
}

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

    // Detect the exact R angle construction before promotion to the requested arithmetic.
    const bool is_backscatter = psms_public_backscatter_angles(
        theta_body, theta_scatter, phi_body, phi_scatter);

    // Extract the surrounding-medium density used in the boundary terms.
    double density_sw = medium["density"];
    
    // Extract the quadrature rule used for overlap integrals in the fluid/gas
    // pathways.
    Rcpp::NumericVector nodes = integration_pts["nodes"];
    Rcpp::NumericVector weights = integration_pts["weights"];
    // Reject unrepresentable quadrature dimensions before native copies or allocations.
    if (nodes.size() >= INT_MAX || weights.size() >= INT_MAX)
        Rcpp::stop("Fluid overlap quadrature exceeds the native integer range.");

    int n_freq = acoustics.nrows();
    std::vector<Rcomplex> f_bs_out(n_freq);

    if (precision == "quad") {
        std::vector<acousticts_quad_t> nodes_q(nodes.size()), weights_q(weights.size());
        for (int i = 0; i < nodes.size(); ++i) nodes_q[i] = static_cast<acousticts_quad_t>(nodes[i]);
        for (int i = 0; i < weights.size(); ++i) weights_q[i] = static_cast<acousticts_quad_t>(weights[i]);
        // Refine known Gauss rules at the same requested order without changing arbitrary rules.
        if (Amn_method == "Amn_fluid" && !adaptive &&
                integration_pts.hasAttribute("acousticTS_gauss_legendre") &&
                Rcpp::as<bool>(integration_pts.attr("acousticTS_gauss_legendre"))) {
            if (nodes.size() < 1 || nodes.size() != weights.size() || nodes.size() >= INT_MAX)
                Rcpp::stop("Prescribed fluid Gauss quadrature has invalid dimensions.");
            psms_automatic::Quadrature<acousticts_quad_t> quadrature;
            quadrature.prepare(static_cast<int>(nodes.size()));
            nodes_q = std::move(quadrature.nodes);
            weights_q = std::move(quadrature.weights);
        }
        for (int f = 0; f < n_freq; ++f) {
            // Each frequency is independent, so the wrapper simply dispatches
            // one retained PSMS solve per frequency and stores the resulting
            // complex backscatter amplitude.
            auto fbs = psms_public_field<acousticts_quad_t>(
                m_max[f], n_max[f],
                static_cast<acousticts_quad_t>(chi_sw[f]),
                static_cast<acousticts_quad_t>(chi_body[f]),
                static_cast<acousticts_quad_t>(xi),
                static_cast<acousticts_quad_t>(theta_body),
                static_cast<acousticts_quad_t>(theta_scatter),
                static_cast<acousticts_quad_t>(phi_body),
                static_cast<acousticts_quad_t>(phi_scatter),
                static_cast<acousticts_quad_t>(density_body),
                static_cast<acousticts_quad_t>(density_sw),
                nodes_q, weights_q,
                Amn_method, is_backscatter, adaptive, vectorized, n_freq > 1
            );
            f_bs_out[f] = to_Rcomplex(fbs);
        }
    } else {
        std::vector<double> nodes_vec = Rcpp::as<std::vector<double>>(nodes);
        std::vector<double> weights_vec = Rcpp::as<std::vector<double>>(weights);
        for (int f = 0; f < n_freq; ++f) {
            auto fbs = psms_public_field<double>(
                m_max[f], n_max[f],
                chi_sw[f], chi_body[f], xi,
                theta_body, theta_scatter, phi_body, phi_scatter,
                density_body, density_sw,
                nodes_vec, weights_vec,
                Amn_method, is_backscatter, adaptive, vectorized
            );
            f_bs_out[f] = to_Rcomplex(fbs);
        }
    }
    return Rcpp::wrap(f_bs_out);
}

// [[Rcpp::export]]
Rcpp::ComplexMatrix prolate_spheroid_rigid_release_batch(
    Rcpp::DataFrame acoustics,
    Rcpp::List body,
    Rcpp::List medium,
    Rcpp::NumericVector theta_body,
    Rcpp::LogicalVector rigid
) {
#if ACOUSTICTS_HAVE_QUADMATH
    const int n_freq = acoustics.nrows();
    const int n_case = theta_body.size();
    if (n_freq < 1 || n_case < 1 || rigid.size() != n_case)
        Rcpp::stop("Fixed-rigid/pressure-release batch dimensions are invalid.");
    const auto chi_sw = Rcpp::as<std::vector<double>>(acoustics["chi_sw"]);
    const auto m_max = Rcpp::as<std::vector<int>>(acoustics["m_max"]);
    const auto n_max = Rcpp::as<std::vector<int>>(acoustics["n_max"]);
    if (chi_sw.size() != size_t(n_freq) || m_max.size() != size_t(n_freq) ||
            n_max.size() != size_t(n_freq))
        Rcpp::stop("Fixed-rigid/pressure-release batch acoustics have inconsistent lengths.");
    const double xi = Rcpp::as<double>(body["xi"]);
    const double phi_body = Rcpp::as<double>(body["phi_body"]);
    const double phi_scatter = Rcpp::as<double>(body["phi_scatter"]);
    const double pi = std::acos(-1.0);
    if (phi_scatter != phi_body + pi && phi_scatter != phi_body - pi)
        Rcpp::stop("Fixed-rigid/pressure-release batching requires monostatic azimuths.");
    if (!medium.containsElementNamed("density"))
        Rcpp::stop("Fixed-rigid/pressure-release batch medium is incomplete.");
    std::vector<SpheroidalDirection<acousticts_quad_t>> incidences;
    std::vector<bool> rigid_flags;
    incidences.reserve(n_case);
    rigid_flags.reserve(n_case);
    for (int case_index = 0; case_index < n_case; ++case_index) {
        if (rigid[case_index] == NA_LOGICAL)
            Rcpp::stop("Fixed-rigid/pressure-release batch labels cannot be missing.");
        incidences.push_back(psms_direction_from_angle(
            static_cast<acousticts_quad_t>(theta_body[case_index])));
        rigid_flags.push_back(rigid[case_index] == TRUE);
    }
    Rcpp::ComplexMatrix output(n_freq, n_case);
    // Each frequency shares native degree blocks only among its own cases.
    for (int f = 0; f < n_freq; ++f) {
        Rcpp::checkUserInterrupt();
        if (m_max[f] < 0 || n_max[f] < m_max[f])
            Rcpp::stop("Invalid spheroidal modal support.");
        const auto fields = psms_adaptive_rigid_release_fields(
            static_cast<acousticts_quad_t>(chi_sw[f]),
            static_cast<acousticts_quad_t>(xi),
            incidences, rigid_flags,
            static_cast<acousticts_quad_t>(phi_scatter) -
                static_cast<acousticts_quad_t>(phi_body), n_freq > 1);
        for (int case_index = 0; case_index < n_case; ++case_index)
            output(f, case_index) = to_Rcomplex(fields[case_index]);
    }
    return output;
#else
    Rcpp::stop("Native quad precision is unavailable for fixed-rigid/pressure-release batching.");
#endif
}

// [[Rcpp::export]]
Rcpp::ComplexMatrix prolate_spheroid_fluid_batch(
    Rcpp::DataFrame acoustics,
    Rcpp::List body,
    Rcpp::List medium,
    Rcpp::IntegerMatrix m_seeds,
    Rcpp::IntegerMatrix n_seeds,
    Rcpp::NumericVector theta_body,
    Rcpp::NumericVector theta_scatter,
    bool gas
) {
#if ACOUSTICTS_HAVE_QUADMATH
    const int n_freq = acoustics.nrows();
    const int n_case = theta_body.size();
    if (n_freq < 1 || n_case < 1 || theta_scatter.size() != n_case ||
            m_seeds.nrow() != n_freq || n_seeds.nrow() != n_freq ||
            m_seeds.ncol() != n_case || n_seeds.ncol() != n_case)
        Rcpp::stop("Liquid/gas batch dimensions are invalid.");
    const auto chi_sw = Rcpp::as<std::vector<double>>(acoustics["chi_sw"]);
    const auto chi_body = Rcpp::as<std::vector<double>>(acoustics["chi_body"]);
    if (chi_sw.size() != size_t(n_freq) || chi_body.size() != size_t(n_freq))
        Rcpp::stop("Liquid/gas batch acoustics have inconsistent lengths.");
    const double xi = Rcpp::as<double>(body["xi"]);
    const double phi_body = Rcpp::as<double>(body["phi_body"]);
    const double phi_scatter = Rcpp::as<double>(body["phi_scatter"]);
    const double density_body = Rcpp::as<double>(body["density"]);
    const double density_sw = Rcpp::as<double>(medium["density"]);
    const auto density_ratio = static_cast<acousticts_quad_t>(density_body) /
        static_cast<acousticts_quad_t>(density_sw);
    const auto phi_difference = static_cast<acousticts_quad_t>(phi_scatter) -
        static_cast<acousticts_quad_t>(phi_body);
    std::vector<SpheroidalDirection<acousticts_quad_t>> incidences, receptions;
    incidences.reserve(n_case);
    receptions.reserve(n_case);
    for (int case_index = 0; case_index < n_case; ++case_index) {
        if (!psms_public_backscatter_angles(theta_body[case_index],
                theta_scatter[case_index], phi_body, phi_scatter))
            Rcpp::stop("Liquid/gas batching requires monostatic directions.");
        incidences.push_back(psms_direction_from_angle(
            static_cast<acousticts_quad_t>(theta_body[case_index])));
        receptions.push_back(psms_direction_from_angle(
            static_cast<acousticts_quad_t>(theta_scatter[case_index])));
    }
    Rcpp::ComplexMatrix output(n_freq, n_case);
    // Each frequency retains its angle-specific modal seeds and stopping tests.
    for (int f = 0; f < n_freq; ++f) {
        Rcpp::checkUserInterrupt();
        std::vector<int> m_limits, n_limits;
        m_limits.reserve(n_case);
        n_limits.reserve(n_case);
        for (int case_index = 0; case_index < n_case; ++case_index) {
            const int m_seed = m_seeds(f, case_index);
            const int n_seed = n_seeds(f, case_index);
            if (m_seed == NA_INTEGER || n_seed == NA_INTEGER ||
                    m_seed < 0 || n_seed < m_seed)
                Rcpp::stop("Invalid liquid/gas modal support.");
            m_limits.push_back(m_seed);
            n_limits.push_back(n_seed);
        }
        std::vector<std::complex<acousticts_quad_t>> fields;
        if (gas) {
            fields = psms_gas_support::fields(
                static_cast<acousticts_quad_t>(chi_sw[f]),
                static_cast<acousticts_quad_t>(chi_body[f]),
                static_cast<acousticts_quad_t>(xi), n_limits,
                density_ratio, incidences, phi_difference);
        } else {
            fields = psms_adaptive_liquid_fields(
                m_limits, n_limits,
                static_cast<acousticts_quad_t>(chi_sw[f]),
                static_cast<acousticts_quad_t>(chi_body[f]),
                static_cast<acousticts_quad_t>(xi), incidences, receptions,
                phi_difference, density_ratio);
        }
        for (int case_index = 0; case_index < n_case; ++case_index)
            output(f, case_index) = to_Rcomplex(fields[case_index]);
    }
    return output;
#else
    Rcpp::stop("Native quad precision is unavailable for liquid/gas batching.");
#endif
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
        std::vector<acousticts_quad_t> nodes_q(nodes.size()), weights_q(weights.size());
        for (int i = 0; i < nodes.size(); ++i) nodes_q[i] = static_cast<acousticts_quad_t>(nodes[i]);
        for (int i = 0; i < weights.size(); ++i) weights_q[i] = static_cast<acousticts_quad_t>(weights[i]);

        for (int f = 0; f < n_freq; ++f) {
            auto t_blocks = psms_tmatrix_blocks<acousticts_quad_t>(
                m_max[f], n_max[f],
                static_cast<acousticts_quad_t>(chi_sw[f]),
                static_cast<acousticts_quad_t>(chi_body[f]),
                static_cast<acousticts_quad_t>(xi),
                static_cast<acousticts_quad_t>(density_body),
                static_cast<acousticts_quad_t>(density_sw),
                nodes_q, weights_q,
                Amn_method
            );
            auto smn_inc = compute_smn_direction_matrix<acousticts_quad_t>(
                m_max[f], n_max[f], static_cast<acousticts_quad_t>(chi_sw[f]),
                psms_direction_from_angle(static_cast<acousticts_quad_t>(theta_body))
            );
            auto smn_scat = compute_smn_direction_matrix<acousticts_quad_t>(
                m_max[f], n_max[f], static_cast<acousticts_quad_t>(chi_sw[f]),
                psms_direction_from_angle(static_cast<acousticts_quad_t>(theta_scatter))
            );
            auto azimuth = compute_azimuth<acousticts_quad_t>(
                m_max[f],
                static_cast<acousticts_quad_t>(phi_body),
                static_cast<acousticts_quad_t>(phi_scatter)
            );

            auto f_val = compute_fbs_from_tmatrix_blocks<acousticts_quad_t>(
                m_max[f], n_max[f], azimuth, smn_inc, smn_scat, t_blocks
            );
            f_scat[f] = to_Rcomplex(f_val);
            t_store[f] = prolate_tmatrix_blocks_to_rcpp<acousticts_quad_t>(
                m_max[f], n_max[f], t_blocks
            );
        }
    } else {
        std::vector<double> nodes_vec = Rcpp::as<std::vector<double>>(nodes);
        std::vector<double> weights_vec = Rcpp::as<std::vector<double>>(weights);

        for (int f = 0; f < n_freq; ++f) {
            auto t_blocks = psms_tmatrix_blocks<double>(
                m_max[f], n_max[f], chi_sw[f], chi_body[f], xi,
                density_body, density_sw, nodes_vec, weights_vec, Amn_method
            );
            auto smn_inc = compute_smn_direction_matrix<double>(
                m_max[f], n_max[f], chi_sw[f], psms_direction_from_angle(theta_body)
            );
            auto smn_scat = compute_smn_direction_matrix<double>(
                m_max[f], n_max[f], chi_sw[f], psms_direction_from_angle(theta_scatter)
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
    auto smn_inc = compute_smn_direction_matrix<T>(m_max, n_max, chi_sw, psms_direction_from_angle(theta_body));
    auto smn_scat = compute_smn_direction_matrix<T>(m_max, n_max, chi_sw, psms_direction_from_angle(theta_scatter));
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
    auto smn_inc = compute_smn_direction_matrix<T>(m_max, n_max, chi_sw, psms_direction_from_angle(theta_body));

    std::vector<SpheroidalDirection<T>> receive_directions(theta_scatter.size());
    for (size_t i = 0; i < theta_scatter.size(); ++i) {
        receive_directions[i] = psms_direction_from_angle(theta_scatter[i]);
    }
    auto smn_scat = compute_smn_direction_matrices<T>(m_max, n_max, chi_sw, receive_directions);

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

template<typename T>
Rcpp::ComplexVector prolate_scattering_points_from_retained_blocks_single(
    const Rcpp::List& blocks,
    T chi_sw,
    T k_sw,
    T theta_body,
    T phi_body,
    const std::vector<T>& theta_scatter,
    const std::vector<T>& phi_scatter
) {
    // Evaluate paired receive directions while reusing the incident angular
    // basis and retained T-matrix blocks for the whole point batch.
    if (theta_scatter.size() != phi_scatter.size()) {
        Rcpp::stop("'theta_scatter' and 'phi_scatter' must have equal length.");
    }

    int m_max = blocks.size() - 1;
    Rcpp::List last_block = blocks[m_max];
    Rcpp::IntegerVector n_seq_last = last_block["n_seq"];
    int n_max = n_seq_last[n_seq_last.size() - 1];

    auto t_blocks = prolate_tmatrix_blocks_from_rcpp<T>(blocks);
    auto smn_inc = compute_smn_direction_matrix<T>(m_max, n_max, chi_sw, psms_direction_from_angle(theta_body));

    // Reuse each order's spheroidal expansion for the complete receive-point batch.
    std::vector<SpheroidalDirection<T>> receive_directions(theta_scatter.size());
    for (size_t i = 0; i < theta_scatter.size(); ++i) {
        receive_directions[i] = psms_direction_from_angle(theta_scatter[i]);
    }
    auto smn_scat = compute_smn_direction_matrices<T>(m_max, n_max, chi_sw, receive_directions);
    Rcpp::ComplexVector f_points(theta_scatter.size());
    for (size_t i = 0; i < theta_scatter.size(); ++i) {
        auto azimuth = compute_azimuth<T>(m_max, phi_body, phi_scatter[i]);
        auto raw_sum = compute_fbs_from_tmatrix_blocks<T>(
            m_max, n_max, azimuth, smn_inc, smn_scat[i], t_blocks
        );
        f_points[i] = to_Rcomplex(std::complex<T>(0, -T(2) / k_sw) * raw_sum);
    }

    return f_points;
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
        for (int f = 0; f < n_freq; ++f) {
            f_scat[f] = to_Rcomplex(
                prolate_scattering_from_retained_blocks_single<acousticts_quad_t>(
                    t_matrix[f],
                    static_cast<acousticts_quad_t>(chi_sw[f]),
                    static_cast<acousticts_quad_t>(k_sw[f]),
                    static_cast<acousticts_quad_t>(theta_body),
                    static_cast<acousticts_quad_t>(phi_body),
                    static_cast<acousticts_quad_t>(theta_scatter),
                    static_cast<acousticts_quad_t>(phi_scatter)
                )
            );
        }
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

// Evaluate one retained prolate T-matrix block set over paired receive points
// that share one incident direction.
// [[Rcpp::export]]
Rcpp::ComplexVector prolate_spheroid_scattering_points_from_tmatrix_cpp(
    Rcpp::DataFrame acoustics,
    Rcpp::List t_matrix,
    double theta_body,
    double phi_body,
    Rcpp::NumericVector theta_scatter,
    Rcpp::NumericVector phi_scatter,
    std::string precision = "double"
) {
    double chi_sw = Rcpp::as<std::vector<double>>(acoustics["chi_sw"])[0];
    double k_sw = Rcpp::as<std::vector<double>>(acoustics["k_sw"])[0];

    if (precision == "quad") {
        std::vector<acousticts_quad_t> theta_q(theta_scatter.size()), phi_q(phi_scatter.size());
        for (int i = 0; i < theta_scatter.size(); ++i) theta_q[i] = static_cast<acousticts_quad_t>(theta_scatter[i]);
        for (int i = 0; i < phi_scatter.size(); ++i) phi_q[i] = static_cast<acousticts_quad_t>(phi_scatter[i]);
        return prolate_scattering_points_from_retained_blocks_single<acousticts_quad_t>(
            t_matrix,
            static_cast<acousticts_quad_t>(chi_sw),
            static_cast<acousticts_quad_t>(k_sw),
            static_cast<acousticts_quad_t>(theta_body),
            static_cast<acousticts_quad_t>(phi_body),
            theta_q,
            phi_q
        );
    }

    return prolate_scattering_points_from_retained_blocks_single<double>(
        t_matrix,
        chi_sw,
        k_sw,
        theta_body,
        phi_body,
        Rcpp::as<std::vector<double>>(theta_scatter),
        Rcpp::as<std::vector<double>>(phi_scatter)
    );
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
        std::vector<acousticts_quad_t> theta_q(theta_scatter.size()), phi_q(phi_scatter.size());
        for (int i = 0; i < theta_scatter.size(); ++i) theta_q[i] = static_cast<acousticts_quad_t>(theta_scatter[i]);
        for (int j = 0; j < phi_scatter.size(); ++j) phi_q[j] = static_cast<acousticts_quad_t>(phi_scatter[j]);
        return prolate_scattering_grid_from_retained_blocks_single<acousticts_quad_t>(
            t_matrix,
            static_cast<acousticts_quad_t>(chi_sw),
            static_cast<acousticts_quad_t>(k_sw),
            static_cast<acousticts_quad_t>(theta_body),
            static_cast<acousticts_quad_t>(phi_body),
            theta_q,
            phi_q
        );
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

#include "../tests/testthat/native-solvers.h"
