#pragma once
// ============================================================================
// INCLUDES
// ============================================================================
// [[Rcpp::depends(BH)]]
#include <algorithm>
#include <array>
#include <cfenv>
#include <climits>
#include <cmath>
#include <complex>
#include <cstdint>
#include <cstring>
#include <limits>
#include <stdexcept>
#include <type_traits>

// Configure defines this only after verifying the complete C++/Fortran
// binary128 backend.  Standalone compilation defaults to the portable path.
#ifndef ACOUSTICTS_HAVE_QUADMATH
#define ACOUSTICTS_HAVE_QUADMATH 0
#endif

#if ACOUSTICTS_HAVE_QUADMATH
#include <quadmath.h>
#endif

// Native quad is enabled only when configure verifies both sides of the
// C++/Fortran interface.  A long-double placeholder keeps the templates
// compilable on unsupported libc++ toolchains, but public callers are prevented
// from selecting it as genuine quad precision.
#if ACOUSTICTS_HAVE_QUADMATH
using acousticts_quad_t = __float128;
#else
using acousticts_quad_t = long double;
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

// numeric_limits<__float128> is unspecialized in some standard libraries.
// Converting a double NaN preserves the sentinel in every supported precision.
template<typename T>
inline T precnan() {
    return static_cast<T>(std::numeric_limits<double>::quiet_NaN());
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
    // Reuse quadmath powers for the repeated integer exponents in modal tables.
    constexpr int limit = 512;
    static thread_local acousticts_quad_t values[2 * limit + 1] = {};
    static thread_local bool ready[2 * limit + 1] = {};
    if (exponent >= -limit && exponent <= limit) {
        const int index = exponent + limit;
        if (!ready[index]) {
            values[index] = powq(acousticts_quad_t(10), acousticts_quad_t(exponent));
            ready[index] = true;
        }
        return values[index];
    }
    return powq(double_to_quad(10.0), double_to_quad(static_cast<double>(exponent)));
#else
    return std::pow(acousticts_quad_t(10), acousticts_quad_t(exponent));
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
    return std::cos(x);
#endif
}

template<typename T>
inline T precpi() {
    return std::acos(T(-1));
}

template<>
inline acousticts_quad_t precpi<acousticts_quad_t>() {
#if ACOUSTICTS_HAVE_QUADMATH
    return acosq(acousticts_quad_t(-1));
#else
    return std::acos(acousticts_quad_t(-1));
#endif
}

template<typename T>
inline T precsqrt(T x) { return std::sqrt(x); }

template<>
inline acousticts_quad_t precsqrt(acousticts_quad_t x) {
#if ACOUSTICTS_HAVE_QUADMATH
    return sqrtq(x);
#else
    return std::sqrt(x);
#endif
}

template<typename T>
inline T precabs(T x) { return std::abs(x); }

template<>
inline acousticts_quad_t precabs(acousticts_quad_t x) {
#if ACOUSTICTS_HAVE_QUADMATH
    return fabsq(x);
#else
    return std::abs(x);
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

// Test finiteness without narrowing values to the public double output range.
template<typename T>
inline bool is_finite_real(T x) {
    return std::isfinite(x);
}

template<>
inline bool is_finite_real<acousticts_quad_t>(acousticts_quad_t x) {
#if ACOUSTICTS_HAVE_QUADMATH
    return finiteq(x);
#else
    return std::isfinite(x);
#endif
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

#if ACOUSTICTS_HAVE_QUADMATH && defined(__FLT128_MANT_DIG__) && \
    __FLT128_MANT_DIG__ == 113 && __FLT128_MAX_EXP__ == 16384 && \
    defined(__BYTE_ORDER__) && __BYTE_ORDER__ == __ORDER_LITTLE_ENDIAN__
#define ACOUSTICTS_HAVE_EXACT_QUAD_DOT 1
#else
#define ACOUSTICTS_HAVE_EXACT_QUAD_DOT 0
#endif

#if ACOUSTICTS_HAVE_EXACT_QUAD_DOT
namespace quad_dot_product {

using Quad = acousticts_quad_t;
static_assert(sizeof(Quad) == 16, "The exact dot product requires binary128 storage.");

struct BinaryParts {
    std::uint64_t low;
    std::uint64_t high;
    unsigned offset;
    bool negative;
    bool finite;
};

inline BinaryParts decode(Quad value) {
    std::uint64_t words[2];
    std::memcpy(words, &value, sizeof(value));
    const unsigned exponent = unsigned((words[1] >> 48) & 0x7fff);
    return {words[0], (words[1] & UINT64_C(0x0000ffffffffffff)) |
        (exponent ? UINT64_C(0x0001000000000000) : 0),
        exponent ? exponent - 1 : 0, (words[1] >> 63) != 0, exponent != 0x7fff};
}

inline Quad encode(std::uint64_t low, std::uint64_t high) {
    const std::uint64_t words[2]{low, high};
    Quad value;
    std::memcpy(&value, words, sizeof(value));
    return value;
}

inline void multiply(const BinaryParts& first, const BinaryParts& second,
    std::array<std::uint32_t, 8>& product) {
#if defined(__SIZEOF_INT128__) && !defined(ACOUSTICTS_FORCE_QUAD_DOT_LIMB32)
    // Four integer products recover all 226 significand bits without rounding.
    using Wide = __uint128_t;
    const Wide low = Wide(first.low) * second.low;
    const Wide cross_first = Wide(first.low) * second.high;
    const Wide cross_second = Wide(first.high) * second.low;
    const Wide middle = (low >> 64) + std::uint64_t(cross_first) +
        Wide(std::uint64_t(cross_second));
    const Wide upper = Wide(first.high) * second.high + (cross_first >> 64) +
        (cross_second >> 64) + (middle >> 64);
    const std::uint64_t words[4]{std::uint64_t(low), std::uint64_t(middle),
        std::uint64_t(upper), std::uint64_t(upper >> 64)};
    for (int i = 0; i < 4; ++i) {
        product[2 * i] = std::uint32_t(words[i]);
        product[2 * i + 1] = std::uint32_t(words[i] >> 32);
    }
#else
    // The portable limb path has the same exact product and bounded uint64 sums.
    const std::uint32_t first_words[4]{std::uint32_t(first.low),
        std::uint32_t(first.low >> 32), std::uint32_t(first.high),
        std::uint32_t(first.high >> 32)};
    const std::uint32_t second_words[4]{std::uint32_t(second.low),
        std::uint32_t(second.low >> 32), std::uint32_t(second.high),
        std::uint32_t(second.high >> 32)};
    product.fill(0);
    for (int i = 0; i < 4; ++i) {
        std::uint64_t carry = 0;
        for (int j = 0; j < 4; ++j) {
            const std::uint64_t next = std::uint64_t(first_words[i]) * second_words[j] +
                product[i + j] + carry;
            product[i + j] = std::uint32_t(next);
            carry = next >> 32;
        }
        product[i + 4] = std::uint32_t(carry);
    }
#endif
}

class Accumulator {
    // Adapt Neal (2015), arXiv:1505.05571, to exact binary128 products.
    // These bins cover every binary128 product exponent plus INT_MAX-term growth.
    static constexpr int bin_count = 2064;
    // A product of two minimum subnormals is exactly 2^(-32988).
    static constexpr int product_exponent = -32988;
    static constexpr std::int64_t radix = INT64_C(1) << 32;
    std::array<std::int64_t, bin_count> bins{};
    int low = bin_count;
    int high = -1;
    unsigned pending = 0;
    unsigned carry_interval;

    void add_bin(int index, std::uint32_t value, bool negative) {
        if (value == 0) return;
        if (index >= bin_count) throw std::overflow_error("Exact accumulator capacity exceeded.");
        low = std::min(low, index);
        high = std::max(high, index);
        bins[index] += negative ? -std::int64_t(value) : std::int64_t(value);
    }

    void normalize() {
        if (high < low) return;
        // Propagate signed carries into unsigned lower limbs and a signed top limb.
        const int last = high;
        std::int64_t carry = 0;
        for (int index = low; index <= last; ++index) {
            const std::int64_t value = bins[index] + carry;
            carry = value / radix;
            std::int64_t remainder = value % radix;
            if (remainder < 0) {
                remainder += radix;
                --carry;
            }
            bins[index] = remainder;
        }
        if (carry != 0) {
            if (last + 1 >= bin_count) throw std::overflow_error("Exact accumulator carry exceeded.");
            bins[++high] = carry;
        }
        // Remove sign extension so repeated carries cannot grow the active range.
        while (high > low && bins[high] == -1) {
            bins[high - 1] -= radix;
            bins[high--] = 0;
        }
        while (high >= low && bins[high] == 0) --high;
        pending = 0;
    }

    std::uint32_t word(int index) const {
        return index < low || index > high ? 0 : std::uint32_t(bins[index]);
    }

    std::uint32_t word_at_bit(int bit) const {
        const int index = bit / 32;
        const int offset = bit % 32;
        const std::uint64_t combined = std::uint64_t(word(index)) |
            (std::uint64_t(word(index + 1)) << 32);
        return std::uint32_t(combined >> offset);
    }

    bool any_below(int bit) const {
        if (bit <= 0) return false;
        const int whole = bit / 32;
        for (int index = low; index < whole && index <= high; ++index) {
            if (bins[index] != 0) return true;
        }
        const int remainder = bit % 32;
        const std::uint32_t mask = remainder ? std::uint32_t((UINT64_C(1) << remainder) - 1) : 0;
        return (word(whole) & mask) != 0;
    }

public:
    explicit Accumulator(unsigned interval = 1u << 20) : carry_interval(interval) {
        if (interval == 0 || interval > (1u << 20)) {
            throw std::invalid_argument("Invalid exact-accumulator carry interval.");
        }
    }

    void reset() {
        // Reuse storage while clearing only bins touched by the preceding dot product.
        if (high >= low) std::fill(bins.begin() + low, bins.begin() + high + 1, 0);
        low = bin_count;
        high = -1;
        pending = 0;
    }

    bool add(Quad first_value, Quad second_value) {
        const auto first = decode(first_value);
        const auto second = decode(second_value);
        if (!first.finite || !second.finite) return false;
        if ((first.low | first.high) == 0 || (second.low | second.high) == 0) return true;
        std::array<std::uint32_t, 8> product;
        multiply(first, second, product);
        const unsigned shift = first.offset + second.offset;
        const int initial = int(shift / 32);
        const unsigned offset = shift % 32;
        const bool negative = first.negative != second.negative;
        std::uint32_t carry = 0;
        // Each shifted exact product contributes to at most nine overlapping bins.
        for (int limb = 0; limb < 8; ++limb) {
            const std::uint64_t shifted = (std::uint64_t(product[limb]) << offset) | carry;
            add_bin(initial + limb, std::uint32_t(shifted), negative);
            carry = std::uint32_t(shifted >> 32);
        }
        add_bin(initial + 8, carry, negative);
        // At most 2^20 unsigned32 additions keep each signed64 bin below 2^53.
        if (++pending >= carry_interval) normalize();
        return true;
    }

    Quad rounded() {
        // Finalization is destructive, so reset before accumulating another dot product.
        normalize();
        if (high < low) return encode(0, 0);
        const bool negative = bins[high] < 0;
        if (negative) {
            for (int index = low; index <= high; ++index) bins[index] = -bins[index];
            normalize();
        }
        std::uint32_t leading = word(high);
        int leading_bit = 0;
        while (leading >>= 1) ++leading_bit;
        const int highest_bit = 32 * high + leading_bit;
        // Binary128 keeps 113 significant bits, or a fixed 2^(-16494) quantum.
        int discarded_bits = std::max(highest_bit - 112, 16494);
        std::array<std::uint32_t, 4> significand;
        for (int i = 0; i < 4; ++i) significand[i] = word_at_bit(discarded_bits + 32 * i);
        // Round once, including subnormal underflow and exact halfway cases.
        const bool guard = (word_at_bit(discarded_bits - 1) & 1) != 0;
        const bool sticky = any_below(discarded_bits - 1);
        if (guard && (sticky || (significand[0] & 1))) {
            for (int i = 0; i < 4; ++i) {
                if (++significand[i] != 0) break;
            }
        }
        if (significand[3] & (UINT32_C(1) << 17)) {
            for (int i = 0; i < 3; ++i) {
                significand[i] = (significand[i] >> 1) | (significand[i + 1] << 31);
            }
            significand[3] >>= 1;
            ++discarded_bits;
        }
        const std::uint64_t sign = negative ? UINT64_C(1) << 63 : 0;
        const bool normal = (significand[3] & (UINT32_C(1) << 16)) != 0;
        const int exponent = normal ? discarded_bits + 112 + product_exponent + 16383 : 0;
        if (exponent >= 0x7fff) return encode(0, sign | UINT64_C(0x7fff000000000000));
        return encode(std::uint64_t(significand[0]) | (std::uint64_t(significand[1]) << 32),
            sign | (std::uint64_t(exponent) << 48) | std::uint64_t(significand[2]) |
            (std::uint64_t(significand[3] & 0xffff) << 32));
    }
};

inline Quad kahan_dot(const Quad* first, const Quad* second, size_t count,
    size_t first_stride = 1, size_t second_stride = 1) {
    Quad sum = 0, correction = 0;
    for (size_t i = 0; i < count; ++i) {
        const Quad product = first[i * first_stride] * second[i * second_stride];
        const Quad next = product - correction;
        const Quad updated = sum + next;
        correction = (updated - sum) - next;
        sum = updated;
    }
    return sum;
}

inline Quad dot(const Quad* first, const Quad* second, size_t count,
    Accumulator& accumulator,
    size_t first_stride = 1, size_t second_stride = 1) {
    if (count > size_t(INT_MAX)) {
        return kahan_dot(first, second, count, first_stride, second_stride);
    }
    if (std::fegetround() != FE_TONEAREST) {
        return kahan_dot(first, second, count, first_stride, second_stride);
    }
    accumulator.reset();
    for (size_t i = 0; i < count; ++i) {
        if (!accumulator.add(first[i * first_stride], second[i * second_stride])) {
            return kahan_dot(first, second, count, first_stride, second_stride);
        }
    }
    return accumulator.rounded();
}

} // namespace quad_dot_product
#endif
