#pragma once

namespace psms_packed {
template<typename T> struct Decimal {
    T mantissa = 0;
    long long exponent = 0;
};

struct ScalingDiagnostics {
    size_t scaled_underflows = 0, incident_underflows = 0;
    long long largest_column_exponent = 0;
};

template<typename T> T scale_decimal(T value, long long exponent) {
    if (value == T(0) || !is_finite_real(value)) return value;
    // Split decimal multipliers so representable subnormal results survive.
    int step = std::numeric_limits<T>::max_exponent10 / 2;
#if ACOUSTICTS_HAVE_QUADMATH
    if constexpr (std::is_same_v<T, acousticts_quad_t>) step = FLT128_MAX_10_EXP / 2;
#endif
    if (step < 1) throw std::runtime_error("Native decimal exponent range is unavailable.");
    while (exponent > step) {
        value *= pow10_typed<T>(step);
        exponent -= step;
        if (!is_finite_real(value)) return value;
    }
    while (exponent < -step) {
        value *= pow10_typed<T>(-step);
        exponent += step;
        if (value == T(0)) return value;
    }
    return value * pow10_typed<T>(static_cast<int>(exponent));
}

template<typename T> std::complex<T> scale_decimal(const std::complex<T>& value,
        long long exponent) {
    return {scale_decimal(value.real(), exponent), scale_decimal(value.imag(), exponent)};
}

template<typename T> Decimal<T> normalize(T mantissa, long long exponent = 0) {
    if (!is_finite_real(mantissa)) throw std::runtime_error("Packed radial mantissa is nonfinite.");
    if (mantissa == T(0)) return {T(0), 0};
    int shift;
#if ACOUSTICTS_HAVE_QUADMATH
    if constexpr (std::is_same_v<T, acousticts_quad_t>) {
        shift = static_cast<int>(floorq(log10q(precabs(mantissa))));
    } else
#endif
    {
        shift = static_cast<int>(std::floor(std::log10(precabs(mantissa))));
    }
    // Correct the selected decade by native comparisons at exact powers of ten.
    mantissa = scale_decimal(mantissa, -static_cast<long long>(shift));
    exponent += shift;
    while (precabs(mantissa) >= T(10)) { mantissa /= T(10); ++exponent; }
    while (precabs(mantissa) < T(1)) { mantissa *= T(10); --exponent; }
    return {mantissa, exponent};
}

template<typename T> Decimal<T> component(const std::vector<T>& values,
        const std::vector<int>& exponents, size_t index) {
    if (index >= values.size() || index >= exponents.size() || is_na_int(exponents[index]))
        throw std::runtime_error("Packed radial value or exponent is missing.");
    return normalize(values[index], exponents[index]);
}

template<typename T> T shifted(const Decimal<T>& value, long long exponent,
        ScalingDiagnostics& diagnostics) {
    const T result = scale_decimal(value.mantissa, value.exponent - exponent);
    if (value.mantissa != T(0) && result == T(0)) ++diagnostics.scaled_underflows;
    return result;
}

template<typename T> bool below_radial_floor(const Decimal<T>& value) {
    if (value.mantissa == T(0)) return true;
    const auto floor = normalize(T(1e-300));
    return value.exponent < floor.exponent || (value.exponent == floor.exponent &&
        precabs(value.mantissa) < floor.mantissa);
}

template<typename T> struct ExteriorPair {
    std::complex<T> value, derivative;
    T incident_value = 0, incident_derivative = 0;
    long long exponent = 0;
};

template<typename T> ExteriorPair<T> exterior_pair(const ProfcnBatchResult<T>& data,
        int local, int degree, ScalingDiagnostics& diagnostics) {
    const size_t index = degree + size_t(data.lnum) * local;
    const auto first = component(data.r1c, data.ir1e, index);
    const auto first_derivative = component(data.r1dc, data.ir1de, index);
    const auto second = component(data.r2c, data.ir2e, index);
    const auto second_derivative = component(data.r2dc, data.ir2de, index);
    if (first.mantissa == T(0) && first_derivative.mantissa == T(0))
        throw std::runtime_error("Packed first-kind radial value and derivative both vanished.");
    // Preserve the existing second-kind tiny-pair guard before changing scales.
    if (below_radial_floor(second) && below_radial_floor(second_derivative))
        throw std::runtime_error("Packed second-kind radial pair failed the native radial guard.");
    long long exponent = std::numeric_limits<long long>::lowest();
    for (const auto& value : {first, first_derivative, second, second_derivative})
        if (value.mantissa != T(0)) exponent = std::max(exponent, value.exponent);
    ExteriorPair<T> result;
    result.exponent = exponent;
    result.value = {shifted(first, exponent, diagnostics), shifted(second, exponent, diagnostics)};
    result.derivative = {shifted(first_derivative, exponent, diagnostics),
        shifted(second_derivative, exponent, diagnostics)};
    // Incident forcing retains its physical first-kind scale independently of R3.
    result.incident_value = scale_decimal(first.mantissa, first.exponent);
    result.incident_derivative = scale_decimal(first_derivative.mantissa, first_derivative.exponent);
    if (!is_finite_real(result.incident_value) || !is_finite_real(result.incident_derivative))
        throw std::runtime_error("Incident radial values exceed the requested native range.");
    if ((first.mantissa != T(0) && result.incident_value == T(0)) ||
            (first_derivative.mantissa != T(0) && result.incident_derivative == T(0)))
        ++diagnostics.incident_underflows;
    diagnostics.largest_column_exponent = std::max(diagnostics.largest_column_exponent, exponent);
    return result;
}

template<typename T> std::pair<T, T> interior_pair(const ProfcnBatchResult<T>& data,
        int local, int degree, const Decimal<T>& density, ScalingDiagnostics& diagnostics) {
    const size_t index = degree + size_t(data.lnum) * local;
    const auto value = component(data.r1c, data.ir1e, index);
    const auto derivative = component(data.r1dc, data.ir1de, index);
    // Include density before selecting the common interior boundary scale.
    const auto pressure = normalize(value.mantissa * density.mantissa,
        value.exponent + density.exponent);
    if (pressure.mantissa == T(0) && derivative.mantissa == T(0))
        throw std::runtime_error("Packed interior radial pair vanished.");
    const long long exponent = pressure.mantissa == T(0) ? derivative.exponent :
        derivative.mantissa == T(0) ? pressure.exponent : std::max(pressure.exponent, derivative.exponent);
    const T scaled_value = shifted(pressure, exponent, diagnostics);
    const T scaled_derivative = shifted(derivative, exponent, diagnostics);
    const T scale = std::max(precabs(scaled_value), precabs(scaled_derivative));
    if (!(scale > T(0))) throw std::runtime_error("Packed interior normalization vanished.");
    return {scaled_value / scale, scaled_derivative / scale};
}
}

// ============================================================================
// PSMS SPHEROIDAL RADIAL HELPERS
// ============================================================================
// RADIAL SPHEROIDAL WAVE FUNCTION OF THE FIRST KIND [Rmn - SCALAR]
// ============================================================================
template<typename T>
struct RmnResult {
    // Real-valued kind-1 and kind-2 radial functions plus validity flags used
    // when building outgoing/incoming combinations.
    T val1;
    T der1;
    T val2;
    T der2;
    bool valid1;
    bool valid2;
};

template<typename T>
struct RadialValue {
    // Packed real and imaginary radial components for one (m, n) entry.
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
    // profcn returns the first- and second-kind radial functions separately.
    // Reconstruct those ingredients here so the caller can assemble kinds 3
    // and 4 as outgoing/incoming combinations.
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

    T R1 = precnan<T>();
    if (result1.r1c.size() > static_cast<size_t>(idx) &&
        !is_na_real(result1.r1c[idx]))
    {
        R1 = result1.r1c[idx] * pw1;
    }

    T Rd1 = precnan<T>();
    if (result1.r1dc.size() > static_cast<size_t>(idx) &&
        !is_na_real(result1.r1dc[idx]))
    {
        Rd1 = result1.r1dc[idx] * pwde1;
    }

    T R2 = precnan<T>();
    if (result2.r2c.size() > static_cast<size_t>(idx) &&
        !is_na_real(result2.r2c[idx]))
    {
        R2 = result2.r2c[idx] * pw2;
    }

    T Rd2 = precnan<T>();
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
std::pair<std::complex<T>, std::complex<T>> combine_higher_radial_kind(
    const RmnResult<T>& result, int kind
) {
    T second = result.val2, second_derivative = result.der2;
    const T tiny = T(1e-300);
    // Apply the same second-kind pair guard used by vector and matrix extraction.
    if (!is_finite_real(second) || !is_finite_real(second_derivative) ||
            (precabs(second) < tiny && precabs(second_derivative) < tiny)) {
        second = second_derivative = precnan<T>();
    }
    if (kind == 4) {
        second = -second;
        second_derivative = -second_derivative;
    }
    return {{result.val1, second}, {result.der1, second_derivative}};
}

template<typename T>
std::pair<std::complex<T>, std::complex<T>> Rmn_scalar(
    int m, int n, T c, T x1, int kind = 1
) {
    int lnum = checked_profcn_integer(static_cast<long long>(n) - m + 1,
        "inclusive radial degree count");
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

        T R = precnan<T>();
        if (val->size() > static_cast<size_t>(idx) &&
            !is_na_real((*val)[idx]))
        {
            R = (*val)[idx] * pw;
        }

        T Rd = precnan<T>();
        if (der->size() > static_cast<size_t>(idx) &&
            !is_na_real((*der)[idx]))
        {
            Rd = (*der)[idx] * pwde;
        }

        // Rmn(2) can underflow to a meaningless all-zero pair in extreme
        // cases. Mark those as invalid instead of letting them propagate.
        if (kind == 2) {
            const T tiny = T(1e-300);
            bool R_bad = !is_finite_real(R) ||
                        !is_finite_real(Rd) ||
                        (precabs(R) < tiny && precabs(Rd) < tiny);
            if (R_bad) {
                R = precnan<T>();
                Rd = precnan<T>();
            }
        }
        return {std::complex<T>(R, 0), std::complex<T>(Rd, 0)};
    }

    // Higher-order conventions use outgoing/incoming combinations built from
    // first- and second-kind radial functions.
    RmnResult<T> result = Rmn_higher_order<T>(m, n, lnum_safe, c, x1_r);
    return combine_higher_radial_kind(result, kind);
}

// ============================================================================
// RADIAL PROLATE SPHEROIDAL WAVE FUNCTIONS [Rmn - VECTOR/MATRIX]
// ------------------------------------------------------------------------
// This outputs vectors of size n x m with scalar xi
// ============================================================================
template<typename T>
struct RmnMatrixResult {
    // value and derivative use the same modal layout as Smn_matrix: rows are
    // modal requests and columns are degree or outer-product entries.
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
        // Real kinds read directly from the profcn storage arrays after the
        // base-10 exponents have been applied.
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
        out.val_real = precnan<T>();
        if (rc->size() > static_cast<size_t>(idx) && !is_na_real((*rc)[idx])) {
            out.val_real = (*rc)[idx] * pw;
        }
        out.der_real = precnan<T>();
        if (rdc->size() > static_cast<size_t>(idx) && !is_na_real((*rdc)[idx])) {
            out.der_real = (*rdc)[idx] * pwde;
        }
        out.val_imag = 0;
        out.der_imag = 0;
        if (kind == 2) {
            const T tiny = T(1e-300);
            bool bad = !is_finite_real(out.val_real) ||
                    !is_finite_real(out.der_real) ||
                    (precabs(out.val_real) < tiny && precabs(out.der_real) < tiny);
            if (bad) {
                out.val_real = precnan<T>();
                out.der_real = precnan<T>();
            }
        }
        return out;
    }

    // Higher orders rebuild the outgoing/incoming combinations from the
    // stored kind-1 and kind-2 components.
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
    out.val_real = precnan<T>();
    if (result1.r1c.size() > static_cast<size_t>(idx) &&
        !is_na_real(result1.r1c[idx])) {
        out.val_real = result1.r1c[idx] * pw1;
    }
    out.der_real = precnan<T>();
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
    out.val_imag = precnan<T>();
    if (result2.r2c.size() > static_cast<size_t>(idx) &&
        !is_na_real(result2.r2c[idx])) {
        out.val_imag = result2.r2c[idx] * pw2;
    }
    out.der_imag = precnan<T>();
    if (result2.r2dc.size() > static_cast<size_t>(idx) &&
        !is_na_real(result2.r2dc[idx])) {
        out.der_imag = result2.r2dc[idx] * pwde2;
    }

    const T tiny = T(1e-300);
    bool bad = !is_finite_real(out.val_imag) ||
            !is_finite_real(out.der_imag) ||
            (precabs(out.val_imag) < tiny && precabs(out.der_imag) < tiny);
    if (bad) {
        out.val_imag = precnan<T>();
        out.der_imag = precnan<T>();
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
            : precnan<T>();
        out.der_real = (offset < rdc.size() && !is_na_real(rdc[offset]))
            ? rdc[offset]
            : precnan<T>();
        out.val_imag = T(0);
        out.der_imag = T(0);
        if (kind == 2) {
            const T tiny = T(1e-300);
            bool bad = !is_finite_real(out.val_real) ||
                !is_finite_real(out.der_real) ||
                (precabs(out.val_real) < tiny && precabs(out.der_real) < tiny);
            if (bad) {
                out.val_real = precnan<T>();
                out.der_real = precnan<T>();
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
        : precnan<T>();
    out.der_real = (offset < result.r1dc.size() && !is_na_real(result.r1dc[offset]))
        ? result.r1dc[offset]
        : precnan<T>();
    out.val_imag = (offset < result.r2c.size() && !is_na_real(result.r2c[offset]))
        ? result.r2c[offset]
        : precnan<T>();
    out.der_imag = (offset < result.r2dc.size() && !is_na_real(result.r2dc[offset]))
        ? result.r2dc[offset]
        : precnan<T>();

    const T tiny = T(1e-300);
    bool bad = !is_finite_real(out.val_imag) ||
        !is_finite_real(out.der_imag) ||
        (precabs(out.val_imag) < tiny && precabs(out.der_imag) < tiny);
    if (bad) {
        out.val_imag = precnan<T>();
        out.der_imag = precnan<T>();
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

    // Scalar case: evaluate one order and one degree.
    if (m_len == 1 && n_len == 1) {
        if (n[0] < m[0])
            throw std::invalid_argument("'n' must be >= 'm'.");
        auto res = Rmn_scalar<T>(m[0], n[0], c, x1, kind);
        out.value[0][0] = res.first;
        out.derivative[0][0] = res.second;
        return out;
    }

    // Shared-order case: one m with multiple degree requests.
    if (m_len == 1 && n_len > 1) {
        int m_val = m[0];
        int n_min = *std::min_element(n.begin(), n.end());
        int n_max = *std::max_element(n.begin(), n.end());
        if (n_min < m_val)
            throw std::invalid_argument("All 'n' values must be >= 'm'.");
        int lnum = checked_profcn_integer(static_cast<long long>(n_max) - m_val + 1,
            "inclusive radial degree count");
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
                out.value[0][i] = std::complex<T>(precnan<T>(), 0);
                out.derivative[0][i] = std::complex<T>(precnan<T>(), 0);
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

    // Pairwise case: each m[i] is matched with n[i]. Group by unique m so the
    // profcn setup can be reused within each order block.
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
            int lnum = checked_profcn_integer(static_cast<long long>(n_max_local) - m_val + 1,
                "inclusive radial degree count");
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

    // Outer-product case: evaluate every requested order against every degree.
    for (int i = 0; i < m_len; ++i) {
        int m_val = m[i];
        int n_max = *std::max_element(n.begin(), n.end());
        int lnum = checked_profcn_integer(static_cast<long long>(n_max) - m_val + 1,
            "inclusive radial degree count");
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
                out.value[i][j] = std::complex<T>(precnan<T>(), 0);
                out.derivative[i][j] = std::complex<T>(precnan<T>(), 0);
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
