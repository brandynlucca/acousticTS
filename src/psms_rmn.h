#pragma once

// ============================================================================
// PSMS SPHEROIDAL RADIAL HELPERS
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
