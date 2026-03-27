#pragma once

// ============================================================================
// PSMS SPHEROIDAL ANGULAR HELPERS
// ============================================================================
// PROLATE SPHEROIDAL WAVE FUNCTION OF THE FIRST KIND [Smn - SCALAR]
// ============================================================================
template<typename T>
struct SmnResult {
    // Angular values S_mn(c, eta) and their eta-derivatives evaluated at the
    // requested arguments.
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
    
    // profcn stores the angular output in degree-major blocks for each eta.
    // Unpack that layout into plain vectors for the single (m, n) request.
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
    // Rows correspond to the requested modal combinations, and each row stores
    // all angular evaluation points for that combination.
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
        if (precabs(e) > T(1)) {
            throw std::invalid_argument("|eta| must be <= 1 for the angular prolate domain.");
        }

    int ioprad = 0;
    int iopnorm = normalize ? 1 : 0;
    int iopang = 2;
    T x1 = static_cast<T>(0);

    SmnMatrixResult<T> out;

    // Scalar case: evaluate one (m, n) pair at every eta value.
    if (m_len == 1 && n_len == 1) {
        if (n[0] < m[0])
            throw std::invalid_argument("'n' must be >= 'm'.");
        SmnResult<T> res = Smn_scalar<T>(m[0], n[0], c, eta, normalize);
        out.value.push_back(res.value);
        out.derivative.push_back(res.derivative);
        return out;
    }

    // Shared-order case: one m with multiple degrees n. A single profcn call
    // can serve the whole degree vector.
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

    // Pairwise case: each m[i] is matched with n[i]. Group by unique m so the
    // expensive profcn setup can be reused within each order block.
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

    // Outer-product case: evaluate every requested order against every
    // requested degree and eta.
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
    // Convert the templated storage layout back to the matrix/vector shapes
    // expected by the public R wrapper.
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
