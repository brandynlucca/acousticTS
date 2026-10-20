#pragma once

#include <boost/math/special_functions/bessel.hpp>

template<typename T> struct SpheroidalDirection {
    T cosine = 1, sine = 0;
    bool exact_axis = true;
};

template<typename T> SpheroidalDirection<T> psms_direction_from_cosine(T cosine) {
    if (!is_finite_real(cosine) || precabs(cosine) > T(1))
        throw std::invalid_argument("A spheroidal direction cosine must be finite and within [-1,1].");
    const T sine = precsqrt(std::max(T(0), (T(1) - cosine) * (T(1) + cosine)));
    return {cosine, sine, sine == T(0)};
}

template<typename T> SpheroidalDirection<T> psms_direction_from_angle(T theta) {
    if (!is_finite_real(theta))
        throw std::invalid_argument("A spheroidal direction angle must be finite.");
    const T public_pi = T(std::acos(-1.0)), native_pi = precpi<T>();
    // Recognize exact public cardinal angles without merging adjacent representable directions.
    if (theta == T(0)) return {T(1), T(0), true};
    if (precabs(theta) == public_pi || precabs(theta) == native_pi) return {T(-1), T(0), true};
    if (precabs(theta) == public_pi / T(2) || precabs(theta) == native_pi / T(2))
        return {T(0), T(1), false};
    T sine;
#if ACOUSTICTS_HAVE_QUADMATH
    if constexpr (std::is_same_v<T, acousticts_quad_t>) sine = precabs(sinq(theta));
    else
#endif
        sine = precabs(std::sin(theta));
    return {preccos(theta), sine, sine == T(0)};
}

// ============================================================================
// PSMS BACKSCATTER AND KERNEL HELPERS
// ============================================================================
// RADIAL PROLATE SPHEROIDAL WAVE FUNCTION MATRICES
// ----------------------------------------------------------------------------
template<typename T>
struct RadialMatrixResult {
    // Real radial functions and their xi-derivatives on the retained
    // spheroidal modal grid.
    std::vector<std::vector<T>> value;
    std::vector<std::vector<T>> derivative;
};

template<typename T>
struct RadialComplexMatrixResult {
    // Complex outgoing/incoming radial functions and derivatives on the same
    // modal grid.
    std::vector<std::vector<std::complex<T>>> value;
    std::vector<std::vector<std::complex<T>>> derivative;
};

template<typename T>
struct ExternalRadialResult {
    // Exterior incident and scattered radial ingredients evaluated on the body
    // boundary.
    RadialMatrixResult<T> incident;
    RadialComplexMatrixResult<T> scattering;
};

template<typename T>
struct ExternalBoundaryResult {
    // Exterior angular basis values together with the matching radial terms.
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
    out.value.assign(m_max + 1, std::vector<T>(n_max + 1, precnan<T>()));
    out.derivative.assign(
        m_max + 1, std::vector<T>(n_max + 1, precnan<T>())
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
                precnan<T>(), precnan<T>()
            )
        )
    );
    out.derivative.assign(
        m_max + 1, 
        std::vector<std::complex<T>>(
            n_max + 1, 
            std::complex<T>(
                precnan<T>(), precnan<T>()
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
        m_max + 1, std::vector<T>(n_max + 1, precnan<T>())
    );
    out.derivative.assign(
        m_max + 1, std::vector<T>(n_max + 1, precnan<T>())
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
    // Assemble E^{m(1)} by combining the exterior incident field with the
    // interior regular field at the boundary.
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
    // Assemble E^{m(3)} by combining the exterior scattered field with the
    // interior regular field at the boundary.
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
    // The simplified path discards degree coupling and keeps only the diagonal
    // entries of the incident coupling matrix.
    std::vector<std::vector<std::complex<T>>> lower(m_max + 1);
    for (int m = 0; m <= m_max; ++m) {
        int size = n_max - m + 1;
        lower[m].resize(size, std::complex<T>(precnan<T>(), 0));
        for (int li = 0; li < size; ++li) {
            int ell = m + li;
            T R_ml_1_val = i1_val_mat[m][ell];
            T R_ml_1_d   = i1_der_mat[m][ell];
            int n = ell;
            T Rmn_i1_val = e1_val_mat[m][n];
            T Rmn_i1_d   = e1_der_mat[m][n];
            if (is_na_real(Rmn_i1_val)) {
                lower[m][li] = std::complex<T>(precnan<T>(), 0);
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
            n_max + 1, std::complex<T>(precnan<T>(), 0)
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
    // The simplified path discards degree coupling and keeps only the diagonal
    // entries of the scattered coupling matrix.
    std::vector<std::vector<std::complex<T>>> lower(m_max + 1);
    for (int m = 0; m <= m_max; ++m) {
        int size = n_max - m + 1;
        lower[m].resize(
            size, std::complex<T>(precnan<T>(),
            precnan<T>())
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
            std::complex<T>(precnan<T>(),
            precnan<T>())
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
    // Collapse the incident kernel into the right-hand side of the modal
    // system for each retained order m.
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
std::vector<std::vector<std::complex<T>>> kernel_incident_matrix(
    int m_max,
    int n_max,
    const std::vector<std::vector<std::vector<T>>>& expansion_matrix,
    const std::vector<std::vector<std::complex<T>>>& E1_coupling
) {
    std::vector<std::vector<std::complex<T>>> result(m_max + 1);
    // Build the dense incident-side kernel block for each retained order.
    for (int m = 0; m <= m_max; ++m) {
        int size = n_max - m + 1;
        std::vector<std::complex<T>> mat(size * size, std::complex<T>(0, 0));
        for (int li = 0; li < size; ++li) {
            for (int ni = 0; ni < size; ++ni) {
                int n = m + ni;
                T alpha = expansion_matrix[m][li][ni];
                std::complex<T> E1_nl = E1_coupling[m][li * size + ni];
                mat[li * size + ni] = -imaginary_unit_power<T>(n) * alpha * E1_nl;
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
    // Build the dense scattered-side kernel block for each retained order.
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
std::vector<std::vector<std::vector<T>>> compute_smn_matrices(
    int m_max,
    int n_max,
    T chi_sw,
    const std::vector<T>& eta,
    bool normalize = true
) {
    // Retain the triangular modal layout separately for each angular argument.
    std::vector<std::vector<std::vector<T>>> result(
        eta.size(), std::vector<std::vector<T>>(
            m_max + 1, std::vector<T>(n_max + 1, precnan<T>()))
    );
    if (eta.empty()) return result;

    for (int m = 0; m <= m_max; ++m) {
        int size = n_max - m + 1;
        if (size <= 0) continue;
        int ioprad = 0;
        int iopnorm = normalize ? 1 : 0;
        // Share the native expansion setup across angles and omit unused derivatives.
        int iopang = 1;
        auto profcn_result = cprofcn<T>(
            chi_sw, m, size, eta, ioprad, iopnorm, iopang, T(0)
        );

        for (size_t angle = 0; angle < eta.size(); ++angle) {
            for (int i = 0; i < size; ++i) {
                result[angle][m][m + i] = extract_angular_value_from_batch<T>(
                    profcn_result, i, size, angle
                );
            }
        }
    }
    return result;
}

template<typename T>
std::vector<std::vector<T>> compute_smn_matrix(
    int m_max, int n_max, T chi_sw, T eta_scalar, bool normalize = true
) {
    return compute_smn_matrices<T>(m_max, n_max, chi_sw, {eta_scalar}, normalize)[0];
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

        // Build the overlap matrix that maps interior angular modes onto the
        // exterior basis at the selected quadrature rule.
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
    // Full set of precomputed ingredients for the batched fluid/gas solve.
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
        m_max + 1, std::vector<T>(n_max + 1, precnan<T>())
    );
    out.radial_external.incident.value.assign(
        m_max + 1, std::vector<T>(n_max + 1, precnan<T>())
    );
    out.radial_external.incident.derivative.assign(
        m_max + 1, std::vector<T>(n_max + 1, precnan<T>())
    );
    out.radial_external.scattering.value.assign(
        m_max + 1,
        std::vector<std::complex<T>>(
            n_max + 1,
            std::complex<T>(
                precnan<T>(), precnan<T>()
            )
        )
    );
    out.radial_external.scattering.derivative.assign(
        m_max + 1,
        std::vector<std::complex<T>>(
            n_max + 1,
            std::complex<T>(
                precnan<T>(), precnan<T>()
            )
        )
    );
    out.radial_internal.value.assign(
        m_max + 1, std::vector<T>(n_max + 1, precnan<T>())
    );
    out.radial_internal.derivative.assign(
        m_max + 1, std::vector<T>(n_max + 1, precnan<T>())
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
            chi_sw, m_start, m_count, lnum, ext_args, 2, iopnorm, 1, x1_r
        );
        auto int_result_batch = cprofcn_mblock<T>(
            chi_body, m_start, m_count, lnum, nodes, 1, iopnorm, 1, x1_r
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
    // Right-hand sides and scattered kernels for the per-order dense solves.
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

KernelResult<double> fluid_boundary_kernels_vectorized_double(
    int m_max,
    int n_max,
    double chi_sw,
    double chi_body,
    double xi,
    double density_body,
    double density_sw,
    const std::vector<double>& nodes,
    const std::vector<double>& weights,
    const std::vector<std::vector<double>>& smn_matrix,
    const ExternalRadialResult<double>& radial_external
) {
    auto expansion_matrix = alpha_integration_coefficient<double>(
        m_max, n_max, chi_sw, chi_body, nodes, weights
    );

    auto radial_internal_kind1 = radial_internal_incident_matrix<double>(
        m_max, n_max, chi_body, xi
    );

    KernelResult<double> out;
    out.rhs.resize(m_max + 1);
    out.K3_kernel.resize(m_max + 1);

    for (int m = 0; m <= m_max; ++m) {
        int size = n_max - m + 1;
        if (size <= 0) {
            out.rhs[m] = std::vector<std::complex<double>>();
            out.K3_kernel[m] = std::vector<std::complex<double>>();
            continue;
        }

        arma::Mat<std::complex<double>> alpha(size, size, arma::fill::zeros);
        arma::Mat<std::complex<double>> E1(size, size, arma::fill::zeros);
        arma::Mat<std::complex<double>> E3(size, size, arma::fill::zeros);
        arma::Row<std::complex<double>> phase(size, arma::fill::zeros);

        for (int ni = 0; ni < size; ++ni) {
            int n = m + ni;
            double snm_val = smn_matrix[m][n];
            if (!std::isnan(snm_val) && snm_val != 0.0) {
                phase(ni) = imaginary_unit_power<double>(n) * snm_val;
            }
        }

        for (int li = 0; li < size; ++li) {
            int ell = m + li;
            double R_ml_1_val = radial_internal_kind1.value[m][ell];
            double R_ml_1_d = radial_internal_kind1.derivative[m][ell];
            double factor = (R_ml_1_d != 0.0)
                ? (density_body * R_ml_1_val) / (density_sw * R_ml_1_d)
                : 0.0;

            for (int ni = 0; ni < size; ++ni) {
                alpha(li, ni) = expansion_matrix[m][li][ni];
                E1(li, ni) = std::complex<double>(
                    radial_external.incident.value[m][m + ni] -
                        factor * radial_external.incident.derivative[m][m + ni],
                    0.0
                );
                E3(li, ni) =
                    radial_external.scattering.value[m][m + ni] -
                    factor * radial_external.scattering.derivative[m][m + ni];
            }
        }

        arma::Mat<std::complex<double>> K3 = alpha % E3;
        K3.each_row() %= phase;
        arma::Mat<std::complex<double>> rhs_terms = alpha % E1;
        rhs_terms.each_row() %= phase;
        arma::Col<std::complex<double>> rhs_col = -arma::sum(rhs_terms, 1);

        std::vector<std::complex<double>> rhs(size);
        std::vector<std::complex<double>> K3_vec(static_cast<size_t>(size) * static_cast<size_t>(size));
        for (int i = 0; i < size; ++i) {
            rhs[i] = rhs_col(i);
            for (int j = 0; j < size; ++j) {
                K3_vec[i * size + j] = K3(i, j);
            }
        }

        out.rhs[m] = std::move(rhs);
        out.K3_kernel[m] = std::move(K3_vec);
    }

    return out;
}

// Internal PSMS dense-system layer: backscatter-specific reflection helpers,
// LUP/SVD solve paths, and shared boundary-matrix assembly utilities.
// -----------------------------------------------------------
#include "psms_linear_systems.h"

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
    const std::vector<T>& weights,
    bool gas_like = false
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
    out.Amn_tri = solve_fluid_Amn<T>(kernels.rhs, kernels.K3_kernel, gas_like);
    return out;
}

FluidAmnResult<double> fluid_Amn_triangular_vectorized_double(
    int m_max,
    int n_max,
    double chi_sw,
    double chi_body,
    double xi,
    double eta_body,
    double density_body,
    double density_sw,
    const std::vector<double>& nodes,
    const std::vector<double>& weights,
    bool gas_like = false
) {
    auto external_boundary = external_boundary_matrices<double>(
        m_max, n_max, chi_sw, xi, eta_body, true
    );
    auto kernels = fluid_boundary_kernels_vectorized_double(
        m_max, n_max, chi_sw, chi_body, xi, density_body, density_sw,
        nodes, weights, external_boundary.smn, external_boundary.radial
    );

    FluidAmnResult<double> out;
    out.smn_body = std::move(external_boundary.smn);
    out.Amn_tri = solve_fluid_Amn<double>(kernels.rhs, kernels.K3_kernel, gas_like);
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
    const std::vector<T>& weights,
    bool gas_like = false
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
    auto Amn = solve_fluid_Amn<T>(rhs, K3_kernel, gas_like);

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
    const std::vector<T>& weights,
    bool vectorized = false,
    bool gas_like = false
) {
    if constexpr (is_quad_precision_v<T>) {
        // Quad precision benefits the most from the batched pathway because the
        // spheroidal backend is the dominant cost.
        return fluid_Amn_triangular_batched<T>(
            m_max, n_max, chi_sw, chi_body, xi, eta_body, density_body, density_sw,
            nodes, weights, gas_like
        );
    }
    if constexpr (std::is_same_v<T, double>) {
        if (vectorized) {
            return fluid_Amn_triangular_vectorized_double(
                m_max, n_max, chi_sw, chi_body, xi, eta_body, density_body,
                density_sw, nodes, weights, gas_like
            );
        }
    }
    return fluid_Amn_triangular_from_boundary<T>(
        m_max, n_max, chi_sw, chi_body, xi, eta_body, density_body, density_sw,
        nodes, weights, gas_like
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
                precnan<T>(),
                precnan<T>()
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
                    precnan<T>(), precnan<T>()
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

template<typename T>
std::vector<std::vector<std::complex<T>>> diagonal_tmatrix_blocks_from_expansion(
    int m_max,
    int n_max,
    const std::vector<std::vector<std::complex<T>>>& Amn
) {
    std::vector<std::vector<std::complex<T>>> blocks(m_max + 1);
    // Convert the rectangular A_mn storage into diagonal T-matrix blocks,
    // since the pressure-release and rigid pathways decouple by degree.
    for (int m = 0; m <= m_max; ++m) {
        int size = n_max - m + 1;
        blocks[m].assign(size * size, std::complex<T>(0, 0));
        for (int n = m; n <= n_max; ++n) {
            std::complex<T> val = Amn[m][n];
            if (is_na_real(val.real()) || is_na_real(val.imag())) {
                continue;
            }
            int idx = n - m;
            blocks[m][idx * size + idx] = val;
        }
    }
    return blocks;
}

template<typename T>
std::vector<std::vector<std::complex<T>>> diagonal_tmatrix_blocks_from_triangular(
    int m_max,
    int n_max,
    const std::vector<std::vector<std::complex<T>>>& Amn_tri
) {
    std::vector<std::vector<std::complex<T>>> blocks(m_max + 1);
    // Convert the lower-triangular storage into square diagonal blocks.
    for (int m = 0; m <= m_max; ++m) {
        int size = n_max - m + 1;
        blocks[m].assign(size * size, std::complex<T>(0, 0));
        for (int i = 0; i < size; ++i) {
            std::complex<T> val = Amn_tri[m][i];
            if (is_na_real(val.real()) || is_na_real(val.imag())) {
                continue;
            }
            blocks[m][i * size + i] = val;
        }
    }
    return blocks;
}

template<typename T>
std::vector<std::vector<std::complex<T>>> psms_tmatrix_blocks(
    int m_max,
    int n_max,
    T chi_sw,
    T chi_body,
    T xi,
    T density_body,
    T density_sw,
    const std::vector<T>& nodes,
    const std::vector<T>& weights,
    const std::string& Amn_method
) {
    if (Amn_method == "Amn_fixed_rigid") {
        return diagonal_tmatrix_blocks_from_expansion<T>(
            m_max, n_max,
            fixed_rigid_Amn_expansion_matrix<T>(m_max, n_max, chi_sw, xi)
        );
    }

    if (Amn_method == "Amn_pressure_release") {
        return diagonal_tmatrix_blocks_from_expansion<T>(
            m_max, n_max,
            pressure_release_Amn_expansion_matrix<T>(m_max, n_max, chi_sw, xi)
        );
    }

    if (Amn_method == "Amn_fluid_simplify") {
        // For the full fluid and simplified-fluid retained-block routes, keep
        // all retained diagonal modes active before building the kernels.
        std::vector<std::vector<T>> smn_active(
            m_max + 1, std::vector<T>(n_max + 1, precnan<T>())
        );
        for (int m = 0; m <= m_max; ++m) {
            for (int n = m; n <= n_max; ++n) {
                smn_active[m][n] = T(1);
            }
        }
        auto radial_external = radial_external_matrices<T>(
            m_max, n_max, chi_sw, xi, smn_active
        );
        return diagonal_tmatrix_blocks_from_triangular<T>(
            m_max, n_max,
            simplified_fluid_Amn_triangular<T>(
                m_max, n_max, chi_sw, chi_body, xi, density_body, density_sw,
                radial_external
            )
        );
    }

    std::vector<std::vector<T>> smn_active(
        m_max + 1, std::vector<T>(n_max + 1, precnan<T>())
    );
    for (int m = 0; m <= m_max; ++m) {
        for (int n = m; n <= n_max; ++n) {
            smn_active[m][n] = T(1);
        }
    }
    auto radial_external = radial_external_matrices<T>(
        m_max, n_max, chi_sw, xi, smn_active
    );
    auto expansion_matrix = alpha_integration_coefficient<T>(
        m_max, n_max, chi_sw, chi_body, nodes, weights
    );
    auto radial_internal = radial_internal_incident_matrix<T>(
        m_max, n_max, chi_body, xi
    );
    auto E1_coupling = boundary_coupling_incident_matrix<T>(
        m_max, n_max, density_body, density_sw,
        radial_external.incident.value, radial_external.incident.derivative,
        radial_internal.value, radial_internal.derivative,
        smn_active
    );
    auto E3_coupling = boundary_coupling_scattering_matrix<T>(
        m_max, n_max, density_body, density_sw,
        radial_external.scattering.value, radial_external.scattering.derivative,
        radial_internal.value, radial_internal.derivative,
        smn_active
    );
    auto K1_kernel = kernel_incident_matrix<T>(
        m_max, n_max, expansion_matrix, E1_coupling
    );
    auto K3_kernel = kernel_scattering_matrix<T>(
        m_max, n_max, smn_active, expansion_matrix, E3_coupling
    );

    return solve_fluid_t_blocks<T>(K1_kernel, K3_kernel);
}

template<typename T>
std::complex<T> compute_fbs_from_tmatrix_blocks(
    int m_max,
    int n_max,
    const std::vector<T>& azimuth,
    const std::vector<std::vector<T>>& smn_body_matrix,
    const std::vector<std::vector<T>>& smn_scatter_matrix,
    const std::vector<std::vector<std::complex<T>>>& T_blocks
) {
    std::complex<T> fbs_sum(0, 0), c(0, 0);

    // Multiply each retained T-matrix block by the incident angular vector,
    // then project the result onto the scattered angular basis.
    for (int m = 0; m <= m_max; ++m) {
        int size = n_max - m + 1;
        T nu_m = (m == 0) ? T(1) : T(2);

        for (int li = 0; li < size; ++li) {
            int ell = m + li;
            T smn_scatter = smn_scatter_matrix[m][ell];
            if (is_na_real(smn_scatter)) {
                continue;
            }

            std::complex<T> coeff(0, 0);
            for (int ni = 0; ni < size; ++ni) {
                int n = m + ni;
                T smn_incident = smn_body_matrix[m][n];
                std::complex<T> T_val = T_blocks[m][li * size + ni];
                if (is_na_real(smn_incident) ||
                    is_na_real(T_val.real()) || is_na_real(T_val.imag())) {
                    continue;
                }
                coeff += T_val * smn_incident;
            }

            std::complex<T> term = nu_m * azimuth[m] * smn_scatter * coeff;
            std::complex<T> y = term - c;
            std::complex<T> t = fbs_sum + y;
            c = (t - fbs_sum) - y;
            fbs_sum = t;
        }
    }

    return fbs_sum;
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
    // Relative tolerance used by the adaptive modal stopping rule.
    return static_cast<T>(1e-8L);
}

template<>
inline acousticts_quad_t psms_modal_rel_tol<acousticts_quad_t>() {
    return static_cast<acousticts_quad_t>(1e-14L);
}

template<typename T>
inline T psms_modal_abs_tol() {
    // Absolute floor used by the adaptive modal stopping rule.
    return static_cast<T>(1e-12L);
}

template<>
inline acousticts_quad_t psms_modal_abs_tol<acousticts_quad_t>() {
    return static_cast<acousticts_quad_t>(1e-18L);
}

template<typename T>
inline bool psms_modal_band_is_negligible(
    int m,
    int m_max,
    T band_abs,
    T prev_band_abs,
    T max_band
) {
    // Treat a whole m-band as negligible only after the band has entered the
    // small-signal regime and stopped changing appreciably.
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
    // Apply the same idea within a single m-band to the degree-by-degree tail.
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
            chi_sw, m_start, m_count, lnum, {eta_body}, 2, 1, 1, x1_r
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
            chi_sw, m_start, m_count, lnum, {eta_body}, 2, 1, 1, x1_r
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
std::vector<std::complex<T>> solve_fluid_mode_system_svd(
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
    arma::Col<std::complex<double>> A = solve_fluid_system_svd(K3_arma, b);

    std::vector<std::complex<T>> out(size);
    for (int i = 0; i < size; ++i) {
        out[i] = std::complex<T>(
            static_cast<T>(A(i).real()),
            static_cast<T>(A(i).imag())
        );
    }
    return out;
}

template<typename T>
std::vector<std::complex<T>> solve_fluid_mode_system_default(
    const std::vector<std::complex<T>>& K3_kernel,
    const std::vector<std::complex<T>>& rhs,
    int size
) {
    return solve_fluid_mode_system_svd<T>(K3_kernel, rhs, size);
}

template<>
std::vector<std::complex<acousticts_quad_t>> solve_fluid_mode_system_default<acousticts_quad_t>(
    const std::vector<std::complex<acousticts_quad_t>>& K3_kernel,
    const std::vector<std::complex<acousticts_quad_t>>& rhs,
    int size
) {
#if ACOUSTICTS_HAVE_QUADMATH
    return solve_linear_system_lup_refined<acousticts_quad_t>(K3_kernel, rhs, size);
#else
    return solve_fluid_mode_system_svd<acousticts_quad_t>(K3_kernel, rhs, size);
#endif
}

template<typename T>
std::vector<std::complex<T>> solve_fluid_mode_system_gas(
    const std::vector<std::complex<T>>& K3_kernel,
    const std::vector<std::complex<T>>& rhs,
    int size
) {
    // Gas-filled kernel blocks become much more ill-conditioned than the
    // liquid-filled case because the interior sound-speed contrast is extreme.
    // Use the pseudoinverse solve here so the gas branch can remain stable
    // without perturbing the default liquid-filled pathway.
    return solve_fluid_mode_system_svd<T>(K3_kernel, rhs, size);
}

template<>
std::vector<std::complex<acousticts_quad_t>> solve_fluid_mode_system_gas<acousticts_quad_t>(
    const std::vector<std::complex<acousticts_quad_t>>& K3_kernel,
    const std::vector<std::complex<acousticts_quad_t>>& rhs,
    int size
) {
    return solve_fluid_mode_system_svd<acousticts_quad_t>(K3_kernel, rhs, size);
}

template<typename T>
std::vector<std::complex<T>> solve_fluid_mode_system(
    const std::vector<std::complex<T>>& K3_kernel,
    const std::vector<std::complex<T>>& rhs,
    int size,
    bool gas_like = false
) {
    if (gas_like) {
        return solve_fluid_mode_system_gas<T>(K3_kernel, rhs, size);
    }
    return solve_fluid_mode_system_default<T>(K3_kernel, rhs, size);
}

template<typename T>
int psms_gas_degree_limit_for_m(
    int m,
    int n_max_fortran,
    T chi_body
) {
    T n_limit = static_cast<T>(n_max_fortran);
    if (chi_body > T(60)) {
        n_limit = n_limit -
            T(0.031) * static_cast<T>(m * m) -
            (T(80) / std::max(T(1), static_cast<T>(n_max_fortran))) * static_cast<T>(m);
    } else {
        n_limit = n_limit - static_cast<T>(m);
    }

    int n_end = static_cast<int>(std::ceil(static_cast<double>(n_limit)));
    return std::max(m, n_end);
}

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
    bool adaptive = false,
    bool gas_like = false
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
            chi_sw, m_start, m_count, lnum, ext_args, 2, iopnorm, 1, x1_r
        );
        auto int_result_batch = cprofcn_mblock<T>(
            chi_body, m_start, m_count, lnum, nodes, 1, iopnorm, 1, x1_r
        );

        for (int local_m = 0; local_m < m_count; ++local_m) {
            int m = m_start + local_m;
            int n_end = gas_like
                ? psms_gas_degree_limit_for_m<T>(m, n_max, chi_body)
                : n_max;
            int size = n_end - m + 1;
            if (size <= 0) {
                continue;
            }

            if (adaptive) {
                // Use a cheap proxy to decide whether the modal tail can be
                // truncated before the dense kernel system is assembled.
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

            auto Amn = solve_fluid_mode_system<T>(K3, rhs, size, gas_like);
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

template<typename T>
struct SpheroidalCoupledBlock {
    std::vector<std::complex<T>> matrix;
    std::vector<T> incident_kernel, column_scales, row_scales;
    LUPFactorization<T> factors;
};

template<typename T>
struct SpheroidalCoupledCache {
    bool ready = false;
    bool even_only = false;
    bool coefficient_overlap = false, symmetric_quadrature = false;
    int n_max = -1, m_max = -1;
    T chi_sw = 0, chi_body = 0, xi = 0, density_ratio = 0;
    std::vector<T> nodes, weights;
    std::vector<SpheroidalCoupledBlock<T>> blocks;
};

template<typename T>
T psms_direction_cosine(T theta) {
    return psms_direction_from_angle(theta).cosine;
}

// Solve for scattered coefficients without dividing by incident angular zeros.
template<typename T>
std::complex<T> psms_full_coupled_field(
    int m_max, int n_max, T chi_sw, T chi_body, T xi,
    T eta_body, T eta_scatter, T phi_difference,
    T density_ratio, const std::vector<T>& nodes,
    const std::vector<T>& weights, bool backscatter,
    bool coefficient_overlap = false, bool adaptive = false
) {
    if (m_max < 0 || n_max < m_max) {
        throw std::invalid_argument("Spheroidal modal limits require 0 <= m_max <= n_max.");
    }
    // Reflect symmetric quadrature onto the positive half of the angular domain.
    std::vector<T> qnodes, qweights;
    bool symmetric = true;
    for (size_t k = 0; !coefficient_overlap && k < nodes.size(); ++k) {
        const size_t j = nodes.size() - 1 - k;
        if (nodes[k] != -nodes[j] || weights[k] != weights[j]) {
            symmetric = false;
        }
    }
    for (size_t k = 0; !coefficient_overlap && k < nodes.size(); ++k) {
        if (symmetric && nodes[k] < T(0)) continue;
        qnodes.push_back(nodes[k]);
        const T multiplier = symmetric && nodes[k] > T(0) ? T(2) : T(1);
        qweights.push_back(weights[k] * multiplier);
    }
    // Backscatter uses angular parity, so its receive angle needs no evaluation.
    std::vector<T> args = {eta_body};
    if (!backscatter) args.push_back(eta_scatter);
    const size_t quadrature_offset = args.size();
    const std::vector<T> direction_arguments = args;
    args.insert(args.end(), qnodes.begin(), qnodes.end());
    std::complex<T> total(0, 0), correction(0, 0);
    const bool truncate_tails = adaptive && backscatter && !coefficient_overlap;
    T max_band = T(0), previous_band = T(-1);
    int small_bands = 0;
    const int chunk = psms_profcn_mblock_chunk_size<T>();
    if (precabs(eta_body) == T(1)) m_max = 0;
    const bool even_only = eta_body == T(0);
    // Retain one exact physical system within a 64 MiB allocation budget.
    static thread_local SpheroidalCoupledCache<T> cache;
    const bool cache_supported = coefficient_overlap || is_quad_precision_v<T>;
    size_t cache_bytes = sizeof(SpheroidalCoupledCache<T>) +
        2 * (m_max + 1) * sizeof(SpheroidalCoupledBlock<T>);
    if (!coefficient_overlap) cache_bytes += (nodes.size() + weights.size()) * sizeof(T);
    for (int m = 0; cache_supported && m <= m_max; ++m) {
        for (int parity = 0; parity < 2; ++parity) {
            if (even_only && parity == 1) continue;
            const size_t dim = (n_max - m - parity + 2) / 2;
            cache_bytes += dim * dim * (2 * sizeof(std::complex<T>) + sizeof(T)) +
                dim * (2 * sizeof(T) + 2 * sizeof(int));
        }
    }
    const bool retain_system = cache_supported && cache_bytes <= psms_retained_cache_budget_bytes;
    const bool reuse_system = retain_system && cache.ready &&
        cache.coefficient_overlap == coefficient_overlap &&
        cache.symmetric_quadrature == symmetric &&
        (coefficient_overlap || (cache.nodes == nodes && cache.weights == weights)) &&
        cache.even_only == even_only && cache.n_max == n_max &&
        cache.m_max == m_max && cache.chi_sw == chi_sw && cache.chi_body == chi_body &&
        cache.xi == xi && cache.density_ratio == density_ratio;
    if (!reuse_system) {
        // Release old keys and block capacity before applying the retained-memory budget.
        cache = SpheroidalCoupledCache<T>();
        if (retain_system) cache.blocks.resize(2 * (m_max + 1));
    }
    for (int start = 0; start <= m_max; start += chunk) {
        // Share special-function setup across neighboring azimuthal orders.
        Rcpp::checkUserInterrupt();
        const int count = std::min(chunk, m_max - start + 1);
        const int lnum = n_max - start + 1;
        ProfcnBatchResult<T> external, internal;
        ProfcnCoefficientBatchResult<T> external_coefficients, internal_coefficients;
        if (reuse_system) {
            // New directions need only the unchanged native angular evaluator.
            external = cprofcn_mblock<T>(
                chi_sw, start, count, lnum, direction_arguments, 0, 1, 1, xi - T(1)
            );
        } else if (coefficient_overlap) {
            // Use normalized Legendre expansions for the full fluid angular overlaps.
            external_coefficients = cprofcn_mblock_coefficients<T>(
                chi_sw, start, count, lnum, args, 2, 1, xi - T(1)
            );
            internal_coefficients = cprofcn_mblock_coefficients<T>(
                chi_body, start, count, lnum, {}, 1, 0, xi - T(1)
            );
        } else {
            external = cprofcn_mblock<T>(
                chi_sw, start, count, lnum, args, 2, 1, 1, xi - T(1)
            );
            internal = cprofcn_mblock<T>(
                chi_body, start, count, lnum, qnodes, 1, 1, 1, xi - T(1)
            );
        }
        const auto& external_functions = coefficient_overlap && !reuse_system ?
            external_coefficients.functions : external;
        const auto& internal_functions = coefficient_overlap ?
            internal_coefficients.functions : internal;
        for (int local = 0; local < count; ++local) {
            const int m = start + local, size = n_max - m + 1;
            const T azimuth = backscatter ? (m % 2 ? T(-1) : T(1)) :
                preccos(T(m) * phi_difference);
            const T multiplicity = m == 0 ? T(1) : T(2);
            std::vector<std::complex<T>> outgoing_terms(truncate_tails ? size : 0);
            // Solve even and odd degree blocks without discarding coupled modes.
            for (int parity = 0; parity < 2; ++parity) {
                // Odd angular modes have zero broadside forcing (DLMF 30.4.3).
                if (even_only && parity == 1) continue;
                const int dim = (size - parity + 1) / 2;
                if (dim <= 0) continue;
                std::vector<T> incident_angular(dim), receive_angular(dim);
                std::vector<T> interior_value(dim), interior_derivative(dim);
                std::vector<std::complex<T>> external_value(dim), external_derivative(dim);
                std::vector<std::complex<T>> incident_coefficients(dim);
                std::vector<T> external_angular(reuse_system ? 0 : dim * qnodes.size());
                std::vector<T> weighted_internal_angular(reuse_system ? 0 : dim * qnodes.size());
                for (int i = 0; i < dim; ++i) {
                    // Cache angular samples and radial boundary values for reuse.
                    const int index = parity + 2 * i, degree = m + index;
                    incident_angular[i] = extract_angular_value_from_mblock<T>(
                        external_functions, local, index, 0
                    );
                    receive_angular[i] = backscatter ?
                        (index % 2 ? -incident_angular[i] : incident_angular[i]) :
                        extract_angular_value_from_mblock<T>(external_functions, local, index, 1);
                    incident_coefficients[i] = imaginary_unit_power<T>(degree) * incident_angular[i];
                    if (reuse_system) continue;
                    auto external_radial = extract_radial_from_mblock<T>(external_functions, local, index, 3);
                    auto internal_radial = extract_radial_from_mblock<T>(internal_functions, local, index, 1);
                    external_value[i] = {external_radial.val_real, external_radial.val_imag};
                    external_derivative[i] = {external_radial.der_real, external_radial.der_imag};
                    // Scale the interior pair before multiplying radial functions.
                    const T scale = std::max(
                        precabs(internal_radial.der_real),
                        precabs(density_ratio * internal_radial.val_real)
                    );
                    if (!(scale > T(0))) {
                        throw std::runtime_error(
                            "Interior spheroidal radial values underflowed at the retained modal order."
                        );
                    }
                    interior_value[i] = density_ratio * internal_radial.val_real / scale;
                    interior_derivative[i] = internal_radial.der_real / scale;
                    for (size_t k = 0; k < qnodes.size(); ++k) {
                        external_angular[i * qnodes.size() + k] =
                            extract_angular_value_from_mblock<T>(external, local, index, k + quadrature_offset);
                        weighted_internal_angular[i * qnodes.size() + k] =
                            extract_angular_value_from_mblock<T>(internal, local, index, k) * qweights[k];
                    }
                }
                std::vector<std::complex<T>> matrix(reuse_system ? 0 : dim * dim), rhs(dim, {0, 0});
                std::vector<T> incident_kernel(retain_system && !reuse_system ? dim * dim : 0);
                std::vector<T> scales(dim), row_scales(dim);
                std::vector<std::complex<T>> scattered_coefficients;
                if (reuse_system) {
                    // Rebuild the angle-dependent forcing with its original summation order.
                    const auto& block = cache.blocks[2 * m + parity];
                    for (int row = 0; row < dim; ++row) {
                        for (int col = 0; col < dim; ++col) {
                            rhs[row] -= block.incident_kernel[row * dim + col] *
                                incident_coefficients[col];
                        }
                        rhs[row] /= block.row_scales[row];
                    }
                    scales = block.column_scales;
                    scattered_coefficients = solve_factored_system_lup_refined<T>(
                        block.matrix, rhs, block.factors, dim, 2
                    );
                } else {
                    // Match pressure and normal velocity in the interior angular basis.
                    for (int row = 0; row < dim; ++row) {
                        for (int col = 0; col < dim; ++col) {
                            T overlap = coefficient_overlap ? spheroidal_coefficient_overlap<T>(
                                external_coefficients, internal_coefficients, local,
                                parity + 2 * col, parity + 2 * row
                            ) : T(0);
                            for (size_t k = 0; k < qnodes.size(); ++k) {
                                overlap += external_angular[col * qnodes.size() + k] *
                                    weighted_internal_angular[row * qnodes.size() + k];
                            }
                            matrix[row * dim + col] = overlap * (
                                external_value[col] * interior_derivative[row] -
                                interior_value[row] * external_derivative[col]
                            );
                            // Real scalar products give the identical real boundary expression.
                            const T forcing = matrix[row * dim + col].real();
                            rhs[row] -= forcing * incident_coefficients[col];
                            if (retain_system) incident_kernel[row * dim + col] = forcing;
                        }
                    }
                    // Balance radial column magnitudes, then boundary-equation rows.
                    for (int col = 0; col < dim; ++col) {
                        T scale = T(0);
                        for (int row = 0; row < dim; ++row) {
                            const auto& value = matrix[row * dim + col];
                            scale = std::max(scale, precabs(value.real()) + precabs(value.imag()));
                        }
                        if (!(scale > T(0))) {
                            throw std::runtime_error(
                                "Spheroidal coupling column is numerically zero."
                            );
                        }
                        scales[col] = scale;
                        for (int row = 0; row < dim; ++row) matrix[row * dim + col] /= scale;
                    }
                    for (int row = 0; row < dim; ++row) {
                        T scale = T(0);
                        for (int col = 0; col < dim; ++col) {
                            const auto& value = matrix[row * dim + col];
                            scale = std::max(scale, precabs(value.real()) + precabs(value.imag()));
                        }
                        if (!(scale > T(0))) {
                            throw std::runtime_error(
                                "Spheroidal coupling row is numerically zero. Check modal convergence and precision."
                            );
                        }
                        for (int col = 0; col < dim; ++col) matrix[row * dim + col] /= scale;
                        rhs[row] /= scale;
                        row_scales[row] = scale;
                    }
                    auto factors = lup_decompose<T>(matrix, dim);
                    scattered_coefficients = solve_factored_system_lup_refined<T>(matrix, rhs, factors, dim, 2);
                    if (retain_system) {
                        auto& block = cache.blocks[2 * m + parity];
                        block.matrix = std::move(matrix);
                        block.incident_kernel = std::move(incident_kernel);
                        block.column_scales = scales;
                        block.row_scales = std::move(row_scales);
                        block.factors = std::move(factors);
                    }
                }
                // Undo column scaling and sum the receive field with compensation.
                for (int i = 0; i < dim; ++i) {
                    const int degree = m + parity + 2 * i;
                    const auto term = multiplicity * azimuth * receive_angular[i] *
                        imaginary_unit_power<T>(-degree) * (scattered_coefficients[i] / scales[i]);
                    if (truncate_tails) {
                        outgoing_terms[parity + 2 * i] = term;
                        continue;
                    }
                    const auto y = term - correction;
                    const auto updated = total + y;
                    correction = (updated - total) - y;
                    total = updated;
                }
            }
            if (truncate_tails) {
                // Apply the existing tail rule after solving every retained coupled degree.
                std::complex<T> band_sum(0, 0), band_correction(0, 0);
                T max_term = T(0), previous_term = T(-1);
                int small_terms = 0;
                for (int i = 0; i < size; ++i) {
                    const auto& term = outgoing_terms[i];
                    const T term_abs = precsqrt(complex_abs_sq(term));
                    max_term = std::max(max_term, term_abs);
                    if (psms_modal_term_is_negligible<T>(
                            i, size - 1, term_abs, previous_term, max_term)) {
                        ++small_terms;
                    } else {
                        small_terms = 0;
                    }
                    previous_term = term_abs;
                    const auto adjusted = term - band_correction;
                    const auto updated = band_sum + adjusted;
                    band_correction = (updated - band_sum) - adjusted;
                    band_sum = updated;
                    if (small_terms >= 2) break;
                }
                const auto adjusted = band_sum - correction;
                const auto updated = total + adjusted;
                correction = (updated - total) - adjusted;
                total = updated;
                const T band_abs = precsqrt(complex_abs_sq(band_sum));
                max_band = std::max(max_band, band_abs);
                if (psms_modal_band_is_negligible<T>(
                        m, m_max, band_abs, previous_band, max_band)) {
                    ++small_bands;
                } else {
                    small_bands = 0;
                }
                previous_band = band_abs;
                // A newly populated partial cache stays unready when later orders are omitted.
                if (small_bands >= 2) return total;
            }
        }
    }
    if (retain_system && !reuse_system) {
        cache.coefficient_overlap = coefficient_overlap;
        cache.symmetric_quadrature = symmetric;
        if (!coefficient_overlap) {
            cache.nodes = nodes;
            cache.weights = weights;
        }
        cache.n_max = n_max;
        cache.even_only = even_only;
        cache.m_max = m_max;
        cache.chi_sw = chi_sw;
        cache.chi_body = chi_body;
        cache.xi = xi;
        cache.density_ratio = density_ratio;
        cache.ready = true;
    }
    return total;
}

namespace psms_automatic {
template<typename T> T sine(T value) {
#if ACOUSTICTS_HAVE_QUADMATH
    if constexpr (std::is_same_v<T, acousticts_quad_t>) return sinq(value);
    else
#endif
    return std::sin(value);
}

template<typename T> T logarithm(T value) {
#if ACOUSTICTS_HAVE_QUADMATH
    if constexpr (std::is_same_v<T, acousticts_quad_t>) return logq(value);
    else
#endif
    return std::log(value);
}

template<typename T> T exponential(T value) {
#if ACOUSTICTS_HAVE_QUADMATH
    if constexpr (std::is_same_v<T, acousticts_quad_t>) return expq(value);
    else
#endif
    return std::exp(value);
}

template<typename T> T bessel_j(int order, T argument) {
#if ACOUSTICTS_HAVE_QUADMATH
    if constexpr (std::is_same_v<T, acousticts_quad_t>) return jnq(order, argument);
    else
#endif
    return boost::math::cyl_bessel_j(order, argument);
}

template<typename T> void add(T& sum, T& correction, const T& value) {
    const T adjusted = value - correction;
    const T updated = sum + adjusted;
    correction = (updated - sum) - adjusted;
    sum = updated;
}

template<typename T> void finite(const std::complex<T>& value) {
    if (!is_finite_real(value.real()) || !is_finite_real(value.imag()))
        throw std::runtime_error("Full-fluid boundary arithmetic is nonfinite.");
}

template<typename T> T finite_norm(T value) {
    if (!is_finite_real(value) || value < T(0))
        throw std::runtime_error("An automatic spheroidal norm or error budget exceeds the native numerical range.");
    return value;
}

template<typename T> T numerical_budget(T reference) {
    return finite_norm(psms_modal_abs_tol<T>() + psms_modal_rel_tol<T>() * finite_norm(reference));
}

template<typename T> void add_field(std::complex<T>& sum, std::complex<T>& correction,
        const std::complex<T>& value) {
    finite(value);
    add(sum, correction, value);
    // Finite terms can still overflow their compensated accumulated field.
    finite(sum);
    finite(correction);
}

template<typename T> T norm(const std::vector<std::complex<T>>& values) {
    // Scaling the norm avoids squaring subnormal or very large coefficients.
    T scale = 0, sum = 1;
    for (const auto& value : values) {
        finite(value);
        for (T component : {precabs(value.real()), precabs(value.imag())}) {
            if (component == T(0)) continue;
            if (scale < component) {
                const T ratio = scale / component;
                sum = T(1) + sum * ratio * ratio;
                scale = component;
            } else {
                const T ratio = component / scale;
                sum += ratio * ratio;
            }
        }
    }
    return finite_norm(scale == T(0) ? T(0) : scale * precsqrt(sum));
}

template<typename T> T coefficient(const ProfcnCoefficientBatchResult<T>& batch,
        int degree, int legendre) {
    if (legendre >= batch.lengths[degree]) return T(0);
    return batch.coefficients[legendre + size_t(batch.capacity) * degree];
}

inline int grow_degree(int degree, int order) {
    // A geometric search schedules retries without deciding modal acceptance.
    const long long increment = std::max(2LL, (static_cast<long long>(degree) - order + 1) / 4);
    if (degree > INT_MAX - increment)
        throw std::length_error("Modal refinement exceeds the native degree range.");
    checked_profcn_integer(2 * (static_cast<long long>(degree) + increment) + 1,
        "twice the refined degree plus one");
    return degree + static_cast<int>(increment);
}

inline std::pair<int, int> reuse_interior_bracket(int degree, int order, int lower, int upper) {
    if (order < 0 || degree < order || lower < order || upper <= lower)
        throw std::invalid_argument("An interior convergence bracket needs distinct ordered supports.");
    const int first = std::max(degree, lower);
    const int second = std::max(degree, upper);
    return {first, second > first ? second : grow_degree(first, order)};
}

template<typename T> struct BoundaryCheck {
    T pressure_error = 0, derivative_error = 0;
    T pressure_reference = 0, derivative_reference = 0;
    T incident_pressure_error = 0, incident_derivative_error = 0;
    T quadrature_pressure_change = 0, quadrature_derivative_change = 0;
    T pressure_tail_bound = 0, derivative_tail_bound = 0;
    T incident_pressure_tail_bound = 0, incident_derivative_tail_bound = 0;
    T incident_pressure_recurrence_error = 0, incident_derivative_recurrence_error = 0;
    int incident_degree = -1, incident_recurrence_start = -1;
    std::string residual_method = "unassessed", analytic_failure;
    int quadrature_order = 0, trace_prefix = 0;
    bool quadrature_converged = false;
};

template<typename T> void validate_boundary_errors(const BoundaryCheck<T>& check) {
    finite_norm(check.pressure_error);
    finite_norm(check.derivative_error);
    finite_norm(check.incident_pressure_error);
    finite_norm(check.incident_derivative_error);
}

template<typename T> struct OrderSolution {
    int m = 0, degree = 0, interior_degree = 0, attempts = 0;
    int interior_comparison_lower = -1, interior_comparison_upper = -1;
    int overlap_quadrature_order = 0;
    std::string overlap_method = "coefficients";
    bool robin = false, singular = false;
    T gram_frobenius = 0, coefficient_change = 0;
    int minimum_radial_accuracy = INT_MAX, minimum_angular_accuracy = INT_MAX;
    std::complex<T> field = {T(0), T(0)};
    T embedded_roundoff_scale = 0;
    std::array<std::complex<T>, 2> embedded_fields{};
    std::array<bool, 2> embedded_valid{{true, true}};
    std::vector<std::complex<T>> coefficients, scaled_coefficients, outgoing_terms;
    std::vector<long long> radial_exponents;
    std::array<std::vector<std::complex<T>>, 6> traces;
    BoundaryCheck<T> boundary;
    psms_packed::ScalingDiagnostics scaling;
    size_t working_bytes = 0;
    size_t overlap_entries_reused = 0, overlap_entries_computed = 0;
};

template<typename T> struct Result {
    bool converged = false, matched_medium = false;
    bool modal_tail_truncated = false;
    T modal_tail_change = 0;
    std::complex<T> field = {precnan<T>(), precnan<T>()};
    int last_order = -1, last_degree = -1, last_interior_degree = -1;
    T pressure_relative = 0, derivative_relative = 0;
    T omitted_pressure_bound = 0, omitted_derivative_bound = 0;
    T added_order_norm = 0;
    // C++ allocation estimates exclude private workspaces inside the Fortran evaluator.
    size_t retained_bytes = 0, peak_working_bytes = 0;
    BoundaryCheck<T> last_boundary;
    std::string detail;
    std::vector<OrderSolution<T>> orders;
};

template<typename T> struct LegendreWorkspace {
    int m = -1;
    size_t support = 0;
    int scale_step = 0;
    T upper_scale = 0, lower_scale = 0;
    std::vector<T> diagonal_factors, first_factors, second_factors;

    void prepare(int order, size_t count) {
        if (order == m && count == support) return;
        m = order;
        support = count;
        diagonal_factors.resize(m);
        first_factors.resize(count);
        second_factors.resize(count);
        // Recurrence factors depend on degree and order, never on the quadrature node.
        for (int index = 1; index <= m; ++index)
            diagonal_factors[index - 1] = precsqrt((T(2) * index + T(1)) / (T(2) * index));
        for (size_t index = 0; index < count; ++index) {
            const T degree = T(m) + T(index) + T(1);
            const T denominator = degree * degree - T(m) * T(m);
            first_factors[index] = precsqrt((T(4) * degree * degree - T(1)) / denominator);
            second_factors[index] = index == 0 ? T(0) : precsqrt(((T(2) * degree + T(1)) *
                ((degree - T(1)) * (degree - T(1)) - T(m) * T(m))) /
                ((T(2) * degree - T(3)) * denominator));
        }
        scale_step = std::numeric_limits<T>::max_exponent10 / 2;
#if ACOUSTICTS_HAVE_QUADMATH
        if constexpr (std::is_same_v<T, acousticts_quad_t>) scale_step = FLT128_MAX_10_EXP / 2;
#endif
        upper_scale = pow10_typed<T>(scale_step);
        lower_scale = pow10_typed<T>(-scale_step);
    }
};

template<typename T> std::vector<T> normalized_ferrers_values(int m,
        const SpheroidalDirection<T>& direction, size_t count, LegendreWorkspace<T>& workspace) {
    if (m < 0 || !is_finite_real(direction.cosine) || precabs(direction.cosine) > T(1) ||
            !is_finite_real(direction.sine) || direction.sine < T(0) || direction.sine > T(1) ||
            (direction.exact_axis && (direction.sine != T(0) || precabs(direction.cosine) != T(1))))
        throw std::invalid_argument("The normalized Ferrers direction is outside its native domain.");
    workspace.prepare(m, count);
    std::vector<T> values(count);
    if (direction.exact_axis) {
        // Exact Ferrers endpoint values avoid an unnecessary long recurrence at a true axis.
        if (m == 0)
            for (size_t degree = 0; degree < count; ++degree)
                values[degree] = (direction.cosine < T(0) && degree % 2 ? T(-1) : T(1)) *
                    precsqrt((T(2) * degree + T(1)) / T(2));
        return values;
    }
    T diagonal = T(1) / precsqrt(T(2));
    long long exponent = 0;
    // The positive unit-norm Ferrers diagonal uses the original sine, including near an axis.
    for (int order = 1; order <= m; ++order) {
        diagonal *= direction.sine * workspace.diagonal_factors[order - 1];
        if (diagonal != T(0) && precabs(diagonal) < workspace.lower_scale) {
            diagonal *= workspace.upper_scale;
            exponent -= workspace.scale_step;
        }
    }
    T previous = 0, current = diagonal;
    for (size_t index = 0; index < count; ++index) {
        values[index] = psms_packed::scale_decimal(current, exponent);
        if (!is_finite_real(values[index]))
            throw std::runtime_error("A normalized Ferrers value exceeds the native numerical range.");
        const T next = workspace.first_factors[index] * direction.cosine * current -
            workspace.second_factors[index] * previous;
        previous = current;
        current = next;
        const T magnitude = std::max(precabs(previous), precabs(current));
        if (magnitude > workspace.upper_scale) {
            previous *= workspace.lower_scale;
            current *= workspace.lower_scale;
            exponent += workspace.scale_step;
        } else if (magnitude != T(0) && magnitude < workspace.lower_scale) {
            previous *= workspace.upper_scale;
            current *= workspace.upper_scale;
            exponent -= workspace.scale_step;
        }
    }
    return values;
}

template<typename T> void evaluate_coefficient_angular_values(ProfcnCoefficientBatchResult<T>& batch,
        int m, const std::vector<SpheroidalDirection<T>>& directions, LegendreWorkspace<T>& workspace) {
    auto& functions = batch.functions;
    if (functions.m_count != 1 || directions.size() != size_t(functions.narg))
        throw std::invalid_argument("An angular coefficient block requires one order and its complete direction list.");
    // All normalized expansion coefficients contribute to each angular value without an angle cutoff.
    for (size_t angle = 0; angle < directions.size(); ++angle) {
        const auto basis = normalized_ferrers_values(m, directions[angle], batch.capacity, workspace);
        for (int degree = 0; degree < functions.lnum; ++degree) {
            const size_t offset = degree + size_t(functions.lnum) * angle;
            const T* coefficients = batch.coefficients.data() + size_t(batch.capacity) * degree;
            T value = 0;
#if ACOUSTICTS_HAVE_EXACT_QUAD_DOT
            if constexpr (std::is_same_v<T, acousticts_quad_t>) {
                static thread_local quad_dot_product::Accumulator accumulator;
                value = quad_dot_product::dot(coefficients, basis.data(), batch.lengths[degree], accumulator);
            }
            else
#endif
            {
                T correction = 0;
                for (int index = 0; index < batch.lengths[degree]; ++index)
                    add(value, correction, coefficients[index] * basis[index]);
            }
            if (!is_finite_real(value))
                throw std::runtime_error("A coefficient angular expansion is nonfinite in native precision.");
            functions.s1c[offset] = value;
            functions.is1e[offset] = 0;
            // The disabled Fortran angular evaluator supplies no accuracy estimate for this sum.
            functions.naccs[offset] = NA_INTEGER;
        }
    }
}

template<typename T> struct IncidentSequence {
    T argument = T(-1);
    int degree = -1, start = 0;
    T pressure_change = 0, radial_change = 0, angular_change = 0;
    std::vector<T> values, derivatives, logarithmic_bounds;

    void reset(T value) {
        if (argument == value) return;
        argument = value;
        degree = -1;
        start = 0;
        values.clear();
        derivatives.clear();
        logarithmic_bounds.assign(1, T(0));
    }

    size_t bytes() const {
        return (values.capacity() + derivatives.capacity() + logarithmic_bounds.capacity()) * sizeof(T);
    }
};

template<typename T> struct Quadrature {
    int order = 0;
    int previous_order = 0;
    std::vector<T> nodes, weights;
    std::vector<T> previous_nodes, previous_weights;
    LegendreWorkspace<T> legendre;
    IncidentSequence<T> incident;

    size_t bytes() const {
        return (nodes.capacity() + weights.capacity() + previous_nodes.capacity() +
            previous_weights.capacity() + legendre.diagonal_factors.capacity() +
            legendre.first_factors.capacity() + legendre.second_factors.capacity()) * sizeof(T) + incident.bytes();
    }

    void prepare(int count) {
        if (count == order) return;
        if (count < 1) throw std::length_error("Boundary quadrature order is invalid.");
        // Keep the two exact rules used by successive residual refinements.
        if (count == previous_order) {
            std::swap(order, previous_order);
            nodes.swap(previous_nodes);
            weights.swap(previous_weights);
            return;
        }
        previous_order = order;
        previous_nodes.swap(nodes);
        previous_weights.swap(weights);
        order = 0;
        nodes.resize(count);
        weights.resize(count);
        const T epsilon = preceps<T>();
        const int iterations = is_quad_precision_v<T> ? 113 : std::numeric_limits<T>::digits;
        for (int index = 0; index < count / 2 + count % 2; ++index) {
            // Newton refinement preserves native precision and exact node symmetry.
            T root = 2 * index + 1 == count ? T(0) :
                preccos(precpi<T>() * (T(index) + T(3) / T(4)) / (T(count) + T(1) / T(2)));
#if ACOUSTICTS_HAVE_QUADMATH
            if constexpr (std::is_same_v<T, acousticts_quad_t>) {
                // Refine a cheap initial guess before applying the full native Newton iteration.
                if (root != T(0)) {
                    double guess = static_cast<double>(root), previous_guess = NAN;
                    for (int iteration = 0; iteration < std::numeric_limits<double>::digits; ++iteration) {
                        double lower = 1, current = guess;
                        for (int degree = 2; degree <= count; ++degree) {
                            const double next = ((2.0 * degree - 1) * guess * current -
                                (degree - 1.0) * lower) / degree;
                            lower = current;
                            current = next;
                        }
                        const double derivative = count * (guess * current - lower) / (guess * guess - 1);
                        const double updated = guess - current / derivative;
                        if (!std::isfinite(updated) || !(updated > 0 && updated < 1)) break;
                        if (updated == guess || updated == previous_guess ||
                                std::abs(updated - guess) <= std::numeric_limits<double>::epsilon()) {
                            root = T(updated);
                            break;
                        }
                        previous_guess = guess;
                        guess = updated;
                    }
                }
            }
#endif
            T previous = precnan<T>(), derivative = 0;
            bool converged = false;
            for (int iteration = 0; iteration < iterations; ++iteration) {
                T lower = 1, current = root;
                for (int degree = 2; degree <= count; ++degree) {
                    const T next = ((T(2) * degree - T(1)) * root * current -
                        (T(degree) - T(1)) * lower) / T(degree);
                    lower = current;
                    current = next;
                }
                derivative = T(count) * (root * current - lower) / (root * root - T(1));
                const T updated = root - current / derivative;
                if (updated == root || updated == previous || precabs(updated - root) <= epsilon) {
                    root = updated;
                    converged = true;
                    break;
                }
                previous = root;
                root = updated;
            }
            if (!converged || !is_finite_real(root))
                throw std::runtime_error("Native boundary quadrature roots did not converge.");
            T lower = 1, current = root;
            for (int degree = 2; degree <= count; ++degree) {
                const T next = ((T(2) * degree - T(1)) * root * current -
                    (T(degree) - T(1)) * lower) / T(degree);
                lower = current;
                current = next;
            }
            derivative = T(count) * (root * current - lower) / (root * root - T(1));
            nodes[index] = -root;
            nodes[count - index - 1] = root;
            const T weight = T(2) / ((T(1) - root * root) * derivative * derivative);
            weights[index] = weight;
            weights[count - index - 1] = weight;
        }
        for (int index = 0; index < count; ++index) {
            if (!(weights[index] > T(0)) || !is_finite_real(weights[index]) ||
                    !(nodes[index] > T(-1) && nodes[index] < T(1)) ||
                    (index > 0 && !(nodes[index] > nodes[index - 1])))
                throw std::runtime_error("Native boundary quadrature lost a distinct node or positive weight.");
        }
        order = count;
    }
};

template<typename T> struct OverlapQuadrature {
    const std::vector<T>* supplied_nodes = nullptr;
    const std::vector<T>* supplied_weights = nullptr;
    bool adaptive = false;
    int exact_order = 0;

    bool active() const { return supplied_nodes != nullptr; }

    void validate() const {
        if ((supplied_nodes == nullptr) != (supplied_weights == nullptr))
            throw std::invalid_argument("Fluid overlap quadrature requires both nodes and weights.");
        if (!active() || adaptive) return;
        if (supplied_nodes->empty() || supplied_nodes->size() != supplied_weights->size() ||
                supplied_nodes->size() >= size_t(INT_MAX))
            throw std::invalid_argument("Prescribed fluid overlap quadrature has invalid dimensions.");
        for (size_t index = 0; index < supplied_nodes->size(); ++index)
            if (!is_finite_real((*supplied_nodes)[index]) || precabs((*supplied_nodes)[index]) > T(1) ||
                    !is_finite_real((*supplied_weights)[index]))
                throw std::invalid_argument("Prescribed fluid overlap nodes and weights must be finite on [-1,1].");
    }

    const std::vector<T>& nodes() const { return *supplied_nodes; }
    const std::vector<T>& weights() const { return *supplied_weights; }

    void prepare(int m, int external_length, int internal_length) {
        if (!adaptive) return;
        if (m < 0 || external_length < 1 || internal_length < 1)
            throw std::invalid_argument("Finite angular expansion dimensions are invalid.");
        // Products of Ferrers expansions are polynomials of degree at most L_ext + L_int.
        const long long product_degree = 2LL * m + external_length + internal_length - 2;
        const long long count = (product_degree + 2) / 2;
        if (count < 1 || count >= INT_MAX)
            throw std::length_error("Degree-exact fluid overlap quadrature exceeds the native index range.");
        // DLMF 3.5(v) and 14.17.6 reduce this exact finite-expansion Gauss rule to a coefficient dot.
        exact_order = static_cast<int>(count);
    }

    size_t rank_bound(int m, int parity) const {
        std::vector<T> arguments;
        arguments.reserve(nodes().size());
        for (size_t index = 0; index < nodes().size(); ++index) {
            const T argument = precabs(nodes()[index]);
            if (weights()[index] == T(0) || (parity && argument == T(0)) ||
                    (m > 0 && argument == T(1))) continue;
            arguments.push_back(argument);
        }
        // Opposite nodes give proportional samples within a fixed angular parity.
        std::sort(arguments.begin(), arguments.end());
        return size_t(std::unique(arguments.begin(), arguments.end()) - arguments.begin());
    }

    size_t bytes() const { return 0; }
};

template<typename T> T expansion_at_node(const ProfcnCoefficientBatchResult<T>& batch,
        int mode, const std::vector<T>& basis) {
    const int parity = mode % 2, length = batch.lengths[mode];
    const T* coefficients = batch.coefficients.data() + size_t(batch.capacity) * mode;
#if ACOUSTICTS_HAVE_EXACT_QUAD_DOT
    if constexpr (std::is_same_v<T, acousticts_quad_t>) {
        static thread_local quad_dot_product::Accumulator accumulator;
        const size_t count = length > parity ? size_t(length - parity + 1) / 2 : 0;
        return quad_dot_product::dot(coefficients + parity, basis.data() + parity,
            count, accumulator, 2, 2);
    }
#endif
    T value = 0, correction = 0;
    for (int index = parity; index < length; index += 2)
        add(value, correction, coefficients[index] * basis[index]);
    return value;
}

template<typename T> struct OverlapSamples {
    size_t count = 0;
    std::vector<T> exterior, weighted_interior;

    T integral(int external_mode, int internal_mode) const {
        const T* first = exterior.data() + size_t(external_mode) * count;
        const T* second = weighted_interior.data() + size_t(internal_mode) * count;
#if ACOUSTICTS_HAVE_EXACT_QUAD_DOT
        if constexpr (std::is_same_v<T, acousticts_quad_t>) {
            static thread_local quad_dot_product::Accumulator accumulator;
            return quad_dot_product::dot(first, second, count, accumulator);
        }
#endif
        T value = 0, correction = 0;
        for (size_t node = 0; node < count; ++node) add(value, correction, first[node] * second[node]);
        return value;
    }

    size_t bytes() const { return (exterior.capacity() + weighted_interior.capacity()) * sizeof(T); }
};

template<typename T> OverlapSamples<T> sample_overlaps(int m,
        const ProfcnCoefficientBatchResult<T>& exterior,
        const ProfcnCoefficientBatchResult<T>& interior,
        OverlapQuadrature<T>& quadrature, bool even_only) {
    const int external_length = *std::max_element(exterior.lengths.begin(), exterior.lengths.end());
    const int internal_length = *std::max_element(interior.lengths.begin(), interior.lengths.end());
    quadrature.prepare(m, external_length, internal_length);
    OverlapSamples<T> samples;
    samples.count = quadrature.nodes().size();
    const size_t external_size = size_t(exterior.functions.lnum);
    const size_t internal_size = size_t(interior.functions.lnum);
    if (samples.count == 0 || external_size > samples.exterior.max_size() / samples.count ||
            internal_size > samples.weighted_interior.max_size() / samples.count)
        throw std::length_error("Fluid overlap quadrature samples exceed native allocation dimensions.");
    samples.exterior.resize(external_size * samples.count);
    samples.weighted_interior.resize(internal_size * samples.count);
    LegendreWorkspace<T> workspace;
    const size_t support = size_t(std::max(external_length, internal_length));
    // Evaluate the complete normalized angular expansions at the actual integration points.
    for (size_t node = 0; node < samples.count; ++node) {
        const auto basis = normalized_ferrers_values(m,
            psms_direction_from_cosine(quadrature.nodes()[node]), support, workspace);
        for (int mode = 0; mode < exterior.functions.lnum; ++mode) {
            if (even_only && mode % 2) continue;
            const T value = expansion_at_node(exterior, mode, basis);
            if (!is_finite_real(value)) throw std::runtime_error("Exterior quadrature angular value is nonfinite.");
            samples.exterior[size_t(mode) * samples.count + node] = value;
        }
        for (int mode = 0; mode < interior.functions.lnum; ++mode) {
            if (even_only && mode % 2) continue;
            const T value = expansion_at_node(interior, mode, basis) * quadrature.weights()[node];
            if (!is_finite_real(value)) throw std::runtime_error("Weighted interior quadrature angular value is nonfinite.");
            samples.weighted_interior[size_t(mode) * samples.count + node] = value;
        }
    }
    return samples;
}

template<typename T> std::array<std::complex<T>, 6> evaluate_traces(
        int m, T eta, const std::array<std::vector<std::complex<T>>, 6>& traces,
        const LegendreWorkspace<T>& workspace) {
    std::array<std::complex<T>, 6> sums{}, corrections{};
    if (workspace.support == 0) return sums;
    // The native Ferrers basis has positive P_m^m and unit integral norm.
    T diagonal = T(1) / precsqrt(T(2));
    const T transverse = precsqrt(std::max(T(0), (T(1) - eta) * (T(1) + eta)));
    const int scale_step = workspace.scale_step;
    const T upper_scale = workspace.upper_scale, lower_scale = workspace.lower_scale;
    long long basis_exponent = 0;
    for (int order = 1; order <= m; ++order) {
        diagonal *= transverse * workspace.diagonal_factors[order - 1];
        if (diagonal != T(0) && precabs(diagonal) < lower_scale) {
            diagonal *= upper_scale;
            basis_exponent -= scale_step;
        }
    }
    T previous = 0, current = diagonal;
    for (size_t index = 0; index < workspace.support; ++index) {
        const T value = basis_exponent == 0 ? current : psms_packed::scale_decimal(current, basis_exponent);
        for (size_t trace = 0; trace < traces.size(); ++trace)
            add(sums[trace], corrections[trace], value * traces[trace][index]);
        const T first = workspace.first_factors[index], second = workspace.second_factors[index];
        const T next = first * eta * current - second * previous;
        previous = current;
        current = next;
        // Preserve a recessive initial value until the degree recurrence reaches its support.
        const T magnitude = std::max(precabs(previous), precabs(current));
        if (magnitude > upper_scale) {
            previous *= lower_scale;
            current *= lower_scale;
            basis_exponent += scale_step;
        } else if (magnitude != T(0) && magnitude < lower_scale) {
            previous *= upper_scale;
            current *= upper_scale;
            basis_exponent -= scale_step;
        }
    }
    return sums;
}

template<typename T> std::array<T, 10> integrate_boundary(const OrderSolution<T>& solution,
        T chi, T xi, const SpheroidalDirection<T>& direction, Quadrature<T>& quadrature, size_t trace_prefix) {
    std::array<T, 10> sums{}, corrections{};
    const T incident_cosine = direction.cosine, incident_sine = direction.sine;
    const T focal_radius = precsqrt((xi - T(1)) * (xi + T(1)));
    const T transverse_phase = chi * focal_radius * incident_sine;
    const T transverse_derivative = chi * xi * incident_sine / focal_radius;
    const auto azimuth_phase = imaginary_unit_power<T>(solution.m) / T(2);
    quadrature.legendre.prepare(solution.m, trace_prefix);
    for (size_t index = 0; index < quadrature.nodes.size(); ++index) {
        const T eta = quadrature.nodes[index], weight = quadrature.weights[index];
        const T transverse = precsqrt(std::max(T(0), (T(1) - eta) * (T(1) + eta)));
        const T argument = transverse_phase * transverse;
        const T bessel = bessel_j(solution.m, argument);
        const T bessel_derivative = solution.m == 0 ? -bessel_j(1, argument) :
            (bessel_j(solution.m - 1, argument) - bessel_j(solution.m + 1, argument)) / T(2);
        const T phase = chi * xi * incident_cosine * eta;
        const auto phase_factor = azimuth_phase * std::complex<T>(preccos(phase), sine(phase));
        // The physical plane wave is twice the reduced cosine-series trace.
        const auto exact_pressure = phase_factor * bessel;
        const auto exact_derivative = phase_factor * std::complex<T>(
            transverse_derivative * transverse * bessel_derivative,
            chi * incident_cosine * eta * bessel);
        const auto traces = evaluate_traces(solution.m, eta, solution.traces, quadrature.legendre);
        const auto incident_pressure_error = exact_pressure - traces[4];
        const auto incident_derivative_error = exact_derivative - traces[5];
        const auto exterior_pressure = traces[0] + incident_pressure_error;
        const auto exterior_derivative = traces[1] + incident_derivative_error;
        const std::array<std::complex<T>, 10> values = {
            exterior_pressure - traces[2], exterior_derivative - traces[3],
            incident_pressure_error, incident_derivative_error, exterior_pressure,
            exterior_derivative, traces[2], traces[3], exact_pressure, exact_derivative};
        for (size_t component = 0; component < values.size(); ++component) {
            finite(values[component]);
            const T square = complex_abs_sq(values[component]);
            if (!is_finite_real(square))
                throw std::runtime_error("Boundary norm exceeds the requested native range.");
            add(sums[component], corrections[component], weight * square);
        }
    }
    for (T value : sums) finite_norm(value);
    return sums;
}

template<typename T> std::pair<T, T> boundary_budgets(const BoundaryCheck<T>& check,
        int m, T pressure_reference, T derivative_reference) {
    finite_norm(pressure_reference);
    finite_norm(derivative_reference);
    finite_norm(check.pressure_reference);
    finite_norm(check.derivative_reference);
    const T multiplicity = m == 0 ? T(1) : T(2);
    // Telescoping weights allocate numerical error across arbitrarily many orders.
    const T allocation = T(1) / precsqrt(T(2) * multiplicity * (T(m) + T(1)) * (T(m) + T(2)));
    const T pressure_scale = std::max(pressure_reference,
        precsqrt(multiplicity) * check.pressure_reference);
    const T derivative_scale = std::max(derivative_reference,
        precsqrt(multiplicity) * check.derivative_reference);
    return {finite_norm(allocation * numerical_budget(pressure_scale)),
        finite_norm(allocation * numerical_budget(derivative_scale))};
}

template<typename T> bool boundary_passes(const BoundaryCheck<T>& check,
        const std::pair<T, T>& budgets) {
    validate_boundary_errors(check);
    finite_norm(budgets.first);
    finite_norm(budgets.second);
    return check.pressure_error <= budgets.first && check.derivative_error <= budgets.second &&
        check.incident_pressure_error <= budgets.first && check.incident_derivative_error <= budgets.second;
}

template<typename T> struct TailNorm {
    T scale = 0, squares = 1;

    void append(const std::complex<T>& value) {
        finite(value);
        for (T component : {precabs(value.real()), precabs(value.imag())}) {
            if (component == T(0)) continue;
            if (scale < component) {
                const T ratio = scale / component;
                squares = T(1) + squares * ratio * ratio;
                scale = component;
            } else {
                const T ratio = component / scale;
                squares += ratio * ratio;
            }
        }
    }

    T value() const { return finite_norm(scale == T(0) ? T(0) : scale * precsqrt(squares)); }
};

template<typename T> struct TracePrefix {
    size_t length = 0;
    std::array<T, 6> tails{};
};

template<typename T> TracePrefix<T> select_trace_prefix(const OrderSolution<T>& solution,
        const std::pair<T, T>& budgets, bool constrain_pressure, bool constrain_derivative) {
    TracePrefix<T> result;
    result.length = solution.traces[0].size();
    // Unit-norm Ferrers orthogonality (DLMF 14.17.6) identifies coefficient and L2 tail norms.
    std::array<TailNorm<T>, 6> accumulated;
    // Reserve a quarter of each numerical budget for exact Parseval tail bounds.
    const T pressure_budget = budgets.first / T(4), derivative_budget = budgets.second / T(4);
    for (size_t index = result.length; index-- > 0;) {
        const auto exterior_pressure = solution.traces[0][index] - solution.traces[4][index];
        const auto exterior_derivative = solution.traces[1][index] - solution.traces[5][index];
        const std::array<std::complex<T>, 6> components = {
            exterior_pressure - solution.traces[2][index],
            exterior_derivative - solution.traces[3][index],
            solution.traces[4][index], solution.traces[5][index],
            exterior_pressure, exterior_derivative};
        std::array<T, 6> tails;
        for (size_t trace = 0; trace < components.size(); ++trace) {
            accumulated[trace].append(components[trace]);
            tails[trace] = accumulated[trace].value();
        }
        if ((constrain_pressure && (tails[0] > pressure_budget || tails[2] > pressure_budget)) ||
                (constrain_derivative && (tails[1] > derivative_budget || tails[3] > derivative_budget)))
            break;
        result.length = index;
        result.tails = tails;
    }
    return result;
}

template<typename T> void check_boundary_quadrature(OrderSolution<T>& solution, T chi, T xi,
        const SpheroidalDirection<T>& direction, T pressure_reference, T derivative_reference, Quadrature<T>& quadrature,
        bool constrain_pressure = true, bool constrain_derivative = true) {
    if (solution.traces[0].size() > size_t(INT_MAX - solution.m))
        throw std::length_error("Boundary reconstruction exceeds the native degree range.");
    const auto initial_budgets = boundary_budgets(BoundaryCheck<T>(), solution.m,
        pressure_reference, derivative_reference);
    const auto prefix = select_trace_prefix(solution, initial_budgets, constrain_pressure, constrain_derivative);
    const T interior_pressure_norm = norm(solution.traces[2]);
    const T interior_derivative_norm = norm(solution.traces[3]);
    // This rule exactly integrates every polynomial norm retained for evaluation.
    int count = std::max(1, solution.m + int(prefix.length));
    auto& check = solution.boundary;
    check.residual_method = "quadrature";
    check.trace_prefix = int(prefix.length);
    check.pressure_tail_bound = prefix.tails[0];
    check.derivative_tail_bound = prefix.tails[1];
    check.incident_pressure_tail_bound = prefix.tails[2];
    check.incident_derivative_tail_bound = prefix.tails[3];
    std::array<T, 10> previous{};
    bool have_previous = false;
    for (;;) {
        Rcpp::checkUserInterrupt();
        quadrature.prepare(count);
        const auto values = integrate_boundary(solution, chi, xi, direction, quadrature, prefix.length);
        check.pressure_error = precsqrt(values[0]) + prefix.tails[0];
        check.derivative_error = precsqrt(values[1]) + prefix.tails[1];
        check.incident_pressure_error = precsqrt(values[2]) + prefix.tails[2];
        check.incident_derivative_error = precsqrt(values[3]) + prefix.tails[3];
        validate_boundary_errors(check);
        // Full interior norms and lower exterior norm bounds prevent relaxed reference scales.
        check.pressure_reference = std::max({std::max(T(0), precsqrt(values[4]) - prefix.tails[4]),
            interior_pressure_norm, precsqrt(values[8])});
        check.derivative_reference = std::max({std::max(T(0), precsqrt(values[5]) - prefix.tails[5]),
            interior_derivative_norm, precsqrt(values[9])});
        check.quadrature_order = count;
        const auto budgets = boundary_budgets(check, solution.m, pressure_reference, derivative_reference);
        const bool pressure_ok = check.pressure_error <= budgets.first &&
            check.incident_pressure_error <= budgets.first;
        const bool derivative_ok = check.derivative_error <= budgets.second &&
            check.incident_derivative_error <= budgets.second;
        if ((constrain_pressure && !pressure_ok) || (constrain_derivative && !derivative_ok)) return;
        if (have_previous) {
            check.quadrature_pressure_change = precsqrt(std::max(
                precabs(values[0] - previous[0]), precabs(values[2] - previous[2])));
            check.quadrature_derivative_change = precsqrt(std::max(
                precabs(values[1] - previous[1]), precabs(values[3] - previous[3])));
            check.quadrature_converged = (!constrain_pressure || check.quadrature_pressure_change <= budgets.first) &&
                (!constrain_derivative || check.quadrature_derivative_change <= budgets.second);
            if (check.quadrature_converged) {
                // Add Parseval tails after accounting for the observed integration changes.
                const T pressure_change = check.quadrature_pressure_change * check.quadrature_pressure_change;
                const T derivative_change = check.quadrature_derivative_change * check.quadrature_derivative_change;
                check.pressure_error = precsqrt(values[0] + pressure_change) + prefix.tails[0];
                check.derivative_error = precsqrt(values[1] + derivative_change) + prefix.tails[1];
                check.incident_pressure_error = precsqrt(values[2] + pressure_change) + prefix.tails[2];
                check.incident_derivative_error = precsqrt(values[3] + derivative_change) + prefix.tails[3];
                validate_boundary_errors(check);
                check.quadrature_converged = (!constrain_pressure ||
                    (check.pressure_error <= budgets.first && check.incident_pressure_error <= budgets.first)) &&
                    (!constrain_derivative ||
                    (check.derivative_error <= budgets.second && check.incident_derivative_error <= budgets.second));
                if (check.quadrature_converged) return;
            }
        }
        previous = values;
        have_previous = true;
        if (count > INT_MAX / 2)
            throw std::length_error("Boundary quadrature refinement exceeds the native index range.");
        count *= 2;
    }
}

template<typename T> void check_boundary_quadrature(OrderSolution<T>& solution, T chi, T xi,
        T incident_cosine, T pressure_reference, T derivative_reference, Quadrature<T>& quadrature,
        bool constrain_pressure = true, bool constrain_derivative = true) {
    check_boundary_quadrature(solution, chi, xi, psms_direction_from_cosine(incident_cosine),
        pressure_reference, derivative_reference, quadrature, constrain_pressure, constrain_derivative);
}
}

namespace psms_automatic {
template<typename T> std::pair<std::vector<T>, std::vector<T>> incident_miller_sequence(
        T argument, int degree, int start) {
    if (!is_finite_real(argument) || !(argument > T(0)) || degree < 0 ||
            start <= degree || start > INT_MAX - 2)
        throw std::runtime_error("The incident spherical-Bessel recurrence has invalid native support.");
    std::vector<T> work(size_t(start) + 2, T(0));
    work[start] = T(1);
    int exponent = std::numeric_limits<T>::max_exponent10 / 4;
#if ACOUSTICTS_HAVE_QUADMATH
    if constexpr (std::is_same_v<T, acousticts_quad_t>) exponent = FLT128_MAX_10_EXP / 4;
#endif
    const T upper = pow10_typed<T>(exponent), lower = T(1) / upper;
    // Miller's backward recurrence selects the minimal solution (DLMF 3.6(iii)).
    for (int n = start; n >= 1; --n) {
        const T factor = (T(2) * n + T(1)) / argument;
        if (!is_finite_real(factor))
            throw std::runtime_error("The incident spherical-Bessel recurrence exceeds the native range.");
        while (precabs(work[n]) > upper / std::max(T(1), factor) || precabs(work[n + 1]) > upper) {
            for (int index = n; index <= start + 1; ++index) work[index] *= lower;
        }
        work[n - 1] = factor * work[n] - work[n + 1];
    }
    const T first = sine(argument) / argument;
    const T second = (first - preccos(argument)) / argument;
    const int anchor = precabs(first) >= precabs(second) ? 0 : 1;
    const T scale = (anchor == 0 ? first : second) / work[anchor];
    if (!is_finite_real(scale) || scale == T(0))
        throw std::runtime_error("The incident spherical-Bessel normalization is unresolved.");
    for (int index = 0; index <= degree + 1; ++index) work[index] *= scale;
    std::vector<T> values(degree + 1), derivatives(degree + 1);
    for (int index = 0; index <= degree; ++index) {
        values[index] = work[index];
        // The adjacent-order derivative avoids division by a small argument (DLMF 10.51.1).
        derivatives[index] = index == 0 ? -work[1] :
            (T(index) * work[index - 1] - T(index + 1) * work[index + 1]) / (T(2) * index + T(1));
        if (!is_finite_real(values[index]) || !is_finite_real(derivatives[index]))
            throw std::runtime_error("The incident spherical-Bessel sequence is nonfinite.");
    }
    return {std::move(values), std::move(derivatives)};
}

template<typename T> T incident_degree_tail(IncidentSequence<T>& sequence, int first_degree) {
    if (first_degree < 0 || !(sequence.argument > T(0)))
        throw std::runtime_error("The incident degree-tail bound has invalid support.");
    const T argument = sequence.argument;
    const T ratio = (argument / (T(2) * first_degree + T(1))) *
        (argument / (T(2) * first_degree + T(3)));
    if (!(ratio < T(1))) return T(std::numeric_limits<double>::infinity());
    const T log_argument = logarithm(argument);
    while (sequence.logarithmic_bounds.size() <= size_t(first_degree)) {
        const T degree = T(sequence.logarithmic_bounds.size());
        sequence.logarithmic_bounds.push_back(sequence.logarithmic_bounds.back() +
            log_argument - logarithm(T(2) * degree + T(1)));
    }
    // Rayleigh/addition and DLMF 10.14.4 bound the entire spherical degree tail.
    const T log_bound = sequence.logarithmic_bounds[first_degree] +
        (logarithm((T(2) * first_degree + T(1)) / T(2)) - logarithm(T(1) - ratio)) / T(2);
    T bound = exponential(log_bound);
    if (bound == T(0)) {
        // A representable positive upper allowance avoids silently dropping an underflowed tail.
#if ACOUSTICTS_HAVE_QUADMATH
        if constexpr (std::is_same_v<T, acousticts_quad_t>) bound = nextafterq(T(0), T(1));
        else
#endif
        bound = std::nextafter(T(0), T(1));
    }
    return bound;
}

template<typename T> void prepare_incident_sequence(IncidentSequence<T>& sequence,
        int degree, T radial_rate, T angular_rate, T pressure_budget, T derivative_budget,
        bool constrain_pressure, bool constrain_derivative) {
    finite_norm(pressure_budget);
    finite_norm(derivative_budget);
    auto derivative_change = [&]() {
        return finite_norm(precabs(radial_rate) * sequence.radial_change +
            precabs(angular_rate) * sequence.angular_change);
    };
    if (sequence.degree >= degree &&
            (!constrain_pressure || sequence.pressure_change <= pressure_budget) &&
            (!constrain_derivative || derivative_change() <= derivative_budget)) return;
    if (sequence.degree < degree) {
        if (degree > INT_MAX - 2)
            throw std::length_error("Incident coefficient support exceeds the native index range.");
        sequence.start = std::max(sequence.start, degree + 1);
        auto initial = incident_miller_sequence(sequence.argument, degree, sequence.start);
        sequence.values = std::move(initial.first);
        sequence.derivatives = std::move(initial.second);
        sequence.degree = degree;
    }
    T previous_error = T(std::numeric_limits<double>::infinity());
    for (;;) {
        Rcpp::checkUserInterrupt();
        sequence.start = grow_degree(sequence.start, 0);
        auto refined = incident_miller_sequence(sequence.argument, sequence.degree, sequence.start);
        TailNorm<T> pressure, radial, angular;
        for (int index = 0; index <= sequence.degree; ++index) {
            const T weight = precsqrt((T(2) * index + T(1)) / T(2));
            const T difference = refined.first[index] - sequence.values[index];
            pressure.append({weight * difference, T(0)});
            radial.append({weight * (refined.second[index] - sequence.derivatives[index]), T(0)});
            // Differentiating the addition theorem gives the polar derivative's squared norm.
            angular.append({weight * precsqrt(T(index) * (T(index) + T(1)) / T(2)) * difference, T(0)});
        }
        sequence.pressure_change = pressure.value();
        sequence.radial_change = radial.value();
        sequence.angular_change = angular.value();
        sequence.values = std::move(refined.first);
        sequence.derivatives = std::move(refined.second);
        const T derivative_error = derivative_change();
        const T error = std::max(constrain_pressure ? sequence.pressure_change / pressure_budget : T(0),
            constrain_derivative ? derivative_error / derivative_budget : T(0));
        finite_norm(error);
        if (error <= T(1)) return;
        if (error >= previous_error)
            throw std::runtime_error("Incident recurrence refinement did not resolve its allocated numerical allowance.");
        previous_error = error;
    }
}

template<typename T> struct IncidentCoefficients {
    int degree = -1;
    T pressure_tail = 0, derivative_tail = 0;
    T pressure_recurrence = 0, derivative_recurrence = 0;
    std::vector<std::complex<T>> pressure, derivative;
};

template<typename T> IncidentCoefficients<T> incident_coefficients(int m, T chi, T xi,
        const SpheroidalDirection<T>& direction, T pressure_budget, T derivative_budget, Quadrature<T>& workspace,
        bool constrain_pressure = true, bool constrain_derivative = true) {
    const T incident_cosine = direction.cosine, incident_sine = direction.sine;
    if (!(chi > T(0)) || !(xi > T(1)) || m < 0 || precabs(incident_cosine) > T(1))
        throw std::runtime_error("Analytic incident coefficients require positive finite spheroidal coordinates.");
    const T transverse = precsqrt((xi - T(1)) * (xi + T(1)));
    // Rayleigh geometry retains both components instead of recovering a small sine from its cosine.
    const T axial_component = xi * incident_cosine, transverse_component = transverse * incident_sine;
    const T distance_squared = axial_component * axial_component + transverse_component * transverse_component;
    const T distance = precsqrt(distance_squared), argument = chi * distance;
    const T cosine = direction.exact_axis ? incident_cosine : axial_component / distance;
    const T sine_angle = direction.exact_axis ? T(0) : transverse_component / distance;
    const T radial_rate = chi * xi * (incident_cosine * incident_cosine + incident_sine * incident_sine) / distance;
    const T angular_rate = incident_cosine * incident_sine / (transverse * distance_squared);
    const T derivative_magnitude = chi * precsqrt(incident_cosine * incident_cosine +
        (xi * incident_sine / transverse) * (xi * incident_sine / transverse));
    if (!is_finite_real(argument) || !(argument > T(0)) ||
            !is_finite_real(radial_rate) || !is_finite_real(angular_rate))
        throw std::runtime_error("Analytic incident geometry exceeds the native numerical range.");
    finite_norm(derivative_magnitude);
    auto& sequence = workspace.incident;
    sequence.reset(argument);
    const T multiplicity_root = m == 0 ? T(1) : precsqrt(T(2));
    const T pressure_allowance = pressure_budget / T(8);
    const T derivative_allowance = derivative_budget / T(8);
    if (!(pressure_allowance > T(0)) || !(derivative_allowance > T(0)))
        throw std::runtime_error("Analytic incident error allowances are below the native range.");
    IncidentCoefficients<T> result;
    result.degree = std::max(m, 1);
    for (;;) {
        Rcpp::checkUserInterrupt();
        if (result.degree > INT_MAX - 2)
            throw std::length_error("Incident degree-tail refinement exceeds the native index range.");
        result.pressure_tail = incident_degree_tail(sequence, result.degree + 1) / multiplicity_root;
        // Multiplication by the linear phase derivative raises spherical degree by at most one.
        result.derivative_tail = derivative_magnitude * incident_degree_tail(sequence, result.degree) / multiplicity_root;
        if (is_finite_real(result.pressure_tail) && is_finite_real(result.derivative_tail) &&
                (!constrain_pressure || result.pressure_tail <= pressure_allowance) &&
                (!constrain_derivative || result.derivative_tail <= derivative_allowance)) break;
        result.degree = grow_degree(result.degree, 0);
    }
    prepare_incident_sequence(sequence, result.degree, radial_rate, angular_rate,
        pressure_allowance * multiplicity_root, derivative_allowance * multiplicity_root,
        constrain_pressure, constrain_derivative);
    result.pressure_recurrence = sequence.pressure_change / multiplicity_root;
    result.derivative_recurrence = (precabs(radial_rate) * sequence.radial_change +
        precabs(angular_rate) * sequence.angular_change) / multiplicity_root;
    const size_t count = size_t(result.degree - m) + 1;
    result.pressure.resize(count);
    result.derivative.resize(count);
    workspace.legendre.prepare(m, count);
    T diagonal = T(1) / precsqrt(T(2)), diagonal_derivative = 0;
    long long exponent = 0;
    auto scale_pair = [&](T& first, T& second) {
        const T largest = std::max(precabs(first), precabs(second));
        if (largest != T(0) && largest < workspace.legendre.lower_scale) {
            first *= workspace.legendre.upper_scale;
            second *= workspace.legendre.upper_scale;
            exponent -= workspace.legendre.scale_step;
        }
    };
    // Differentiate in polar angle without dividing by an endpoint sine.
    for (int order = 1; order <= m; ++order) {
        const T factor = workspace.legendre.diagonal_factors[order - 1];
        const T next_derivative = factor * (cosine * diagonal + sine_angle * diagonal_derivative);
        diagonal *= factor * sine_angle;
        diagonal_derivative = next_derivative;
        scale_pair(diagonal, diagonal_derivative);
    }
    T previous = 0, current = diagonal, previous_derivative = 0, current_derivative = diagonal_derivative;
    for (size_t index = 0; index < count; ++index) {
        const int degree = m + int(index);
        const T value = psms_packed::scale_decimal(current, exponent);
        const T derivative = psms_packed::scale_decimal(current_derivative, exponent);
        const auto phase = imaginary_unit_power<T>(degree);
        // Rayleigh and spherical addition use the same unit-Ferrers normalization (DLMF 10.60.7, 14.30.9).
        result.pressure[index] = phase * (sequence.values[degree] * value);
        result.derivative[index] = phase * (radial_rate * sequence.derivatives[degree] * value +
            angular_rate * sequence.values[degree] * derivative);
        finite(result.pressure[index]);
        finite(result.derivative[index]);
        const T first = workspace.legendre.first_factors[index], second = workspace.legendre.second_factors[index];
        const T next = first * cosine * current - second * previous;
        const T next_derivative = first * (cosine * current_derivative - sine_angle * current) -
            second * previous_derivative;
        previous = current;
        current = next;
        previous_derivative = current_derivative;
        current_derivative = next_derivative;
        const T largest = std::max({precabs(previous), precabs(current),
            precabs(previous_derivative), precabs(current_derivative)});
        if (largest > workspace.legendre.upper_scale) {
            previous *= workspace.legendre.lower_scale;
            current *= workspace.legendre.lower_scale;
            previous_derivative *= workspace.legendre.lower_scale;
            current_derivative *= workspace.legendre.lower_scale;
            exponent += workspace.legendre.scale_step;
        } else if (largest != T(0) && largest < workspace.legendre.lower_scale) {
            previous *= workspace.legendre.upper_scale;
            current *= workspace.legendre.upper_scale;
            previous_derivative *= workspace.legendre.upper_scale;
            current_derivative *= workspace.legendre.upper_scale;
            exponent -= workspace.legendre.scale_step;
        }
    }
    return result;
}

template<typename T> IncidentCoefficients<T> incident_coefficients(int m, T chi, T xi,
        T incident_cosine, T pressure_budget, T derivative_budget, Quadrature<T>& workspace,
        bool constrain_pressure = true, bool constrain_derivative = true) {
    return incident_coefficients(m, chi, xi, psms_direction_from_cosine(incident_cosine),
        pressure_budget, derivative_budget, workspace, constrain_pressure, constrain_derivative);
}

template<typename T> void check_boundary_coefficients(OrderSolution<T>& solution, T chi, T xi,
        const SpheroidalDirection<T>& direction, T pressure_reference, T derivative_reference, Quadrature<T>& workspace,
        bool constrain_pressure = true, bool constrain_derivative = true) {
    const auto initial = boundary_budgets(BoundaryCheck<T>(), solution.m,
        pressure_reference, derivative_reference);
    const auto incident = incident_coefficients(solution.m, chi, xi, direction,
        initial.first, initial.second, workspace, constrain_pressure, constrain_derivative);
    const T pressure_allowance = finite_norm(incident.pressure_tail + incident.pressure_recurrence);
    const T derivative_allowance = finite_norm(incident.derivative_tail + incident.derivative_recurrence);
    size_t count = incident.pressure.size();
    for (const auto& trace : solution.traces) count = std::max(count, trace.size());
    std::array<TailNorm<T>, 10> norms;
    auto at = [](const std::vector<std::complex<T>>& values, size_t index) {
        return index < values.size() ? values[index] : std::complex<T>(0, 0);
    };
    // Parseval includes every exported solver coefficient and the complete incident-tail allowance.
    for (size_t index = 0; index < count; ++index) {
        const auto exact_pressure = at(incident.pressure, index), exact_derivative = at(incident.derivative, index);
        const auto incident_pressure_error = exact_pressure - at(solution.traces[4], index);
        const auto incident_derivative_error = exact_derivative - at(solution.traces[5], index);
        const auto exterior_pressure = at(solution.traces[0], index) + incident_pressure_error;
        const auto exterior_derivative = at(solution.traces[1], index) + incident_derivative_error;
        const auto interior_pressure = at(solution.traces[2], index), interior_derivative = at(solution.traces[3], index);
        const std::array<std::complex<T>, 10> values = {
            exterior_pressure - interior_pressure, exterior_derivative - interior_derivative,
            incident_pressure_error, incident_derivative_error, exterior_pressure, exterior_derivative,
            interior_pressure, interior_derivative, exact_pressure, exact_derivative};
        for (size_t component = 0; component < values.size(); ++component) norms[component].append(values[component]);
    }
    auto& check = solution.boundary;
    check = BoundaryCheck<T>();
    check.residual_method = "analytic_coefficients";
    check.incident_degree = incident.degree;
    check.incident_recurrence_start = workspace.incident.start;
    check.incident_pressure_tail_bound = incident.pressure_tail;
    check.incident_derivative_tail_bound = incident.derivative_tail;
    check.incident_pressure_recurrence_error = incident.pressure_recurrence;
    check.incident_derivative_recurrence_error = incident.derivative_recurrence;
    check.pressure_error = finite_norm(norms[0].value() + pressure_allowance);
    check.derivative_error = finite_norm(norms[1].value() + derivative_allowance);
    check.incident_pressure_error = finite_norm(norms[2].value() + pressure_allowance);
    check.incident_derivative_error = finite_norm(norms[3].value() + derivative_allowance);
    check.pressure_reference = std::max({std::max(T(0), norms[4].value() - pressure_allowance),
        norms[6].value(), std::max(T(0), norms[8].value() - pressure_allowance)});
    check.derivative_reference = std::max({std::max(T(0), norms[5].value() - derivative_allowance),
        norms[7].value(), std::max(T(0), norms[9].value() - derivative_allowance)});
    const auto budgets = boundary_budgets(check, solution.m, pressure_reference, derivative_reference);
    // This legacy readiness flag covers either method, while residual_method records the actual check.
    check.quadrature_converged = (!constrain_pressure ||
        (check.pressure_error <= budgets.first && check.incident_pressure_error <= budgets.first)) &&
        (!constrain_derivative ||
        (check.derivative_error <= budgets.second && check.incident_derivative_error <= budgets.second));
}

template<typename T> void check_boundary_coefficients(OrderSolution<T>& solution, T chi, T xi,
        T incident_cosine, T pressure_reference, T derivative_reference, Quadrature<T>& workspace,
        bool constrain_pressure = true, bool constrain_derivative = true) {
    check_boundary_coefficients(solution, chi, xi, psms_direction_from_cosine(incident_cosine),
        pressure_reference, derivative_reference, workspace, constrain_pressure, constrain_derivative);
}

template<typename T> void check_boundary(OrderSolution<T>& solution, T chi, T xi,
        const SpheroidalDirection<T>& direction, T pressure_reference, T derivative_reference, Quadrature<T>& workspace,
        bool constrain_pressure = true, bool constrain_derivative = true) {
    try {
        check_boundary_coefficients(solution, chi, xi, direction, pressure_reference,
            derivative_reference, workspace, constrain_pressure, constrain_derivative);
    } catch (const std::runtime_error& error) {
        // A failed numerical incident expansion uses the independent refined-quadrature checker.
        solution.boundary = BoundaryCheck<T>();
        solution.boundary.analytic_failure = error.what();
        check_boundary_quadrature(solution, chi, xi, direction, pressure_reference,
            derivative_reference, workspace, constrain_pressure, constrain_derivative);
    }
}

template<typename T> void check_boundary(OrderSolution<T>& solution, T chi, T xi,
        T incident_cosine, T pressure_reference, T derivative_reference, Quadrature<T>& workspace,
        bool constrain_pressure = true, bool constrain_derivative = true) {
    check_boundary(solution, chi, xi, psms_direction_from_cosine(incident_cosine),
        pressure_reference, derivative_reference, workspace, constrain_pressure, constrain_derivative);
}

template<typename T> size_t batch_bytes(const ProfcnCoefficientBatchResult<T>& batch) {
    const auto& functions = batch.functions;
    return (batch.coefficients.capacity() + functions.r1c.capacity() +
        functions.r1dc.capacity() + functions.r2c.capacity() + functions.r2dc.capacity() +
        functions.s1c.capacity() + functions.s1dc.capacity()) * sizeof(T) +
        (batch.lengths.capacity() + functions.ir1e.capacity() + functions.ir1de.capacity() +
        functions.ir2e.capacity() + functions.ir2de.capacity() + functions.naccr.capacity() +
        functions.is1e.capacity() + functions.is1de.capacity() + functions.naccs.capacity()) * sizeof(int);
}

template<typename T> OrderSolution<T> solve_order(int m, int degree, int interior_degree,
        T chi_body, T xi, T incident_cosine, T phi_difference, T density_ratio,
        bool backscatter, bool robin, const ProfcnCoefficientBatchResult<T>& exterior,
        OverlapQuadrature<T>* overlap_quadrature = nullptr,
        SpheroidalCoefficientOverlapCache<T>* overlap_cache = nullptr,
        bool export_traces = true, bool retain_outgoing_terms = false,
        bool checked_direct = false, bool matrix_free_robin = false,
        bool check_degree_support = false) {
    using Complex = std::complex<T>;
    OrderSolution<T> result;
    result.m = m;
    result.degree = degree;
    result.interior_degree = interior_degree;
    result.robin = robin;
    const int size = checked_profcn_integer(static_cast<long long>(degree) - m + 1,
        "inclusive exterior degree count");
    const auto density = psms_packed::normalize(density_ratio);
    // Each retry allocates one azimuthal order and preserves packed radial outputs.
    const int interior_size = checked_profcn_integer(static_cast<long long>(interior_degree) - m + 1,
        "inclusive interior degree count");
    const auto interior = cprofcn_mblock_coefficients<T>(chi_body, m, 1,
        interior_size, {}, 1, 0, xi - T(1), true);
    OverlapSamples<T> quadrature_samples;
    const bool use_quadrature = overlap_quadrature != nullptr && overlap_quadrature->active();
    const bool sample_quadrature = use_quadrature && !overlap_quadrature->adaptive;
    if (use_quadrature) {
        if (overlap_quadrature->adaptive) {
            const int external_length = *std::max_element(exterior.lengths.begin(), exterior.lengths.end());
            const int internal_length = *std::max_element(interior.lengths.begin(), interior.lengths.end());
            overlap_quadrature->prepare(m, external_length, internal_length);
            result.overlap_quadrature_order = overlap_quadrature->exact_order;
            result.overlap_method = "gauss_exactness_coefficient_identity";
        } else {
            quadrature_samples = sample_overlaps(m, exterior, interior, *overlap_quadrature,
                incident_cosine == T(0));
            result.overlap_quadrature_order = int(quadrature_samples.count);
            result.overlap_method = "prescribed_quadrature";
        }
        if (overlap_quadrature->adaptive && backscatter) result.outgoing_terms.resize(size);
    }
    if (retain_outgoing_terms) result.outgoing_terms.resize(size);
    const int support = std::max(exterior.capacity, interior.capacity);
    result.coefficients.resize(size);
    result.scaled_coefficients.resize(size);
    result.radial_exponents.resize(size);
    std::array<std::vector<Complex>, 6> trace_corrections;
    for (size_t trace = 0; export_traces && trace < result.traces.size(); ++trace) {
        result.traces[trace].resize(support);
        trace_corrections[trace].resize(support);
    }
    const size_t base_bytes = batch_bytes(exterior) + batch_bytes(interior) + quadrature_samples.bytes() +
        result.outgoing_terms.capacity() * sizeof(Complex) +
        (export_traces ? size_t(support) * 12 * sizeof(Complex) : 0) +
        size_t(size) * (2 * sizeof(Complex) + sizeof(long long));
    result.working_bytes = base_bytes + (overlap_cache ? overlap_cache->bytes() : 0);
    const T multiplicity = m == 0 ? T(1) : T(2);
    const T azimuth = backscatter ? (m % 2 ? T(-1) : T(1)) : preccos(T(m) * phi_difference);
    Complex field_correction(0, 0);
    std::array<Complex, 2> embedded_corrections{};
    const bool embedded = check_degree_support && !robin;
    for (int parity = 0; parity < 2; ++parity) {
        if (incident_cosine == T(0) && parity == 1) continue;
        const int dim = (size - parity + 1) / 2;
        const int internal_dim = (interior_size - parity + 1) / 2;
        if (dim <= 0 || internal_dim <= 0) continue;
        if (use_quadrature && !overlap_quadrature->adaptive &&
                size_t(dim) > 2 * overlap_quadrature->rank_bound(m, parity))
            throw std::runtime_error("Prescribed fluid overlap quadrature is insufficient: the retained parity block exceeds its exact quadrature rank bound. Increase n_integration or enable adaptive quadrature.");
        if (size_t(dim) > size_t(INT_MAX) / size_t(std::max(dim, internal_dim)))
            throw std::length_error("Full-fluid block exceeds the native matrix index range.");
        std::vector<T> incident_value(dim), incident_derivative(dim), receive_angular(dim);
        std::vector<T> interior_value(internal_dim), interior_derivative(internal_dim);
        std::vector<T> column_scales(dim), overlaps(size_t(internal_dim) * dim);
        std::vector<Complex> external_value(dim), external_derivative(dim), incident_coefficients(dim);
        std::vector<Complex> robin_pressure(robin || export_traces ? internal_dim : 0);
        std::vector<Complex> robin_derivative(export_traces ? internal_dim : 0);
        std::vector<long long> exponents(dim);
        for (int col = 0; col < dim; ++col) {
            const int index = parity + 2 * col;
            const T angular = extract_angular_value_from_mblock<T>(exterior.functions, 0, index, 0);
            receive_angular[col] = backscatter ? (parity ? -angular : angular) :
                extract_angular_value_from_mblock<T>(exterior.functions, 0, index, 1);
            if (!is_finite_real(angular) || !is_finite_real(receive_angular[col]))
                throw std::runtime_error("Full-fluid angular function is nonfinite.");
            incident_coefficients[col] = imaginary_unit_power<T>(m + index) * angular;
            const auto radial = psms_packed::exterior_pair(exterior.functions, 0, index, result.scaling);
            external_value[col] = radial.value;
            external_derivative[col] = radial.derivative;
            incident_value[col] = radial.incident_value;
            incident_derivative[col] = radial.incident_derivative;
            exponents[col] = radial.exponent;
            const int radial_accuracy = exterior.functions.naccr[index];
            const int angular_accuracy = exterior.functions.naccs[index];
            if (!is_na_int(radial_accuracy))
                result.minimum_radial_accuracy = std::min(result.minimum_radial_accuracy, radial_accuracy);
            if (angular != T(0) && !is_na_int(angular_accuracy))
                result.minimum_angular_accuracy = std::min(result.minimum_angular_accuracy, angular_accuracy);
        }
        const bool reuse_coefficients = overlap_cache != nullptr && !sample_quadrature;
        if (reuse_coefficients) {
            // Reuse only unchanged complete coefficient columns within this azimuthal order.
            overlap_cache->fill(m, parity, exterior, interior, overlaps);
            result.overlap_entries_reused += overlap_cache->reused_entries;
            result.overlap_entries_computed += overlap_cache->computed_entries;
            const size_t vectors_bytes = (size_t(4) * dim + size_t(2) * internal_dim +
                size_t(internal_dim) * dim) * sizeof(T) +
                (size_t(3) * dim + size_t(2) * internal_dim) * sizeof(Complex) +
                size_t(dim) * sizeof(long long);
            result.working_bytes = std::max(result.working_bytes,
                base_bytes + vectors_bytes + overlap_cache->peak_bytes);
        }
        for (int row = 0; row < internal_dim; ++row) {
            const auto radial = psms_packed::interior_pair(interior.functions, 0,
                parity + 2 * row, density, result.scaling);
            interior_value[row] = radial.first;
            interior_derivative[row] = radial.second;
            if (robin || export_traces) {
                const Complex denominator(radial.second, radial.first);
                robin_pressure[row] = radial.first / denominator;
                if (export_traces) robin_derivative[row] = radial.second / denominator;
            }
            if (reuse_coefficients) continue;
            for (int col = 0; col < dim; ++col) {
                // Preserve supplied samples for fixed rules and use Ferrers orthogonality otherwise.
                const T overlap = sample_quadrature ?
                    quadrature_samples.integral(parity + 2 * col, parity + 2 * row) :
                    spheroidal_coefficient_overlap(exterior, interior,
                        0, parity + 2 * col, parity + 2 * row);
                if (!is_finite_real(overlap)) throw std::runtime_error("Angular overlap is nonfinite.");
                overlaps[size_t(row) * dim + col] = overlap;
                ++result.overlap_entries_computed;
            }
        }

        std::vector<Complex> scattered;
        bool matrix_free_solved = false;
        if (robin && matrix_free_robin && !sample_quadrature) {
            // Outgoing Robin traces remove radial scales from the coupled operator.
            std::vector<Complex> outgoing_robin(dim), exterior_ratio(dim), forcing(dim);
            for (int col = 0; col < dim; ++col) {
                outgoing_robin[col] = external_derivative[col] + Complex(0, 1) * external_value[col];
                if (!psms_robin_solver::divide(external_value[col], outgoing_robin[col], exterior_ratio[col]))
                    throw std::runtime_error("Invalid outgoing Robin trace divisor.");
                forcing[col] = Complex(incident_derivative[col], incident_value[col]) * incident_coefficients[col];
            }
            auto rhs = psms_robin_solver::apply_map(overlaps, robin_pressure, forcing, internal_dim, dim);
            for (int col = 0; col < dim; ++col) rhs[col] -= incident_value[col] * incident_coefficients[col];
            // Bound retained arrays and restarted Arnoldi storage for workspace reporting.
            const size_t columns = size_t(dim), rows = size_t(internal_dim);
            const size_t restart = size_t(std::min(dim, psms_robin_solver::Options<T>().restart));
            const size_t arnoldi_complex = (2 * restart + 8) * columns +
                (restart + 1) * restart + 3 * restart + 1 + 2 * rows;
            const size_t retained_complex = 9 * columns + 2 * rows;
            const size_t retained_real = 4 * columns + 3 * rows + rows * columns + restart;
            result.working_bytes = std::max(result.working_bytes, base_bytes +
                (arnoldi_complex + retained_complex) * sizeof(Complex) + retained_real * sizeof(T) +
                columns * sizeof(long long) + (2 * restart + 1) * sizeof(std::vector<Complex>) +
                (overlap_cache ? overlap_cache->bytes() : 0));
            const auto solution = psms_robin_solver::solve(overlaps, robin_pressure, exterior_ratio,
                rhs, internal_dim, dim);
            if (solution.success) {
                scattered.resize(dim);
                for (int col = 0; col < dim; ++col) {
                    if (!psms_robin_solver::divide(solution.solution[col], outgoing_robin[col], scattered[col]))
                        throw std::runtime_error("Invalid recovered outgoing coefficient.");
                    column_scales[col] = T(1);
                }
                matrix_free_solved = true;
            } else {
                // Retry dense Robin only when the iterative algorithm cannot complete.
                const auto& status = solution.status;
                const bool dense_fallback = status == "diagonal_preconditioner_unavailable" ||
                    status == "projected_pivot_guard" || status == "iteration_budget_exhausted" ||
                    status == "native_iteration_stagnation" || status == "arnoldi_breakdown_without_convergence";
                if (!dense_fallback) throw std::runtime_error("Matrix-free Robin: " + status);
            }
        }
        if (!matrix_free_solved) {
            std::vector<Complex> matrix(size_t(dim) * dim), rhs(dim), boundary_map;
            std::array<std::vector<Complex>, 2> reduced_rhs;
            if (embedded) for (auto& values : reduced_rhs) values.resize(dim);
            if (robin) {
                // Compress the pole-free Robin map without dividing by either radial zero.
                boundary_map.resize(size_t(dim) * dim);
                std::vector<T> weighted_real(internal_dim), weighted_imaginary(internal_dim);
                T gram_squared = 0, gram_correction = 0;
                for (int row = 0; row < dim; ++row) {
                    for (int internal_mode = 0; internal_mode < internal_dim; ++internal_mode) {
                        const T overlap = overlaps[size_t(internal_mode) * dim + row];
                        weighted_real[internal_mode] = overlap * robin_pressure[internal_mode].real();
                        weighted_imaginary[internal_mode] = overlap * robin_pressure[internal_mode].imag();
                    }
                    for (int col = row; col < dim; ++col) {
                        Complex value(0, 0), correction(0, 0);
                        T gram = 0, correction_gram = 0;
#if ACOUSTICTS_HAVE_EXACT_QUAD_DOT
                        if constexpr (std::is_same_v<T, acousticts_quad_t>) {
                            static thread_local quad_dot_product::Accumulator accumulator;
                            value = Complex(quad_dot_product::dot(weighted_real.data(),
                                overlaps.data() + col, internal_dim, accumulator, 1, dim),
                                quad_dot_product::dot(weighted_imaginary.data(),
                                overlaps.data() + col, internal_dim, accumulator, 1, dim));
                            gram = quad_dot_product::dot(overlaps.data() + row,
                                overlaps.data() + col, internal_dim, accumulator, dim, dim);
                        } else
#endif
                        {
                            for (int internal_mode = 0; internal_mode < internal_dim; ++internal_mode) {
                                const T product = overlaps[size_t(internal_mode) * dim + row] *
                                    overlaps[size_t(internal_mode) * dim + col];
                                add(value, correction, product * robin_pressure[internal_mode]);
                                add(gram, correction_gram, product);
                            }
                        }
                        boundary_map[size_t(row) * dim + col] = value;
                        boundary_map[size_t(col) * dim + row] = value;
                        const T defect = (row == col ? T(1) : T(0)) - gram;
                        const T contribution = (row == col ? T(1) : T(2)) * defect * defect;
                        add(gram_squared, gram_correction, contribution);
                    }
                }
                result.gram_frobenius = std::max(result.gram_frobenius, precsqrt(gram_squared));
            }
            result.working_bytes = std::max(result.working_bytes, base_bytes +
                size_t(dim) * dim * sizeof(Complex) * (robin ? 3 : 2) +
                size_t(internal_dim) * dim * sizeof(T) +
                size_t(dim + internal_dim) * (12 * sizeof(Complex) + 8 * sizeof(T)) +
                (overlap_cache ? overlap_cache->bytes() : 0));
            for (int row = 0; row < dim; ++row) {
                Complex correction(0, 0);
                std::array<Complex, 2> reduced_corrections{};
                for (int col = 0; col < dim; ++col) {
                    Complex forcing;
                    if (robin) {
                        const Complex map = boundary_map[size_t(row) * dim + col];
                        const T identity = row == col ? T(1) : T(0);
                        matrix[size_t(row) * dim + col] = identity * external_value[col] -
                            map * (external_derivative[col] + Complex(0, 1) * external_value[col]);
                        forcing = identity * incident_value[col] -
                            map * Complex(incident_derivative[col], incident_value[col]);
                    } else {
                        const T overlap = overlaps[size_t(row) * dim + col];
                        matrix[size_t(row) * dim + col] = overlap *
                            (interior_derivative[row] * external_value[col] -
                            interior_value[row] * external_derivative[col]);
                        forcing = overlap * (interior_derivative[row] * incident_value[col] -
                            interior_value[row] * incident_derivative[col]);
                    }
                    add(rhs[row], correction, -forcing * incident_coefficients[col]);
                    // Omit incident forcing from the deleted degrees as well.
                    if (embedded) for (int reduction = 0; reduction < 2; ++reduction)
                        if (row < dim - reduction - 1 && col < dim - reduction - 1)
                            add(reduced_rhs[reduction][row], reduced_corrections[reduction],
                                -forcing * incident_coefficients[col]);
                }
            }
            // Keep complete pivoting, native guards and the existing residual refinement.
            for (int col = 0; col < dim; ++col) {
                T scale = 0;
                for (int row = 0; row < dim; ++row) {
                    const auto value = matrix[size_t(row) * dim + col];
                    finite(value);
                    scale = std::max(scale, precabs(value.real()) + precabs(value.imag()));
                }
                if (!(scale > T(0))) { result.singular = true; return result; }
                column_scales[col] = scale;
                for (int row = 0; row < dim; ++row) matrix[size_t(row) * dim + col] /= scale;
            }
            for (int row = 0; row < dim; ++row) {
                T scale = 0;
                for (int col = 0; col < dim; ++col) {
                    const auto value = matrix[size_t(row) * dim + col];
                    scale = std::max(scale, precabs(value.real()) + precabs(value.imag()));
                }
                if (!(scale > T(0))) { result.singular = true; return result; }
                for (int col = 0; col < dim; ++col) matrix[size_t(row) * dim + col] /= scale;
                rhs[row] /= scale;
                if (embedded) for (auto& values : reduced_rhs) values[row] /= scale;
            }
            const auto factors = lup_decompose<T>(matrix, dim);
            if (factors.singular) { result.singular = true; return result; }
            scattered = solve_factored_system_lup_refined<T>(matrix, rhs, factors, dim, 2);
            if (embedded) {
                const auto reduced = psms_embedded_support::reduced_solutions<T>(matrix, factors, dim, reduced_rhs);
                // Restore physical radial scaling before adding each reduced far field.
                for (int reduction = 0; reduction < 2; ++reduction) {
                    result.embedded_valid[reduction] = result.embedded_valid[reduction] && reduced.valid[reduction];
                    if (!reduced.valid[reduction]) continue;
                    for (int col = 0; col < dim - reduction - 1; ++col) {
                        const int index = parity + 2 * col;
                        const auto coefficient = reduced.solutions[reduction][col] / column_scales[col];
                        const auto term = psms_packed::scale_decimal(multiplicity * azimuth *
                            receive_angular[col] * imaginary_unit_power<T>(-m - index) * coefficient, -exponents[col]);
                        finite(term);
                        result.embedded_roundoff_scale += precabs(term.real()) + precabs(term.imag());
                        add_field(result.embedded_fields[reduction], embedded_corrections[reduction], term);
                    }
                }
            }
            if (checked_direct && !robin &&
                    !psms_direct_refinement::refine(matrix, rhs, scattered, factors, dim,
                        psms_modal_rel_tol<T>()).accepted) {
                result.singular = true;
                return result;
            }
        }
        std::vector<Complex> exterior_pressure(export_traces ? dim : 0);
        std::vector<Complex> exterior_derivative(export_traces ? dim : 0), exterior_robin(export_traces ? dim : 0);
        std::vector<Complex> internal_pressure(export_traces ? internal_dim : 0);
        std::vector<Complex> internal_derivative(export_traces ? internal_dim : 0);
        for (int col = 0; col < dim; ++col) {
            scattered[col] /= column_scales[col];
            finite(scattered[col]);
            if (export_traces) {
                exterior_pressure[col] = incident_value[col] * incident_coefficients[col] +
                    external_value[col] * scattered[col];
                exterior_derivative[col] = incident_derivative[col] * incident_coefficients[col] +
                    external_derivative[col] * scattered[col];
                exterior_robin[col] = exterior_derivative[col] + Complex(0, 1) * exterior_pressure[col];
            }
            const int index = parity + 2 * col;
            result.scaled_coefficients[index] = scattered[col];
            result.radial_exponents[index] = exponents[col];
            result.coefficients[index] = psms_packed::scale_decimal(scattered[col], -exponents[col]);
            const auto term = psms_packed::scale_decimal(multiplicity * azimuth *
                receive_angular[col] * imaginary_unit_power<T>(-m - index) * scattered[col], -exponents[col]);
            finite(term);
            add_field(result.field, field_correction, term);
            if (embedded) result.embedded_roundoff_scale += precabs(term.real()) + precabs(term.imag());
            if (!result.outgoing_terms.empty()) result.outgoing_terms[index] = term;
        }
        // Fixed-support evaluation needs the coupled field but no boundary trace reconstruction.
        if (!export_traces) continue;
        for (int row = 0; row < internal_dim; ++row) {
            Complex trace(0, 0), correction(0, 0);
            for (int col = 0; col < dim; ++col)
                add(trace, correction, overlaps[size_t(row) * dim + col] * exterior_robin[col]);
            internal_pressure[row] = robin_pressure[row] * trace;
            internal_derivative[row] = robin_derivative[row] * trace;
        }
        // Export complete Legendre traces for an independent exact-incident check.
        for (int legendre = parity; legendre < support; legendre += 2) {
            for (int col = 0; col < dim; ++col) {
                const T basis = coefficient(exterior, parity + 2 * col, legendre);
                add(result.traces[0][legendre], trace_corrections[0][legendre], basis * exterior_pressure[col]);
                add(result.traces[1][legendre], trace_corrections[1][legendre], basis * exterior_derivative[col]);
                add(result.traces[4][legendre], trace_corrections[4][legendre],
                    basis * incident_value[col] * incident_coefficients[col]);
                add(result.traces[5][legendre], trace_corrections[5][legendre],
                    basis * incident_derivative[col] * incident_coefficients[col]);
            }
            for (int row = 0; row < internal_dim; ++row) {
                const T basis = coefficient(interior, parity + 2 * row, legendre);
                add(result.traces[2][legendre], trace_corrections[2][legendre], basis * internal_pressure[row]);
                add(result.traces[3][legendre], trace_corrections[3][legendre], basis * internal_derivative[row]);
            }
        }
    }
    if (result.minimum_angular_accuracy == INT_MAX) result.minimum_angular_accuracy = NA_INTEGER;
    return result;
}
}

namespace psms_automatic {
template<typename T> T combine_norms(T first, T second) {
    finite_norm(first);
    finite_norm(second);
    const T scale = std::max(first, second);
    if (scale == T(0)) return T(0);
    first /= scale;
    second /= scale;
    return finite_norm(scale * precsqrt(first * first + second * second));
}

template<typename T> T combine_omitted_bound(T error, T omitted) {
    finite_norm(error);
    if (is_na_real(omitted) || omitted < T(0))
        throw std::runtime_error("The omitted spheroidal forcing bound is invalid.");
    // Infinity is intentional until the azimuthal factorial bound becomes useful.
    return is_finite_real(omitted) ? combine_norms(error, omitted) : omitted;
}

template<typename T> T coefficient_difference(const std::vector<std::complex<T>>& first,
        const std::vector<std::complex<T>>& second) {
    std::vector<std::complex<T>> difference(std::max(first.size(), second.size()));
    for (size_t index = 0; index < difference.size(); ++index)
        difference[index] = (index < first.size() ? first[index] : std::complex<T>(0, 0)) -
            (index < second.size() ? second[index] : std::complex<T>(0, 0));
    return norm(difference);
}

template<typename T> T coefficient_budget(int m, T total_norm, T current_norm) {
    finite_norm(total_norm);
    finite_norm(current_norm);
    const T multiplicity = m == 0 ? T(1) : T(2);
    const T scale = std::max(total_norm, precsqrt(multiplicity) * current_norm);
    return finite_norm(numerical_budget(scale) /
        precsqrt(T(2) * multiplicity * (T(m) + T(1)) * (T(m) + T(2))));
}

template<typename T> std::pair<T, T> omitted_incident_bound(int maximum_order,
        T chi, T xi, const SpheroidalDirection<T>& direction) {
    const T incident_cosine = direction.cosine, sine_incidence = direction.sine;
    if (direction.exact_axis) return {T(0), T(0)};
    const T radius = precsqrt((xi - T(1)) * (xi + T(1)));
    const T transverse = chi * radius * sine_incidence;
    const T ratio = transverse / (T(2) * (T(maximum_order) + T(1)));
    const T infinity = T(std::numeric_limits<double>::infinity());
    if (!(ratio < T(1))) return {infinity, infinity};
    // DLMF 10.14.4 and the derivative recurrence bound the entire omitted m tail.
    const T logarithmic_argument = logarithm(transverse / T(2));
    T logarithmic_factorial = 0;
    for (int order = 2; order <= maximum_order; ++order)
        logarithmic_factorial += logarithm(T(order));
    const T first = exponential(T(maximum_order) * logarithmic_argument - logarithmic_factorial);
    const T second = first * transverse / (T(2) * (T(maximum_order) + T(1)));
    const T third = second * transverse / (T(2) * (T(maximum_order) + T(2)));
    const T denominator = precsqrt(T(1) - ratio * ratio);
    const T axial_derivative = precabs(chi * incident_cosine);
    const T transverse_derivative = chi * xi * sine_incidence / radius;
    return {second / denominator,
        (axial_derivative * second + transverse_derivative * (first + third) / T(2)) / denominator};
}

template<typename T> std::pair<T, T> omitted_incident_bound(int maximum_order,
        T chi, T xi, T incident_cosine) {
    return omitted_incident_bound(maximum_order, chi, xi, psms_direction_from_cosine(incident_cosine));
}

template<typename T> size_t retained_bytes(const OrderSolution<T>& order) {
    size_t bytes = sizeof(OrderSolution<T>) +
        (order.coefficients.capacity() + order.scaled_coefficients.capacity() +
            order.outgoing_terms.capacity()) * sizeof(std::complex<T>) +
        order.radial_exponents.capacity() * sizeof(long long);
    for (const auto& trace : order.traces) bytes += trace.capacity() * sizeof(std::complex<T>);
    return bytes;
}

template<typename T> std::complex<T> adaptive_outgoing_band(const OrderSolution<T>& order) {
    finite(order.field);
    if (order.outgoing_terms.empty()) return order.field;
    std::complex<T> band(0, 0), correction(0, 0);
    T max_term = 0, previous_term = -1;
    int small_terms = 0;
    // Treat the existing stopping rule as a proposal after computing every retained degree.
    for (size_t index = 0; index < order.outgoing_terms.size(); ++index) {
        const auto& term = order.outgoing_terms[index];
        const T magnitude = combine_norms(precabs(term.real()), precabs(term.imag()));
        max_term = std::max(max_term, magnitude);
        if (psms_modal_term_is_negligible<T>(int(index), int(order.outgoing_terms.size()) - 1,
                magnitude, previous_term, max_term)) ++small_terms;
        else small_terms = 0;
        previous_term = magnitude;
        add_field(band, correction, term);
        if (small_terms >= 2) break;
    }
    // Small leading terms cannot justify omitting a later contribution to the complete band.
    const auto difference = band - order.field;
    const T change = combine_norms(precabs(difference.real()), precabs(difference.imag()));
    const T reference = combine_norms(precabs(order.field.real()), precabs(order.field.imag()));
    return change <= numerical_budget(reference) ? band : order.field;
}

template<typename T> std::complex<T> adaptive_outgoing_sum(const std::vector<OrderSolution<T>>& orders) {
    if (orders.empty()) return std::complex<T>(0, 0);
    std::complex<T> complete(0, 0), complete_correction(0, 0);
    for (const auto& order : orders) add_field(complete, complete_correction, order.field);
    std::complex<T> total(0, 0), correction(0, 0);
    T max_band = 0, previous_band = -1;
    int small_bands = 0;
    const int last_order = orders.back().m;
    // Candidate stopping never substitutes for evaluating uncomputed azimuthal orders.
    for (const auto& order : orders) {
        const auto band = adaptive_outgoing_band(order);
        add_field(total, correction, band);
        const T magnitude = combine_norms(precabs(band.real()), precabs(band.imag()));
        max_band = std::max(max_band, magnitude);
        if (psms_modal_band_is_negligible<T>(order.m, last_order, magnitude, previous_band, max_band))
            ++small_bands;
        else small_bands = 0;
        previous_band = magnitude;
        if (small_bands >= 2) break;
    }
    // Check accumulated omissions against the complete field, including cancellation across bands.
    const auto difference = total - complete;
    const T change = combine_norms(precabs(difference.real()), precabs(difference.imag()));
    const T reference = combine_norms(precabs(complete.real()), precabs(complete.imag()));
    return change <= numerical_budget(reference) ? total : complete;
}
}

template<typename T> psms_automatic::Result<T> psms_automatic_full_fluid_field(
        int m_seed, int n_seed, T chi_sw, T chi_body, T xi,
        const SpheroidalDirection<T>& incidence, const SpheroidalDirection<T>& reception,
        T phi_difference, T density_ratio, bool backscatter,
        const std::vector<T>* overlap_nodes = nullptr,
        const std::vector<T>* overlap_weights = nullptr, bool adaptive = false) {
    using namespace psms_automatic;
    using Complex = std::complex<T>;
    Result<T> result;
    try {
        OverlapQuadrature<T> overlap_quadrature;
        overlap_quadrature.supplied_nodes = overlap_nodes;
        overlap_quadrature.supplied_weights = overlap_weights;
        overlap_quadrature.adaptive = adaptive;
        overlap_quadrature.validate();
        const T eta_body = incidence.cosine, eta_scatter = reception.cosine;
        if (m_seed < 0 || n_seed < m_seed || !is_finite_real(chi_sw) ||
                !is_finite_real(chi_body) || !is_finite_real(xi) ||
                !is_finite_real(density_ratio) || !(chi_sw > T(0)) ||
                !(chi_body > T(0)) || !(xi > T(1)) || !(density_ratio > T(0)) ||
                !is_finite_real(eta_body) || precabs(eta_body) > T(1) ||
                !is_finite_real(eta_scatter) || precabs(eta_scatter) > T(1) ||
                !is_finite_real(phi_difference) || !is_finite_real(incidence.sine) ||
                !is_finite_real(reception.sine) || incidence.sine < T(0) || incidence.sine > T(1) ||
                reception.sine < T(0) || reception.sine > T(1))
            throw std::invalid_argument("Automatic full-fluid inputs are outside their native physical domain.");
        if (chi_sw == chi_body && density_ratio == T(1)) {
            // An exactly matched medium has no boundary contrast or scattered field.
            result.converged = true;
            result.matched_medium = true;
            result.field = Complex(0, 0);
            result.detail = "Exactly matched medium.";
            return result;
        }
        const bool axial = incidence.exact_axis;
        if (axial) m_seed = 0;
        const T sine_squared = incidence.sine * incidence.sine;
        const T exact_incident_pressure_norm = T(1) / precsqrt(T(2));
        const T exact_incident_derivative_norm = chi_sw * precsqrt((eta_body * eta_body +
            xi * xi * sine_squared / ((xi - T(1)) * (xi + T(1)))) / T(6));
        finite_norm(exact_incident_pressure_norm);
        finite_norm(exact_incident_derivative_norm);
        T pressure_reference_norm = exact_incident_pressure_norm;
        T derivative_reference_norm = exact_incident_derivative_norm;
        T accumulated_pressure_reference = 0, accumulated_derivative_reference = 0;
        T pressure_error_norm = 0, derivative_error_norm = 0, outgoing_norm = 0;
        T added_order_norm = 0;
        Complex field(0, 0), field_correction(0, 0);
        Quadrature<T> quadrature;
        int shell_end = m_seed;
        const int shell_size = psms_profcn_mblock_chunk_size<T>();
        for (int m = 0;; ++m) {
            Rcpp::checkUserInterrupt();
            result.last_order = m;
            int degree = std::max(m, n_seed), interior_degree = degree;
            int interior_bracket_lower = -1, interior_bracket_upper = -1;
            SpheroidalCoefficientOverlapCache<T> overlap_cache;
            bool robin = false, have_previous = false;
            std::vector<Complex> previous_coefficients;
            BoundaryCheck<T> previous_boundary;
            OrderSolution<T> accepted;
            int attempts = 0;
            for (;;) {
                result.last_degree = degree;
                result.last_interior_degree = interior_degree;
                std::vector<T> arguments = {eta_body};
                if (!backscatter) arguments.push_back(eta_scatter);
                auto exterior = cprofcn_mblock_coefficients<T>(chi_sw, m, 1,
                    checked_profcn_integer(static_cast<long long>(degree) - m + 1,
                        "inclusive exterior degree count"), arguments, 2, 0, xi - T(1), true);
                std::vector<SpheroidalDirection<T>> directions = {incidence};
                if (!backscatter) directions.push_back(reception);
                evaluate_coefficient_angular_values(exterior, m, directions, quadrature.legendre);
                std::vector<Complex> previous_internal;
                std::array<std::vector<Complex>, 4> previous_internal_traces;
                bool have_previous_internal = false, increase_exterior = false;
                int previous_internal_degree = -1, scheduled_interior_upper = -1;
                if (robin && interior_bracket_lower >= m) {
                    // Re-solve both successful supports for the new exterior basis before growing either.
                    const auto bracket = reuse_interior_bracket(degree, m,
                        interior_bracket_lower, interior_bracket_upper);
                    interior_degree = bracket.first;
                    scheduled_interior_upper = bracket.second;
                }
                OrderSolution<T> trial;
                for (;;) {
                    Rcpp::checkUserInterrupt();
                    result.last_interior_degree = interior_degree;
                    trial = solve_order(m, degree, interior_degree, chi_body, xi, eta_body,
                        phi_difference, density_ratio, backscatter, robin, exterior,
                        overlap_quadrature.active() ? &overlap_quadrature : nullptr,
                        overlap_quadrature.active() && !adaptive ? nullptr : &overlap_cache);
                    ++attempts;
                    size_t retained_trace_bytes = 0;
                    for (const auto& trace : previous_internal_traces)
                        retained_trace_bytes += trace.capacity() * sizeof(Complex);
                    result.peak_working_bytes = std::max(result.peak_working_bytes,
                        result.retained_bytes + trial.working_bytes + retained_trace_bytes +
                        (previous_coefficients.capacity() + previous_internal.capacity()) * sizeof(Complex) +
                        quadrature.bytes() + overlap_quadrature.bytes());
                    if (trial.singular) {
                        if (!robin) {
                            // Change to equivalent Robin equations instead of weakening the LU guard.
                            robin = true;
                            interior_degree = degree;
                            previous_internal.clear();
                            have_previous_internal = false;
                            previous_internal_degree = -1;
                            scheduled_interior_upper = -1;
                            continue;
                        }
                        increase_exterior = true;
                        break;
                    }
                    check_boundary(trial, chi_sw, xi, incidence, pressure_reference_norm,
                        derivative_reference_norm, quadrature);
                    result.last_boundary = trial.boundary;
                    if (!robin) break;
                    // Interior and exterior refinement are tested independently for the Robin system.
                    const T budget = coefficient_budget(m, outgoing_norm, norm(trial.coefficients));
                    const T change = have_previous_internal ?
                        coefficient_difference(previous_internal, trial.coefficients) : T(0);
                    const auto boundary_limits = boundary_budgets(trial.boundary, m,
                        pressure_reference_norm, derivative_reference_norm);
                    bool traces_stable = have_previous_internal;
                    for (size_t trace = 0; traces_stable && trace < previous_internal_traces.size(); ++trace) {
                        const T trace_change = coefficient_difference(previous_internal_traces[trace], trial.traces[trace]);
                        traces_stable = trace_change <= (trace % 2 ? boundary_limits.second : boundary_limits.first);
                    }
                    if (have_previous_internal && change <= budget && traces_stable) {
                        interior_bracket_lower = previous_internal_degree;
                        interior_bracket_upper = interior_degree;
                        trial.interior_comparison_lower = interior_bracket_lower;
                        trial.interior_comparison_upper = interior_bracket_upper;
                        break;
                    }
                    previous_internal = trial.coefficients;
                    for (size_t trace = 0; trace < previous_internal_traces.size(); ++trace)
                        previous_internal_traces[trace] = trial.traces[trace];
                    have_previous_internal = true;
                    previous_internal_degree = interior_degree;
                    interior_degree = scheduled_interior_upper > interior_degree ?
                        scheduled_interior_upper : grow_degree(interior_degree, m);
                    scheduled_interior_upper = -1;
                }
                if (!increase_exterior) {
                    const auto budgets = boundary_budgets(trial.boundary, m,
                        pressure_reference_norm, derivative_reference_norm);
                    const bool boundary_ok = trial.boundary.quadrature_converged &&
                        boundary_passes(trial.boundary, budgets);
                    trial.coefficient_change = have_previous ?
                        coefficient_difference(previous_coefficients, trial.coefficients) : T(0);
                    const T budget = coefficient_budget(m, outgoing_norm, norm(trial.coefficients));
                    if (boundary_ok && have_previous && trial.coefficient_change <= budget) {
                        trial.attempts = attempts;
                        accepted = std::move(trial);
                        break;
                    }
                    if (!boundary_ok && have_previous && trial.coefficient_change == T(0) &&
                            trial.boundary.pressure_error == previous_boundary.pressure_error &&
                            trial.boundary.derivative_error == previous_boundary.derivative_error &&
                            trial.boundary.incident_pressure_error == previous_boundary.incident_pressure_error &&
                            trial.boundary.incident_derivative_error == previous_boundary.incident_derivative_error)
                        throw std::runtime_error("Native coefficients and boundary residuals stagnated before the requested tolerance.");
                    previous_coefficients = trial.coefficients;
                    previous_boundary = trial.boundary;
                    have_previous = true;
                }
                degree = grow_degree(degree, m);
                interior_degree = robin ? std::max(interior_degree, degree) : degree;
            }
            const T multiplicity_root = m == 0 ? T(1) : precsqrt(T(2));
            pressure_error_norm = combine_norms(pressure_error_norm,
                multiplicity_root * accepted.boundary.pressure_error);
            derivative_error_norm = combine_norms(derivative_error_norm,
                multiplicity_root * accepted.boundary.derivative_error);
            accumulated_pressure_reference = combine_norms(accumulated_pressure_reference,
                multiplicity_root * accepted.boundary.pressure_reference);
            accumulated_derivative_reference = combine_norms(accumulated_derivative_reference,
                multiplicity_root * accepted.boundary.derivative_reference);
            pressure_reference_norm = std::max(exact_incident_pressure_norm, accumulated_pressure_reference);
            derivative_reference_norm = std::max(exact_incident_derivative_norm, accumulated_derivative_reference);
            const T order_norm = multiplicity_root * norm(accepted.coefficients);
            outgoing_norm = combine_norms(outgoing_norm, order_norm);
            if (m > m_seed) added_order_norm = combine_norms(added_order_norm, order_norm);
            add_field(field, field_correction, accepted.field);
            // Completed orders retain coefficients and diagnostics, never direction-stale LU factors.
            for (auto& trace : accepted.traces) std::vector<Complex>().swap(trace);
            result.retained_bytes += retained_bytes(accepted);
            result.orders.push_back(std::move(accepted));
            if (m < shell_end) continue;
            const auto omitted = omitted_incident_bound(m, chi_sw, xi, incidence);
            result.omitted_pressure_bound = omitted.first;
            result.omitted_derivative_bound = omitted.second;
            const T pressure_budget = numerical_budget(pressure_reference_norm);
            const T derivative_budget = numerical_budget(derivative_reference_norm);
            const bool azimuth_ok = omitted.first <= pressure_budget / precsqrt(T(2)) &&
                omitted.second <= derivative_budget / precsqrt(T(2));
            const bool outgoing_ok = axial || (m > m_seed && added_order_norm <=
                numerical_budget(outgoing_norm));
            result.added_order_norm = added_order_norm;
            result.pressure_relative = combine_omitted_bound(pressure_error_norm, omitted.first) / pressure_reference_norm;
            result.derivative_relative = combine_omitted_bound(derivative_error_norm, omitted.second) / derivative_reference_norm;
            if (azimuth_ok && outgoing_ok) {
                if (combine_norms(pressure_error_norm, omitted.first) > pressure_budget ||
                        combine_norms(derivative_error_norm, omitted.second) > derivative_budget)
                    throw std::runtime_error("The accumulated boundary error exceeds its allocated numerical tolerance.");
                finite(field);
                result.converged = true;
                result.field = field;
                result.detail = "Full incident boundary traces and nested modal coefficients satisfy the numerical tolerances. Each order records analytic incident-tail and recurrence allowances or native quadrature refinement for its boundary residual. The DLMF bound applies to omitted incident forcing. Observed outgoing-shell convergence does not exclude an isolated higher-order resonance.";
                if (overlap_quadrature.active() && adaptive && backscatter) {
                    const auto truncated = adaptive_outgoing_sum(result.orders);
                    const T change = norm(std::vector<Complex>{truncated - field});
                    // Retain the complete field if accumulated omitted terms exceed its numerical budget.
                    if (change <= numerical_budget(norm(std::vector<Complex>{field}))) {
                        result.field = truncated;
                        result.modal_tail_change = change;
                        result.modal_tail_truncated = truncated != field;
                    }
                }
                return result;
            }
            added_order_norm = T(0);
            if (axial) throw std::runtime_error("The analytic axial azimuthal check failed.");
            if (m > INT_MAX - shell_size)
                throw std::length_error("Azimuthal refinement exceeds the native index range.");
            shell_end = m + shell_size;
        }
    } catch (const std::exception& error) {
        result.detail = std::string(error.what()) + " [m=" + std::to_string(result.last_order) +
            ", exterior degree=" + std::to_string(result.last_degree) +
            ", interior degree=" + std::to_string(result.last_interior_degree) + "]";
        if (overlap_nodes != nullptr && !adaptive)
            result.detail += " [prescribed overlap quadrature points=" + std::to_string(overlap_nodes->size()) +
                "; quadrature was not enlarged or replaced by coefficient overlaps]";
        return result;
    }
}

template<typename T> psms_automatic::Result<T> psms_automatic_full_fluid_field(
        int m_seed, int n_seed, T chi_sw, T chi_body, T xi, T eta_body,
        T eta_scatter, T phi_difference, T density_ratio, bool backscatter,
        const std::vector<T>* overlap_nodes = nullptr,
        const std::vector<T>* overlap_weights = nullptr, bool adaptive = false) {
    return psms_automatic_full_fluid_field(m_seed, n_seed, chi_sw, chi_body, xi,
        psms_direction_from_cosine(eta_body), psms_direction_from_cosine(eta_scatter),
        phi_difference, density_ratio, backscatter, overlap_nodes, overlap_weights, adaptive);
}

namespace psms_automatic {
template<typename T> OrderSolution<T> solve_hard_order(int m, int degree,
        T incident_cosine, T phi_difference, bool backscatter, bool rigid,
        const ProfcnCoefficientBatchResult<T>& exterior, bool export_traces = true,
        bool retain_outgoing_terms = false) {
    using Complex = std::complex<T>;
    OrderSolution<T> result;
    result.m = m;
    result.degree = degree;
    result.interior_degree = -1;
    const int size = checked_profcn_integer(static_cast<long long>(degree) - m + 1,
        "inclusive exterior degree count"), support = exterior.capacity;
    result.coefficients.resize(size);
    result.scaled_coefficients.resize(size);
    result.radial_exponents.resize(size);
    if (retain_outgoing_terms) result.outgoing_terms.resize(size);
    std::array<std::vector<Complex>, 6> trace_corrections;
    for (size_t trace = 0; export_traces && trace < result.traces.size(); ++trace) {
        result.traces[trace].resize(support);
        trace_corrections[trace].resize(support);
    }
    result.working_bytes = batch_bytes(exterior) +
        (export_traces ? size_t(support) * 12 * sizeof(Complex) : 0) +
        result.outgoing_terms.capacity() * sizeof(Complex) +
        size_t(size) * (2 * sizeof(Complex) + sizeof(long long));
    const T multiplicity = m == 0 ? T(1) : T(2);
    const T azimuth = backscatter ? (m % 2 ? T(-1) : T(1)) : preccos(T(m) * phi_difference);
    Complex field_correction(0, 0);
    for (int index = 0; index < size; ++index) {
        const int parity = index % 2;
        if (incident_cosine == T(0) && parity == 1) continue;
        const T angular = extract_angular_value_from_mblock<T>(exterior.functions, 0, index, 0);
        const T receive_angular = backscatter ? (parity ? -angular : angular) :
            extract_angular_value_from_mblock<T>(exterior.functions, 0, index, 1);
        if (!is_finite_real(angular) || !is_finite_real(receive_angular))
            throw std::runtime_error("Hard-boundary angular function is nonfinite.");
        const Complex incident = imaginary_unit_power<T>(m + index) * angular;
        const auto radial = psms_packed::exterior_pair(exterior.functions, 0, index, result.scaling);
        const Complex divisor = rigid ? radial.derivative : radial.value;
        if (divisor == Complex(0, 0))
            throw std::runtime_error("The outgoing hard-boundary radial divisor vanished in native precision.");
        // Each separated boundary equation is diagonal and keeps its outgoing exponent packed.
        const T forcing = rigid ? radial.incident_derivative : radial.incident_value;
        const Complex scattered = -incident * (forcing / divisor);
        finite(scattered);
        Complex incident_pressure, incident_derivative, pressure, derivative;
        if (export_traces) {
            incident_pressure = radial.incident_value * incident;
            incident_derivative = radial.incident_derivative * incident;
            pressure = incident_pressure + radial.value * scattered;
            derivative = incident_derivative + radial.derivative * scattered;
        }
        result.scaled_coefficients[index] = scattered;
        result.radial_exponents[index] = radial.exponent;
        result.coefficients[index] = psms_packed::scale_decimal(scattered, -radial.exponent);
        const Complex term = psms_packed::scale_decimal(multiplicity * azimuth * receive_angular *
            imaginary_unit_power<T>(-m - index) * scattered, -radial.exponent);
        finite(term);
        add_field(result.field, field_correction, term);
        if (retain_outgoing_terms) result.outgoing_terms[index] = term;
        const int radial_accuracy = exterior.functions.naccr[index];
        const int angular_accuracy = exterior.functions.naccs[index];
        if (!is_na_int(radial_accuracy))
            result.minimum_radial_accuracy = std::min(result.minimum_radial_accuracy, radial_accuracy);
        if (angular != T(0) && !is_na_int(angular_accuracy))
            result.minimum_angular_accuracy = std::min(result.minimum_angular_accuracy, angular_accuracy);
        // Fixed-support evaluation omits traces used only by adaptive boundary checks.
        if (!export_traces) continue;
        // Retained traces are compared with the complete plane wave on the actual boundary.
        for (int legendre = parity; legendre < exterior.lengths[index]; legendre += 2) {
            const T basis = coefficient(exterior, index, legendre);
            add(result.traces[0][legendre], trace_corrections[0][legendre], basis * pressure);
            add(result.traces[1][legendre], trace_corrections[1][legendre], basis * derivative);
            add(result.traces[4][legendre], trace_corrections[4][legendre], basis * incident_pressure);
            add(result.traces[5][legendre], trace_corrections[5][legendre], basis * incident_derivative);
        }
    }
    if (result.minimum_angular_accuracy == INT_MAX) result.minimum_angular_accuracy = NA_INTEGER;
    return result;
}
}

template<typename T> psms_automatic::Result<T> psms_automatic_hard_boundary_field(
        int m_seed, int n_seed, T chi_sw, T xi,
        const SpheroidalDirection<T>& incidence, const SpheroidalDirection<T>& reception,
        T phi_difference, bool backscatter, bool rigid) {
    using namespace psms_automatic;
    using Complex = std::complex<T>;
    Result<T> result;
    // The unconstrained trace is an output field, not a boundary error.
    if (rigid) result.pressure_relative = precnan<T>();
    else result.derivative_relative = precnan<T>();
    try {
        const T eta_body = incidence.cosine, eta_scatter = reception.cosine;
        if (m_seed < 0 || n_seed < m_seed || !is_finite_real(chi_sw) ||
                !is_finite_real(xi) || !(chi_sw > T(0)) || !(xi > T(1)) ||
                !is_finite_real(eta_body) || precabs(eta_body) > T(1) ||
                !is_finite_real(eta_scatter) || precabs(eta_scatter) > T(1) ||
                !is_finite_real(phi_difference) || !is_finite_real(incidence.sine) ||
                !is_finite_real(reception.sine) || incidence.sine < T(0) || incidence.sine > T(1) ||
                reception.sine < T(0) || reception.sine > T(1))
            throw std::invalid_argument("Automatic hard-boundary inputs are outside their native physical domain.");
        const bool axial = incidence.exact_axis;
        if (axial) m_seed = 0;
        const T sine_squared = incidence.sine * incidence.sine;
        const T pressure_reference = T(1) / precsqrt(T(2));
        const T derivative_reference = chi_sw * precsqrt((eta_body * eta_body +
            xi * xi * sine_squared / ((xi - T(1)) * (xi + T(1)))) / T(6));
        finite_norm(pressure_reference);
        finite_norm(derivative_reference);
        const T incident_reference = rigid ? derivative_reference : pressure_reference;
        const T boundary_budget = numerical_budget(incident_reference);
        T boundary_error_norm = 0, outgoing_norm = 0, added_order_norm = 0;
        Complex field(0, 0), field_correction(0, 0);
        Quadrature<T> quadrature;
        int shell_end = m_seed;
        const int shell_size = psms_profcn_mblock_chunk_size<T>();
        for (int m = 0;; ++m) {
            Rcpp::checkUserInterrupt();
            result.last_order = m;
            int degree = std::max(m, n_seed), attempts = 0;
            bool have_previous = false;
            std::vector<Complex> previous_coefficients;
            T previous_error = 0, previous_incident_error = 0;
            OrderSolution<T> accepted;
            for (;;) {
                result.last_degree = degree;
                std::vector<T> arguments = {eta_body};
                if (!backscatter) arguments.push_back(eta_scatter);
                auto exterior = cprofcn_mblock_coefficients<T>(chi_sw, m, 1,
                    checked_profcn_integer(static_cast<long long>(degree) - m + 1,
                        "inclusive exterior degree count"), arguments, 2, 0, xi - T(1), true);
                std::vector<SpheroidalDirection<T>> directions = {incidence};
                if (!backscatter) directions.push_back(reception);
                evaluate_coefficient_angular_values(exterior, m, directions, quadrature.legendre);
                auto trial = solve_hard_order(m, degree, eta_body, phi_difference,
                    backscatter, rigid, exterior);
                ++attempts;
                // Only the prescribed Dirichlet or Neumann trace controls boundary acceptance.
                check_boundary(trial, chi_sw, xi, incidence, pressure_reference,
                    derivative_reference, quadrature, !rigid, rigid);
                result.last_boundary = trial.boundary;
                result.peak_working_bytes = std::max(result.peak_working_bytes,
                    result.retained_bytes + trial.working_bytes +
                    previous_coefficients.capacity() * sizeof(Complex) + quadrature.bytes());
                const auto limits = boundary_budgets(trial.boundary, m, pressure_reference, derivative_reference);
                const T error = rigid ? trial.boundary.derivative_error : trial.boundary.pressure_error;
                const T incident_error = rigid ? trial.boundary.incident_derivative_error : trial.boundary.incident_pressure_error;
                const T limit = rigid ? limits.second : limits.first;
                const bool boundary_ok = trial.boundary.quadrature_converged &&
                    error <= limit && incident_error <= limit;
                trial.coefficient_change = have_previous ?
                    coefficient_difference(previous_coefficients, trial.coefficients) : T(0);
                const T budget = coefficient_budget(m, outgoing_norm, norm(trial.coefficients));
                if (boundary_ok && have_previous && trial.coefficient_change <= budget) {
                    trial.attempts = attempts;
                    accepted = std::move(trial);
                    break;
                }
                if (!boundary_ok && have_previous && trial.coefficient_change == T(0) &&
                        error == previous_error && incident_error == previous_incident_error)
                    throw std::runtime_error("Native hard-boundary coefficients and residuals stagnated before the requested tolerance.");
                previous_coefficients = trial.coefficients;
                previous_error = error;
                previous_incident_error = incident_error;
                have_previous = true;
                degree = grow_degree(degree, m);
            }
            const T multiplicity_root = m == 0 ? T(1) : precsqrt(T(2));
            const T error = rigid ? accepted.boundary.derivative_error : accepted.boundary.pressure_error;
            boundary_error_norm = combine_norms(boundary_error_norm, multiplicity_root * error);
            const T order_norm = multiplicity_root * norm(accepted.coefficients);
            outgoing_norm = combine_norms(outgoing_norm, order_norm);
            if (m > m_seed) added_order_norm = combine_norms(added_order_norm, order_norm);
            add_field(field, field_correction, accepted.field);
            for (auto& trace : accepted.traces) std::vector<Complex>().swap(trace);
            result.retained_bytes += retained_bytes(accepted);
            result.orders.push_back(std::move(accepted));
            if (m < shell_end) continue;
            const auto omitted = omitted_incident_bound(m, chi_sw, xi, incidence);
            result.omitted_pressure_bound = omitted.first;
            result.omitted_derivative_bound = omitted.second;
            const T omitted_forcing = rigid ? omitted.second : omitted.first;
            const bool azimuth_ok = omitted_forcing <= boundary_budget / precsqrt(T(2));
            const bool outgoing_ok = axial || (m > m_seed && added_order_norm <=
                numerical_budget(outgoing_norm));
            result.added_order_norm = added_order_norm;
            const T relative = combine_omitted_bound(boundary_error_norm, omitted_forcing) / incident_reference;
            if (rigid) result.derivative_relative = relative;
            else result.pressure_relative = relative;
            if (azimuth_ok && outgoing_ok) {
                if (combine_norms(boundary_error_norm, omitted_forcing) > boundary_budget)
                    throw std::runtime_error("The accumulated hard-boundary error exceeds its allocated numerical tolerance.");
                finite(field);
                result.converged = true;
                result.field = field;
                result.detail = std::string(rigid ? "The complete zero-normal-derivative boundary condition" :
                    "The complete zero-pressure boundary condition") +
                    " and nested outgoing coefficients satisfy the numerical tolerances. Each order records analytic incident-tail and recurrence allowances or native quadrature refinement for its boundary residual. The DLMF bound covers omitted incident forcing. Outgoing-shell convergence is an observed numerical check.";
                return result;
            }
            added_order_norm = T(0);
            if (axial) throw std::runtime_error("The analytic axial hard-boundary azimuthal check failed.");
            if (m > INT_MAX - shell_size)
                throw std::length_error("Hard-boundary azimuthal refinement exceeds the native index range.");
            shell_end = m + shell_size;
        }
    } catch (const std::exception& error) {
        result.detail = std::string(error.what()) + " [m=" + std::to_string(result.last_order) +
            ", exterior degree=" + std::to_string(result.last_degree) + "]";
        return result;
    }
}

template<typename T> psms_automatic::Result<T> psms_automatic_hard_boundary_field(
        int m_seed, int n_seed, T chi_sw, T xi, T eta_body, T eta_scatter,
        T phi_difference, bool backscatter, bool rigid) {
    return psms_automatic_hard_boundary_field(m_seed, n_seed, chi_sw, xi,
        psms_direction_from_cosine(eta_body), psms_direction_from_cosine(eta_scatter),
        phi_difference, backscatter, rigid);
}


#if ACOUSTICTS_HAVE_QUADMATH
inline std::complex<acousticts_quad_t> psms_adaptive_liquid_field(
        int m_seed, int n_seed, acousticts_quad_t chi_sw, acousticts_quad_t chi_body,
        acousticts_quad_t xi, const SpheroidalDirection<acousticts_quad_t>& incidence,
        const SpheroidalDirection<acousticts_quad_t>& reception,
        acousticts_quad_t phi_difference, acousticts_quad_t density_ratio) {
    using T = acousticts_quad_t;
    using Complex = std::complex<T>;
    using namespace psms_automatic;
    // Preserve the native physical-domain checks and the exact matched-medium solution.
    if (m_seed < 0 || n_seed < m_seed || !is_finite_real(chi_sw) || !(chi_sw > T(0)) ||
            !is_finite_real(xi) || !(xi > T(1)) || !is_finite_real(phi_difference) ||
            !is_finite_real(incidence.cosine) || precabs(incidence.cosine) > T(1) ||
            !is_finite_real(reception.cosine) || precabs(reception.cosine) > T(1) ||
            !is_finite_real(incidence.sine) || incidence.sine < T(0) || incidence.sine > T(1) ||
            !is_finite_real(reception.sine) || reception.sine < T(0) || reception.sine > T(1) ||
            !is_finite_real(chi_body) || !(chi_body > T(0)) ||
            !is_finite_real(density_ratio) || !(density_ratio > T(0)))
        throw std::invalid_argument("Adaptive liquid spheroidal inputs are outside their native physical domain.");
    checked_profcn_integer(static_cast<long long>(n_seed) + 1, "inclusive maximum degree count");
    if (chi_sw == chi_body && density_ratio == T(1)) return Complex(0, 0);
    const std::vector<T> arguments{incidence.cosine};
    const std::vector<SpheroidalDirection<T>> directions{incidence};
    LegendreWorkspace<T> workspace;
    auto field_norm = [](const Complex& value) {
        return combine_norms(precabs(value.real()), precabs(value.imag()));
    };
    // Require both the TS target and the native complex-amplitude tolerance.
    const T target_strength_tolerance_db = T(1e-8);
    const T ts_amplitude_allowance = -expm1q(-logq(T(10)) * target_strength_tolerance_db / T(20));
    const T amplitude_allowance = std::min(ts_amplitude_allowance, psms_modal_rel_tol<T>());
    Complex field(0, 0), field_correction(0, 0);
    T previous_order = -1, penultimate_order = -1, total_degree_estimate = 0;
    int small_orders = 0;
    // Refine only the unresolved order while retaining completed contributions.
    for (int m = 0; ; ++m) {
        Rcpp::checkUserInterrupt();
        // Removing two coefficients per parity requires at least six retained degrees.
        int degree = std::max(n_seed, checked_profcn_integer(
            static_cast<long long>(m) + 5, "embedded degree support"));
        SpheroidalCoefficientOverlapCache<T> overlap_cache;
        OrderSolution<T> order;
        T estimate = 0, allowance = 0;
        bool degree_pass = false;
        for (int refinement = 0; refinement < 12; ++refinement) {
            auto exterior = cprofcn_mblock_coefficients<T>(chi_sw, m, 1,
                degree - m + 1, arguments, 2, 0, xi - T(1), true);
            evaluate_coefficient_angular_values(exterior, m, directions, workspace);
            // Complete coefficient overlaps provide automatic finite-expansion integration.
            order = solve_order<T>(m, degree, degree, chi_body, xi, incidence.cosine,
                phi_difference, density_ratio, true, false, exterior, nullptr,
                &overlap_cache, false, false, false, false, true);
            if (order.singular)
                throw std::runtime_error("Adaptive liquid coupled system is singular.");
            if (!order.embedded_valid[0] || !order.embedded_valid[1])
                throw std::runtime_error("Adaptive liquid degree comparison failed its arithmetic check.");
            const T scale = field_norm(field + order.field);
            allowance = amplitude_allowance * scale / (T(2) * T(m + 1) * T(m + 2));
            const T last = field_norm(order.field - order.embedded_fields[0]);
            const T preceding = field_norm(order.embedded_fields[0] - order.embedded_fields[1]);
            // Use modal magnitudes to account for cancellation in the rounding scale.
            const T arithmetic_floor = T(64) * preceps<T>() * order.embedded_roundoff_scale;
            degree_pass = false;
            // A geometric remainder estimates continued contraction of the degree increments.
            if (last <= arithmetic_floor && preceding <= arithmetic_floor) {
                estimate = arithmetic_floor;
                degree_pass = estimate <= allowance;
            } else if (preceding > T(0) && last < preceding) {
                const T ratio = last / preceding;
                estimate = last / (T(1) - ratio);
                degree_pass = estimate <= allowance && preceding <= allowance;
            }
            if (degree_pass) break;
            degree = checked_profcn_integer((5LL * (degree + 1LL) + 3LL) / 4LL - 1LL,
                "adaptive liquid degree refinement");
        }
        if (!degree_pass)
            throw std::runtime_error("Adaptive liquid degree tail did not converge within twelve refinements.");
        total_degree_estimate += estimate;
        add_field(field, field_correction, order.field);
        const T magnitude = field_norm(order.field);
        const T order_allowance = amplitude_allowance * field_norm(field) / T(2);
        small_orders = m >= 2 && magnitude <= order_allowance ? small_orders + 1 : 0;
        bool order_pass = incidence.exact_axis;
        if (small_orders >= 2 && previous_order > T(0) && penultimate_order > T(0)) {
            const T ratio = std::max(magnitude / previous_order, previous_order / penultimate_order);
            order_pass = ratio < T(1) && magnitude * ratio / (T(1) - ratio) <= order_allowance;
        }
        if (order_pass) {
            // Recheck the total degree allowance against the final accumulated field.
            if (total_degree_estimate > amplitude_allowance * field_norm(field) / T(2))
                throw std::runtime_error("Cancellation invalidated the adaptive liquid degree allowance.");
            finite(field);
            return field;
        }
        penultimate_order = previous_order;
        previous_order = magnitude;
    }
}
#endif

template<typename T> std::complex<T> psms_fixed_support_field(
        int m_max, int n_max, T chi_sw, T chi_body, T xi,
        const SpheroidalDirection<T>& incidence, const SpheroidalDirection<T>& reception,
        T phi_difference, T density_ratio, bool backscatter, bool hard, bool rigid, bool adaptive,
        const std::vector<T>* overlap_nodes = nullptr,
        const std::vector<T>* overlap_weights = nullptr, bool gas_recovery = false) {
    using namespace psms_automatic;
    using Complex = std::complex<T>;
    // Validate the prescribed support and preserve the supplied liquid quadrature rule.
    if (m_max < 0 || n_max < m_max || !is_finite_real(chi_sw) || !(chi_sw > T(0)) ||
            !is_finite_real(xi) || !(xi > T(1)) || !is_finite_real(phi_difference) ||
            !is_finite_real(incidence.cosine) || precabs(incidence.cosine) > T(1) ||
            !is_finite_real(reception.cosine) || precabs(reception.cosine) > T(1) ||
            !is_finite_real(incidence.sine) || incidence.sine < T(0) || incidence.sine > T(1) ||
            !is_finite_real(reception.sine) || reception.sine < T(0) || reception.sine > T(1) ||
            (!hard && (!is_finite_real(chi_body) || !(chi_body > T(0)) ||
                !is_finite_real(density_ratio) || !(density_ratio > T(0)))))
        throw std::invalid_argument("Fixed-support spheroidal inputs are outside their native physical domain.");
    checked_profcn_integer(static_cast<long long>(n_max) + 1, "inclusive maximum degree count");
    OverlapQuadrature<T> overlap_quadrature;
    overlap_quadrature.supplied_nodes = overlap_nodes;
    overlap_quadrature.supplied_weights = overlap_weights;
    overlap_quadrature.adaptive = adaptive;
    overlap_quadrature.validate();
    if (!hard && chi_sw == chi_body && density_ratio == T(1)) return Complex(0, 0);

    // Exact axial incidence has only m=0, while all other directions retain the requested orders.
    if (incidence.exact_axis) m_max = 0;
    std::vector<T> arguments = {incidence.cosine};
    std::vector<SpheroidalDirection<T>> directions = {incidence};
    if (!backscatter) {
        arguments.push_back(reception.cosine);
        directions.push_back(reception);
    }
    LegendreWorkspace<T> workspace;
    Quadrature<T> recovery_workspace;
    Complex field(0, 0), correction(0, 0);
    Complex truncated(0, 0), truncated_correction(0, 0);
    const bool truncate_tails = adaptive && backscatter;
    bool candidate_complete = false;
    T max_band = 0, previous_band = -1;
    int small_bands = 0;
    const bool recovery_enabled = gas_recovery && is_quad_precision_v<T> &&
        !hard && !overlap_quadrature.active();
    bool recovery_required = false;
    int previous_recovery_degree = n_max;
    for (int stage = 0; stage < (recovery_enabled ? 2 : 1); ++stage) {
        const bool recovering = stage == 1;
        if (recovering) {
            // Discard the entire direct prefix and its adaptive stopping state.
            field = correction = truncated = truncated_correction = Complex(0, 0);
            candidate_complete = false;
            max_band = 0;
            previous_band = -1;
            small_bands = 0;
        }
        for (int m = 0; m <= m_max; ++m) {
            Rcpp::checkUserInterrupt();
            int degree = recovering ? previous_recovery_degree : n_max;
            SpheroidalCoefficientOverlapCache<T> overlap_cache;
            // Share the current normalized angular expansion and packed radial boundary values.
            auto exterior = cprofcn_mblock_coefficients<T>(chi_sw, m, 1,
                degree - m + 1, arguments, 2, 0, xi - T(1), true);
            evaluate_coefficient_angular_values(exterior, m, directions, workspace);
            auto order = hard ?
                solve_hard_order(m, n_max, incidence.cosine, phi_difference,
                    backscatter, rigid, exterior, false, truncate_tails) :
                solve_order<T>(m, degree, degree, chi_body, xi, incidence.cosine,
                    phi_difference, density_ratio, backscatter, recovering, exterior,
                    overlap_quadrature.active() ? &overlap_quadrature : nullptr,
                    recovering ? &overlap_cache : nullptr, recovering, truncate_tails,
                    recovery_enabled && !recovering, recovering);

            // A rejected direct block requires physical recovery of every retained order.
            if (recovery_enabled && !recovering && order.singular) {
                recovery_required = true;
                break;
            }
            if (recovering) {
                const T pressure_reference = T(1) / precsqrt(T(2));
                const T derivative_reference = chi_sw * precsqrt((incidence.cosine * incidence.cosine +
                    xi * xi * incidence.sine * incidence.sine / ((xi - T(1)) * (xi + T(1)))) / T(6));
                std::vector<Complex> previous_coefficients;
                BoundaryCheck<T> previous_boundary;
                bool have_previous = false;
                for (;;) {
                    Rcpp::checkUserInterrupt();
                    if (!order.singular) {
                        check_boundary(order, chi_sw, xi, incidence, pressure_reference,
                            derivative_reference, recovery_workspace);
                        const auto budgets = boundary_budgets(order.boundary, m,
                            pressure_reference, derivative_reference);
                        if (order.boundary.quadrature_converged && boundary_passes(order.boundary, budgets)) break;
                        // Stop only when native coefficients and all four residual errors are unchanged.
                        if (have_previous && coefficient_difference(previous_coefficients, order.coefficients) == T(0) &&
                                order.boundary.pressure_error == previous_boundary.pressure_error &&
                                order.boundary.derivative_error == previous_boundary.derivative_error &&
                                order.boundary.incident_pressure_error == previous_boundary.incident_pressure_error &&
                                order.boundary.incident_derivative_error == previous_boundary.incident_derivative_error)
                            throw std::runtime_error("Native coefficients and boundary residuals stagnated before the requested tolerance.");
                        previous_coefficients = order.coefficients;
                        previous_boundary = order.boundary;
                        have_previous = true;
                    }
                    // Enlarge only the unresolved order and reuse unchanged complete overlap columns.
                    degree = grow_degree(degree, m);
                    exterior = cprofcn_mblock_coefficients<T>(chi_sw, m, 1,
                        degree - m + 1, arguments, 2, 0, xi - T(1), true);
                    evaluate_coefficient_angular_values(exterior, m, directions, workspace);
                    order = solve_order<T>(m, degree, degree, chi_body, xi, incidence.cosine,
                        phi_difference, density_ratio, backscatter, true, exterior, nullptr,
                        &overlap_cache, true, truncate_tails, false, true);
                }
                previous_recovery_degree = degree;
            }
            if (order.singular && !hard && (!overlap_quadrature.active() || overlap_quadrature.adaptive)) {
                // Retry only the failed coefficient-overlap order without enlarging either degree support.
                order = solve_order<T>(m, n_max, n_max, chi_body, xi, incidence.cosine,
                    phi_difference, density_ratio, backscatter, true, exterior,
                    overlap_quadrature.active() ? &overlap_quadrature : nullptr,
                    nullptr, true, truncate_tails);
                if (order.singular)
                    throw std::runtime_error("PSMS fixed-support direct and Robin systems are singular at m=" +
                        std::to_string(m) + ", degree=" + std::to_string(n_max) + ".");
                // Finite Robin projection differs from direct elimination and requires physical trace checks.
                const T pressure_reference = T(1) / precsqrt(T(2));
                const T derivative_reference = chi_sw * precsqrt((incidence.cosine * incidence.cosine +
                    xi * xi * incidence.sine * incidence.sine / ((xi - T(1)) * (xi + T(1)))) / T(6));
                check_boundary(order, chi_sw, xi, incidence, pressure_reference,
                    derivative_reference, recovery_workspace);
                const auto budgets = boundary_budgets(order.boundary, m,
                    pressure_reference, derivative_reference);
                if (!order.boundary.quadrature_converged || !boundary_passes(order.boundary, budgets))
                    throw std::runtime_error("PSMS same-support Robin recovery failed the physical boundary checks at m=" +
                        std::to_string(m) + ", degree=" + std::to_string(n_max) + ".");
            }
            // Explicit liquid quadrature retains its original equations and singularity rejection.
            if (order.singular)
                throw std::runtime_error("PSMS fixed-support full-fluid system is singular at m=" +
                    std::to_string(m) + ", degree=" + std::to_string(n_max) +
                    ". Check the prescribed modal support and precision.");
            // Solve every prescribed order once so a small prefix cannot hide a later contributing band.
            add_field(field, correction, order.field);
            if (truncate_tails && !candidate_complete) {
                const auto band = adaptive_outgoing_band(order);
                add_field(truncated, truncated_correction, band);
                const T magnitude = combine_norms(precabs(band.real()), precabs(band.imag()));
                max_band = std::max(max_band, magnitude);
                if (psms_modal_band_is_negligible<T>(m, m_max, magnitude, previous_band, max_band))
                    ++small_bands;
                else small_bands = 0;
                previous_band = magnitude;
                if (small_bands >= 2) candidate_complete = true;
            }
        }
        if (!recovery_required || recovering) break;
    }
    finite(field);
    if (truncate_tails) {
        // Accept truncation only within the existing budget for the complete prescribed-support field.
        const auto difference = truncated - field;
        const T change = combine_norms(precabs(difference.real()), precabs(difference.imag()));
        const T reference = combine_norms(precabs(field.real()), precabs(field.imag()));
        if (change <= numerical_budget(reference)) return truncated;
    }
    return field;
}

template<typename T> std::vector<std::vector<std::vector<T>>> compute_smn_direction_matrices(
        int m_max, int n_max, T chi, const std::vector<SpheroidalDirection<T>>& directions) {
    if (m_max < 0 || n_max < m_max)
        throw std::invalid_argument("Spheroidal angular matrix dimensions are invalid.");
    checked_profcn_integer(2LL * n_max + 1, "twice the maximum angular degree plus one");
    std::vector<std::vector<std::vector<T>>> result(directions.size(),
        std::vector<std::vector<T>>(m_max + 1, std::vector<T>(n_max + 1, precnan<T>())));
    if (directions.empty()) return result;
    std::vector<T> arguments;
    for (const auto& direction : directions) arguments.push_back(direction.cosine);
    psms_automatic::LegendreWorkspace<T> workspace;
    // Reuse each normalized expansion across all requested receive directions.
    for (int m = 0; m <= m_max; ++m) {
        auto batch = cprofcn_mblock_coefficients<T>(chi, m, 1, n_max - m + 1,
            arguments, 0, 0, T(0));
        psms_automatic::evaluate_coefficient_angular_values(batch, m, directions, workspace);
        for (size_t angle = 0; angle < directions.size(); ++angle)
            for (int degree = m; degree <= n_max; ++degree)
                result[angle][m][degree] = extract_angular_value_from_mblock(batch.functions,
                    0, degree - m, int(angle));
    }
    return result;
}

template<typename T> std::vector<std::vector<T>> compute_smn_direction_matrix(
        int m_max, int n_max, T chi, const SpheroidalDirection<T>& direction) {
    return compute_smn_direction_matrices<T>(m_max, n_max, chi, {direction})[0];
}

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
    const std::string& Amn_method,
    bool is_backscatter,
    bool adaptive = false,
    bool vectorized = false
) {
    T eta_body = preccos(theta_body);
    T eta_scatter = preccos(theta_scatter);
    if (Amn_method == "Amn_fluid_gas" || Amn_method == "Amn_fluid_coefficients" ||
        (is_quad_precision_v<T> && Amn_method == "Amn_fluid")) {
        eta_body = psms_direction_cosine(theta_body);
        eta_scatter = psms_direction_cosine(theta_scatter);
        return psms_full_coupled_field<T>(
            m_max, n_max, chi_sw, chi_body, xi, eta_body, eta_scatter,
            phi_body - phi_scatter, density_body / density_sw,
            nodes, weights, is_backscatter, Amn_method != "Amn_fluid",
            adaptive && Amn_method == "Amn_fluid"
        );
    }
    if (is_backscatter && Amn_method == "Amn_fixed_rigid") {
        return psms_backscatter_fixed_rigid<T>(
            m_max, n_max, chi_sw, xi, eta_body, adaptive
        );
    }

    if (is_backscatter && Amn_method == "Amn_pressure_release") {
        return psms_backscatter_pressure_release<T>(
            m_max, n_max, chi_sw, xi, eta_body, adaptive
        );
    }

#if ACOUSTICTS_HAVE_QUADMATH
    if constexpr (is_quad_precision_v<T>) {
        // The gas shortcut also changes the retained degree limit for each m.
        // Keep gas solves on the full retained system below so recognizing
        // double-origin backscatter angles does not change their truncation.
        if (is_backscatter && Amn_method == "Amn_fluid") {
            return psms_backscatter_fluid_adaptive<T>(
                m_max, n_max, chi_sw, chi_body, xi, eta_body,
                density_body, density_sw, nodes, weights, true, adaptive, false
            );
        }
    }
#endif

    // Determine appropriate Amn expansion matrix computation
    std::vector<std::vector<std::complex<T>>> Amn;
    std::vector<std::vector<T>> smn_body_matrix;

    // Fluid-filled
    if (Amn_method == "Amn_fluid" || Amn_method == "Amn_fluid_gas") {
        auto fluid_result = fluid_Amn_triangular<T>(
            m_max, n_max, chi_sw, chi_body, xi,
            eta_body, density_body, density_sw, nodes, weights, vectorized,
            Amn_method == "Amn_fluid_gas"
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
    auto smn_scatter_matrix = compute_smn_matrix<T>(
        m_max, n_max, chi_sw, eta_scatter, true
    );

    // Compute linear scattering coefficient
    return compute_fbs<T>(
        m_max, n_max, azimuth, smn_body_matrix, smn_scatter_matrix, Amn
    );
}

// ============================================================================
