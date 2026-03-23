#pragma once

// ============================================================================
// PSMS BACKSCATTER AND KERNEL HELPERS
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
    const std::vector<double>& weights
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
    out.Amn_tri = solve_fluid_Amn<double>(kernels.rhs, kernels.K3_kernel);
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
    const std::vector<T>& weights,
    bool vectorized = false
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
    if constexpr (std::is_same_v<T, double>) {
        if (vectorized) {
            return fluid_Amn_triangular_vectorized_double(
                m_max, n_max, chi_sw, chi_body, xi, eta_body, density_body,
                density_sw, nodes, weights
            );
        }
    }
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
    bool adaptive = false,
    bool vectorized = false
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
            eta_body, density_body, density_sw, nodes, weights, vectorized
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
