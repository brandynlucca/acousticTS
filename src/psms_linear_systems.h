#pragma once

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
    // Combined LU factors together with the row/column permutations used by
    // the complete-pivoting decomposition.
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
    // Apply the stored permutations, perform forward substitution through L,
    // then backward substitution through U.
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
    // Dense matrix-vector product on the flattened row-major kernel storage.
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
std::vector<std::vector<std::complex<T>>> solve_fluid_Amn_default(
    const std::vector<std::vector<std::complex<T>>>& rhs,
    const std::vector<std::vector<std::complex<T>>>& K3_kernel
) {
    return solve_fluid_Amn_divide_and_conquer<T>(rhs, K3_kernel);
}

template<typename T>
std::vector<std::vector<std::complex<T>>> solve_fluid_Amn_gas(
    const std::vector<std::vector<std::complex<T>>>& rhs,
    const std::vector<std::vector<std::complex<T>>>& K3_kernel
) {
    // Keep the gas-filled branch fully separate from the liquid-filled path so
    // future stabilisation work cannot spill back into liquid_filled PSMS.
    return solve_fluid_Amn_divide_and_conquer<T>(rhs, K3_kernel);
}

template<typename T>
std::vector<std::vector<std::complex<T>>> solve_fluid_Amn(
    const std::vector<std::vector<std::complex<T>>>& rhs,
    const std::vector<std::vector<std::complex<T>>>& K3_kernel,
    bool gas_like = false
) {
    if (gas_like) {
        return solve_fluid_Amn_gas<T>(rhs, K3_kernel);
    }
    return solve_fluid_Amn_default<T>(rhs, K3_kernel);
}

template<typename T>
std::vector<std::vector<std::complex<T>>> solve_fluid_t_blocks_native(
    const std::vector<std::vector<std::complex<T>>>& K1_kernel,
    const std::vector<std::vector<std::complex<T>>>& K3_kernel
) {
    int m_max = static_cast<int>(K3_kernel.size()) - 1;
    std::vector<std::vector<std::complex<T>>> result(m_max + 1);

    // Each retained order m produces one dense block solve for the T-matrix
    // block associated with that order.
    for (int m = 0; m <= m_max; ++m) {
        int size = static_cast<int>(std::sqrt(K3_kernel[m].size()));
        auto fac = lup_decompose<T>(K3_kernel[m], size);
        if (fac.singular) {
            throw std::runtime_error("Fluid kernel system is singular or ill-conditioned");
        }

        std::vector<std::complex<T>> block(size * size, std::complex<T>(0, 0));
        for (int col = 0; col < size; ++col) {
            std::vector<std::complex<T>> rhs_col(size, std::complex<T>(0, 0));
            for (int row = 0; row < size; ++row) {
                rhs_col[row] = K1_kernel[m][row * size + col];
            }
            auto sol_col = lup_solve(fac, rhs_col);
            for (int row = 0; row < size; ++row) {
                block[row * size + col] = sol_col[row];
            }
        }

        result[m] = block;
    }

    return result;
}

template<>
inline std::vector<std::vector<std::complex<double>>> solve_fluid_t_blocks_native<double>(
    const std::vector<std::vector<std::complex<double>>>& K1_kernel,
    const std::vector<std::vector<std::complex<double>>>& K3_kernel
) {
    int m_max = static_cast<int>(K3_kernel.size()) - 1;
    std::vector<std::vector<std::complex<double>>> result(m_max + 1);

    // The double-precision specialization uses Armadillo first, then falls
    // back to an SVD pseudoinverse if the direct solve is unstable.
    for (int m = 0; m <= m_max; ++m) {
        int size = static_cast<int>(std::sqrt(K3_kernel[m].size()));
        arma::Mat<std::complex<double>> K3_arma(size, size);
        arma::Mat<std::complex<double>> K1_arma(size, size);
        for (int row = 0; row < size; ++row) {
            for (int col = 0; col < size; ++col) {
                K3_arma(row, col) = K3_kernel[m][row * size + col];
                K1_arma(row, col) = K1_kernel[m][row * size + col];
            }
        }

        arma::Mat<std::complex<double>> block;
        bool ok = arma::solve(block, K3_arma, K1_arma, arma::solve_opts::fast);
        if (!ok || !block.is_finite()) {
            ok = arma::solve(block, K3_arma, K1_arma);
        }
        if (!ok || !block.is_finite()) {
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
            arma::Mat<std::complex<double>> diag_dinv = arma::diagmat(
                arma::conv_to<arma::Col<std::complex<double>>>::from(d_inv)
            );
            block = V * diag_dinv * U.t() * K1_arma;
        }

        std::vector<std::complex<double>> out_block(size * size, std::complex<double>(0, 0));
        for (int row = 0; row < size; ++row) {
            for (int col = 0; col < size; ++col) {
                out_block[row * size + col] = block(row, col);
            }
        }
        result[m] = std::move(out_block);
    }

    return result;
}

template<typename T>
std::vector<std::vector<std::complex<T>>> solve_fluid_t_blocks(
    const std::vector<std::vector<std::complex<T>>>& K1_kernel,
    const std::vector<std::vector<std::complex<T>>>& K3_kernel
) {
    return solve_fluid_t_blocks_native<T>(K1_kernel, K3_kernel);
}

template<>
std::vector<std::vector<std::complex<acousticts_quad_t>>>
solve_fluid_Amn_default<acousticts_quad_t>(
    const std::vector<std::vector<std::complex<acousticts_quad_t>>>& rhs,
    const std::vector<std::vector<std::complex<acousticts_quad_t>>>& K3_kernel
) {
#if ACOUSTICTS_HAVE_QUADMATH
    return solve_fluid_Amn_native<acousticts_quad_t>(rhs, K3_kernel);
#else
    return solve_fluid_Amn_divide_and_conquer<acousticts_quad_t>(
        rhs, K3_kernel
    );
#endif
}

template<>
std::vector<std::vector<std::complex<acousticts_quad_t>>>
solve_fluid_Amn_gas<acousticts_quad_t>(
    const std::vector<std::vector<std::complex<acousticts_quad_t>>>& rhs,
    const std::vector<std::vector<std::complex<acousticts_quad_t>>>& K3_kernel
) {
    // The gas-filled branch is intentionally solved with the pseudoinverse
    // backend even in quad mode. That keeps the extreme fluid-fluid contrast
    // from tripping the complete-pivoting solve while leaving liquid-filled
    // PSMS on the faster native quad path.
    return solve_fluid_Amn_divide_and_conquer<acousticts_quad_t>(
        rhs, K3_kernel
    );
}

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
