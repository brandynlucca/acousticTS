#pragma once

// EXPANSION COEFFICIENT MATRIX SOLVERS
// -----------------------------------------------------------
// Double precision uses column-scaled Armadillo SVD for the dense modal kernel
// systems. Quad precision uses the native pivoted solver below.
inline arma::cx_mat solve_fluid_system_svd(
    const arma::cx_mat& matrix,
    const arma::cx_mat& rhs
) {
    // Radial functions at consecutive degrees can differ by many orders of
    // magnitude. Solve for column-scaled coefficients, then undo the scaling.
    arma::rowvec scale = arma::max(arma::abs(matrix), 0);
    arma::cx_mat balanced = matrix;
    for (arma::uword j = 0; j < matrix.n_cols; ++j) {
        if (scale[j] == 0.0) scale[j] = 1.0;
        balanced.col(j) /= scale[j];
    }

    arma::cx_mat U, V;
    arma::vec singular;
    if (!arma::svd(U, singular, V, balanced)) {
        throw std::runtime_error("SVD failed");
    }
    double tolerance = static_cast<double>(matrix.n_cols) * singular.max() *
        std::numeric_limits<double>::epsilon();
    bool full_rank = arma::all(singular > tolerance);

    // Column scaling changes the minimum-norm convention for deficient
    // systems. Retain the original unscaled pseudoinverse for those cases.
    if (!full_rank) {
        if (!arma::svd(U, singular, V, matrix)) {
            throw std::runtime_error("SVD failed");
        }
        tolerance = static_cast<double>(matrix.n_cols) * singular.max() *
            std::numeric_limits<double>::epsilon();
    }

    arma::cx_mat projected = U.t() * rhs;
    for (arma::uword i = 0; i < singular.n_elem; ++i) {
        if (singular[i] > tolerance) projected.row(i) /= singular[i];
        else projected.row(i).zeros();
    }
    arma::cx_mat result = V * projected;
    if (full_rank) {
        for (arma::uword i = 0; i < result.n_rows; ++i) result.row(i) /= scale[i];
    }
    return result;
}

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
        arma::Col<std::complex<double>> A = solve_fluid_system_svd(K3_arma, b);
        
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
                precnan<T>(),
                precnan<T>()
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
                    precnan<T>(),
                    precnan<T>()
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
class LUPPivotBounds {
    T half_maximum;

public:
    explicit LUPPivotBounds(const std::complex<T>& pivot) { update(pivot); }

    void update(const std::complex<T>& pivot) {
        half_maximum = std::max(precabs(pivot.real()), precabs(pivot.imag())) * T(0.5);
    }

    bool rejects(const std::complex<T>& candidate) const {
        if constexpr (is_quad_precision_v<T>) {
            return precabs(candidate.real()) <= half_maximum &&
                precabs(candidate.imag()) <= half_maximum;
        }
        return false;
    }
};

#if ACOUSTICTS_HAVE_EXACT_QUAD_DOT
template<>
class LUPPivotBounds<acousticts_quad_t> {
    using Quad = acousticts_quad_t;
    struct Magnitude {
        std::uint64_t low, high;

        explicit Magnitude(Quad value) {
            std::uint64_t words[2];
            std::memcpy(words, &value, sizeof(value));
            low = words[0];
            high = words[1] & UINT64_C(0x7fffffffffffffff);
        }

        bool finite() const { return (high >> 48) != 0x7fff; }

        bool operator<=(const Magnitude& other) const {
            return high < other.high || (high == other.high && low <= other.low);
        }
    };
    Quad half_maximum;
    Magnitude minimum{Quad(0)}, maximum{Quad(0)}, half{Quad(0)};
    bool finite_pivot;

public:
    explicit LUPPivotBounds(const std::complex<Quad>& pivot) { update(pivot); }

    void update(const std::complex<Quad>& pivot) {
        const Quad real = precabs(pivot.real()), imag = precabs(pivot.imag());
        minimum = Magnitude(std::min(real, imag));
        maximum = Magnitude(std::max(real, imag));
        half_maximum = std::max(real, imag) * Quad(0.5);
        half = Magnitude(half_maximum);
        finite_pivot = Magnitude(real).finite() && Magnitude(imag).finite();
    }

    bool rejects(const std::complex<Quad>& candidate) const {
        const Magnitude real(candidate.real()), imag(candidate.imag());
        if (finite_pivot && real.finite() && imag.finite()) {
            // Sign-cleared binary128 words preserve finite magnitude ordering.
            if (real <= half && imag <= half) return true;
            // Sorted component dominance cannot increase the rounded squared norm.
            return real <= maximum && imag <= maximum &&
                (real <= minimum || imag <= minimum);
        }
        // Keep the original floating comparisons for nonfinite components.
        return precabs(candidate.real()) <= half_maximum &&
            precabs(candidate.imag()) <= half_maximum;
    }
};
#endif

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
        LUPPivotBounds<T> pivot_bounds(fac.lu[k * n + k]);
        for (int i = k + 1; i < n; ++i) {
            for (int j = k; j < n; ++j) {
                // Screen only candidates that cannot exceed the current pivot norm.
                if (pivot_bounds.rejects(fac.lu[i * n + j])) continue;
                T cand_abs_sq = complex_abs_sq(fac.lu[i * n + j]);
                if (cand_abs_sq > pivot_abs_sq) {
                    pivot_row = i;
                    pivot_col = j;
                    pivot_abs_sq = cand_abs_sq;
                    pivot_bounds.update(fac.lu[i * n + j]);
                }
            }
        }
        for (int j = k; j < n; ++j) {
            if (pivot_bounds.rejects(fac.lu[k * n + j])) continue;
            T cand_abs_sq = complex_abs_sq(fac.lu[k * n + j]);
            if (cand_abs_sq > pivot_abs_sq) {
                pivot_row = k;
                pivot_col = j;
                pivot_abs_sq = cand_abs_sq;
                pivot_bounds.update(fac.lu[k * n + j]);
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
        m_max + 1, std::vector<T>(n_max + 1, precnan<T>())
    );
    out.incident.derivative.assign(
        m_max + 1, std::vector<T>(n_max + 1, precnan<T>())
    );
    out.scattering.value.assign(
        m_max + 1,
        std::vector<std::complex<T>>(
            n_max + 1,
            std::complex<T>(
                precnan<T>(), precnan<T>()
            )
        )
    );
    out.scattering.derivative.assign(
        m_max + 1,
        std::vector<std::complex<T>>(
            n_max + 1,
            std::complex<T>(
                precnan<T>(), precnan<T>()
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
        m_max + 1, std::vector<T>(n_max + 1, precnan<T>())
    );
    out.radial.incident.value.assign(
        m_max + 1, std::vector<T>(n_max + 1, precnan<T>())
    );
    out.radial.incident.derivative.assign(
        m_max + 1, std::vector<T>(n_max + 1, precnan<T>())
    );
    out.radial.scattering.value.assign(
        m_max + 1,
        std::vector<std::complex<T>>(
            n_max + 1,
            std::complex<T>(
                precnan<T>(), precnan<T>()
            )
        )
    );
    out.radial.scattering.derivative.assign(
        m_max + 1,
        std::vector<std::complex<T>>(
            n_max + 1,
            std::complex<T>(
                precnan<T>(), precnan<T>()
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

        auto result_both = cprofcn<T>(chi_sw, m, size, {eta_body}, 2, iopnorm, 1, x1_r);

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
std::vector<std::complex<T>> solve_factored_system_lup_refined(
    const std::vector<std::complex<T>>& A,
    const std::vector<std::complex<T>>& b,
    const LUPFactorization<T>& fac,
    int n,
    int n_refine = 1
) {
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

#include <algorithm>
#include <cfenv>
#include <complex>
#include <cstddef>
#include <limits>
#include <new>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

// Native right-preconditioned GMRES with a factorized Robin operator.
namespace psms_robin_solver {
template<typename T> using Complex = std::complex<T>;

template<typename T> bool finite(const Complex<T>& value) {
    return is_finite_real(value.real()) && is_finite_real(value.imag());
}

template<typename T> void add(T value, T& sum, T& correction) {
    const T adjusted = value - correction;
    const T next = sum + adjusted;
    correction = (next - sum) - adjusted;
    sum = next;
}

template<typename T> void add(const Complex<T>& value, Complex<T>& sum, Complex<T>& correction) {
    T real = sum.real(), imag = sum.imag(), cr = correction.real(), ci = correction.imag();
    add(value.real(), real, cr);
    add(value.imag(), imag, ci);
    sum = {real, imag};
    correction = {cr, ci};
}

template<typename T> T component_max(const Complex<T>& value) {
    return std::max(precabs(value.real()), precabs(value.imag()));
}

template<typename T> T magnitude(const Complex<T>& value) {
    const T scale = component_max(value);
    if (scale == T(0)) return T(0);
    const T real = value.real() / scale, imag = value.imag() / scale;
    return scale * precsqrt(real * real + imag * imag);
}

template<typename T> T norm_infinity(const std::vector<Complex<T>>& values) {
    T result = 0;
    for (const auto& value : values) {
        if (!finite(value)) return precnan<T>();
        const T absolute = precabs(value.real()) + precabs(value.imag());
        if (!is_finite_real(absolute)) return precnan<T>();
        result = std::max(result, absolute);
    }
    return result;
}

template<typename T> T norm_two(const std::vector<Complex<T>>& values) {
    T scale = 0;
    for (const auto& value : values) {
        if (!finite(value)) return precnan<T>();
        scale = std::max(scale, component_max(value));
    }
    if (scale == T(0)) return T(0);
    T sum = 0, correction = 0;
    for (const auto& value : values) {
        const T real = value.real() / scale, imag = value.imag() / scale;
        add(real * real, sum, correction);
        add(imag * imag, sum, correction);
    }
    return scale * precsqrt(sum);
}

template<typename T> bool divide(const Complex<T>& numerator, const Complex<T>& denominator,
        Complex<T>& result) {
    const T scale = component_max(denominator);
    if (!finite(numerator) || !finite(denominator) || !(scale > T(0))) return false;
    auto first = numerator / scale;
    const auto second = denominator / scale;
    if (!finite(first)) return false;
    const T square = second.real() * second.real() + second.imag() * second.imag();
    first /= square;
    result = {first.real() * second.real() + first.imag() * second.imag(),
        first.imag() * second.real() - first.real() * second.imag()};
    return finite(result);
}

template<typename T> Complex<T> inner_product(const std::vector<Complex<T>>& first,
        const std::vector<Complex<T>>& second) {
#if ACOUSTICTS_HAVE_EXACT_QUAD_DOT
    if constexpr (std::is_same_v<T, acousticts_quad_t>) {
        if (std::fegetround() == FE_TONEAREST) {
            static thread_local quad_dot_product::Accumulator accumulator;
            accumulator.reset();
            bool valid = true;
            for (std::size_t index = 0; index < first.size(); ++index) {
                valid = accumulator.add(first[index].real(), second[index].real()) && valid;
                valid = accumulator.add(first[index].imag(), second[index].imag()) && valid;
            }
            const T real = valid ? accumulator.rounded() : precnan<T>();
            accumulator.reset();
            for (std::size_t index = 0; index < first.size(); ++index) {
                valid = accumulator.add(first[index].real(), second[index].imag()) && valid;
                valid = accumulator.add(-first[index].imag(), second[index].real()) && valid;
            }
            if (valid) return {real, accumulator.rounded()};
        }
    }
#endif
    T real = 0, imag = 0, cr = 0, ci = 0;
    for (std::size_t index = 0; index < first.size(); ++index) {
        add(first[index].real() * second[index].real(), real, cr);
        add(first[index].imag() * second[index].imag(), real, cr);
        add(first[index].real() * second[index].imag(), imag, ci);
        add(-first[index].imag() * second[index].real(), imag, ci);
    }
    return {real, imag};
}

template<typename T> struct Options {
    int restart = 32, maximum_iterations = 256;
    T relative_tolerance = T(64) * preceps<T>();
    T absolute_tolerance = 0;
};

template<typename T> bool rotation(const Complex<T>& first, const Complex<T>& second,
        T& cosine, Complex<T>& sine, Complex<T>& diagonal) {
    const T a = magnitude(first), b = magnitude(second), scale = std::max(a, b);
    if (!is_finite_real(scale)) return false;
    if (b == T(0)) { cosine = 1; sine = 0; diagonal = first; return true; }
    if (a == T(0)) { cosine = 0; sine = std::conj(second) / b; diagonal = b; return true; }
    const T length = scale * precsqrt((a / scale) * (a / scale) + (b / scale) * (b / scale));
    if (!is_finite_real(length)) return false;
    const auto phase = first / a;
    cosine = a / length;
    sine = phase * (std::conj(second) / length);
    diagonal = phase * length;
    return finite(sine) && finite(diagonal);
}



template<typename T> struct Result {
    bool success = false, needs_lu_fallback = true;
    std::string status = "not_solved";
    int iterations = 0, restarts = 0, true_residual_checks = 0;
    T residual_norm = 0, backward_error = 0, estimated_residual = 0, operator_norm = 0;
    T rhs_relative_residual = 0;
    bool operator_norm_is_bound = true;
    std::vector<Complex<T>> solution;
};


// Accumulate real-coefficient projections in the requested native precision.
template<typename T> Complex<T> projection(const T* coefficients, std::size_t stride,
        const std::vector<Complex<T>>& input) {
#if ACOUSTICTS_HAVE_EXACT_QUAD_DOT
    if constexpr (std::is_same_v<T, acousticts_quad_t>) {
        if (std::fegetround() == FE_TONEAREST) {
            static thread_local quad_dot_product::Accumulator accumulator;
            accumulator.reset();
            bool valid = true;
            for (std::size_t index = 0; index < input.size(); ++index)
                valid = accumulator.add(coefficients[index * stride], input[index].real()) && valid;
            const T real = valid ? accumulator.rounded() : precnan<T>();
            accumulator.reset();
            for (std::size_t index = 0; index < input.size(); ++index)
                valid = accumulator.add(coefficients[index * stride], input[index].imag()) && valid;
            if (valid) return {real, accumulator.rounded()};
        }
    }
#endif
    Complex<T> sum(0, 0), correction(0, 0);
    for (std::size_t index = 0; index < input.size(); ++index)
        add(coefficients[index * stride] * input[index], sum, correction);
    return sum;
}

inline bool valid_dimensions(std::size_t count, int rows, int columns) {
    return rows >= 0 && columns >= 0 && (columns == 0 || std::size_t(rows) <=
        std::numeric_limits<std::size_t>::max() / std::size_t(columns)) &&
        count == std::size_t(rows) * std::size_t(columns);
}

template<typename T> std::vector<Complex<T>> weighted_projection(const std::vector<T>& overlap,
        const std::vector<Complex<T>>& interior_response, const std::vector<Complex<T>>& input,
        int rows, int columns) {
    std::vector<Complex<T>> weighted(rows);
    for (int row = 0; row < rows; ++row) {
        weighted[row] = interior_response[row] * projection(overlap.data() + std::size_t(row) * columns, 1, input);
        if (!finite(weighted[row])) throw std::runtime_error("Factorized Robin projection range.");
    }
    return weighted;
}

template<typename T> std::vector<Complex<T>> apply_map(const std::vector<T>& overlap,
        const std::vector<Complex<T>>& interior_response, const std::vector<Complex<T>>& input,
        int rows, int columns) {
    if (!valid_dimensions(overlap.size(), rows, columns) ||
            interior_response.size() != std::size_t(rows) || input.size() != std::size_t(columns))
        throw std::runtime_error("Factorized Robin map dimension mismatch.");
    if (columns == 0) return {};
    for (const auto& value : overlap)
        if (!is_finite_real(value)) throw std::runtime_error("Nonfinite Robin overlap.");
    for (const auto& value : interior_response)
        if (!finite(value)) throw std::runtime_error("Nonfinite Robin response.");
    for (const auto& value : input)
        if (!finite(value)) throw std::runtime_error("Nonfinite Robin input.");
    const auto weighted = weighted_projection(overlap, interior_response, input, rows, columns);
    std::vector<Complex<T>> output(columns);
    for (int column = 0; column < columns; ++column) {
        output[column] = rows == 0 ? Complex<T>(0, 0) :
            projection(overlap.data() + column, columns, weighted);
        if (!finite(output[column])) throw std::runtime_error("Factorized Robin map range.");
    }
    return output;
}

template<typename T> struct Operator {
    int rows, size;
    const std::vector<T>& overlap;
    const std::vector<Complex<T>>& interior_response;
    const std::vector<Complex<T>>& exterior_ratio;
    std::vector<Complex<T>> diagonal;
    T norm_bound = 0;

    bool initialize() {
        std::vector<T> weighted_row_sums(rows);
        for (int row = 0; row < rows; ++row) {
            T sum = 0, correction = 0;
            for (int column = 0; column < size; ++column)
                add(precabs(overlap[std::size_t(row) * size + column]), sum, correction);
            weighted_row_sums[row] = (precabs(interior_response[row].real()) + precabs(interior_response[row].imag())) * sum;
            if (!is_finite_real(weighted_row_sums[row])) return false;
        }
        diagonal.resize(size);
        std::vector<Complex<T>> weighted_column(rows);
        for (int column = 0; column < size; ++column) {
            T bound = precabs(exterior_ratio[column].real()) + precabs(exterior_ratio[column].imag()), correction = 0;
            for (int row = 0; row < rows; ++row) {
                const T coefficient = overlap[std::size_t(row) * size + column];
                weighted_column[row] = coefficient * interior_response[row];
                add(precabs(coefficient) * weighted_row_sums[row], bound, correction);
            }
            diagonal[column] = exterior_ratio[column] - (rows == 0 ? Complex<T>(0, 0) :
                projection(overlap.data() + column, size, weighted_column));
            if (!finite(diagonal[column]) || !is_finite_real(bound)) return false;
            norm_bound = std::max(norm_bound, bound);
        }
        return true;
    }

    // Evaluate b-Ax afresh, retaining the affine terms in the native exact sum.
    bool action(const std::vector<Complex<T>>& input, std::vector<Complex<T>>& output,
            const std::vector<Complex<T>>* rhs = nullptr) const {
        if (input.size() != std::size_t(size)) return false;
        for (const auto& value : input) if (!finite(value)) return false;
        const auto weighted = weighted_projection(overlap, interior_response, input, rows, size);
        output.assign(size, Complex<T>(0, 0));
        for (int column = 0; column < size; ++column) {
#if ACOUSTICTS_HAVE_EXACT_QUAD_DOT
            if constexpr (std::is_same_v<T, acousticts_quad_t>) {
                if (std::fegetround() == FE_TONEAREST) {
                    static thread_local quad_dot_product::Accumulator accumulator;
                    const T sign = rhs ? T(-1) : T(1);
                    accumulator.reset();
                    bool valid = !rhs || accumulator.add((*rhs)[column].real(), T(1));
                    valid = accumulator.add(sign * exterior_ratio[column].real(), input[column].real()) && valid;
                    valid = accumulator.add(-sign * exterior_ratio[column].imag(), input[column].imag()) && valid;
                    for (int row = 0; row < rows; ++row)
                        valid = accumulator.add(-sign * overlap[std::size_t(row) * size + column],
                            weighted[row].real()) && valid;
                    const T real = valid ? accumulator.rounded() : precnan<T>();
                    accumulator.reset();
                    valid = (!rhs || accumulator.add((*rhs)[column].imag(), T(1))) && valid;
                    valid = accumulator.add(sign * exterior_ratio[column].real(), input[column].imag()) && valid;
                    valid = accumulator.add(sign * exterior_ratio[column].imag(), input[column].real()) && valid;
                    for (int row = 0; row < rows; ++row)
                        valid = accumulator.add(-sign * overlap[std::size_t(row) * size + column],
                            weighted[row].imag()) && valid;
                    if (!valid) return false;
                    output[column] = {real, accumulator.rounded()};
                    if (!finite(output[column])) return false;
                    continue;
                }
            }
#endif
            Complex<T> correction(0, 0);
            const T sign = rhs ? T(-1) : T(1);
            output[column] = rhs ? (*rhs)[column] : Complex<T>(0, 0);
            add(sign * exterior_ratio[column] * input[column], output[column], correction);
            for (int row = 0; row < rows; ++row)
                add(-sign * overlap[std::size_t(row) * size + column] * weighted[row],
                    output[column], correction);
            if (!finite(output[column])) return false;
        }
        return true;
    }

    bool apply(const std::vector<Complex<T>>& input, std::vector<Complex<T>>& output) const {
        return action(input, output);
    }
};

template<typename T> Result<T> solve_impl(const std::vector<T>& overlap, const std::vector<Complex<T>>& interior_response,
        const std::vector<Complex<T>>& exterior_ratio, const std::vector<Complex<T>>& rhs,
        int rows, int columns, const Options<T>& options = Options<T>()) {
    Result<T> result;
    if (rhs.size() > std::size_t(std::numeric_limits<int>::max())) {
        result.status = "dimension_range";
        return result;
    }
    const int size = int(rhs.size());
    if (!valid_dimensions(overlap.size(), rows, columns) || columns != size ||
            interior_response.size() != std::size_t(rows) || exterior_ratio.size() != rhs.size()) {
        result.status = "dimension_mismatch"; return result;
    }
    if (options.restart <= 0 || options.maximum_iterations < 0 ||
            !is_finite_real(options.relative_tolerance) || !(options.relative_tolerance > T(0)) ||
            !is_finite_real(options.absolute_tolerance) || options.absolute_tolerance < T(0)) {
        result.status = "invalid_controls";
        return result;
    }
    for (std::size_t index = 0; index < rhs.size(); ++index) {
        if (!finite(rhs[index])) {
            result.status = "nonfinite_input";
            return result;
        }
    }
    result.solution.assign(size, Complex<T>(0, 0));


    for (const auto& value : overlap)
        if (!is_finite_real(value)) { result.status = "nonfinite_overlap"; return result; }
    for (const auto& value : interior_response)
        if (!finite(value)) { result.status = "nonfinite_response"; return result; }
    for (const auto& value : exterior_ratio)
        if (!finite(value)) { result.status = "nonfinite_diagonal"; return result; }
    Operator<T> operation{rows, size, overlap, interior_response, exterior_ratio};
    if (!operation.initialize()) { result.status = "operator_bound_range"; return result; }
    result.operator_norm = operation.norm_bound;
    const T diagonal_guard = T(64) * preceps<T>() * std::max(T(1), result.operator_norm);
    for (const auto& value : operation.diagonal) {
        if (!(magnitude(value) > diagonal_guard)) {
            result.status = "diagonal_preconditioner_unavailable";
            return result;
        }
    }
    std::vector<Complex<T>> residual(rhs), candidate;
    const T rhs_norm = norm_infinity(rhs);
    if (!is_finite_real(rhs_norm)) { result.status = "rhs_norm_range"; return result; }
    auto inspect = [&](const std::vector<Complex<T>>& value) {
        ++result.true_residual_checks;
        if (!operation.action(value, residual, &rhs)) { result.status = "true_residual_range"; return -1; }
        result.residual_norm = norm_infinity(residual);
        const T reference = result.operator_norm * norm_infinity(value) + rhs_norm;
        const T allowance = options.absolute_tolerance + options.relative_tolerance * reference;
        if (!is_finite_real(result.residual_norm) || !is_finite_real(reference) || !is_finite_real(allowance)) {
            result.status = "residual_budget_range";
            return -1;
        }
        result.backward_error = reference == T(0) ? result.residual_norm : result.residual_norm / reference;
        result.rhs_relative_residual = rhs_norm == T(0) ? result.residual_norm : result.residual_norm / rhs_norm;
        const T rhs_allowance = T(64) * preceps<T>() * rhs_norm;
        return result.residual_norm <= allowance && result.residual_norm <= rhs_allowance ? 1 : 0;
    };
    auto accept = [&]() {
        result.success = true;
        result.needs_lu_fallback = false;
        result.status = "converged_true_residual";
    };
    int checked = inspect(result.solution);
    if (checked < 0) return result;
    if (checked > 0) { accept(); return result; }
    const int restart = std::min(size, options.restart);
    // Retain each actual native diagonal-preconditioned Arnoldi vector.
    while (result.iterations < options.maximum_iterations) {
        const int width = std::min(restart, options.maximum_iterations - result.iterations);
        const T beta = norm_two(residual);
        if (!is_finite_real(beta) || !(beta > T(0))) {
            result.status = "residual_norm_range";
            return result;
        }
        std::vector<std::vector<Complex<T>>> basis(std::size_t(width) + 1, std::vector<Complex<T>>(size));
        std::vector<std::vector<Complex<T>>> preconditioned(width, std::vector<Complex<T>>(size));
        std::vector<Complex<T>> hessenberg((std::size_t(width) + 1) * width), transformed(std::size_t(width) + 1);
        std::vector<Complex<T>> sines(width), coefficients;
        std::vector<T> cosines(width);
        const auto base = result.solution;
        for (int row = 0; row < size; ++row) basis[0][row] = residual[row] / beta;
        transformed[0] = beta;
        bool restart_ready = false;
        for (int column = 0; column < width; ++column) {

            for (int row = 0; row < size; ++row) {
                if (!divide(basis[column][row], operation.diagonal[row],
                        preconditioned[column][row])) {
                    result.status = "diagonal_preconditioner_range";
                    return result;
                }
            }
            std::vector<Complex<T>> next;
            if (!operation.apply(preconditioned[column], next)) { result.status = "arnoldi_product_range"; return result; }
            // Two modified Gram-Schmidt passes preserve native complex orthogonality.
            for (int pass = 0; pass < 2; ++pass) {
                for (int row = 0; row <= column; ++row) {
                    const auto projection = inner_product(basis[row], next);
                    if (!finite(projection)) { result.status = "arnoldi_inner_product_range"; return result; }
                    hessenberg[std::size_t(row) * width + column] += projection;
                    for (int index = 0; index < size; ++index) next[index] -= projection * basis[row][index];
                }
            }
            const T next_norm = norm_two(next);
            if (!is_finite_real(next_norm)) { result.status = "arnoldi_norm_range"; return result; }
            hessenberg[std::size_t(column + 1) * width + column] = next_norm;
            if (next_norm > T(0))
                for (int row = 0; row < size; ++row) basis[column + 1][row] = next[row] / next_norm;
            for (int previous = 0; previous < column; ++previous) {
                auto& upper = hessenberg[std::size_t(previous) * width + column];
                auto& lower = hessenberg[std::size_t(previous + 1) * width + column];
                const auto value = cosines[previous] * upper + sines[previous] * lower;
                lower = -std::conj(sines[previous]) * upper + cosines[previous] * lower;
                upper = value;
            }
            auto& diagonal = hessenberg[std::size_t(column) * width + column];
            if (!rotation(diagonal, hessenberg[std::size_t(column + 1) * width + column],
                    cosines[column], sines[column], diagonal)) {
                result.status = "arnoldi_rotation_range";
                return result;
            }
            hessenberg[std::size_t(column + 1) * width + column] = 0;
            transformed[column + 1] = -std::conj(sines[column]) * transformed[column];
            transformed[column] *= cosines[column];
            result.estimated_residual = magnitude(transformed[column + 1]);
            ++result.iterations;
            const T trigger = options.absolute_tolerance + options.relative_tolerance * rhs_norm;
            const bool at_end = column + 1 == width;
            if (!(result.estimated_residual <= trigger || next_norm == T(0) || at_end)) continue;
            coefficients.assign(column + 1, Complex<T>(0, 0));
            T triangular_norm = 0;
            for (int row = 0; row <= column; ++row)
                for (int col = row; col <= column; ++col)
                    triangular_norm = std::max(triangular_norm, magnitude(hessenberg[std::size_t(row) * width + col]));
            const T triangular_guard = T(64) * preceps<T>() * std::max(T(1), triangular_norm);
            for (int row = column; row >= 0; --row) {
                Complex<T> value = transformed[row], correction(0, 0);
                for (int col = row + 1; col <= column; ++col)
                    add(-hessenberg[std::size_t(row) * width + col] * coefficients[col], value, correction);
                const auto pivot = hessenberg[std::size_t(row) * width + row];
                if (!(magnitude(pivot) > triangular_guard) || !divide(value, pivot, coefficients[row])) {
                    result.status = "projected_pivot_guard";
                    return result;
                }
            }
            candidate = base;
            for (int row = 0; row < size; ++row) {
                Complex<T> correction(0, 0);
                for (int col = 0; col <= column; ++col)
                    add(preconditioned[col][row] * coefficients[col], candidate[row], correction);
            }
            checked = inspect(candidate);
            if (checked < 0) return result;
            // Accept only after evaluating the original factorized residual.
            if (checked > 0) {
                result.solution = std::move(candidate);
                accept();
                return result;
            }
            if (next_norm == T(0)) { result.status = "arnoldi_breakdown_without_convergence"; return result; }
            if (at_end) {
                if (candidate == base) { result.status = "native_iteration_stagnation"; return result; }
                result.solution = std::move(candidate);
                restart_ready = true;
            }
        }
        if (!restart_ready) { result.status = "restart_state_invalid"; return result; }
        ++result.restarts;
    }
    result.status = "iteration_budget_exhausted";
    return result;
}

template<typename T> Result<T> solve(const std::vector<T>& overlap, const std::vector<Complex<T>>& interior_response,
        const std::vector<Complex<T>>& exterior_ratio, const std::vector<Complex<T>>& rhs,
        int rows, int columns, const Options<T>& options = Options<T>()) {
    try {
        return solve_impl(overlap, interior_response, exterior_ratio, rhs, rows, columns, options);
    } catch (const std::runtime_error&) {
        Result<T> result;
        result.status = "factorized_operator_range";
        return result;
    } catch (const std::bad_alloc&) {
        Result<T> result;
        result.status = "workspace_allocation_failure";
        return result;
    } catch (const std::length_error&) {
        Result<T> result;
        result.status = "workspace_dimension_range";
        return result;
    }
}
}

namespace psms_direct_refinement {
template<typename T> struct Result {
    bool accepted = false;
    T initial_backward_error = 0, backward_error = 0;
    bool forward_checked = false;
    T forward_error = 0;
    int correction_count = 0;
    std::string status = "not_checked";
};

template<typename T> void add(T& sum, T& correction, const T& value) {
    const T adjusted = value - correction;
    const T updated = sum + adjusted;
    correction = (updated - sum) - adjusted;
    sum = updated;
}

class ForwardArithmeticError : public std::runtime_error {
public:
    ForwardArithmeticError() : std::runtime_error("Nonfinite native forward-error estimate.") {}
};

template<typename T> void require_finite(const std::complex<T>& value) {
    if (!is_finite_real(value.real()) || !is_finite_real(value.imag()))
        throw ForwardArithmeticError();
}

template<typename T> T abs1(const std::complex<T>& value) {
    return precabs(value.real()) + precabs(value.imag());
}

template<typename T> bool finite(const std::complex<T>& value) {
    return is_finite_real(value.real()) && is_finite_real(value.imag());
}

template<typename T> std::complex<T> residual_row(const std::vector<std::complex<T>>& matrix,
        const std::vector<std::complex<T>>& rhs, const std::vector<std::complex<T>>& solution,
        int row, int dimension) {
#if ACOUSTICTS_HAVE_EXACT_QUAD_DOT
    if constexpr (std::is_same_v<T,acousticts_quad_t>) {
        if (std::fegetround() == FE_TONEAREST) {
            static thread_local quad_dot_product::Accumulator real, imaginary;
            real.reset(); imaginary.reset();
            bool valid = real.add(rhs[row].real(),T(1)) && imaginary.add(rhs[row].imag(),T(1));
            for (int column=0; column<dimension; ++column) {
                const auto a=matrix[std::size_t(row)*dimension+column], x=solution[column];
                valid = real.add(-a.real(),x.real()) && valid;
                valid = real.add(a.imag(),x.imag()) && valid;
                valid = imaginary.add(-a.real(),x.imag()) && valid;
                valid = imaginary.add(-a.imag(),x.real()) && valid;
            }
            return valid ? std::complex<T>(real.rounded(),imaginary.rounded()) :
                std::complex<T>(precnan<T>(),precnan<T>());
        }
    }
#endif
    std::complex<T> residual=rhs[row], correction(0,0);
    for (int column=0; column<dimension; ++column)
        add(residual,correction,-matrix[std::size_t(row)*dimension+column]*solution[column]);
    return residual;
}

template<typename T> bool measure(const std::vector<std::complex<T>>& matrix,
        const std::vector<std::complex<T>>& rhs, const std::vector<std::complex<T>>& solution,
        int dimension, T& backward_error, std::vector<std::complex<T>>& residual) {
    backward_error=0;
    bool valid=dimension>=0 && matrix.size()==std::size_t(dimension)*dimension &&
        rhs.size()==std::size_t(dimension) && solution.size()==std::size_t(dimension);
    residual.resize(valid ? dimension : 0);
    std::vector<T> solution_abs(valid ? dimension : 0);
    for (int column=0; valid && column<dimension; ++column) {
        solution_abs[column]=abs1(solution[column]);
        valid=finite(solution[column]) && is_finite_real(solution_abs[column]);
    }
    // Use the componentwise abs1 denominator from LAPACK ZGERFS.
    for (int row=0; valid && row<dimension; ++row) {
        T denominator=abs1(rhs[row]), correction=0;
        valid=finite(rhs[row]) && is_finite_real(denominator);
        for (int column=0; valid && column<dimension; ++column) {
            const auto value=matrix[std::size_t(row)*dimension+column];
            const T term=abs1(value)*solution_abs[column];
            valid=finite(value) && is_finite_real(term);
            if (valid) add(denominator,correction,term);
        }
        if (!valid || !is_finite_real(denominator)) { valid=false; break; }
        residual[row]=residual_row(matrix,rhs,solution,row,dimension);
        const T residual_abs=abs1(residual[row]);
        valid=finite(residual[row]) && is_finite_real(residual_abs);
        if (!valid || (denominator==T(0) && residual_abs!=T(0))) { valid=false; break; }
        const T row_error=denominator>T(0) ? residual_abs/denominator : T(0);
        valid=is_finite_real(row_error);
        if (valid) backward_error=std::max(backward_error,row_error);
    }
    return valid;
}

template<typename T> T modulus(const std::complex<T>& value) {
    require_finite(value);
    const T scale=std::max(precabs(value.real()),precabs(value.imag()));
    if (scale==T(0)) return T(0);
    const T real=value.real()/scale,imag=value.imag()/scale;
    return scale*precsqrt(real*real+imag*imag);
}
template<typename T> T vector_one_norm(const std::vector<std::complex<T>>& values) {
    T sum=0,correction=0;
    for (const auto& value:values) add(sum,correction,modulus(value));
    require_finite(std::complex<T>(sum,0));
    return sum;
}
template<typename T> std::vector<std::complex<T>> adjoint_solve(
        const LUPFactorization<T>& factors,const std::vector<std::complex<T>>& rhs) {
    const int dimension=int(rhs.size());
    if (factors.singular || factors.row_piv.size()!=rhs.size() ||
            factors.col_piv.size()!=rhs.size() || factors.lu.size()!=rhs.size()*rhs.size())
        throw std::invalid_argument("Invalid complete-pivot adjoint factors.");
    std::vector<std::complex<T>> transformed(dimension),solution(dimension);
    // P*A*Q=L*U gives A^H=Q*U^H*L^H*P for the stored permutations.
    for (int row=0; row<dimension; ++row) {
        auto value=rhs[factors.col_piv[row]];
        for (int column=0; column<row; ++column)
            value-=std::conj(factors.lu[size_t(column)*dimension+row])*transformed[column];
        transformed[row]=value/std::conj(factors.lu[size_t(row)*dimension+row]);
        require_finite(transformed[row]);
    }
    for (int row=dimension-1; row>=0; --row) {
        for (int column=row+1; column<dimension; ++column)
            transformed[row]-=std::conj(factors.lu[size_t(column)*dimension+row])*transformed[column];
        require_finite(transformed[row]);
    }
    for (int row=0; row<dimension; ++row) solution[factors.row_piv[row]]=transformed[row];
    return solution;
}
template<typename T> struct ForwardEstimate {
    T weighted_inverse_norm=0,relative_error=0;
};
template<typename T> ForwardEstimate<T> forward_error_estimate(
        const LUPFactorization<T>& factors,const std::vector<T>& weights,T solution_norm) {
    using Complex=std::complex<T>;
    const int dimension=int(weights.size());
    if (dimension<=0 || factors.singular || factors.row_piv.size()!=weights.size())
        throw std::invalid_argument("Invalid forward-error estimate system.");
    require_finite(std::complex<T>(solution_norm, 0));
    if (solution_norm < T(0)) throw std::invalid_argument("Negative solution norm.");
    for (const T weight : weights) {
        require_finite(std::complex<T>(weight, 0));
        if (weight < T(0)) throw std::invalid_argument("Negative forward-error weight.");
    }
    ForwardEstimate<T> result;
    T safe_minimum=std::numeric_limits<T>::min();
#if ACOUSTICTS_HAVE_QUADMATH
    if constexpr (std::is_same_v<T,acousticts_quad_t>) safe_minimum=scalbnq(T(1),-16382);
#endif
    // Estimate the 1-norm of diag(W)*A^-H, equal to the required infinity norm.
    auto action=[&](const std::vector<Complex>& input) {
        auto output=adjoint_solve(factors,input);
        for (int index=0; index<dimension; ++index) {
            output[index]*=weights[index];
            require_finite(output[index]);
        }
        return output;
    };
    auto adjoint_action=[&](std::vector<Complex> input) {
        for (int index=0; index<dimension; ++index) input[index]*=weights[index];
        auto output=lup_solve(factors,input);
        for (const auto& value:output) require_finite(value);
        return output;
    };
    auto phases=[&](std::vector<Complex>& values) {
        for (auto& value:values) {
            const T magnitude=modulus(value);
            value=magnitude>safe_minimum ? value/magnitude : Complex(1,0);
        }
    };
    auto maximum_index=[&](const std::vector<Complex>& values) {
        int selected=0;
        T maximum=modulus(values[0]);
        for (int index=1; index<dimension; ++index) {
            const T magnitude=modulus(values[index]);
            if (magnitude>maximum) { maximum=magnitude; selected=index; }
        }
        return selected;
    };
    // Retain ZLACN2's phase search, cycling stop and final alternating probe.
    auto values=action(std::vector<Complex>(dimension,Complex(T(1)/T(dimension),0)));
    result.weighted_inverse_norm=vector_one_norm(values);
    if (dimension>1) {
        phases(values);
        values=adjoint_action(values);
        int selected=maximum_index(values);
        for (int iteration=2; iteration<=5; ++iteration) {
            std::fill(values.begin(),values.end(),Complex(0,0));
            values[selected]=Complex(1,0);
            values=action(values);
            const T previous=result.weighted_inverse_norm;
            result.weighted_inverse_norm=vector_one_norm(values);
            if (result.weighted_inverse_norm<=previous) break;
            phases(values);
            values=adjoint_action(values);
            const int previous_index=selected;
            selected=maximum_index(values);
            if (modulus(values[previous_index])==modulus(values[selected])) break;
        }
        for (int index=0; index<dimension; ++index)
            values[index]=Complex((index%2 ? T(-1) : T(1))*(T(1)+T(index)/T(dimension-1)),0);
        values=action(values);
        result.weighted_inverse_norm=std::max(result.weighted_inverse_norm,
            T(2)*(vector_one_norm(values)/T(3*dimension)));
    }
    result.relative_error=solution_norm>T(0) ? result.weighted_inverse_norm/solution_norm : result.weighted_inverse_norm;
    require_finite(Complex(result.weighted_inverse_norm,result.relative_error));
    return result;
}

template<typename T> Result<T> refine(const std::vector<std::complex<T>>& matrix,
        const std::vector<std::complex<T>>& rhs, std::vector<std::complex<T>>& solution,
        const LUPFactorization<T>& factors, int dimension, T forward_tolerance = T(0)) {
    Result<T> result;
    if (dimension < 0 || matrix.size() != std::size_t(dimension) * dimension ||
            rhs.size() != std::size_t(dimension) || solution.size() != std::size_t(dimension) ||
            factors.singular || factors.row_piv.size() != std::size_t(dimension) ||
            factors.col_piv.size() != std::size_t(dimension) ||
            factors.lu.size() != std::size_t(dimension) * dimension)
        throw std::invalid_argument("Invalid checked direct system or factorization.");
    if (!is_finite_real(forward_tolerance) || forward_tolerance < T(0))
        throw std::invalid_argument("Invalid forward-error tolerance.");
    std::vector<std::complex<T>> residual;
    bool valid = measure(matrix, rhs, solution, dimension, result.backward_error, residual);
    result.initial_backward_error = valid ? result.backward_error : precnan<T>();
    T previous_backward_error = T(3);
    // ZGERFS permits five corrections while the backward error decreases by a factor of two.
    while (valid && result.backward_error > T(64) * preceps<T>() &&
            T(2) * result.backward_error <= previous_backward_error && result.correction_count < 5) {
        previous_backward_error = result.backward_error;
        const auto correction = lup_solve<T>(factors, residual);
        for (int index = 0; index < dimension; ++index) solution[index] += correction[index];
        ++result.correction_count;
        valid = measure(matrix, rhs, solution, dimension, result.backward_error, residual);
    }
    result.accepted = valid && result.backward_error <= T(64) * preceps<T>();
    result.status = !valid ? "nonfinite" :
        (result.accepted ? (result.correction_count == 0 ? "accepted_initial" : "accepted_refined") :
        (result.correction_count >= 5 ? "correction_limit" : "insufficient_improvement"));
    if (!valid) result.backward_error = precnan<T>();
    if (result.accepted && forward_tolerance > T(0) && dimension > 0) {
        result.forward_checked = true;
        try {
            std::vector<T> weights(dimension);
            T solution_norm = 0, safe_minimum = std::numeric_limits<T>::min();
#if ACOUSTICTS_HAVE_QUADMATH
            if constexpr (std::is_same_v<T, acousticts_quad_t>) safe_minimum = scalbnq(T(1), -16382);
#endif
            const T safe1 = T(dimension + 1) * safe_minimum;
            const T safe2 = safe1 / preceps<T>();
            for (const auto& value : solution) solution_norm = std::max(solution_norm, abs1(value));
            // ZGERFS weights include residual and native rounding error in the equilibrated system.
            for (int row = 0; row < dimension; ++row) {
                T denominator = abs1(rhs[row]), correction = 0;
                for (int column = 0; column < dimension; ++column)
                    add(denominator, correction,
                        abs1(matrix[std::size_t(row) * dimension + column]) * abs1(solution[column]));
                weights[row] = abs1(residual[row]) + T(dimension + 1) * preceps<T>() * denominator +
                    (denominator <= safe2 ? safe1 : T(0));
            }
            result.forward_error = forward_error_estimate(factors, weights, solution_norm).relative_error;
            result.accepted = result.forward_error <= forward_tolerance;
            if (!result.accepted) result.status = "forward_error_budget";
        } catch (const ForwardArithmeticError&) {
            result.forward_error = precnan<T>();
            result.accepted = false;
            result.status = "forward_error_nonfinite";
        }
    }
    return result;
}
}

template<typename T>
std::vector<std::complex<T>> solve_linear_system_lup_refined(
    const std::vector<std::complex<T>>& A,
    const std::vector<std::complex<T>>& b,
    int n,
    int n_refine = 1
) {
    const auto factors = lup_decompose<T>(A, n);
    return solve_factored_system_lup_refined<T>(A, b, factors, n, n_refine);
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
std::vector<std::complex<T>> solve_matrix_system_lup_equilibrated(
    const std::vector<std::complex<T>>& matrix,
    const std::vector<std::complex<T>>& rhs,
    int size
) {
    if (size <= 0 || matrix.size() != size_t(size) * size || rhs.size() != matrix.size())
        throw std::invalid_argument("Fluid T-matrix system dimensions do not agree.");
    const auto require_finite = [](const std::complex<T>& value) {
        if (!is_finite_real(value.real()) || !is_finite_real(value.imag()))
            throw std::runtime_error("Fluid T-matrix system contains a nonfinite value.");
    };
    std::vector<std::complex<T>> balanced = matrix;
    std::vector<T> column_scales(size), row_scales(size);

    // Balance outgoing columns with component maxima that cannot overflow by squaring.
    for (int col = 0; col < size; ++col) {
        T scale = T(0);
        for (int row = 0; row < size; ++row) {
            const auto& value = balanced[size_t(row) * size + col];
            require_finite(value);
            scale = std::max({scale, precabs(value.real()), precabs(value.imag())});
        }
        if (!(scale > T(0)))
            throw std::runtime_error("Fluid kernel system is singular or ill-conditioned");
        column_scales[col] = scale;
        for (int row = 0; row < size; ++row) balanced[size_t(row) * size + col] /= scale;
    }

    // Scale every boundary equation identically across its matrix and incident right-hand sides.
    for (int row = 0; row < size; ++row) {
        T scale = T(0);
        for (int col = 0; col < size; ++col) {
            const auto& value = balanced[size_t(row) * size + col];
            scale = std::max({scale, precabs(value.real()), precabs(value.imag())});
        }
        if (!(scale > T(0)))
            throw std::runtime_error("Fluid kernel system is singular or ill-conditioned");
        row_scales[row] = scale;
        for (int col = 0; col < size; ++col) balanced[size_t(row) * size + col] /= scale;
    }
    const auto factors = lup_decompose<T>(balanced, size);
    if (factors.singular)
        throw std::runtime_error("Fluid kernel system is singular or ill-conditioned");

    // Reuse one native factorization and restore outgoing coefficient rows after each refined solve.
    std::vector<std::complex<T>> result(matrix.size()), rhs_column(size);
    for (int col = 0; col < size; ++col) {
        for (int row = 0; row < size; ++row) {
            rhs_column[row] = rhs[size_t(row) * size + col] / row_scales[row];
            require_finite(rhs_column[row]);
        }
        const auto solution = solve_factored_system_lup_refined<T>(
            balanced, rhs_column, factors, size, 2);
        for (int row = 0; row < size; ++row) {
            result[size_t(row) * size + col] = solution[row] / column_scales[row];
            require_finite(result[size_t(row) * size + col]);
        }
    }
    return result;
}

template<typename T>
std::vector<std::vector<std::complex<T>>> solve_fluid_t_blocks_native(
    const std::vector<std::vector<std::complex<T>>>& K1_kernel,
    const std::vector<std::vector<std::complex<T>>>& K3_kernel
) {
    int m_max = static_cast<int>(K3_kernel.size()) - 1;
    std::vector<std::vector<std::complex<T>>> result(m_max + 1);

    // Each retained order shares one balanced native factorization across all incident columns.
    for (int m = 0; m <= m_max; ++m) {
        int size = static_cast<int>(std::sqrt(K3_kernel[m].size()));
        result[m] = solve_matrix_system_lup_equilibrated<T>(K3_kernel[m], K1_kernel[m], size);
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

namespace psms_embedded_support {

template<typename T> struct Result {
    std::array<std::vector<std::complex<T>>, 2> solutions;
    std::array<bool, 2> valid{{false, false}};
    std::array<T, 2> backward_error{{precnan<T>(), precnan<T>()}};
    std::array<int, 2> correction_count{{0, 0}};
};

template<typename T> Result<T> reduced_solutions(
        const std::vector<std::complex<T>>& matrix,
        const LUPFactorization<T>& factors, int dimension,
        const std::array<std::vector<std::complex<T>>, 2>& reduced_rhs) {
    using Complex = std::complex<T>;
    Result<T> result;
    if (dimension < 2 || factors.singular ||
            matrix.size() != std::size_t(dimension) * dimension ||
            factors.lu.size() != matrix.size() ||
            factors.row_piv.size() != std::size_t(dimension) ||
            factors.col_piv.size() != std::size_t(dimension)) return result;
    const T roundoff_guard = T(64) * preceps<T>();
    std::array<std::vector<Complex>, 2> inverse_columns;
    std::array<bool, 2> inverse_valid{{false, false}};

    // Reuse the existing factorization for the last two inverse columns.
    for (int column = 0; column < 2 && column + 1 < dimension; ++column) {
        try {
            std::vector<Complex> unit_rhs(dimension, Complex(0, 0));
            unit_rhs[dimension - column - 1] = Complex(1, 0);
            inverse_columns[column] = solve_factored_system_lup_refined<T>(
                matrix, unit_rhs, factors, dimension, 2);
            inverse_valid[column] = true;
            for (const auto& value : inverse_columns[column])
                inverse_valid[column] = inverse_valid[column] && psms_direct_refinement::finite(value);
        } catch (const std::exception&) {
            inverse_valid[column] = false;
        }
    }

    for (int comparison = 0; comparison < 2; ++comparison) {
        const int removed = comparison + 1, retained = dimension - removed;
        if (retained <= 0 || !inverse_valid[0] ||
                (removed == 2 && !inverse_valid[1]) ||
                reduced_rhs[comparison].size() != std::size_t(dimension)) continue;
        bool rhs_valid = true;
        for (int row = 0; row < dimension; ++row) {
            const auto& value = reduced_rhs[comparison][row];
            rhs_valid = rhs_valid && psms_direct_refinement::finite(value) &&
                (row < retained || value == Complex(0, 0));
        }
        if (!rhs_valid) continue;
        try {
            const auto forced = solve_factored_system_lup_refined<T>(
                matrix, reduced_rhs[comparison], factors, dimension, 2);
            std::vector<Complex> small_matrix(std::size_t(removed) * removed), small_rhs(removed);
            T small_scale = 0;
            for (int row = 0; row < removed; ++row) {
                small_rhs[row] = forced[retained + row];
                for (int column = 0; column < removed; ++column) {
                    const auto value = inverse_columns[removed - column - 1][retained + row];
                    small_matrix[std::size_t(row) * removed + column] = value;
                    small_scale = std::max(small_scale, psms_direct_refinement::abs1(value));
                }
            }
            if (!(small_scale > T(0)) || !is_finite_real(small_scale)) continue;

            // Normalize the small system before applying the existing complete-pivot guard.
            for (auto& value : small_matrix) value /= small_scale;
            for (auto& value : small_rhs) value /= small_scale;
            bool small_valid = true;
            for (const auto& value : small_matrix)
                small_valid = small_valid && psms_direct_refinement::finite(value);
            for (const auto& value : small_rhs)
                small_valid = small_valid && psms_direct_refinement::finite(value);
            if (!small_valid) continue;
            const auto small_factors = lup_decompose<T>(small_matrix, removed);
            if (small_factors.singular) continue;
            const auto weights = solve_factored_system_lup_refined<T>(
                small_matrix, small_rhs, small_factors, removed, 2);

            // Apply the block inverse identity: https://nhigham.com/2023/06/01/what-is-the-schur-complement-of-a-matrix/.
            std::vector<Complex> reduced(dimension, Complex(0, 0));
            bool finite_solution = true;
            for (int row = 0; row < retained; ++row) {
                reduced[row] = forced[row];
                Complex correction(0, 0);
                for (int column = 0; column < removed; ++column)
                    psms_direct_refinement::add(reduced[row], correction,
                        -inverse_columns[removed - column - 1][row] * weights[column]);
                finite_solution = finite_solution && psms_direct_refinement::finite(reduced[row]);
            }
            if (!finite_solution) continue;

            // Check the retained equations with the existing componentwise arithmetic guard.
            T backward_error = 0;
            std::vector<Complex> residual_rhs(dimension, Complex(0, 0));
            const auto measure_reduced_residual = [&]() {
                backward_error = 0;
                for (int row = 0; row < retained; ++row) {
                    residual_rhs[row] = psms_direct_refinement::residual_row(
                        matrix, reduced_rhs[comparison], reduced, row, dimension);
                    const T residual_magnitude = psms_direct_refinement::abs1(residual_rhs[row]);
                    T denominator = psms_direct_refinement::abs1(reduced_rhs[comparison][row]);
                    T correction = 0;
                    for (int column = 0; column < retained; ++column) {
                        const T contribution = psms_direct_refinement::abs1(
                            matrix[std::size_t(row) * dimension + column]) *
                            psms_direct_refinement::abs1(reduced[column]);
                        if (!is_finite_real(contribution)) return false;
                        psms_direct_refinement::add(denominator, correction, contribution);
                    }
                    if (!psms_direct_refinement::finite(residual_rhs[row]) ||
                            !is_finite_real(residual_magnitude) || !is_finite_real(denominator) ||
                            (denominator == T(0) && residual_magnitude != T(0))) return false;
                    const T row_error = denominator > T(0) ? residual_magnitude / denominator : T(0);
                    if (!is_finite_real(row_error)) return false;
                    backward_error = std::max(backward_error, row_error);
                }
                return true;
            };
            bool residual_valid = measure_reduced_residual();

            // Refine rejected reconstructions through the same block inverse without refactorization.
            while (residual_valid && backward_error > roundoff_guard &&
                    result.correction_count[comparison] < 2) {
                const T previous_error = backward_error;
                const auto forced_correction = solve_factored_system_lup_refined<T>(
                    matrix, residual_rhs, factors, dimension, 2);
                std::vector<Complex> correction_rhs(removed);
                for (int row = 0; row < removed; ++row)
                    correction_rhs[row] = forced_correction[retained + row] / small_scale;
                const auto correction_weights = solve_factored_system_lup_refined<T>(
                    small_matrix, correction_rhs, small_factors, removed, 2);
                for (int row = 0; row < retained; ++row) {
                    Complex change = forced_correction[row], correction(0, 0);
                    for (int column = 0; column < removed; ++column)
                        psms_direct_refinement::add(change, correction,
                            -inverse_columns[removed - column - 1][row] * correction_weights[column]);
                    reduced[row] += change;
                    residual_valid = residual_valid && psms_direct_refinement::finite(reduced[row]);
                }
                ++result.correction_count[comparison];
                if (!residual_valid) break;
                residual_valid = measure_reduced_residual();
                if (residual_valid && backward_error >= previous_error) break;
            }
            result.backward_error[comparison] = residual_valid ? backward_error : precnan<T>();
            if (!residual_valid || backward_error > roundoff_guard) continue;
            reduced.resize(retained);
            result.solutions[comparison] = std::move(reduced);
            result.valid[comparison] = true;
        } catch (const std::exception&) {
            result.valid[comparison] = false;
        }
    }
    return result;
}

}
