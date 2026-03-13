#include <Rcpp.h>
#include <cmath>
#include <complex>
#include <vector>
#include "bessel_helpers.h" 
using namespace Rcpp;

// Constants
static const double pi = M_PI;
static const double tol = 1E-300;
static const std::complex<double> i_unit(0.0, 1.0);

// =============================================================================
// Helper: Compute determinant of a 4x4 complex matrix
// =============================================================================
std::complex<double> det4x4(std::complex<double> A[4][4]) {
    std::complex<double> det(0.0, 0.0);
    
    for (int i = 0; i < 4; ++i) {
        // 3x3 minor
        std::complex<double> minor[3][3];
        for (int j = 1; j < 4; ++j) {
            int col = 0;
            for (int k = 0; k < 4; ++k) {
                if (k != i) {
                    minor[j-1][col++] = A[j][k];
                }
            }
        }
        
        std::complex<double> det3 = 
            minor[0][0] * (minor[1][1] * minor[2][2] - minor[1][2] * minor[2][1]) -
            minor[0][1] * (minor[1][0] * minor[2][2] - minor[1][2] * minor[2][0]) +
            minor[0][2] * (minor[1][0] * minor[2][1] - minor[1][1] * minor[2][0]);
        
        double sign = (i % 2 == 0) ? 1.0 : -1.0;
        det += sign * A[0][i] * det3;
    }
    
    return det;
}

// =============================================================================
// Helper: Compute determinant of a 6x6 complex matrix using LU decomposition
// =============================================================================
std::complex<double> det6x6(std::complex<double> A[6][6]) {
    std::complex<double> LU[6][6];
    for (int i = 0; i < 6; ++i) {
        for (int j = 0; j < 6; ++j) {
            LU[i][j] = A[i][j];
        }
    }
    
    std::complex<double> det(1.0, 0.0);
    int sign = 1;
    
    for (int k = 0; k < 6; ++k) {
        // Find pivot
        int pivot = k;
        double max_val = std::abs(LU[k][k]);
        for (int i = k + 1; i < 6; ++i) {
            if (std::abs(LU[i][k]) > max_val) {
                max_val = std::abs(LU[i][k]);
                pivot = i;
            }
        }
        
        if (pivot != k) {
            for (int j = 0; j < 6; ++j) {
                std::swap(LU[k][j], LU[pivot][j]);
            }
            sign = -sign;
        }
        
        if (std::abs(LU[k][k]) < tol) {
            return std::complex<double>(0.0, 0.0);
        }
        
        det *= LU[k][k];
        
        for (int i = k + 1; i < 6; ++i) {
            std::complex<double> factor = LU[i][k] / LU[k][k];
            for (int j = k + 1; j < 6; ++j) {
                LU[i][j] -= factor * LU[k][j];
            }
        }
    }
    
    return det * std::complex<double>(sign, 0.0);
}

std::complex<double> det6x6_scaled(std::complex<double> A[6][6]) {
    std::complex<double> LU[6][6];
    
    // Scale each row by its max absolute value
    double row_scales[6];
    for (int i = 0; i < 6; ++i) {
        double max_val = 0.0;
        for (int j = 0; j < 6; ++j) {
            if (std::abs(A[i][j]) > max_val) {
                max_val = std::abs(A[i][j]);
            }
        }
        row_scales[i] = (max_val > 0.0) ? max_val : 1.0;
        for (int j = 0; j < 6; ++j) {
            LU[i][j] = A[i][j] / row_scales[i];
        }
    }
    
    std::complex<double> det(1.0, 0.0);
    int sign = 1;
    
    for (int k = 0; k < 6; ++k) {
        int pivot = k;
        double max_val = std::abs(LU[k][k]);
        for (int i = k + 1; i < 6; ++i) {
            if (std::abs(LU[i][k]) > max_val) {
                max_val = std::abs(LU[i][k]);
                pivot = i;
            }
        }
        
        if (pivot != k) {
            for (int j = 0; j < 6; ++j) {
                std::swap(LU[k][j], LU[pivot][j]);
            }
            std::swap(row_scales[k], row_scales[pivot]);
            sign = -sign;
        }
        
        if (std::abs(LU[k][k]) == 0.0) {
            return std::complex<double>(0.0, 0.0);
        }
        
        det *= LU[k][k];
        
        for (int i = k + 1; i < 6; ++i) {
            std::complex<double> factor = LU[i][k] / LU[k][k];
            for (int j = k + 1; j < 6; ++j) {
                LU[i][j] -= factor * LU[k][j];
            }
        }
    }

    // Multiply back the row scales
    for (int i = 0; i < 6; ++i) {
        det *= row_scales[i];
    }
    
    return det * std::complex<double>(sign, 0.0);
}

// =============================================================================
// Compute b_m for a single modal order and frequency
// =============================================================================
std::complex<double> compute_single_bm(
    int m,
    double k1a, double kLa_s, double kTa_s,
    double kLa_f, double kTa_f, double k3a,
    double lambda, double mu,
    double rho_ratio_sw, double rho_ratio_fl
) {
    double m_dbl = static_cast<double>(m);
    double m_m_plus_1 = m_dbl * (m_dbl + 1.0);
    double lambda_plus_2mu = lambda + 2.0 * mu;
    
    // Compute spherical Bessel functions using the single-value helpers
    double js_k1a = js_single_impl(m, k1a);
    double jsd_k1a = js_deriv_single_impl(m, k1a, 1);
    std::complex<double> hs_k1a = hs_single_impl(m, k1a);
    std::complex<double> hsd_k1a = hs_deriv_single_impl(m, k1a, 1);
    
    double js_kLa_s = js_single_impl(m, kLa_s);
    double ys_kLa_s = ys_single_impl(m, kLa_s);
    double jsd_kLa_s = js_deriv_single_impl(m, kLa_s, 1);
    double ysd_kLa_s = ys_deriv_single_impl(m, kLa_s, 1);
    double jsdd_kLa_s = js_deriv_single_impl(m, kLa_s, 2);
    double ysdd_kLa_s = ys_deriv_single_impl(m, kLa_s, 2);
    
    double js_kTa_s = js_single_impl(m, kTa_s);
    double ys_kTa_s = ys_single_impl(m, kTa_s);
    double jsd_kTa_s = js_deriv_single_impl(m, kTa_s, 1);
    double ysd_kTa_s = ys_deriv_single_impl(m, kTa_s, 1);
    double jsdd_kTa_s = js_deriv_single_impl(m, kTa_s, 2);
    double ysdd_kTa_s = ys_deriv_single_impl(m, kTa_s, 2);
    
    double js_kLa_f = js_single_impl(m, kLa_f);
    double ys_kLa_f = ys_single_impl(m, kLa_f);
    double jsd_kLa_f = js_deriv_single_impl(m, kLa_f, 1);
    double ysd_kLa_f = ys_deriv_single_impl(m, kLa_f, 1);
    double jsdd_kLa_f = js_deriv_single_impl(m, kLa_f, 2);
    double ysdd_kLa_f = ys_deriv_single_impl(m, kLa_f, 2);
    
    double js_kTa_f = js_single_impl(m, kTa_f);
    double ys_kTa_f = ys_single_impl(m, kTa_f);
    double jsd_kTa_f = js_deriv_single_impl(m, kTa_f, 1);
    double ysd_kTa_f = ys_deriv_single_impl(m, kTa_f, 1);
    double jsdd_kTa_f = js_deriv_single_impl(m, kTa_f, 2);
    double ysdd_kTa_f = ys_deriv_single_impl(m, kTa_f, 2);
    
    double js_k3a = js_single_impl(m, k3a);
    double jsd_k3a = js_deriv_single_impl(m, k3a, 1);
    
    // Compute alpha coefficients
    double a1 = js_k1a * rho_ratio_sw;
    double a2 = k1a * jsd_k1a;
    std::complex<double> a11 = hs_k1a * rho_ratio_sw;
    std::complex<double> a21 = k1a * hsd_k1a;
    
    double a12 = (lambda * js_kLa_s - 2.0 * mu * jsdd_kLa_s) / lambda_plus_2mu;
    double a22 = kLa_s * jsd_kLa_s;
    double a32 = 2.0 * (kLa_s * jsd_kLa_s - js_kLa_s);
    double a42 = (lambda * js_kLa_f - 2.0 * mu * jsdd_kLa_f) / lambda_plus_2mu;
    double a52 = kLa_f * jsd_kLa_f;
    double a62 = 2.0 * (kLa_f * jsd_kLa_f - js_kLa_f);
    
    double kTa_s_sq = kTa_s * kTa_s;
    double kTa_f_sq = kTa_f * kTa_f;
    
    double a13 = -2.0 * m_m_plus_1 / kTa_s_sq * (kTa_s * jsd_kTa_s - js_kTa_s);
    double a23 = m_m_plus_1 * js_kTa_s;
    double a33 = kTa_s_sq * jsdd_kTa_s + (m_dbl + 2.0) * (m_dbl - 1.0) * js_kTa_s;
    double a43 = -2.0 * m_m_plus_1 / kTa_f_sq * (kTa_f * jsd_kTa_f - js_kTa_f);
    double a53 = m_m_plus_1 * js_kTa_f;
    double a63 = kTa_f_sq * jsdd_kTa_f + (m_dbl + 2.0) * (m_dbl - 1.0) * js_kTa_f;
    
    double a14 = (lambda * ys_kLa_s - 2.0 * mu * ysdd_kLa_s) / lambda_plus_2mu;
    double a24 = kLa_s * ysd_kLa_s;
    double a34 = 2.0 * (kLa_s * ysd_kLa_s - ys_kLa_s);
    double a44 = (lambda * ys_kLa_f - 2.0 * mu * ysdd_kLa_f) / lambda_plus_2mu;
    double a54 = kLa_f * ysd_kLa_f;
    double a64 = 2.0 * (kLa_f * ysd_kLa_f - ys_kLa_f);
    
    double a15 = -2.0 * m_m_plus_1 / kTa_s_sq * (kTa_s * ysd_kTa_s - ys_kTa_s);
    double a25 = m_m_plus_1 * ys_kTa_s;
    double a35 = kTa_s_sq * ysdd_kTa_s + (m_dbl + 2.0) * (m_dbl - 1.0) * ys_kTa_s;
    double a45 = -2.0 * m_m_plus_1 / kTa_f_sq * (kTa_f * ysd_kTa_f - ys_kTa_f);
    double a55 = m_m_plus_1 * ys_kTa_f;
    double a65 = kTa_f_sq * ysdd_kTa_f + (m_dbl + 2.0) * (m_dbl - 1.0) * ys_kTa_f;
    
    double a46 = js_k3a * rho_ratio_fl;
    double a56 = k3a * jsd_k3a;
    
    std::complex<double> b_m;
    
    if (m == 0) {
        // 4x4 matrix for m = 0
        std::complex<double> A_num[4][4] = {
            {a1,  a12, a14, 0.0},
            {a2,  a22, a24, 0.0},
            {0.0, a42, a44, a46},
            {0.0, a52, a54, a56}
        };
        
        std::complex<double> A_den[4][4];
        for (int i = 0; i < 4; ++i) {
            for (int j = 0; j < 4; ++j) {
                A_den[i][j] = A_num[i][j];
            }
        }
        A_den[0][0] = a11;
        A_den[1][0] = a21;
        
        std::complex<double> det_num = det4x4(A_num);
        std::complex<double> det_den = det4x4(A_den);
        
        if (std::abs(det_den) < tol) {
            b_m = std::complex<double>(0.0, 0.0);
        } else {
            b_m = det_num / det_den;
        }
    } else {
        // 6x6 matrix for m > 0
        std::complex<double> A_num[6][6] = {
            {a1,  a12, a13, a14, a15, 0.0},
            {a2,  a22, a23, a24, a25, 0.0},
            {0.0, a32, a33, a34, a35, 0.0},
            {0.0, a42, a43, a44, a45, a46},
            {0.0, a52, a53, a54, a55, a56},
            {0.0, a62, a63, a64, a65, 0.0}
        };
        
        std::complex<double> A_den[6][6];
        for (int i = 0; i < 6; ++i) {
            for (int j = 0; j < 6; ++j) {
                A_den[i][j] = A_num[i][j];
            }
        }
        A_den[0][0] = a11;
        A_den[1][0] = a21;
        
        std::complex<double> det_num = det6x6_scaled(A_num);
        std::complex<double> det_den = det6x6_scaled(A_den);
        
        if (std::abs(det_den) < tol) {
            b_m = std::complex<double>(0.0, 0.0);
        } else {
            b_m = det_num / det_den;
        }
    }
    
    return b_m;
}

// =============================================================================
// Main exported function: compute b_m for all frequencies and modal orders
// Returns a complex matrix (n_freq x n_m) of b_m values
// =============================================================================
// [[Rcpp::export]]
ComplexMatrix elastic_shell_boundary_conditions_old(
    NumericVector k1a,       // k_sw * radius_shell for each frequency
    NumericVector kLa_shell, // kL * radius_shell
    NumericVector kTa_shell, // kT * radius_shell
    NumericVector kLa_fluid, // kL * radius_fluid
    NumericVector kTa_fluid, // kT * radius_fluid
    NumericVector k3a_fluid, // k_fluid * radius_fluid
    IntegerVector m_vec,     // modal orders to compute
    double lambda,
    double mu,
    double rho_ratio_sw,     // density_sw / density_shell
    double rho_ratio_fl      // density_fluid / density_shell
) {
    int n_freq = k1a.size();
    int n_m = m_vec.size();
    
    ComplexMatrix result(n_freq, n_m);
    
    for (int f = 0; f < n_freq; ++f) {
        for (int mi = 0; mi < n_m; ++mi) {
            int m = m_vec[mi];
            std::complex<double> bm = compute_single_bm(
                m,
                k1a[f], kLa_shell[f], kTa_shell[f],
                kLa_fluid[f], kTa_fluid[f], k3a_fluid[f],
                lambda, mu, rho_ratio_sw, rho_ratio_fl
            );
            result(f, mi) = to_Rcomplex(bm);
        }
    }
    
    return result;
}
// [[Rcpp::export]]
ComplexMatrix elastic_shell_boundary_conditions(
    NumericMatrix ka_matrix, // rows: "k1a_shell", "kLa_shell", ..., columns: freq
    IntegerVector m_limit, // length n_freq, each entry is the max m for that freq
    double lambda,
    double mu,
    double rho_ratio_sw,     // density_sw / density_shell
    double rho_ratio_fl      // density_fluid / density_shell
) {
    int n_freq = ka_matrix.ncol();
    int m_max = Rcpp::max(m_limit);
    ComplexMatrix result(n_freq, m_max + 1);

    // Row indices for ka_matrix
    int idx_k1a_shell = 0;
    int idx_kLa_shell = 1;
    int idx_kTa_shell = 2;
    int idx_kLa_fluid = 5;
    int idx_kTa_fluid = 4;
    int idx_k3a_fluid = 6;

    for (int f = 0; f < n_freq; ++f) {
        int m_limit_f = m_limit[f];
        for (int m = 0; m <= m_max; ++m) {
            if (m <= m_limit_f) {
                std::complex<double> bm = compute_single_bm(
                    m,
                    ka_matrix(idx_k1a_shell, f),
                    ka_matrix(idx_kLa_shell, f),
                    ka_matrix(idx_kTa_shell, f),
                    ka_matrix(idx_kLa_fluid, f),
                    ka_matrix(idx_kTa_fluid, f),
                    ka_matrix(idx_k3a_fluid, f),
                    lambda, mu, rho_ratio_sw, rho_ratio_fl
                );
                result(f, m) = Rcpp::ComplexVector::get_na();
                result(f, m) = to_Rcomplex(bm);
            } else {
                result(f, m) = Rcpp::ComplexVector::get_na();
            }
        }
    }
    
    return result;
}