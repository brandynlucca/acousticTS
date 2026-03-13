#include <Rcpp.h>
#include <cmath>
#include <complex>
#include <vector>
#include "bessel_helpers.h" 
using namespace Rcpp;

// Define PI and Tolerance for consistency
const double pi = M_PI; 
const double tol = 1e-300; 
const std::complex<double> i_unit(0.0, 1.0);

// --------------------------------------------------------------------------
// Helper: Calculate i^n
// --------------------------------------------------------------------------
std::complex<double> i_power_n(int n) {
    int rem = n % 4;
    if (rem < 0) rem += 4;
    
    if (rem == 0) return std::complex<double>(1.0, 0.0);
    if (rem == 1) return std::complex<double>(0.0, 1.0);
    if (rem == 2) return std::complex<double>(-1.0, 0.0);
    return -i_unit;
}

// --------------------------------------------------------------------------
// Helper: Determine if we should return matrix
// --------------------------------------------------------------------------
bool should_return_matrix(int size1, int size2) {
    return (size1 > 1 && size2 > 1 && size1 != size2);
}

// --------------------------------------------------------------------------
// Helper: Compute single J_nu(z) value
// --------------------------------------------------------------------------
std::complex<double> jc_single_impl(std::complex<double> zi, double nui) {
    Rcpp::Function R_besselI("besselI");
    
    double zi_real = zi.real();
    double zi_imag = zi.imag();
    
    bool is_integer = std::abs(nui - std::round(nui)) < tol;
    int int_nu = (int)std::round(nui);
    
    // Handle J_nu(0)
    if (std::abs(zi_real) < tol && std::abs(zi_imag) < tol) {
        if (std::abs(nui) < tol) return std::complex<double>(1.0, 0.0);
        else return std::complex<double>(0.0, 0.0);
    }
    
    // Case A1: Purely Real Argument
    if (std::abs(zi_imag) < tol) {
        double xi = zi_real;
        
        if (xi > 0) {
            return std::complex<double>(R::bessel_j(xi, nui), 0.0);
        } else {
            double abs_xi = std::abs(xi);
            double J_abs_x = R::bessel_j(abs_xi, nui);
            
            if (is_integer) {
                double sign_factor = (int_nu % 2 == 0) ? 1.0 : -1.0;
                return std::complex<double>(sign_factor * J_abs_x, 0.0);
            } else {
                std::complex<double> cplx_factor = std::exp(std::complex<double>(0.0, pi * nui));
                return cplx_factor * J_abs_x;
            }
        }
    }
    
    // Case A2: Purely Imaginary Argument
    if (std::abs(zi_real) < tol) {
        double y_val = zi_imag;
        double abs_y = std::abs(y_val);
        
        double I_nu_abs_y = Rcpp::as<double>(R_besselI(abs_y, nui, Named("expon.scaled") = false));
        
        std::complex<double> exp_factor;
        if (is_integer) {
            exp_factor = i_power_n(int_nu);
        } else {
            exp_factor = std::exp((i_unit * pi * nui) / 2.0);
        }
        
        std::complex<double> final_cplx = exp_factor * I_nu_abs_y;
        
        if (y_val < 0) {
            std::complex<double> sign_factor = std::exp(std::complex<double>(0.0, pi * nui));
            final_cplx = sign_factor * final_cplx;
        }
        
        return final_cplx;
    }
    
    // Case A3: General Complex Argument
    Rcpp::stop("Computation failed: J_nu(z) with general complex z (%.2f + %.2fi) is unsupported.", zi_real, zi_imag);
    return std::complex<double>(0.0, 0.0);
}

// --------------------------------------------------------------------------
// Helper: Compute single Y_nu(z) value
// --------------------------------------------------------------------------
std::complex<double> yc_single_impl(std::complex<double> zi, double nui) {
    Rcpp::Function R_besselI("besselI");
    Rcpp::Function R_besselK("besselK");
    
    double zi_real = zi.real();
    double zi_imag = zi.imag();
    
    bool is_integer = std::abs(nui - std::round(nui)) < tol;
    
    // Handle z = 0
    if (std::abs(zi_real) < tol && std::abs(zi_imag) < tol) {
        return std::complex<double>(R_NegInf, 0.0);
    }
    
    // Case A1: Purely Real Argument
    if (std::abs(zi_imag) < tol) {
        double xi = zi_real;
        
        if (xi > 0) {
            return std::complex<double>(R::bessel_y(xi, nui), 0.0);
        } else {
            double abs_xi = std::abs(xi);
            double Y_nu_x = R::bessel_y(abs_xi, nui);
            double J_nu_x = R::bessel_j(abs_xi, nui);
            double cos_nu_pi = std::cos(pi * nui);
            double sin_nu_pi = std::sin(pi * nui);
            return std::complex<double>(cos_nu_pi * Y_nu_x + sin_nu_pi * J_nu_x, 0.0);
        }
    }
    
    // Case A2: Purely Imaginary Argument
    if (std::abs(zi_real) < tol) {
        double y_val = zi_imag;
        double abs_y = std::abs(y_val);
        
        double I_nu_abs_y = Rcpp::as<double>(R_besselI(abs_y, nui, Named("expon.scaled") = false));
        double K_nu_abs_y = Rcpp::as<double>(R_besselK(abs_y, nui, Named("expon.scaled") = false));
        
        std::complex<double> factor1_exp = std::exp(-i_unit * pi * nui / 2.0);
        std::complex<double> term1 = i_unit * factor1_exp * I_nu_abs_y;
        
        std::complex<double> factor2_exp = std::exp(i_unit * pi * nui / 2.0);
        std::complex<double> term2 = -(2.0 / pi) * factor2_exp * K_nu_abs_y;
        
        std::complex<double> final_cplx = term1 + term2;
        
        if (y_val < 0) {
            if (is_integer) {
                int int_nu = (int)std::round(nui);
                std::complex<double> sign_factor = std::complex<double>(std::cos(pi * int_nu), 0.0);
                final_cplx = sign_factor * final_cplx;
            } else {
                final_cplx = std::conj(term1) + std::conj(term2);
            }
        }
        
        return final_cplx;
    }
    
    // Case A3: General Complex Argument
    Rcpp::stop("Computation failed: Y_nu(z) with general complex z (%.2f + %.2fi) is unsupported.", zi_real, zi_imag);
    return std::complex<double>(0.0, 0.0);
}

// --------------------------------------------------------------------------
// Helper: Compute single spherical j_l(z) value
// --------------------------------------------------------------------------
double js_single_impl(int li, double zi) {
    if (std::abs(zi) < tol) {
        return (li == 0) ? 1.0 : 0.0;
    }
    double nu = li + 0.5;
    std::complex<double> J = jc_single_impl(std::complex<double>(zi, 0.0), nu);
    return std::sqrt(pi / (2.0 * zi)) * J.real();
}

// --------------------------------------------------------------------------
// Helper: Compute single spherical y_l(z) value
// --------------------------------------------------------------------------
double ys_single_impl(int li, double zi) {
    if (std::abs(zi) < tol) {
        return R_NegInf;
    }
    double nu = li + 0.5;
    std::complex<double> Y = yc_single_impl(std::complex<double>(std::abs(zi), 0.0), nu);
    return std::sqrt(pi / (2.0 * std::abs(zi))) * Y.real();
}

// --------------------------------------------------------------------------
// Helper: Compute single spherical h_l(z) value
// --------------------------------------------------------------------------
std::complex<double> hs_single_impl(int li, double zi) {
    if (std::abs(zi) < tol) {
        double re = (li == 0) ? 1.0 : 0.0;
        return std::complex<double>(re, R_NegInf);
    }
    double nu = li + 0.5;
    std::complex<double> J = jc_single_impl(std::complex<double>(zi, 0.0), nu);
    std::complex<double> Y = yc_single_impl(std::complex<double>(zi, 0.0), nu);
    return std::sqrt(pi / (2.0 * zi)) * (J + i_unit * Y);
}

// --------------------------------------------------------------------------
// Helper: Compute js derivative for single (li, zi)
// Spherical Bessel functions satisfy:
//   z^2 j''_l + 2z j'_l + [z^2 - l(l+1)] j_l = 0
// So: j''_l = -(2/z) j'_l + [l(l+1)/z^2 - 1] j_l
// 
// For first derivative: j'_l = j_{l-1} - (l+1)/z * j_l
// --------------------------------------------------------------------------
double js_deriv_single_impl(int li, double zi, int k) {
    if (k == 0) {
        return js_single_impl(li, zi);
    }
    
    if (std::abs(zi) < tol) {
        if (k == 1 && li == 1) return 1.0 / 3.0;
        if (k == 2 && li == 0) return -1.0 / 3.0;
        if (k == 2 && li == 2) return 2.0 / 15.0;
        return 0.0;
    }
    
    double z = zi;
    double z2 = z * z;
    double ll1 = (double)li * (li + 1);  // l(l+1)
    
    // Get j_l and j_{l-1} (or j_{l+1} for l=0)
    double j_l = js_single_impl(li, z);
    
    // First derivative: j'_l = j_{l-1} - (l+1)/z * j_l  (for l > 0)
    //                   j'_0 = -j_1
    double jp_l;
    if (li == 0) {
        jp_l = -js_single_impl(1, z);
    } else {
        double j_lm1 = js_single_impl(li - 1, z);
        jp_l = j_lm1 - (double)(li + 1) / z * j_l;
    }
    
    if (k == 1) {
        return jp_l;
    }
    
    // For k >= 2, use the differential equation iteratively
    // j''_l = -(2/z) j'_l + [l(l+1)/z^2 - 1] j_l
    // j'''_l = d/dz[j''_l] = ...
    
    // Store derivatives: D[i] = j^(i)_l(z)
    std::vector<double> D(k + 1);
    D[0] = j_l;
    D[1] = jp_l;

    // General pattern using Leibniz rule on j'' = a*j' + b*j:
    // j^(n) = sum over terms involving j^(i) for i < n
    // Compute iteratively using the recurrence from the ODE
    for (int n = 2; n <= k; ++n) {
        // j'' = -(2/z) j' + [l(l+1)/z^2 - 1] j
        
        if (n == 2) {
            D[2] = -(2.0/z) * D[1] + (ll1/z2 - 1.0) * D[0];
        } else {
            // Differentiate the ODE (n-2) times using Faà di Bruno / Leibniz
            // z^2 D[n] + 2(n-2) z D[n-1] + (n-2)(n-3) D[n-2]
            //   + 2z D[n-1] + 2(n-2) D[n-2] 
            //   + (z^2 - ll1) D[n-2] + 2(n-2) z D[n-3] + (n-2)(n-3) D[n-4] = 0
                        
            double sum = 0.0;
            double binom = 1.0;
            for (int i = 0; i <= n - 2; ++i) {
                if (i > 0) binom = binom * (n - 2 - i + 1) / i;
                
                // d^i/dz^i [-2/z] = (-2) * (-1)^i * i! / z^{i+1} = 2 * (-1)^{i+1} * i! / z^{i+1}
                double a_i = 2.0 * (i % 2 == 0 ? -1.0 : 1.0);
                double fact_i = 1.0;
                for (int f = 1; f <= i; ++f) fact_i *= f;
                a_i *= fact_i / std::pow(z, i + 1);
                
                // b^(i) = d^i/dz^i [ll1/z^2 - 1] = ll1 * d^i/dz^i[z^{-2}]
                // d^i/dz^i[z^{-2}] = (-2)(-3)...(-2-i+1) / z^{2+i} = (-1)^i * (i+1)! / z^{2+i}
                double b_i;
                if (i == 0) {
                    b_i = ll1 / z2 - 1.0;
                } else {
                    double prod = 1.0;
                    for (int j = 0; j < i; ++j) prod *= -(2.0 + j);
                    b_i = ll1 * prod / std::pow(z, 2 + i);
                }
                
                int idx_a = n - 1 - i;
                int idx_b = n - 2 - i;
                
                if (idx_a >= 0 && idx_a < (int)D.size()) {
                    sum += binom * a_i * D[idx_a];
                }
                if (idx_b >= 0 && idx_b < (int)D.size()) {
                    sum += binom * b_i * D[idx_b];
                }
            }
            D[n] = sum;
        }
    }
    
    return D[k];
}

// --------------------------------------------------------------------------
// Helper: Compute ys derivative for single (li, zi) using ODE method
// Spherical Bessel functions satisfy:
//   z^2 y''_l + 2z y'_l + [z^2 - l(l+1)] y_l = 0
// So: y''_l = -(2/z) y'_l + [l(l+1)/z^2 - 1] y_l
// 
// For first derivative: y'_l = y_{l-1} - (l+1)/z * y_l  (for l > 0)
//                       y'_0 = -y_1
// --------------------------------------------------------------------------
double ys_deriv_single_impl(int li, double zi, int k) {
    if (k == 0) {
        return ys_single_impl(li, zi);
    }
    
    if (std::abs(zi) < tol) {
        return R_NegInf;
    }
    
    double z = zi;
    double z2 = z * z;
    double ll1 = (double)li * (li + 1);  // l(l+1)
    
    // Get y_l
    double y_l = ys_single_impl(li, z);
    
    // First derivative: y'_l = y_{l-1} - (l+1)/z * y_l  (for l > 0)
    //                   y'_0 = -y_1
    double yp_l;
    if (li == 0) {
        yp_l = -ys_single_impl(1, z);
    } else {
        double y_lm1 = ys_single_impl(li - 1, z);
        yp_l = y_lm1 - (double)(li + 1) / z * y_l;
    }
    
    if (k == 1) {
        return yp_l;
    }
    
    // For k >= 2, use the differential equation iteratively
    // y''_l = -(2/z) y'_l + [l(l+1)/z^2 - 1] y_l
    
    std::vector<double> D(k + 1);
    D[0] = y_l;
    D[1] = yp_l;
    
    for (int n = 2; n <= k; ++n) {
        if (n == 2) {
            D[2] = -(2.0/z) * D[1] + (ll1/z2 - 1.0) * D[0];
        } else {
            // Use Leibniz rule on y'' = a*y' + b*y where a = -2/z, b = ll1/z^2 - 1
            double sum = 0.0;
            double binom = 1.0;
            for (int i = 0; i <= n - 2; ++i) {
                if (i > 0) binom = binom * (n - 2 - i + 1) / i;
                
                // a^(i) = d^i/dz^i [-2/z] = 2 * (-1)^{i+1} * i! / z^{i+1}
                double a_i = 2.0 * (i % 2 == 0 ? -1.0 : 1.0);
                double fact_i = 1.0;
                for (int f = 1; f <= i; ++f) fact_i *= f;
                a_i *= fact_i / std::pow(z, i + 1);
                
                // b^(i) = d^i/dz^i [ll1/z^2 - 1]
                double b_i;
                if (i == 0) {
                    b_i = ll1 / z2 - 1.0;
                } else {
                    double prod = 1.0;
                    for (int j = 0; j < i; ++j) prod *= -(2.0 + j);
                    b_i = ll1 * prod / std::pow(z, 2 + i);
                }
                
                int idx_a = n - 1 - i;
                int idx_b = n - 2 - i;
                
                if (idx_a >= 0 && idx_a < (int)D.size()) {
                    sum += binom * a_i * D[idx_a];
                }
                if (idx_b >= 0 && idx_b < (int)D.size()) {
                    sum += binom * b_i * D[idx_b];
                }
            }
            D[n] = sum;
        }
    }
    
    return D[k];
}

// --------------------------------------------------------------------------
// Helper: Compute hs derivative for single (li, zi) using recurrence
// h'_l(z) = h_{l-1}(z) - (l+1)/z * h_l(z)  for l > 0
// h'_0(z) = -h_1(z)
// --------------------------------------------------------------------------
std::complex<double> hs_deriv_single_impl(int li, double zi, int k) {
    if (k == 0) {
        return hs_single_impl(li, zi);
    }
    
    if (std::abs(zi) < tol) {
        return std::complex<double>(R_NaN, R_NaN);
    }
    
    double z = zi;
    double z2 = z * z;
    double ll1 = (double)li * (li + 1);  // l(l+1)
    
    // Get h_l
    std::complex<double> h_l = hs_single_impl(li, z);
    
    // First derivative: h'_l = h_{l-1} - (l+1)/z * h_l  (for l > 0)
    //                   h'_0 = -h_1
    std::complex<double> hp_l;
    if (li == 0) {
        hp_l = -hs_single_impl(1, z);
    } else {
        std::complex<double> h_lm1 = hs_single_impl(li - 1, z);
        hp_l = h_lm1 - (double)(li + 1) / z * h_l;
    }
    
    if (k == 1) {
        return hp_l;
    }
    
    // For k >= 2, use the differential equation iteratively
    // h''_l = -(2/z) h'_l + [l(l+1)/z^2 - 1] h_l
    
    std::vector<std::complex<double> > D(k + 1);
    D[0] = h_l;
    D[1] = hp_l;
    
    for (int n = 2; n <= k; ++n) {
        if (n == 2) {
            D[2] = -(2.0/z) * D[1] + (ll1/z2 - 1.0) * D[0];
        } else {
            // Use Leibniz rule on h'' = a*h' + b*h where a = -2/z, b = ll1/z^2 - 1
            std::complex<double> sum(0.0, 0.0);
            double binom = 1.0;
            for (int i = 0; i <= n - 2; ++i) {
                if (i > 0) binom = binom * (n - 2 - i + 1) / i;
                
                // a^(i) = d^i/dz^i [-2/z] = 2 * (-1)^{i+1} * i! / z^{i+1}
                double a_i = 2.0 * (i % 2 == 0 ? -1.0 : 1.0);
                double fact_i = 1.0;
                for (int f = 1; f <= i; ++f) fact_i *= f;
                a_i *= fact_i / std::pow(z, i + 1);
                
                // b^(i) = d^i/dz^i [ll1/z^2 - 1]
                double b_i;
                if (i == 0) {
                    b_i = ll1 / z2 - 1.0;
                } else {
                    double prod = 1.0;
                    for (int j = 0; j < i; ++j) prod *= -(2.0 + j);
                    b_i = ll1 * prod / std::pow(z, 2 + i);
                }
                
                int idx_a = n - 1 - i;
                int idx_b = n - 2 - i;
                
                if (idx_a >= 0 && idx_a < (int)D.size()) {
                    sum += binom * a_i * D[idx_a];
                }
                if (idx_b >= 0 && idx_b < (int)D.size()) {
                    sum += binom * b_i * D[idx_b];
                }
            }
            D[n] = sum;
        }
    }
    
    return D[k];
}

// ==========================================================================
// CYLINDRICAL BESSEL FUNCTIONS
// ==========================================================================

// --------------------------------------------------------------------------
// 1. Bessel Function of the First Kind (J_nu(z))
// --------------------------------------------------------------------------
// [[Rcpp::export]]
SEXP jc_cpp(ComplexVector z, ComplexVector nu) {
    int z_size = z.size();
    int nu_size = nu.size();
    
    // Check for complex order
    for (int i = 0; i < nu_size; ++i) {
        if (std::abs(nu[i].i) > tol) {
            Rcpp::stop("Computation failed: Complex order nu (%.2f + %.2fi) is unsupported.", nu[i].r, nu[i].i);
        }
    }
    
    if (should_return_matrix(nu_size, z_size)) {
        ComplexMatrix result(nu_size, z_size);
        
        for (int i = 0; i < nu_size; ++i) {
            double nui = nu[i].r;
            for (int j = 0; j < z_size; ++j) {
                std::complex<double> zi(z[j].r, z[j].i);
                std::complex<double> val = jc_single_impl(zi, nui);
                result(i, j) = to_Rcomplex(val);
            }
        }
        return result;
    } else {
        int N = std::max(z_size, nu_size);
        ComplexVector result(N);
        
        for (int i = 0; i < N; ++i) {
            int z_idx = i % z_size;
            int nu_idx = i % nu_size;
            std::complex<double> zi(z[z_idx].r, z[z_idx].i);
            double nui = nu[nu_idx].r;
            std::complex<double> val = jc_single_impl(zi, nui);
            result[i] = to_Rcomplex(val);
        }
        return result;
    }
}

// --------------------------------------------------------------------------
// 2. Bessel Function of the Second Kind (Y_nu(z))
// --------------------------------------------------------------------------
// [[Rcpp::export]]
SEXP yc_cpp(ComplexVector z, ComplexVector nu) {
    int z_size = z.size();
    int nu_size = nu.size();
    
    for (int i = 0; i < nu_size; ++i) {
        if (std::abs(nu[i].i) > tol) {
            Rcpp::stop("Computation failed: Complex order nu (%.2f + %.2fi) is unsupported.", nu[i].r, nu[i].i);
        }
    }
    
    if (should_return_matrix(nu_size, z_size)) {
        ComplexMatrix result(nu_size, z_size);
        
        for (int i = 0; i < nu_size; ++i) {
            double nui = nu[i].r;
            for (int j = 0; j < z_size; ++j) {
                std::complex<double> zi(z[j].r, z[j].i);
                std::complex<double> val = yc_single_impl(zi, nui);
                result(i, j) = to_Rcomplex(val);
            }
        }
        return result;
    } else {
        int N = std::max(z_size, nu_size);
        ComplexVector result(N);
        
        for (int i = 0; i < N; ++i) {
            int z_idx = i % z_size;
            int nu_idx = i % nu_size;
            std::complex<double> zi(z[z_idx].r, z[z_idx].i);
            double nui = nu[nu_idx].r;
            std::complex<double> val = yc_single_impl(zi, nui);
            result[i] = to_Rcomplex(val);
        }
        return result;
    }
}

// --------------------------------------------------------------------------
// 3. Hankel Function of the First Kind H^(1)_nu(z)
// --------------------------------------------------------------------------
// [[Rcpp::export]]
SEXP hc_cpp(ComplexVector z, ComplexVector nu) {
    int z_size = z.size();
    int nu_size = nu.size();
    
    for (int i = 0; i < nu_size; ++i) {
        if (std::abs(nu[i].i) > tol) {
            Rcpp::stop("Computation failed: Complex order nu (%.2f + %.2fi) is unsupported.", nu[i].r, nu[i].i);
        }
    }
    
    if (should_return_matrix(nu_size, z_size)) {
        ComplexMatrix result(nu_size, z_size);
        
        for (int i = 0; i < nu_size; ++i) {
            double nui = nu[i].r;
            for (int j = 0; j < z_size; ++j) {
                std::complex<double> zi(z[j].r, z[j].i);
                std::complex<double> J = jc_single_impl(zi, nui);
                std::complex<double> Y = yc_single_impl(zi, nui);
                std::complex<double> H = J + i_unit * Y;
                result(i, j) = to_Rcomplex(H);
            }
        }
        return result;
    } else {
        int N = std::max(z_size, nu_size);
        ComplexVector result(N);
        
        for (int i = 0; i < N; ++i) {
            int z_idx = i % z_size;
            int nu_idx = i % nu_size;
            std::complex<double> zi(z[z_idx].r, z[z_idx].i);
            double nui = nu[nu_idx].r;
            std::complex<double> J = jc_single_impl(zi, nui);
            std::complex<double> Y = yc_single_impl(zi, nui);
            std::complex<double> H = J + i_unit * Y;
            result[i] = to_Rcomplex(H);
        }
        return result;
    }
}

// --------------------------------------------------------------------------
// 4. Cylindrical Bessel J - k-th derivative
// --------------------------------------------------------------------------
// [[Rcpp::export]]
SEXP jc_deriv_cpp(ComplexVector z, ComplexVector nu, int k) {
    if (k < 0) {
        Rcpp::stop("Derivative order k must be non-negative");
    }
    if (k == 0) {
        return jc_cpp(z, nu);
    }
    
    int z_size = z.size();
    int nu_size = nu.size();
    
    // Precompute binomial coefficients
    std::vector<double> binom(k + 1);
    binom[0] = 1.0;
    for (int j = 1; j <= k; ++j) {
        binom[j] = binom[j - 1] * (k - j + 1) / j;
    }
    double scale = std::pow(0.5, k);
    
    if (should_return_matrix(nu_size, z_size)) {
        ComplexMatrix result(nu_size, z_size);
        
        for (int i = 0; i < nu_size; ++i) {
            double nui = nu[i].r;
            for (int j_idx = 0; j_idx < z_size; ++j_idx) {
                std::complex<double> zi(z[j_idx].r, z[j_idx].i);
                std::complex<double> sum(0.0, 0.0);
                
                for (int j = 0; j <= k; ++j) {
                    double nu_shifted = nui - k + 2 * j;
                    std::complex<double> J = jc_single_impl(zi, nu_shifted);
                    double sign = (j % 2 == 0) ? 1.0 : -1.0;
                    sum += sign * binom[j] * J;
                }
                
                std::complex<double> deriv = scale * sum;
                result(i, j_idx) = to_Rcomplex(deriv);
            }
        }
        return result;
    } else {
        int N = std::max(z_size, nu_size);
        ComplexVector result(N);
        
        for (int i = 0; i < N; ++i) {
            int z_idx = i % z_size;
            int nu_idx = i % nu_size;
            std::complex<double> zi(z[z_idx].r, z[z_idx].i);
            double nui = nu[nu_idx].r;
            
            std::complex<double> sum(0.0, 0.0);
            for (int j = 0; j <= k; ++j) {
                double nu_shifted = nui - k + 2 * j;
                std::complex<double> J = jc_single_impl(zi, nu_shifted);
                double sign = (j % 2 == 0) ? 1.0 : -1.0;
                sum += sign * binom[j] * J;
            }
            
            std::complex<double> deriv = scale * sum;
            result[i] = to_Rcomplex(deriv);
        }
        return result;
    }
}

// --------------------------------------------------------------------------
// 5. Cylindrical Bessel Y - k-th derivative
// --------------------------------------------------------------------------
// [[Rcpp::export]]
SEXP yc_deriv_cpp(ComplexVector z, ComplexVector nu, int k) {
    if (k < 0) {
        Rcpp::stop("Derivative order k must be non-negative");
    }
    if (k == 0) {
        return yc_cpp(z, nu);
    }
    
    int z_size = z.size();
    int nu_size = nu.size();
    
    std::vector<double> binom(k + 1);
    binom[0] = 1.0;
    for (int j = 1; j <= k; ++j) {
        binom[j] = binom[j - 1] * (k - j + 1) / j;
    }
    double scale = std::pow(0.5, k);
    
    if (should_return_matrix(nu_size, z_size)) {
        ComplexMatrix result(nu_size, z_size);
        
        for (int i = 0; i < nu_size; ++i) {
            double nui = nu[i].r;
            for (int j_idx = 0; j_idx < z_size; ++j_idx) {
                std::complex<double> zi(z[j_idx].r, z[j_idx].i);
                std::complex<double> sum(0.0, 0.0);
                
                for (int j = 0; j <= k; ++j) {
                    double nu_shifted = nui - k + 2 * j;
                    std::complex<double> Y = yc_single_impl(zi, nu_shifted);
                    double sign = (j % 2 == 0) ? 1.0 : -1.0;
                    sum += sign * binom[j] * Y;
                }
                
                std::complex<double> deriv = scale * sum;
                result(i, j_idx) = to_Rcomplex(deriv);
            }
        }
        return result;
    } else {
        int N = std::max(z_size, nu_size);
        ComplexVector result(N);
        
        for (int i = 0; i < N; ++i) {
            int z_idx = i % z_size;
            int nu_idx = i % nu_size;
            std::complex<double> zi(z[z_idx].r, z[z_idx].i);
            double nui = nu[nu_idx].r;
            
            std::complex<double> sum(0.0, 0.0);
            for (int j = 0; j <= k; ++j) {
                double nu_shifted = nui - k + 2 * j;
                std::complex<double> Y = yc_single_impl(zi, nu_shifted);
                double sign = (j % 2 == 0) ? 1.0 : -1.0;
                sum += sign * binom[j] * Y;
            }
            
            std::complex<double> deriv = scale * sum;
            result[i] = to_Rcomplex(deriv);
        }
        return result;
    }
}

// --------------------------------------------------------------------------
// 6. Cylindrical Hankel H^(1) - k-th derivative
// --------------------------------------------------------------------------
// [[Rcpp::export]]
SEXP hc_deriv_cpp(ComplexVector z, ComplexVector nu, int k) {
    if (k < 0) {
        Rcpp::stop("Derivative order k must be non-negative");
    }
    if (k == 0) {
        return hc_cpp(z, nu);
    }
    
    int z_size = z.size();
    int nu_size = nu.size();
    
    std::vector<double> binom(k + 1);
    binom[0] = 1.0;
    for (int j = 1; j <= k; ++j) {
        binom[j] = binom[j - 1] * (k - j + 1) / j;
    }
    double scale = std::pow(0.5, k);
    
    if (should_return_matrix(nu_size, z_size)) {
        ComplexMatrix result(nu_size, z_size);
        
        for (int i = 0; i < nu_size; ++i) {
            double nui = nu[i].r;
            for (int j_idx = 0; j_idx < z_size; ++j_idx) {
                std::complex<double> zi(z[j_idx].r, z[j_idx].i);
                std::complex<double> sum(0.0, 0.0);
                
                for (int j = 0; j <= k; ++j) {
                    double nu_shifted = nui - k + 2 * j;
                    std::complex<double> J = jc_single_impl(zi, nu_shifted);
                    std::complex<double> Y = yc_single_impl(zi, nu_shifted);
                    std::complex<double> H = J + i_unit * Y;
                    double sign = (j % 2 == 0) ? 1.0 : -1.0;
                    sum += sign * binom[j] * H;
                }
                
                std::complex<double> deriv = scale * sum;
                result(i, j_idx) = to_Rcomplex(deriv);
            }
        }
        return result;
    } else {
        int N = std::max(z_size, nu_size);
        ComplexVector result(N);
        
        for (int i = 0; i < N; ++i) {
            int z_idx = i % z_size;
            int nu_idx = i % nu_size;
            std::complex<double> zi(z[z_idx].r, z[z_idx].i);
            double nui = nu[nu_idx].r;
            
            std::complex<double> sum(0.0, 0.0);
            for (int j = 0; j <= k; ++j) {
                double nu_shifted = nui - k + 2 * j;
                std::complex<double> J = jc_single_impl(zi, nu_shifted);
                std::complex<double> Y = yc_single_impl(zi, nu_shifted);
                std::complex<double> H = J + i_unit * Y;
                double sign = (j % 2 == 0) ? 1.0 : -1.0;
                sum += sign * binom[j] * H;
            }
            
            std::complex<double> deriv = scale * sum;
            result[i] = to_Rcomplex(deriv);
        }
        return result;
    }
}

// ==========================================================================
// SPHERICAL BESSEL FUNCTIONS
// ==========================================================================

// --------------------------------------------------------------------------
// 7. Spherical Bessel j_l(z)
// --------------------------------------------------------------------------
// [[Rcpp::export]]
SEXP js_cpp(IntegerVector l, NumericVector z) {
    int l_size = l.size();
    int z_size = z.size();
    
    if (should_return_matrix(l_size, z_size)) {
        NumericMatrix result(l_size, z_size);
        
        for (int i = 0; i < l_size; ++i) {
            int li = l[i];
            for (int j = 0; j < z_size; ++j) {
                result(i, j) = js_single_impl(li, z[j]);
            }
        }
        return result;
    } else {
        int N = std::max(l_size, z_size);
        NumericVector result(N);
        
        for (int i = 0; i < N; ++i) {
            int li = l[i % l_size];
            double zi = z[i % z_size];
            result[i] = js_single_impl(li, zi);
        }
        return result;
    }
}

// --------------------------------------------------------------------------
// 8. Spherical Bessel y_l(z)
// --------------------------------------------------------------------------
// [[Rcpp::export]]
SEXP ys_cpp(IntegerVector l, NumericVector z) {
    int l_size = l.size();
    int z_size = z.size();
    
    if (should_return_matrix(l_size, z_size)) {
        NumericMatrix result(l_size, z_size);
        
        for (int i = 0; i < l_size; ++i) {
            int li = l[i];
            for (int j = 0; j < z_size; ++j) {
                result(i, j) = ys_single_impl(li, z[j]);
            }
        }
        return result;
    } else {
        int N = std::max(l_size, z_size);
        NumericVector result(N);
        
        for (int i = 0; i < N; ++i) {
            int li = l[i % l_size];
            double zi = z[i % z_size];
            result[i] = ys_single_impl(li, zi);
        }
        return result;
    }
}

// --------------------------------------------------------------------------
// 9. Spherical Hankel h_l(z)
// --------------------------------------------------------------------------
// [[Rcpp::export]]
SEXP hs_cpp(IntegerVector l, NumericVector z) {
    int l_size = l.size();
    int z_size = z.size();
    
    if (should_return_matrix(l_size, z_size)) {
        ComplexMatrix result(l_size, z_size);
        
        for (int i = 0; i < l_size; ++i) {
            int li = l[i];
            for (int j = 0; j < z_size; ++j) {
                std::complex<double> val = hs_single_impl(li, z[j]);
                result(i, j) = to_Rcomplex(val);
            }
        }
        return result;
    } else {
        int N = std::max(l_size, z_size);
        ComplexVector result(N);
        
        for (int i = 0; i < N; ++i) {
            int li = l[i % l_size];
            double zi = z[i % z_size];
            std::complex<double> val = hs_single_impl(li, zi);
            result[i] = to_Rcomplex(val);
        }
        return result;
    }
}

// --------------------------------------------------------------------------
// 10. Spherical Bessel j_l - k-th derivative
// --------------------------------------------------------------------------
// [[Rcpp::export]]
SEXP js_deriv_cpp(IntegerVector l, NumericVector z, int k) {
    if (k < 0) {
        Rcpp::stop("Derivative order k must be non-negative");
    }
    if (k == 0) {
        return js_cpp(l, z);
    }
    
    int l_size = l.size();
    int z_size = z.size();
    
    if (should_return_matrix(l_size, z_size)) {
        NumericMatrix result(l_size, z_size);
        
        for (int i = 0; i < l_size; ++i) {
            int li = l[i];
            for (int j = 0; j < z_size; ++j) {
                result(i, j) = js_deriv_single_impl(li, z[j], k);
            }
        }
        return result;
    } else {
        int N = std::max(l_size, z_size);
        NumericVector result(N);
        
        for (int i = 0; i < N; ++i) {
            int li = l[i % l_size];
            double zi = z[i % z_size];
            result[i] = js_deriv_single_impl(li, zi, k);
        }
        return result;
    }
}

// --------------------------------------------------------------------------
// 11. Spherical Bessel y_l - k-th derivative
// --------------------------------------------------------------------------
// [[Rcpp::export]]
SEXP ys_deriv_cpp(IntegerVector l, NumericVector z, int k) {
    if (k < 0) {
        Rcpp::stop("Derivative order k must be non-negative");
    }
    if (k == 0) {
        return ys_cpp(l, z);
    }
    
    int l_size = l.size();
    int z_size = z.size();
    
    if (should_return_matrix(l_size, z_size)) {
        NumericMatrix result(l_size, z_size);
        
        for (int i = 0; i < l_size; ++i) {
            int li = l[i];
            for (int j = 0; j < z_size; ++j) {
                result(i, j) = ys_deriv_single_impl(li, z[j], k);
            }
        }
        return result;
    } else {
        int N = std::max(l_size, z_size);
        NumericVector result(N);
        
        for (int i = 0; i < N; ++i) {
            int li = l[i % l_size];
            double zi = z[i % z_size];
            result[i] = ys_deriv_single_impl(li, zi, k);
        }
        return result;
    }
}

// --------------------------------------------------------------------------
// 12. Spherical Hankel h_l - k-th derivative
// --------------------------------------------------------------------------
// [[Rcpp::export]]
SEXP hs_deriv_cpp(IntegerVector l, NumericVector z, int k) {
    if (k < 0) {
        Rcpp::stop("Derivative order k must be non-negative");
    }
    if (k == 0) {
        return hs_cpp(l, z);
    }
    
    int l_size = l.size();
    int z_size = z.size();
    
    if (should_return_matrix(l_size, z_size)) {
        ComplexMatrix result(l_size, z_size);
        
        for (int i = 0; i < l_size; ++i) {
            int li = l[i];
            for (int j = 0; j < z_size; ++j) {
                std::complex<double> val = hs_deriv_single_impl(li, z[j], k);
                result(i, j) = to_Rcomplex(val);
            }
        }
        return result;
    } else {
        int N = std::max(l_size, z_size);
        ComplexVector result(N);
        
        for (int i = 0; i < N; ++i) {
            int li = l[i % l_size];
            double zi = z[i % z_size];
            std::complex<double> val = hs_deriv_single_impl(li, zi, k);
            result[i] = to_Rcomplex(val);
        }
        return result;
    }
}