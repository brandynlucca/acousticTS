// [[Rcpp::depends(BH)]]
#include <Rcpp.h>
#include <boost/math/special_functions/legendre.hpp>
#include <complex>
#include <cmath>

using namespace Rcpp;
using boost::math::tgamma;
using boost::math::legendre_q;
using boost::math::legendre_p;
using std::complex;

// Convert the internal std::complex representation used by the special-
// function code into the Rcomplex layout expected by Rcpp.
Rcomplex toRcomplex(const std::complex<double>& z) {
    Rcomplex r;
    r.r = z.real();
    r.i = z.imag();
    return r;
}

// Evaluate the Gauss hypergeometric series directly. The series is only used
// in the fractional-order Ferrers P helper, so the implementation favors
// clarity over a more elaborate continuation strategy.
complex<double> hyper2F1_series(double a, double b, double c, complex<double> z) {
    const int max_iter = 20000;
    const double tol = 1e-14;

    complex<double> term(1.0,0.0);
    complex<double> sum = term;

    for(int n=1; n<max_iter; ++n) {
        term *= ((a + n - 1.0)*(b + n - 1.0)/((c + n - 1.0)*n)) * z;
        sum += term;
        if(std::abs(term) < tol * std::abs(sum)) break;
    }
    return sum;
}

// Ferrers P_nu(x) on [-1, 1], expressed through the hypergeometric form.
double ferrers_P(double nu, double x) {
    complex<double> z((1.0 - x)/2.0, 0.0);
    complex<double> F = hyper2F1_series(-nu, nu+1.0, 1.0, z);
    return std::real(F);
}

// Numerical contour-style integral used for fractional Legendre P when the
// order is non-integer. The same helper is used on and off the canonical
// interval by switching the branch of the square-root term.
double P_fractional_complex_integral(double nu, double x) {
    const int N = 4000; // integration points
    double h = 2.0 * M_PI / N;

    complex<double> I(0.0,1.0);
    complex<double> sum(0.0, 0.0);

    // r handles the "radius" in the complex plane
    complex<double> r = (std::abs(x) <= 1.0) ? complex<double>(std::sqrt(1.0 - x*x), 0.0) : complex<double>(std::sqrt(x*x - 1.0), 0.0);

    for(int k=0; k<N; ++k) {
        double theta = h * (k + 0.5);
        complex<double> val;

        if(std::abs(x) <= 1.0) {
            val = complex<double>(x,0.0) + r * std::cos(theta) * I; // complex branch
        } else {
            val = complex<double>(x,0.0) + r * std::cos(theta);      // real branch
        }

        sum += std::pow(val, nu);
    }

    return std::real(sum * h / (2.0 * M_PI));
}
// [[Rcpp::export]]
NumericMatrix Pn_cpp(NumericVector n, NumericVector x) {
    int N = n.size();
    int X = x.size();
    NumericMatrix result(N, X);

    for(int i=0; i<N; ++i) {
        double ni = n[i];
        bool isInt = std::abs(ni - std::round(ni)) < 1e-14;

        for(int j=0; j<X; ++j) {
            double xj = x[j];

            if(isInt) {
                // Integer degrees use the standard three-term recurrence.
                int ni_int = static_cast<int>(std::round(ni));
                if(ni_int == 0) result(i,j) = 1.0;
                else if(ni_int == 1) result(i,j) = xj;
                else {
                    double Pn_2 = 1.0;
                    double Pn_1 = xj;
                    double Pn_val = 0.0;
                    for(int k=2; k<=ni_int; ++k) {
                        Pn_val = ((2.0*k-1.0)*xj*Pn_1 - (k-1.0)*Pn_2)/k;
                        Pn_2 = Pn_1;
                        Pn_1 = Pn_val;
                    }
                    result(i,j) = Pn_val;
                }
            } else {
                // Fractional degrees fall back to the integral helper.
                result(i,j) = P_fractional_complex_integral(ni, xj);
            }
        }
    }

    return result;
}

// =============================================================================
// k-th Derivative of Legendre Polynomial P_n(x)
// Uses Rodrigues-type formula and recurrence relations
// =============================================================================
// [[Rcpp::export]]
NumericMatrix Pn_deriv_cpp(NumericVector n, NumericVector x, int k) {
    int N = n.size();
    int X = x.size();
    NumericMatrix result(N, X);
    
    if (k < 0) {
        Rcpp::stop("Derivative order k must be non-negative");
    }
    
    if (k == 0) {
        return Pn_cpp(n, x);
    }
    
    for (int i = 0; i < N; ++i) {
        double ni = n[i];
        int ni_int = static_cast<int>(std::round(ni));
        bool is_integer = std::abs(ni - ni_int) < 1e-10;
        
        for (int j = 0; j < X; ++j) {
            double xj = x[j];
            
            // For integer n, use the relation:
            // d^k/dx^k P_n(x) is related to associated Legendre P_n^k(x)           
            if (is_integer && k > ni_int) {
                // k-th derivative of polynomial of degree n is 0 if k > n
                result(i, j) = 0.0;
            } else if (is_integer) {
                // Using Boost associated Legendre:
                // d^k P_n / dx^k at interior points can be computed as:
                double val;
                if (std::abs(xj * xj - 1.0) < 1e-10) {
                    // General: product_{j=0}^{k-1} (n-j)(n+1+j) / (2^k * k!)
                    double prod = 1.0;
                    for (int m = 0; m < k; ++m) {
                        prod *= (ni - m) * (ni + 1 + m);
                    }
                    double factorial_k = 1.0;
                    for (int m = 1; m <= k; ++m) factorial_k *= m;
                    val = prod / (std::pow(2.0, k) * factorial_k);
                    
                    if (xj < 0) {
                        // P_n^(k)(-1) = (-1)^(n+k) * P_n^(k)(1)
                        val *= ((ni_int + k) % 2 == 0) ? 1.0 : -1.0;
                    }
                } else {
                    // Interior point: use associated Legendre
                    // d^m P_n/dx^m = (-1)^m (1-x^2)^(-m/2) P_n^m(x)
                    
                    double assoc_legendre = boost::math::legendre_p(ni_int, k, xj);
                    double sign = (k % 2 == 0) ? 1.0 : -1.0;
                    double factor = std::pow(1.0 - xj * xj, -k / 2.0);
                    val = sign * factor * assoc_legendre;
                }
                result(i, j) = val;
            } else {
                // Fractional-order derivatives are approximated numerically.
                double h = 1e-6;
                if (k == 1) {
                    NumericVector x_plus = NumericVector::create(xj + h);
                    NumericVector x_minus = NumericVector::create(xj - h);
                    NumericVector ni_vec = NumericVector::create(ni);
                    double fp = Pn_cpp(ni_vec, x_plus)(0, 0);
                    double fm = Pn_cpp(ni_vec, x_minus)(0, 0);
                    result(i, j) = (fp - fm) / (2.0 * h);
                } else {
                    // Higher derivatives use a simple finite-difference stencil.
                    Rcpp::warning("Fractional order derivatives k>1 use finite differences");
                    NumericVector ni_vec = NumericVector::create(ni);
                    double sum = 0.0;
                    for (int m = 0; m <= k; ++m) {
                        double binom = 1.0;
                        for (int b = 0; b < m; ++b) {
                            binom *= (double)(k - b) / (b + 1);
                        }
                        double sign = (m % 2 == 0) ? 1.0 : -1.0;
                        NumericVector x_pt = NumericVector::create(xj + (k/2.0 - m) * h);
                        double f_pt = Pn_cpp(ni_vec, x_pt)(0, 0);
                        sum += sign * binom * f_pt;
                    }
                    result(i, j) = sum / std::pow(h, k);
                }
            }
        }
    }
    
    return result;
}

// Fractional-order helper for Q_nu(x). The implementation mirrors the P_nu
// helper but switches formulas according to whether x lies inside or outside
// the canonical interval.
double Q_fractional_complex_integral(double nu, double x) {
    const int N = 4000;
    
    if(std::abs(x) <= 1.0) {
        // Original code for |x| <= 1
        double h = 2.0 * M_PI / N;
        complex<double> I(0.0,1.0);
        complex<double> sum(0.0, 0.0);
        complex<double> r(std::sqrt(1.0 - x*x), 0.0);

        for(int k=0; k<N; ++k) {
            double theta = h * (k + 0.5);
            complex<double> val = complex<double>(x,0.0) + r * std::cos(theta) * I;
            sum += std::pow(val, nu);
        }
        return std::real(sum * h / (2.0 * M_PI));
    } else {
        // For |x| > 1: P_nu(x) = (1/pi) * integral_0^pi [(x + sqrt(x^2-1)*cos(t))^nu] dt
        double h = M_PI / N;
        double sum = 0.0;
        double sqrt_term = std::sqrt(x*x - 1.0);
        
        for(int k=0; k<N; ++k) {
            double t = h * (k + 0.5);
            double arg = x + sqrt_term * std::cos(t);
            sum += std::pow(arg, nu);
        }
        
        return sum * h / M_PI;
    }
}

// [[Rcpp::export]]
Rcpp::ComplexMatrix Qn_cpp(NumericVector nu, NumericVector x) {
    int nrow = nu.size();
    int ncol = x.size();
    Rcpp::ComplexMatrix result(nrow, ncol);

    for(int i = 0; i < nrow; ++i) {
        double nu_val = nu[i];
        
        for(int j = 0; j < ncol; ++j) {
            double x_val = x[j];
            
            // Singularity at x = ±1
            if(std::abs(std::abs(x_val) - 1.0) < 1e-14) {
                result(i,j).r = R_PosInf;
                result(i,j).i = R_PosInf;
                continue;
            }
            
            if(std::abs(x_val) < 1.0) {
                // On the canonical interval, Q_nu is real-valued.
                double nint = std::round(nu_val);
                bool is_integer = std::abs(nu_val - nint) < 1e-14;
                
                double Q;
                if(is_integer) {
                    Q = legendre_q(static_cast<unsigned>(nint), x_val);
                } else {
                    // Fractional degrees use the Ferrers identity relating
                    // Q_nu to P_nu(x) and P_nu(-x).
                    double P_x  = ferrers_P(nu_val, x_val);
                    double P_nx = ferrers_P(nu_val, -x_val);
                    double sinpi = std::sin(M_PI * nu_val);
                    double cospi = std::cos(M_PI * nu_val);
                    Q = (M_PI / (2.0 * sinpi)) * (cospi * P_x - P_nx);
                }
                
                result(i,j).r = Q;
                result(i,j).i = 0.0;
            } else {
                // Outside [-1, 1], Q_nu generally becomes complex-valued.
                
                // The real part comes from the decaying integral form.
                const int N_integ = 10000;
                double t_max = 30.0;
                double h = t_max / N_integ;
                double sum = 0.0;
                double sqrt_term = std::sqrt(x_val * x_val - 1.0);
                
                for(int k = 0; k < N_integ; ++k) {
                    double t = h * (k + 0.5);
                    double denominator = x_val + sqrt_term * std::cosh(t);
                    sum += 1.0 / std::pow(denominator, nu_val + 1.0);
                }
                
                double Q_real = sum * h;
                
                // The imaginary part follows from the branch-cut jump and is
                // proportional to P_nu(x).
                NumericVector nu_single = NumericVector::create(nu_val);
                NumericVector x_single = NumericVector::create(x_val);
                double P_nu_x = Pn_cpp(nu_single, x_single)(0, 0);
                
                double Q_imag = -M_PI * P_nu_x / 2.0;
                
                result(i,j).r = Q_real;
                result(i,j).i = Q_imag;
            }
        }
    }
    
    return result;
}

// =============================================================================
// k-th Derivative of Legendre Function Q_n(x)
// =============================================================================
// [[Rcpp::export]]
ComplexMatrix Qn_deriv_cpp(NumericVector n, NumericVector x, int k) {
    int N = n.size();
    int X = x.size();
    ComplexMatrix result(N, X);
    
    if (k < 0) {
        Rcpp::stop("Derivative order k must be non-negative");
    }
    
    if (k == 0) {
        return Qn_cpp(n, x);
    }
    
    // Q_n derivatives are evaluated numerically here because the branch
    // structure makes the closed-form relations less convenient to maintain.
    double h = 1e-6;
    
    for (int i = 0; i < N; ++i) {
        NumericVector ni_vec = NumericVector::create(n[i]);
        
        for (int j = 0; j < X; ++j) {
            double xj = x[j];
            
            if (k == 1) {
                // First derivative via a central finite difference.
                NumericVector x_plus = NumericVector::create(xj + h);
                NumericVector x_minus = NumericVector::create(xj - h);
                
                ComplexMatrix Qp = Qn_cpp(ni_vec, x_plus);
                ComplexMatrix Qm = Qn_cpp(ni_vec, x_minus);
                
                std::complex<double> Qp_val(Qp(0, 0).r, Qp(0, 0).i);
                std::complex<double> Qm_val(Qm(0, 0).r, Qm(0, 0).i);
                
                std::complex<double> deriv = (Qp_val - Qm_val) / (2.0 * h);
                result(i, j).r = deriv.real();
                result(i, j).i = deriv.imag();
            } else {
                // Higher derivatives via a symmetric finite-difference stencil.
                std::complex<double> sum(0.0, 0.0);
                for (int m = 0; m <= k; ++m) {
                    double binom = 1.0;
                    for (int b = 0; b < m; ++b) {
                        binom *= (double)(k - b) / (b + 1);
                    }
                    double sign = (m % 2 == 0) ? 1.0 : -1.0;
                    NumericVector x_pt = NumericVector::create(xj + (k/2.0 - m) * h);
                    ComplexMatrix Q_pt = Qn_cpp(ni_vec, x_pt);
                    std::complex<double> Q_val(Q_pt(0, 0).r, Q_pt(0, 0).i);
                    sum += sign * binom * Q_val;
                }
                std::complex<double> deriv = sum / std::pow(h, k);
                result(i, j).r = deriv.real();
                result(i, j).i = deriv.imag();
            }
        }
    }
    
    return result;
}
