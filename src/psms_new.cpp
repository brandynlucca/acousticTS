#include <Rcpp.h>
#include <limits>
#include <cmath>
#include <complex>

extern "C" {
  void cprofcn_interface_new(double* c, int* m, int* lnum, int* ioprad, double* x1, int* iopang, int* iopnorm, int* narg, double* arg,
                          double* r1c, int* ir1e, double* r1dc, int* ir1de, double* r2c, int* ir2e, double* r2dc, int* ir2de, int* naccr,
                          double* s1c, int* is1e, double* s1dc, int* is1de, int* naccs);
}

// Helper: get reduced epsilon (next lower power of ten)
double reduced_epsilon() {
    double eps = std::numeric_limits<double>::epsilon();
    int exp = static_cast<int>(std::floor(std::log10(eps)));
    return std::pow(10.0, exp);
}

// Helper: set the effective epsilon for compatibility with R
double effective_epsilon() {
    double eps = reduced_epsilon();
    // If reduced epsilon is less than 1e-15, use 1e-15 as minimum
    if (eps < 1e-15) return 1e-15;
    return eps;
}

// Use reduced epsilon for rounding
double round_to_epsilon(double x) {
    double eps = effective_epsilon();
    return std::round(x / eps) * eps;
}

// [[Rcpp::export]]
Rcpp::List safe_cprofcn(
    double c,
    int m,
    int lnum,
    Rcpp::NumericVector arg,
    int ioprad,
    int iopnorm,
    int iopang,
    double x1
) {
    int narg = std::max<int>(arg.size(), 1);

    // only padding to guard Fortran access
    // int lnum_safe = std::max(lnum, 2);
    // int narg_safe = std::max(narg, 1);

    std::vector<double> r1c(lnum, 0.0), r1dc(lnum, 0.0),
        r2c(lnum, 0.0), r2dc(lnum, 0.0);
    std::vector<int> ir1e(lnum, 0), ir1de(lnum, 0),
        ir2e(lnum, 0), ir2de(lnum, 0), naccr(lnum, 0);

    std::vector<double> s1c(lnum * narg, 0.0),
        s1dc(lnum * narg, 0.0);
    std::vector<int> is1e(lnum * narg, 0),
        is1de(lnum * narg, 0), naccs(lnum * narg, 0);

    std::vector<double> arg_c(narg);
    for (int i = 0; i < narg; ++i) arg_c[i] = arg[i];

    // *** keep lnum and narg exactly as-is ***
    cprofcn_interface_new(
        &c, &m, &lnum, &ioprad, &x1, &iopang, &iopnorm, &narg, arg_c.data(),
        r1c.data(), ir1e.data(), r1dc.data(), ir1de.data(),
        r2c.data(), ir2e.data(), r2dc.data(), ir2de.data(), naccr.data(),
        s1c.data(), is1e.data(), s1dc.data(), is1de.data(), naccs.data()
    );

    // Return results as a list
    return Rcpp::List::create(
        Rcpp::Named("r1c") = r1c,
        Rcpp::Named("ir1e") = ir1e,
        Rcpp::Named("r1dc") = r1dc,
        Rcpp::Named("ir1de") = ir1de,
        Rcpp::Named("r2c") = r2c,
        Rcpp::Named("ir2e") = ir2e,
        Rcpp::Named("r2dc") = r2dc,
        Rcpp::Named("ir2de") = ir2de,
        Rcpp::Named("naccr") = naccr,
        Rcpp::Named("s1c") = s1c,
        Rcpp::Named("is1e") = is1e,
        Rcpp::Named("s1dc") = s1dc,
        Rcpp::Named("is1de") = is1de,
        Rcpp::Named("naccs") = naccs
    );

    // trim to logical sizes — nothing else
    // return Rcpp::List::create(
    //     Rcpp::Named("r1c") = Rcpp::NumericVector(r1c),
    //     Rcpp::Named("ir1e") = Rcpp::IntegerVector(ir1e),
    //     Rcpp::Named("r1dc") = Rcpp::NumericVector(r1dc),
    //     Rcpp::Named("ir1de") = Rcpp::IntegerVector(ir1de),
    //     Rcpp::Named("r2c") = Rcpp::NumericVector(r2c),
    //     Rcpp::Named("ir2e") = Rcpp::IntegerVector(ir2e),
    //     Rcpp::Named("r2dc") = Rcpp::NumericVector(r2dc),
    //     Rcpp::Named("ir2de") = Rcpp::IntegerVector(ir2de),
    //     Rcpp::Named("naccr") = Rcpp::IntegerVector(naccr),
    //     Rcpp::Named("s1c") = Rcpp::NumericVector(s1c),
    //     Rcpp::Named("is1e") = Rcpp::IntegerVector(is1e),
    //     Rcpp::Named("s1dc") = Rcpp::NumericVector(s1dc),
    //     Rcpp::Named("is1de") = Rcpp::IntegerVector(is1de),
    //     Rcpp::Named("naccs") = Rcpp::IntegerVector(naccs)
    // );
// }
    // return Rcpp::List::create(
    //     Rcpp::Named("r1c")   = Rcpp::NumericVector(r1c.begin(),   r1c.begin()   + lnum),
    //     Rcpp::Named("ir1e")  = Rcpp::IntegerVector(ir1e.begin(),  ir1e.begin()  + lnum),
    //     Rcpp::Named("r1dc")  = Rcpp::NumericVector(r1dc.begin(),  r1dc.begin()  + lnum),
    //     Rcpp::Named("ir1de") = Rcpp::IntegerVector(ir1de.begin(), ir1de.begin() + lnum),
    //     Rcpp::Named("r2c")   = Rcpp::NumericVector(r2c.begin(),   r2c.begin()   + lnum),
    //     Rcpp::Named("ir2e")  = Rcpp::IntegerVector(ir2e.begin(),  ir2e.begin()  + lnum),
    //     Rcpp::Named("r2dc")  = Rcpp::NumericVector(r2dc.begin(),  r2dc.begin()  + lnum),
    //     Rcpp::Named("ir2de") = Rcpp::IntegerVector(ir2de.begin(), ir2de.begin() + lnum),
    //     Rcpp::Named("naccr") = Rcpp::IntegerVector(naccr.begin(), naccr.begin() + lnum),

    //     // keep full flattened Fortran order (lnum * narg)
    //     Rcpp::Named("s1c")   = Rcpp::NumericVector(s1c.begin(),   s1c.begin()   + lnum * narg),
    //     Rcpp::Named("is1e")  = Rcpp::IntegerVector(is1e.begin(),  is1e.begin()  + lnum * narg),
    //     Rcpp::Named("s1dc")  = Rcpp::NumericVector(s1dc.begin(),  s1dc.begin()  + lnum * narg),
    //     Rcpp::Named("is1de") = Rcpp::IntegerVector(is1de.begin(), is1de.begin() + lnum * narg),
    //     Rcpp::Named("naccs") = Rcpp::IntegerVector(naccs.begin(), naccs.begin() + lnum * narg)
    // );
}



// Rcpp::List safe_cprofcn(
//     double c,
//     int m,
//     int lnum,
//     Rcpp::NumericVector arg,
//     int ioprad,
//     int iopnorm,
//     int iopang,
//     double x1
// ) {
//     int narg = arg.size();

//     // Ensure minimum array length of 2 to prevent Fortran bounds errors
//     int lnum_safe = std::max(lnum, 2);
//     int narg_safe = std::max(narg, 1);

//     std::vector<double> r1c(lnum_safe), r1dc(lnum_safe), r2c(lnum_safe), r2dc(lnum_safe);
//     std::vector<int> ir1e(lnum_safe), ir1de(lnum_safe), ir2e(lnum_safe), ir2de(lnum_safe), naccr(lnum_safe);
//     std::vector<double> s1c(lnum_safe * narg_safe), s1dc(lnum_safe * narg_safe);
//     std::vector<int> is1e(lnum_safe * narg_safe), is1de(lnum_safe * narg_safe), naccs(lnum_safe * narg_safe);

//     // Call Fortran routine with "safe" array lengths
//     cprofcn_interface_new(
//         &c, &m, &lnum_safe, &ioprad, &x1, &iopang, &iopnorm, &narg_safe, arg.begin(),
//         r1c.data(), ir1e.data(), r1dc.data(), ir1de.data(),
//         r2c.data(), ir2e.data(), r2dc.data(), ir2de.data(), naccr.data(),
//         s1c.data(), is1e.data(), s1dc.data(), is1de.data(), naccs.data()
//     );

//     // Return only the first lnum entries (trim the extra allocation)
//     return Rcpp::List::create(
//         Rcpp::Named("r1c") = Rcpp::NumericVector(r1c.begin(), r1c.begin() + lnum),
//         Rcpp::Named("ir1e") = Rcpp::IntegerVector(ir1e.begin(), ir1e.begin() + lnum),
//         Rcpp::Named("r1dc") = Rcpp::NumericVector(r1dc.begin(), r1dc.begin() + lnum),
//         Rcpp::Named("ir1de") = Rcpp::IntegerVector(ir1de.begin(), ir1de.begin() + lnum),
//         Rcpp::Named("r2c") = Rcpp::NumericVector(r2c.begin(), r2c.begin() + lnum),
//         Rcpp::Named("ir2e") = Rcpp::IntegerVector(ir2e.begin(), ir2e.begin() + lnum),
//         Rcpp::Named("r2dc") = Rcpp::NumericVector(r2dc.begin(), r2dc.begin() + lnum),
//         Rcpp::Named("ir2de") = Rcpp::IntegerVector(ir2de.begin(), ir2de.begin() + lnum),
//         Rcpp::Named("naccr") = Rcpp::IntegerVector(naccr.begin(), naccr.begin() + lnum),
//         Rcpp::Named("s1c") = Rcpp::NumericVector(s1c.begin(), s1c.begin() + lnum * narg),
//         Rcpp::Named("is1e") = Rcpp::IntegerVector(is1e.begin(), is1e.begin() + lnum * narg),
//         Rcpp::Named("s1dc") = Rcpp::NumericVector(s1dc.begin(), s1dc.begin() + lnum * narg),
//         Rcpp::Named("is1de") = Rcpp::IntegerVector(is1de.begin(), is1de.begin() + lnum * narg),
//         Rcpp::Named("naccs") = Rcpp::IntegerVector(naccs.begin(), naccs.begin() + lnum * narg)
//     );
// }

Rcpp::List cprofcn(
    double c,
    int m,
    int lnum,
    Rcpp::NumericVector arg,
    int ioprad,
    int iopnorm,
    int iopang,
    double x1
) {
    // Validation
    if (m < 0) {
        Rcpp::stop("'m' must be greater than or equal to 0.");
    }
    if (m != static_cast<int>(m)) {
        Rcpp::stop("'m' must be an integer.");
    }

    // Ensure minimum array length of 2 to prevent Fortran bounds errors
    int narg = arg.size();

    // Round arg to 15 decimal places for numerical stability
    Rcpp::NumericVector arg_rounded(narg);
    for (int i = 0; i < narg; ++i) {
        arg_rounded[i] = round_to_epsilon(arg[i]);
    }

    // Default output array sizes
    std::vector<double> r1c(lnum), r1dc(lnum), r2c(lnum), r2dc(lnum);
    std::vector<int> ir1e(lnum), ir1de(lnum), ir2e(lnum), ir2de(lnum), naccr(lnum);
    std::vector<double> s1c(lnum * narg), s1dc(lnum * narg);
    std::vector<int> is1e(lnum * narg), is1de(lnum * narg), naccs(lnum * narg);

    // Call Fortran routine
    cprofcn_interface_new(
        &c, &m, &lnum, &ioprad, &x1, &iopang, &iopnorm, &narg, arg_rounded.begin(),
        r1c.data(), ir1e.data(), r1dc.data(), ir1de.data(), r2c.data(), ir2e.data(), r2dc.data(), ir2de.data(), naccr.data(),
        s1c.data(), is1e.data(), s1dc.data(), is1de.data(), naccs.data()
    );

    // Return results as a list
    return Rcpp::List::create(
        Rcpp::Named("r1c") = r1c,
        Rcpp::Named("ir1e") = ir1e,
        Rcpp::Named("r1dc") = r1dc,
        Rcpp::Named("ir1de") = ir1de,
        Rcpp::Named("r2c") = r2c,
        Rcpp::Named("ir2e") = ir2e,
        Rcpp::Named("r2dc") = r2dc,
        Rcpp::Named("ir2de") = ir2de,
        Rcpp::Named("naccr") = naccr,
        Rcpp::Named("s1c") = s1c,
        Rcpp::Named("is1e") = is1e,
        Rcpp::Named("s1dc") = s1dc,
        Rcpp::Named("is1de") = is1de,
        Rcpp::Named("naccs") = naccs
    );
}



// [[Rcpp::export]]
Rcpp::List Smn_cpp(int m, int n, double c, Rcpp::NumericVector arg, bool normalize = false) {
    // Input validation
    if (n < 0) {
        Rcpp::stop("'n' must be greater than or equal to 0.");
    }
    if (n < m) {
        Rcpp::stop("'n' must be greater than or equal to argument 'm'.");
    }
    if (n != static_cast<int>(n)) {
        Rcpp::stop("'n' must be an integer.");
    }
    // Analytical limit of eta at eta == 1
    if (!std::all_of(arg.begin(), arg.end(), [](double eta) { return std::abs(eta) <= 1.0; })) {
        Rcpp::stop("|eta| must be less than or equal to 1 for the angular prolate domain.");
    }

    // Initialize constants and arguments for input into lower Fortran code
    int ioprad = 0;
    int iopnorm = normalize ? 1 : 0;
    int iopang = 2;
    double x1 = 0.0;
    int lnum = n - m + 1;
    int narg = arg.size();

    // Run C++-wrapper
    Rcpp::List result = cprofcn(c, m, lnum, arg, ioprad, iopnorm, iopang, x1);

    // Zero-based index for C++
    int idx = n - m;

    // Power correction/offsets
    Rcpp::IntegerVector is1e = result["is1e"];
    Rcpp::IntegerVector is1de = result["is1de"];
    Rcpp::NumericVector s1c = result["s1c"];
    Rcpp::NumericVector s1dc = result["s1dc"];

    Rcpp::NumericVector value(narg), derivative(narg);

    // Check for valid exponent and not NA/NaN
    for (int i = 0; i < narg; ++i) {
        int offset = idx + i * lnum;
        double pw = 1.0, pwde = 1.0;
        if (is1e.size() > offset && !Rcpp::IntegerVector::is_na(is1e[offset])) {
        pw = std::pow(10.0, is1e[offset]);
        }
        if (is1de.size() > offset && !Rcpp::IntegerVector::is_na(is1de[offset])) {
        pwde = std::pow(10.0, is1de[offset]);
        }
        value[i] = (s1c.size() > offset && !Rcpp::NumericVector::is_na(s1c[offset])) ? s1c[offset] * pw : NA_REAL;
        derivative[i] = (s1dc.size() > offset && !Rcpp::NumericVector::is_na(s1dc[offset])) ? s1dc[offset] * pwde : NA_REAL;
    }

    return Rcpp::List::create(
        Rcpp::Named("value") = value,
        Rcpp::Named("derivative") = derivative
    );
}

// [[Rcpp::export]]
Rcpp::List Smn_cpp_safe(int m, int n, double c, Rcpp::NumericVector arg, bool normalize = false) {
    // Input validation
    if (n < 0) {
        Rcpp::stop("'n' must be greater than or equal to 0.");
    }
    if (n < m) {
        Rcpp::stop("'n' must be greater than or equal to argument 'm'.");
    }
    if (n != static_cast<int>(n)) {
        Rcpp::stop("'n' must be an integer.");
    }
    // Analytical limit of eta at eta == 1
    if (!std::all_of(arg.begin(), arg.end(), [](double eta) { return std::abs(eta) <= 1.0; })) {
        Rcpp::stop("|eta| must be less than or equal to 1 for the angular prolate domain.");
    }

    // Initialize constants and arguments for input into lower Fortran code
    int ioprad = 0;
    int iopnorm = normalize ? 1 : 0;
    int iopang = 2;
    double x1 = 0.0;
    int lnum = n - m + 1;
    int narg = arg.size();

    // Run C++-wrapper
    Rcpp::List result = safe_cprofcn(c, m, lnum, arg, ioprad, iopnorm, iopang, x1);

    // Zero-based index for C++
    int idx = n - m;

    // Power correction/offsets
    Rcpp::IntegerVector is1e = result["is1e"];
    Rcpp::IntegerVector is1de = result["is1de"];
    Rcpp::NumericVector s1c = result["s1c"];
    Rcpp::NumericVector s1dc = result["s1dc"];

    Rcpp::NumericVector value(narg), derivative(narg);

    // Check for valid exponent and not NA/NaN
    for (int i = 0; i < narg; ++i) {
        int offset = idx + i * lnum;
        double pw = 1.0, pwde = 1.0;
        if (is1e.size() > offset && !Rcpp::IntegerVector::is_na(is1e[offset])) {
        pw = std::pow(10.0, is1e[offset]);
        }
        if (is1de.size() > offset && !Rcpp::IntegerVector::is_na(is1de[offset])) {
        pwde = std::pow(10.0, is1de[offset]);
        }
        value[i] = (s1c.size() > offset && !Rcpp::NumericVector::is_na(s1c[offset])) ? s1c[offset] * pw : NA_REAL;
        derivative[i] = (s1dc.size() > offset && !Rcpp::NumericVector::is_na(s1dc[offset])) ? s1dc[offset] * pwde : NA_REAL;
    }

    return Rcpp::List::create(
        Rcpp::Named("value") = value,
        Rcpp::Named("derivative") = derivative
    );
}

// Helper to avoid recursive Fortran calls for `kind=3` and `kind=4`
struct RmnResult {
    double val1;
    double der1;
    double val2;
    double der2;
    bool valid1;
    bool valid2;
};

RmnResult Rmn_higher_order(int m, int n, int lnum, double c, double x1) {
    int iopnorm = 0;
    int iopang = 2;
    Rcpp::NumericVector arg = Rcpp::NumericVector::create(1.0);

    // Note: 'x1' had the radial shift already applied in the parent function
    // Call Fortran for kind 1
    Rcpp::List result1 = safe_cprofcn(c, m, lnum, arg, 1, iopnorm, iopang, x1);
    // Call Fortran for kind 2
    Rcpp::List result2 = safe_cprofcn(c, m, lnum, arg, 2, iopnorm, iopang, x1);

    int idx = n - m;
    double pw1 = 1.0, pwde1 = 1.0, pw2 = 1.0, pwde2 = 1.0;

    Rcpp::IntegerVector ie1 = result1["ir1e"];
    Rcpp::IntegerVector ide1 = result1["ir1de"];
    Rcpp::NumericVector val1 = result1["r1c"];
    Rcpp::NumericVector der1 = result1["r1dc"];

    Rcpp::IntegerVector ie2 = result2["ir2e"];
    Rcpp::IntegerVector ide2 = result2["ir2de"];
    Rcpp::NumericVector val2 = result2["r2c"];
    Rcpp::NumericVector der2 = result2["r2dc"];

    if (ie1.size() > idx && !Rcpp::IntegerVector::is_na(ie1[idx])) pw1 = std::pow(10.0, ie1[idx]);
    if (ide1.size() > idx && !Rcpp::IntegerVector::is_na(ide1[idx])) pwde1 = std::pow(10.0, ide1[idx]);
    if (ie2.size() > idx && !Rcpp::IntegerVector::is_na(ie2[idx])) pw2 = std::pow(10.0, ie2[idx]);
    if (ide2.size() > idx && !Rcpp::IntegerVector::is_na(ide2[idx])) pwde2 = std::pow(10.0, ide2[idx]);

    double R1 = (val1.size() > idx && !Rcpp::NumericVector::is_na(val1[idx])) ? val1[idx] * pw1 : NA_REAL;
    double Rd1 = (der1.size() > idx && !Rcpp::NumericVector::is_na(der1[idx])) ? der1[idx] * pwde1 : NA_REAL;
    double R2 = (val2.size() > idx && !Rcpp::NumericVector::is_na(val2[idx])) ? val2[idx] * pw2 : NA_REAL;
    double Rd2 = (der2.size() > idx && !Rcpp::NumericVector::is_na(der2[idx])) ? der2[idx] * pwde2 : NA_REAL;

    bool valid1 = !(val1.size() <= idx || Rcpp::NumericVector::is_na(val1[idx]));
    bool valid2 = !(val2.size() <= idx || Rcpp::NumericVector::is_na(val2[idx]));

    return {R1, Rd1, R2, Rd2, valid1, valid2};
}

// Rcpp::List Rmn_cpp(int m, int n, double c, double x1, int kind = 1) {
//     // Validation 
//     if (n < 0) {
//         Rcpp::stop("'n' must be greater than or equal to 0.");
//     }
//     if (n < m) {
//         Rcpp::stop("'n' must be greater than or equal to argument 'm'.");
//     }
//     if (n != static_cast<int>(n)) {
//         Rcpp::stop("'n' must be an integer.");
//     }
//     if (std::abs(x1) < 1.0) {
//         Rcpp::stop("|xi| must be greater than or equal to 1 for the radial prolate domain.");
//     }
//     if (kind < 1 || kind > 4) {
//         Rcpp::stop("'kind' must be either 1, 2, 3, or 4.");
//     }

//     // Initialize constants and arguments for input into lower Fortran code
//     int lnum = n - m + 1;
//     double x1_r = x1 - 1.0; // Radial shift
//     Rcpp::List result_r1, result_r2;
//     // Zero-based index for C++
//     int idx = n - m;

//   // Kind 1 and 2: direct call
//     if (kind == 1 || kind == 2) {
//         int ioprad = kind;    
//         Rcpp::List result = cprofcn(c, m, lnum, Rcpp::NumericVector::create(1.0), ioprad, 0, 0, x1_r);

//         // Power correction/offsets
//         Rcpp::IntegerVector ie = (kind == 1) ? result["ir1e"] : result["ir2e"];
//         Rcpp::IntegerVector ide = (kind == 1) ? result["ir1de"] : result["ir2de"];
//         Rcpp::NumericVector val = (kind == 1) ? result["r1c"] : result["r2c"];
//         Rcpp::NumericVector der = (kind == 1) ? result["r1dc"] : result["r2dc"];

//         double pw = (ie.size() > idx && !Rcpp::IntegerVector::is_na(ie[idx])) ? std::pow(10.0, ie[idx]) : 1.0;
//         double pwde = (ide.size() > idx && !Rcpp::IntegerVector::is_na(ide[idx])) ? std::pow(10.0, ide[idx]) : 1.0;

//         double R = (val.size() > idx && !Rcpp::NumericVector::is_na(val[idx])) ? val[idx] * pw : NA_REAL;
//         double Rd = (der.size() > idx && !Rcpp::NumericVector::is_na(der[idx])) ? der[idx] * pwde : NA_REAL;

//         // Additional check for kind == 2
//         if (kind == 2 && std::abs(R) < 1e-12) {
//             Rcpp::stop("Computed prolate spheroidal radial function of the second kind is numerically zero or unreliable for the given parameters (possibly due to singularity at xi = 1 or unsupported parameter regime).");
//         }

//         return Rcpp::List::create(
//             Rcpp::Named("value") = R,
//             Rcpp::Named("derivative") = Rd
//         );
//     }

//     // Kind 3 and 4: combine kind 1 and 2
//     RmnResult result = Rmn_higher_order(m, n, c, x1_r);

//     std::complex<double> val, der;
//     if (kind == 3) {
//         val = std::complex<double>(result.val1, result.val2);
//         der = std::complex<double>(result.der1, result.der2);
//     } else {
//         val = std::complex<double>(result.val1, -result.val2);
//         der = std::complex<double>(result.der1, -result.der2);
//     }

//     return Rcpp::List::create(
//         Rcpp::Named("value") = val,
//         Rcpp::Named("derivative") = der
//     );
// }

// [[Rcpp::export]]
Rcpp::List Rmn_cpp(int m, int n, double c, double x1, int kind = 1) {
    // Validation 
    if (n < 0) {
        Rcpp::stop("'n' must be greater than or equal to 0.");
    }
    if (n < m) {
        Rcpp::stop("'n' must be greater than or equal to argument 'm'.");
    }
    if (n != static_cast<int>(n)) {
        Rcpp::stop("'n' must be an integer.");
    }
    if (std::abs(x1) < 1.0) {
        Rcpp::stop("|xi| must be greater than or equal to 1 for the radial prolate domain.");
    }
    if (kind < 1 || kind > 4) {
        Rcpp::stop("'kind' must be either 1, 2, 3, or 4.");
    }

    // Initialize constants and arguments for input into lower Fortran code
    int lnum = n - m + 1;
    int lnum_safe = std::max(lnum, 2);
    double x1_r = x1 - 1.0; // Radial shift
    Rcpp::List result_r1, result_r2;
    // Zero-based index for C++
    int idx = n - m;

    // Kind 1 and 2: direct call
    if (kind == 1 || kind == 2) {
        int ioprad = kind;    
        Rcpp::List result = safe_cprofcn(c, m, lnum_safe, Rcpp::NumericVector::create(1.0), ioprad, 0, 0, x1_r);

        // Power correction/offsets
        Rcpp::IntegerVector ie = (kind == 1) ? result["ir1e"] : result["ir2e"];
        Rcpp::IntegerVector ide = (kind == 1) ? result["ir1de"] : result["ir2de"];
        Rcpp::NumericVector val = (kind == 1) ? result["r1c"] : result["r2c"];
        Rcpp::NumericVector der = (kind == 1) ? result["r1dc"] : result["r2dc"];

        double pw = (ie.size() > idx && !Rcpp::IntegerVector::is_na(ie[idx])) ? std::pow(10.0, ie[idx]) : 1.0;
        double pwde = (ide.size() > idx && !Rcpp::IntegerVector::is_na(ide[idx])) ? std::pow(10.0, ide[idx]) : 1.0;

        double R = (val.size() > idx && !Rcpp::NumericVector::is_na(val[idx])) ? val[idx] * pw : NA_REAL;
        double Rd = (der.size() > idx && !Rcpp::NumericVector::is_na(der[idx])) ? der[idx] * pwde : NA_REAL;

        // Additional check for kind == 2
        if (kind == 2 && std::abs(R) < 1e-12) {
            Rcpp::stop("Computed prolate spheroidal radial function of the second kind is numerically zero or unreliable for the given parameters (possibly due to singularity at xi = 1 or unsupported parameter regime).");
        }

        return Rcpp::List::create(
            Rcpp::Named("value") = R,
            Rcpp::Named("derivative") = Rd
        );
    }

    // Kind 3 and 4: combine kind 1 and 2
    RmnResult result = Rmn_higher_order(m, n, lnum_safe, c, x1_r);

    std::complex<double> val, der;
    if (kind == 3) {
        val = std::complex<double>(result.val1, result.val2);
        der = std::complex<double>(result.der1, result.der2);
    } else {
        val = std::complex<double>(result.val1, -result.val2);
        der = std::complex<double>(result.der1, -result.der2);
    }

    return Rcpp::List::create(
        Rcpp::Named("value") = val,
        Rcpp::Named("derivative") = der
    );
}

// [[Rcpp::export]]
Rcpp::List liquid_spheroidal_expansion_matrix(
    int m_max,
    int n_max,
    double chi_sw,
    double chi_body,
    Rcpp::NumericVector nodes,
    Rcpp::NumericVector weights
) {
    // Equivalent to \alpha_{n,l}^{m} related to Eq. 4 from Furusawa (1988)
    Rcpp::List coef_list(m_max + 1);
    
    // Cache spheroidal functions - call ONCE per (m, n) with ALL nodes
    std::vector<std::vector<Rcpp::NumericVector>> smn_ext_cache(m_max + 1);
    std::vector<std::vector<Rcpp::NumericVector>> smn_int_cache(m_max + 1);

    for (int m = 0; m <= m_max; ++m) {
        int size = n_max - m + 1;
        smn_ext_cache[m].resize(size);
        smn_int_cache[m].resize(size);

        for (int i = 0; i < size; ++i) {
            int n = m + i;
            // Call Smn_cpp ONCE with entire nodes vector, always normalize=true
            smn_ext_cache[m][i] = Rcpp::as<Rcpp::NumericVector>(Smn_cpp(m, n, chi_sw, nodes, true)["value"]);
            smn_int_cache[m][i] = Rcpp::as<Rcpp::NumericVector>(Smn_cpp(m, n, chi_body, nodes, true)["value"]);
        }
    }

    // Fill matrices
    for (int m = 0; m <= m_max; ++m) {
        int size = n_max - m + 1;
        Rcpp::NumericMatrix coef_mat(size, size);

        for (int li = 0; li < size; ++li) {
            for (int ni = 0; ni < size; ++ni) {
                // Get cached vectors
                Rcpp::NumericVector smn_ext_val = smn_ext_cache[m][ni];
                Rcpp::NumericVector smn_int_val = smn_int_cache[m][li];

                // Element-wise multiplication and weighted sum
                // Store as (li, ni) to match how it's read
                coef_mat(li, ni) = Rcpp::sum(smn_ext_val * smn_int_val * weights);
            }
        }
        coef_list[m] = coef_mat;
    }
    return coef_list;
}

// [[Rcpp::export]]
Rcpp::List safe_liquid_spheroidal_expansion_matrix(
    int m_max,
    int n_max,
    double chi_sw,
    double chi_body,
    Rcpp::NumericVector nodes,
    Rcpp::NumericVector weights
) {
    // Equivalent to \alpha_{n,l}^{m} related to Eq. 4 from Furusawa (1988)
    Rcpp::List coef_list(m_max + 1);
    
    // Cache spheroidal functions - call ONCE per (m, n) with ALL nodes
    std::vector<std::vector<Rcpp::NumericVector>> smn_ext_cache(m_max + 1);
    std::vector<std::vector<Rcpp::NumericVector>> smn_int_cache(m_max + 1);

    for (int m = 0; m <= m_max; ++m) {
        int size = n_max - m + 1;
        smn_ext_cache[m].resize(size);
        smn_int_cache[m].resize(size);

        for (int i = 0; i < size; ++i) {
            int n = m + i;
            // Call Smn_cpp ONCE with entire nodes vector, always normalize=true
            smn_ext_cache[m][i] = Rcpp::as<Rcpp::NumericVector>(Smn_cpp_safe(m, n, chi_sw, nodes, true)["value"]);
            smn_int_cache[m][i] = Rcpp::as<Rcpp::NumericVector>(Smn_cpp_safe(m, n, chi_body, nodes, true)["value"]);
        }
    }

    // Fill matrices
    for (int m = 0; m <= m_max; ++m) {
        int size = n_max - m + 1;
        Rcpp::NumericMatrix coef_mat(size, size);

        for (int li = 0; li < size; ++li) {
            for (int ni = 0; ni < size; ++ni) {
                // Get cached vectors
                Rcpp::NumericVector smn_ext_val = smn_ext_cache[m][ni];
                Rcpp::NumericVector smn_int_val = smn_int_cache[m][li];

                // Element-wise multiplication and weighted sum
                // Store as (li, ni) to match how it's read
                coef_mat(li, ni) = Rcpp::sum(smn_ext_val * smn_int_val * weights);
            }
        }
        coef_list[m] = coef_mat;
    }
    return coef_list;
}

Rcpp::List radial_external_incoming_matrices(
    int m_max,
    int n_max,
    double chi_sw,
    double xi,
    Rcpp::NumericMatrix smn_matrix
) {
    // Equivalent to R_{n,m}^{1} related to Eq. 4 from Furusawa (1988)
    Rcpp::List result(m_max + 1);

    for (int m = 0; m <= m_max; ++m) {
        int inner_len = n_max - m + 1;
        Rcpp::List m_list(inner_len);

        for (int n = m; n <= n_max; ++n) {
            int ni = n - m;
            double snm_val = smn_matrix(m, n);
            
            if (snm_val == 0.0 || std::isnan(snm_val)) {
                m_list[ni] = Rcpp::List::create(
                    Rcpp::Named("value") = 0.0,
                    Rcpp::Named("derivative") = 0.0
                );
            } else {
                Rcpp::List rmn = Rmn_cpp(m, n, chi_sw, xi, 1);
                m_list[ni] = Rcpp::List::create(
                    Rcpp::Named("value") = rmn["value"],
                    Rcpp::Named("derivative") = rmn["derivative"]
                );
            }
        }
        result[m] = m_list;
    }
    return result;
}

Rcpp::List radial_external_scattering_matrices(
    int m_max,
    int n_max,
    double chi_sw,
    double xi,
    Rcpp::NumericMatrix smn_matrix
) {
    // Equivalent to R_{n,m}^{3} related to Eq. 4 from Furusawa (1988)
    Rcpp::List result(m_max + 1);

    for (int m = 0; m <= m_max; ++m) {
        int inner_len = n_max - m + 1;
        Rcpp::List m_list(inner_len);

        for (int n = m; n <= n_max; ++n) {
            int ni = n - m;
            double snm_val = smn_matrix(m, n);
            
            if (snm_val == 0.0 || std::isnan(snm_val)) {
                m_list[ni] = Rcpp::List::create(
                    Rcpp::Named("value") = std::complex<double>(0.0, 0.0),
                    Rcpp::Named("derivative") = std::complex<double>(0.0, 0.0)
                );
            } else {
                Rcpp::List rmn = Rmn_cpp(m, n, chi_sw, xi, 3);
                m_list[ni] = Rcpp::List::create(
                    Rcpp::Named("value") = rmn["value"],
                    Rcpp::Named("derivative") = rmn["derivative"]
                );
            }
        }
        result[m] = m_list;
    }
    return result;
}

Rcpp::List radial_internal_incoming_matrices(
    int m_max,
    int n_max,
    double chi_body,
    double xi
) {
    // Equivalent to R_{n,m}^{1} related to Eq. 4 from Furusawa (1988)
    Rcpp::List result(m_max + 1);

    for (int m = 0; m <= m_max; ++m) {
        int inner_len = n_max - m + 1;
        Rcpp::List m_list(inner_len);

        // ALWAYS compute for all l from m to n_max (no Smn matrix check here!)
        for (int ell = m; ell <= n_max; ++ell) {
            int li = ell - m;
            Rcpp::List rmn = Rmn_cpp(m, ell, chi_body, xi, 1);
            m_list[li] = Rcpp::List::create(
                Rcpp::Named("value") = rmn["value"],
                Rcpp::Named("derivative") = rmn["derivative"]
            );
        }
        result[m] = m_list;
    }
    return result;
}

// [[Rcpp::export]]
Rcpp::List radial_external_incoming_matrix(
    int m_max,
    int n_max,
    double chi_sw,
    double xi
) {
    Rcpp::NumericMatrix value_mat(m_max + 1, n_max + 1);
    Rcpp::NumericMatrix deriv_mat(m_max + 1, n_max + 1);

    // initialize to NA
    for (int i = 0; i <= m_max; ++i)
        for (int j = 0; j <= n_max; ++j) {
            value_mat(i, j) = NA_REAL;
            deriv_mat(i, j) = NA_REAL;
        }

    for (int m = 0; m <= m_max; ++m) {
        for (int n = m; n <= n_max; ++n) {
            Rcpp::List rmn = Rmn_cpp(m, n, chi_sw, xi, 1);
            double v = Rcpp::as<double>(rmn["value"]);
            double d = Rcpp::as<double>(rmn["derivative"]);
            value_mat(m, n) = v;
            deriv_mat(m, n) = d;
        }
    }

    return Rcpp::List::create(
        Rcpp::Named("value") = value_mat,
        Rcpp::Named("derivative") = deriv_mat
    );
}

// [[Rcpp::export]]
Rcpp::List radial_external_scattering_matrix(
    int m_max,
    int n_max,
    double chi_sw,
    double xi
) {
    Rcpp::ComplexMatrix value_mat(m_max + 1, n_max + 1);
    Rcpp::ComplexMatrix deriv_mat(m_max + 1, n_max + 1);

    // initialize to NA components
    for (int i = 0; i <= m_max; ++i)
        for (int j = 0; j <= n_max; ++j) {
            value_mat(i, j).r = NA_REAL;
            value_mat(i, j).i = NA_REAL;
            deriv_mat(i, j).r = NA_REAL;
            deriv_mat(i, j).i = NA_REAL;
        }

    for (int m = 0; m <= m_max; ++m) {
        for (int n = m; n <= n_max; ++n) {
            Rcpp::List rmn = Rmn_cpp(m, n, chi_sw, xi, 3);
            std::complex<double> v = Rcpp::as<std::complex<double>>(rmn["value"]);
            std::complex<double> d = Rcpp::as<std::complex<double>>(rmn["derivative"]);
            value_mat(m, n).r = v.real();
            value_mat(m, n).i = v.imag();
            deriv_mat(m, n).r = d.real();
            deriv_mat(m, n).i = d.imag();
        }
    }

    return Rcpp::List::create(
        Rcpp::Named("value") = value_mat,
        Rcpp::Named("derivative") = deriv_mat
    );
}

// [[Rcpp::export]]
Rcpp::List radial_internal_incoming_matrix(
    int m_max,
    int n_max,
    double chi_body,
    double xi
) {
    Rcpp::NumericMatrix value_mat(m_max + 1, n_max + 1);
    Rcpp::NumericMatrix deriv_mat(m_max + 1, n_max + 1);

    // initialize to NA
    for (int i = 0; i <= m_max; ++i)
        for (int j = 0; j <= n_max; ++j) {
            value_mat(i, j) = NA_REAL;
            deriv_mat(i, j) = NA_REAL;
        }

    for (int m = 0; m <= m_max; ++m) {
        for (int ell = m; ell <= n_max; ++ell) {
            Rcpp::List rmn = Rmn_cpp(m, ell, chi_body, xi, 1);
            double v = Rcpp::as<double>(rmn["value"]);
            double d = Rcpp::as<double>(rmn["derivative"]);
            value_mat(m, ell) = v;
            deriv_mat(m, ell) = d;
        }
    }

    return Rcpp::List::create(
        Rcpp::Named("value") = value_mat,
        Rcpp::Named("derivative") = deriv_mat
    );
}

Rcpp::List boundary_coupling_regular_matrix(
    int m_max,
    int n_max,
    double density_body,
    double density_sw,
    Rcpp::List Rmn_e1_matrix,
    Rcpp::List Rmn_i1_matrix,
    Rcpp::NumericMatrix smn_matrix
) {
    // Equivalent to E_{n,l}^{m(1)} related to Eq. 4 from Furusawa (1988)
    Rcpp::List result(m_max + 1);
    for (int m = 0; m <= m_max; ++m) {
        int size = n_max - m + 1;
        Rcpp::ComplexMatrix mat(size, size);

        // Initialize to zero
        for (int li = 0; li < size; ++li) {
            for (int ni = 0; ni < size; ++ni) {
                mat(li, ni).r = 0.0;
                mat(li, ni).i = 0.0;
            }
        }

        for (int li = 0; li < size; ++li) {
            int ell = m + li;
            
            // Get R_ml values (internal, kind 1) - ONCE per row
            Rcpp::List Rml_i1 = Rcpp::as<Rcpp::List>(Rcpp::as<Rcpp::List>(Rmn_i1_matrix[m])[li]);
            double R_ml_1_val = Rcpp::as<double>(Rml_i1["value"]);
            double R_ml_1_d   = Rcpp::as<double>(Rml_i1["derivative"]);
            
            for (int ni = 0; ni < size; ++ni) {
                int n = m + ni;
                double snm_n = smn_matrix(m, n);
                
                if (snm_n == 0.0 || std::isnan(snm_n)) {
                    continue; // Skip this element
                }

                // Get R_mn values (external, kind 1)
                Rcpp::List Rmn_e1 = Rcpp::as<Rcpp::List>(Rcpp::as<Rcpp::List>(Rmn_e1_matrix[m])[ni]);
                double Rmn_i1_val = Rcpp::as<double>(Rmn_e1["value"]);
                double Rmn_i1_d   = Rcpp::as<double>(Rmn_e1["derivative"]);

                double factor = (density_body * R_ml_1_val) / (density_sw * R_ml_1_d);
                double val = Rmn_i1_val - factor * Rmn_i1_d;

                mat(li, ni).r = val;
                mat(li, ni).i = 0.0;
            }
        }
        result[m] = mat;
    }
    return result;
}

Rcpp::List boundary_coupling_scattering_matrix(
    int m_max,
    int n_max,
    double density_body,
    double density_sw,
    Rcpp::List Rmn_e3_matrix,
    Rcpp::List Rmn_i1_matrix,
    Rcpp::NumericMatrix smn_matrix
) {
    // Equivalent to E_{n,l}^{m(3)} related to Eq. 4 from Furusawa (1988)
    Rcpp::List result(m_max + 1);
    for (int m = 0; m <= m_max; ++m) {
        int size = n_max - m + 1;
        Rcpp::ComplexMatrix mat(size, size);
        
        // Initialize to zero
        for (int li = 0; li < size; ++li) {
            for (int ni = 0; ni < size; ++ni) {
                mat(li, ni).r = 0.0;
                mat(li, ni).i = 0.0;
            }
        }
        
        for (int li = 0; li < size; ++li) {
            int ell = m + li;
            
            // Get R_ml values (internal, kind 1) - ONCE per row
            Rcpp::List Rml_i1 = Rcpp::as<Rcpp::List>(Rcpp::as<Rcpp::List>(Rmn_i1_matrix[m])[li]);
            double R_ml_1_val = Rcpp::as<double>(Rml_i1["value"]);
            double R_ml_1_d   = Rcpp::as<double>(Rml_i1["derivative"]);
            
            for (int ni = 0; ni < size; ++ni) {
                int n = m + ni;
                double snm_n = smn_matrix(m, n);
                
                if (snm_n == 0.0 || std::isnan(snm_n)) {
                    continue; // Skip this element (already zero)
                }

                // Get R_mn values (external, kind 3)
                Rcpp::List Rmn_e3 = Rcpp::as<Rcpp::List>(Rcpp::as<Rcpp::List>(Rmn_e3_matrix[m])[ni]);
                std::complex<double> Rmn_i3_val = Rcpp::as<std::complex<double>>(Rmn_e3["value"]);
                std::complex<double> Rmn_i3_d   = Rcpp::as<std::complex<double>>(Rmn_e3["derivative"]);
                
                double factor = (density_body * R_ml_1_val) / (density_sw * R_ml_1_d);
                std::complex<double> val = Rmn_i3_val - factor * Rmn_i3_d;
                
                mat(li, ni).r = val.real();
                mat(li, ni).i = val.imag();
            }
        }
        result[m] = mat;
    }
    return result;
}

// [[Rcpp::export]]
Rcpp::List boundary_coupling_regular_diag(
    int m_max,
    int n_max,
    double density_body,
    double density_sw,
    Rcpp::List Rmn_e1_matrix,
    Rcpp::List Rmn_i1_matrix
) {
    // Expect Rmn_*_matrix to be Lists with "value" and "derivative" NumericMatrix of shape (m_max+1) x (n_max+1)
    Rcpp::NumericMatrix e1_val_mat = Rcpp::as<Rcpp::NumericMatrix>(Rmn_e1_matrix["value"]);
    Rcpp::NumericMatrix e1_der_mat = Rcpp::as<Rcpp::NumericMatrix>(Rmn_e1_matrix["derivative"]);
    Rcpp::NumericMatrix i1_val_mat = Rcpp::as<Rcpp::NumericMatrix>(Rmn_i1_matrix["value"]);
    Rcpp::NumericMatrix i1_der_mat = Rcpp::as<Rcpp::NumericMatrix>(Rmn_i1_matrix["derivative"]);

    Rcpp::List result(m_max + 1);
    for (int m = 0; m <= m_max; ++m) {
        int size = n_max - m + 1;
        Rcpp::NumericVector diag_vec(size, NA_REAL);

        for (int li = 0; li < size; ++li) {
            int ell = m + li;
            // internal (row) values - R_ml from i1 matrices at (m, ell)
            double R_ml_1_val = i1_val_mat(m, ell);
            double R_ml_1_d   = i1_der_mat(m, ell);

            // diagonal index n == ell
            int n = ell;
            double Rmn_i1_val = e1_val_mat(m, n);
            double Rmn_i1_d   = e1_der_mat(m, n);

            // If external value is NA -> leave NA
            if (Rcpp::NumericVector::is_na(Rmn_i1_val)) {
                diag_vec[li] = NA_REAL;
                continue;
            }

            // Compute factor safely
            double factor = 0.0;
            if (!Rcpp::NumericVector::is_na(R_ml_1_d) && R_ml_1_d != 0.0) {
                factor = (density_body * R_ml_1_val) / (density_sw * R_ml_1_d);
            }

            double val = Rmn_i1_val - factor * Rmn_i1_d;
            diag_vec[li] = val;
        }
        result[m] = diag_vec;
    }
    return result;
}

// [[Rcpp::export]]
Rcpp::List boundary_coupling_scattering_diag(
    int m_max,
    int n_max,
    double density_body,
    double density_sw,
    Rcpp::List Rmn_e3_matrix,
    Rcpp::List Rmn_i1_matrix
) {
    // Expect Rmn_e3_matrix to be a List with "value" and "derivative" ComplexMatrix,
    // and Rmn_i1_matrix as List with "value"/"derivative" NumericMatrix.
    Rcpp::ComplexMatrix e3_val_mat = Rcpp::as<Rcpp::ComplexMatrix>(Rmn_e3_matrix["value"]);
    Rcpp::ComplexMatrix e3_der_mat = Rcpp::as<Rcpp::ComplexMatrix>(Rmn_e3_matrix["derivative"]);
    Rcpp::NumericMatrix i1_val_mat = Rcpp::as<Rcpp::NumericMatrix>(Rmn_i1_matrix["value"]);
    Rcpp::NumericMatrix i1_der_mat = Rcpp::as<Rcpp::NumericMatrix>(Rmn_i1_matrix["derivative"]);

    Rcpp::List result(m_max + 1);
    for (int m = 0; m <= m_max; ++m) {
        int size = n_max - m + 1;
        Rcpp::ComplexVector diag_vec(size);

        for (int li = 0; li < size; ++li) {
            int ell = m + li;
            // internal (row) values - R_ml from i1 matrices at (m, ell)
            double R_ml_1_val = i1_val_mat(m, ell);
            double R_ml_1_d   = i1_der_mat(m, ell);

            // diagonal index n == ell
            int n = ell;
            std::complex<double> Rmn_i3_val = Rcpp::as<std::complex<double>>(Rcpp::wrap(e3_val_mat(m, n)));
            std::complex<double> Rmn_i3_d = Rcpp::as<std::complex<double>>(Rcpp::wrap(e3_der_mat(m, n)));

            // Compute factor safely
            double factor = 0.0;
            if (!Rcpp::NumericVector::is_na(R_ml_1_d) && R_ml_1_d != 0.0) {
                factor = (density_body * R_ml_1_val) / (density_sw * R_ml_1_d);
            }

            std::complex<double> val = Rmn_i3_val - factor * Rmn_i3_d;

            diag_vec[li].r = val.real();
            diag_vec[li].i = val.imag();
        }
        result[m] = diag_vec;
    }
    return result;
}

Rcpp::List field_kernel_regular_matrix(
    int m_max,
    int n_max,
    Rcpp::NumericMatrix smn_matrix,
    Rcpp::List expansion_matrix,
    Rcpp::List E1_coupling
) {
    // Equivalent to K_{n,l}^{m(1)} from Furusawa (1988)
    Rcpp::List result(m_max + 1);
    
    for (int m = 0; m <= m_max; ++m) {
        int size = n_max - m + 1;
        Rcpp::ComplexMatrix mat(size, size);
        
        for (int li = 0; li < size; ++li) {
            int l = m + li;
            
            for (int ni = 0; ni < size; ++ni) {
                int n = m + ni;
                double snm_val = smn_matrix(m, n);
                
                if (snm_val == 0.0 || std::isnan(snm_val)) {
                    mat(li, ni).r = 0.0;
                    mat(li, ni).i = 0.0;
                    continue;
                }
                
                // Get values from matrices
                Rcpp::NumericMatrix alpha_mat = Rcpp::as<Rcpp::NumericMatrix>(expansion_matrix[m]);
                double alpha = alpha_mat(li, ni);
                
                Rcpp::ComplexMatrix E1_mat = Rcpp::as<Rcpp::ComplexMatrix>(E1_coupling[m]);
                std::complex<double> E1_nl(E1_mat(li, ni).r, E1_mat(li, ni).i);
                
                // K1 = j^n * S_mn * alpha * E1
                std::complex<double> jn = std::pow(std::complex<double>(0.0, 1.0), n);
                std::complex<double> K1_val = jn * snm_val * alpha * E1_nl;
                
                mat(li, ni).r = K1_val.real();
                mat(li, ni).i = K1_val.imag();
            }
        }
        result[m] = mat;
    }
    return result;
}

Rcpp::List field_kernel_scattering_matrix(
    int m_max,
    int n_max,
    Rcpp::NumericMatrix smn_matrix,
    Rcpp::List expansion_matrix,
    Rcpp::List E3_coupling
) {
    // Equivalent to K_{n,l}^{m(3)} from Furusawa (1988)
    Rcpp::List result(m_max + 1);
    
    for (int m = 0; m <= m_max; ++m) {
        int size = n_max - m + 1;
        Rcpp::ComplexMatrix mat(size, size);

        // Initialize entire matrix to zero
        for (int li = 0; li < size; ++li) {
            for (int ni = 0; ni < size; ++ni) {
                mat(li, ni).r = 0.0;
                mat(li, ni).i = 0.0;
            }
        }
        
        for (int ni = 0; ni < size; ++ni) {
            int n = m + ni;
            double snm_val = smn_matrix(m, n);
                       
            for (int li = 0; li < size; ++li) {
                int l = m + li;
                
                Rcpp::NumericMatrix alpha_mat = Rcpp::as<Rcpp::NumericMatrix>(expansion_matrix[m]);
                double alpha = alpha_mat(li, ni);
                
                Rcpp::ComplexMatrix E3_mat = Rcpp::as<Rcpp::ComplexMatrix>(E3_coupling[m]);
                std::complex<double> E3_nl(E3_mat(li, ni).r, E3_mat(li, ni).i);
                
                std::complex<double> jn = std::pow(std::complex<double>(0.0, 1.0), n);
                std::complex<double> K3_val = jn * snm_val * alpha * E3_nl;
                
                mat(li, ni).r = K3_val.real();
                mat(li, ni).i = K3_val.imag();
            }
        }
        result[m] = mat;
    }
    return result;
}

// [[Rcpp::export]]
Rcpp::List liquid_spheroidal_kernel_matrices(
    int m_max,
    int n_max,
    double chi_sw,
    double chi_body,
    double xi,
    double density_body,
    double density_sw,
    Rcpp::NumericVector nodes,
    Rcpp::NumericVector weights,
    Rcpp::NumericMatrix smn_matrix
) {
    // Spheroidal expansion kernel (no snm check - computes all)
    Rcpp::List expansion_matrix = liquid_spheroidal_expansion_matrix(
        m_max, n_max, chi_sw, chi_body, nodes, weights
    );

    // Radial function matrices (skip if snm = 0)
    Rcpp::List radial_external_kind1 = radial_external_incoming_matrices(m_max, n_max, chi_sw, xi, smn_matrix);
    Rcpp::List radial_external_kind3 = radial_external_scattering_matrices(m_max, n_max, chi_sw, xi, smn_matrix);
    Rcpp::List radial_internal_kind1 = radial_internal_incoming_matrices(m_max, n_max, chi_body, xi);

    // Coupling matrices (skip rows if snm = 0)
    Rcpp::List E1_coupling = boundary_coupling_regular_matrix(
        m_max, n_max, density_body, density_sw, radial_external_kind1, radial_internal_kind1, smn_matrix
    );
    Rcpp::List E3_coupling = boundary_coupling_scattering_matrix(
        m_max, n_max, density_body, density_sw, radial_external_kind3, radial_internal_kind1, smn_matrix
    );

    // Kernel matrices (skip elements if snm = 0)
    Rcpp::List K1_kernel = field_kernel_regular_matrix(
        m_max, n_max, smn_matrix, expansion_matrix, E1_coupling
    );
    Rcpp::List K3_kernel = field_kernel_scattering_matrix(
        m_max, n_max, smn_matrix, expansion_matrix, E3_coupling
    );

    return Rcpp::List::create(
        Rcpp::Named("K1_kernel") = K1_kernel,
        Rcpp::Named("K3_kernel") = K3_kernel
    );
}

// [[Rcpp::export]]
Rcpp::List liquid_spheroidal_simplified_expansion(
    int m_max,
    int n_max,
    double chi_sw,
    double chi_body,
    double xi,
    double density_sw,
    double density_body
) {

    // Compute radial lists
    Rcpp::List radial_external_kind1 = radial_external_incoming_matrix(m_max, n_max, chi_sw, xi);
    Rcpp::List radial_external_kind3 = radial_external_scattering_matrix(m_max, n_max, chi_sw, xi);
    Rcpp::List radial_internal_kind1 = radial_internal_incoming_matrix(m_max, n_max, chi_body, xi);

    // Compute diagonal E1 (real) and E3 (complex)
    Rcpp::List E1_diag = boundary_coupling_regular_diag(
        m_max, n_max, density_body, density_sw, radial_external_kind1, radial_internal_kind1
    );
    Rcpp::List E3_diag = boundary_coupling_scattering_diag(
        m_max, n_max, density_body, density_sw, radial_external_kind3, radial_internal_kind1
    );

    // Build Amn = -E3 / E1 as a full (m_max+1) x (n_max+1) complex matrix
    Rcpp::ComplexMatrix amn_mat(m_max + 1, n_max + 1);

    // initialize to NA for both real and imag parts
    for (int mm = 0; mm <= m_max; ++mm) {
        for (int nn = 0; nn <= n_max; ++nn) {
            amn_mat(mm, nn).r = NA_REAL;
            amn_mat(mm, nn).i = NA_REAL;
        }
    }

    for (int m = 0; m <= m_max; ++m) {
        // E1_diag and E3_diag return per-m diagonal vectors (length = n_max - m + 1)
        Rcpp::NumericVector e1_vec = Rcpp::as<Rcpp::NumericVector>(E1_diag[m]);
        Rcpp::ComplexVector e3_vec = Rcpp::as<Rcpp::ComplexVector>(E3_diag[m]);

        int size = std::max(0, static_cast<int>(e1_vec.size()));

        for (int i = 0; i < size; ++i) {
            int n = m + i;
            double e1v = e1_vec[i];

            // If e1 is NA or zero -> leave NA (already set)
            if (Rcpp::NumericVector::is_na(e1v) || e1v == 0.0) {
                continue;
            }

            // Extract real/imag components from ComplexVector element and build std::complex
            std::complex<double> e3v = Rcpp::as<std::complex<double>>(Rcpp::wrap(e3_vec[i]));

            std::complex<double> res = - e1v / e3v;
            amn_mat(m, n).r = res.real();
            amn_mat(m, n).i = res.imag();
        }
    }

    return Rcpp::List::create(Rcpp::Named("amn") = amn_mat);
}