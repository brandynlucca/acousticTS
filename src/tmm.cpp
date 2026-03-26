// [[Rcpp::depends(RcppArmadillo, BH)]]
#include <RcppArmadillo.h>
#include <boost/math/special_functions/legendre.hpp>
#include <vector>
#include <complex>
#include <cmath>
#include "bessel_helpers.h"

using namespace Rcpp;

void gauss_legendre(int n, std::vector<double>& nodes, std::vector<double>& weights);

static const std::complex<double> tmm_i(0.0, 1.0);

// Integer powers of i appear repeatedly in the modal algebra, so evaluate them
// exactly instead of repeatedly calling std::pow on complex values.
static inline std::complex<double> tmm_pow_i(int n) {
    switch (((n % 4) + 4) % 4) {
    case 0: return std::complex<double>(1.0, 0.0);
    case 1: return std::complex<double>(0.0, 1.0);
    case 2: return std::complex<double>(-1.0, 0.0);
    default: return std::complex<double>(0.0, -1.0);
    }
}

// Stable odd double-factorial used in the associated Legendre seed term.
static double tmm_double_factorial_odd(int m) {
    if (m <= 0) {
        return 1.0;
    }
    return std::exp(std::lgamma(2.0 * m + 1.0) -
        m * std::log(2.0) -
        std::lgamma(m + 1.0));
}

// Build P_n^m(mu) for one azimuthal order m across all retained degrees and
// quadrature nodes.
static arma::mat tmm_assoc_legendre_table(int m, int n_max, const arma::vec& mu) {
    int n_terms = n_max - m + 1;
    arma::mat out(mu.n_elem, n_terms, arma::fill::zeros);

    arma::vec p_mm(mu.n_elem, arma::fill::zeros);
    if (m == 0) {
        p_mm.fill(1.0);
    } else {
        double coeff = std::pow(-1.0, m) * tmm_double_factorial_odd(m);
        p_mm = coeff * arma::pow(1.0 - arma::square(mu), 0.5 * static_cast<double>(m));
    }
    out.col(0) = p_mm;

    if (n_terms > 1) {
        out.col(1) = (2.0 * m + 1.0) * mu % p_mm;
    }

    if (n_terms > 2) {
        for (int n = m + 2; n <= n_max; ++n) {
            int j = n - m;
            out.col(j) = (
                (2.0 * n - 1.0) * mu % out.col(j - 1) -
                    (n + m - 1.0) * out.col(j - 2)
            ) / static_cast<double>(n - m);
        }
    }

    return out;
}

// Convert the associated Legendre table into theta-derivatives needed by the
// curved-surface normal derivative.
static arma::mat tmm_assoc_legendre_theta_derivative(int m,
                                                     const arma::ivec& n_seq,
                                                     const arma::vec& mu,
                                                     const arma::mat& p_mat) {
    arma::vec sin_theta = arma::sqrt(arma::clamp(1.0 - arma::square(mu), 0.0, arma::datum::inf));
    arma::mat deriv(p_mat.n_rows, p_mat.n_cols, arma::fill::zeros);

    for (arma::uword j = 0; j < p_mat.n_cols; ++j) {
        int n = n_seq[j];
        arma::vec p_nm = p_mat.col(j);
        arma::vec p_nm1 = (n == m) ? arma::vec(mu.n_elem, arma::fill::zeros) : p_mat.col(j - 1);

        for (arma::uword i = 0; i < mu.n_elem; ++i) {
            if (sin_theta[i] <= 0.0) {
                deriv(i, j) = 0.0;
            } else {
                deriv(i, j) = (n * mu[i] * p_nm[i] - (n + m) * p_nm1[i]) / sin_theta[i];
            }
        }
    }

    return deriv;
}

// Evaluate r(theta) and dr/dtheta for the geometries supported by the
// spherical-coordinate backend.
static void tmm_surface_radius(const std::string& shape,
                               const arma::vec& theta,
                               const NumericVector& shape_values,
                               arma::vec& radius,
                               arma::vec& radius_derivative) {
    if (shape == "Sphere") {
        radius = arma::vec(theta.n_elem, arma::fill::ones) * shape_values[0];
        radius_derivative = arma::vec(theta.n_elem, arma::fill::zeros);
        return;
    }

    if (shape == "ProlateSpheroid") {
        double a = shape_values[0];
        double b = shape_values[1];
        arma::vec cos_theta = arma::cos(theta);
        arma::vec sin_theta = arma::sin(theta);
        arma::vec denom = arma::square(cos_theta) / (a * a) + arma::square(sin_theta) / (b * b);
        radius = 1.0 / arma::sqrt(denom);
        radius_derivative = -sin_theta % cos_theta *
            (1.0 / (b * b) - 1.0 / (a * a)) /
            arma::pow(denom, 1.5);
        return;
    }

    if (shape == "OblateSpheroid") {
        double c_axial = shape_values[0];
        double a_equatorial = shape_values[1];
        arma::vec cos_theta = arma::cos(theta);
        arma::vec sin_theta = arma::sin(theta);
        arma::vec denom = arma::square(cos_theta) / (c_axial * c_axial) +
            arma::square(sin_theta) / (a_equatorial * a_equatorial);
        radius = 1.0 / arma::sqrt(denom);
        radius_derivative = -sin_theta % cos_theta *
            (1.0 / (a_equatorial * a_equatorial) - 1.0 / (c_axial * c_axial)) /
            arma::pow(denom, 1.5);
        return;
    }

    if (shape == "Cylinder") {
        double half_length = shape_values[0];
        double cyl_radius = shape_values[1];
        arma::vec cos_theta = arma::cos(theta);
        arma::vec sin_theta = arma::sin(theta);
        arma::vec eps_vec(theta.n_elem);
        eps_vec.fill(std::numeric_limits<double>::epsilon());
        arma::vec r_end = half_length / arma::max(arma::abs(cos_theta), eps_vec);
        arma::vec r_side = cyl_radius / arma::max(arma::abs(sin_theta), eps_vec);
        radius = arma::min(r_end, r_side);
        radius_derivative = arma::vec(theta.n_elem, arma::fill::zeros);
        if (theta.n_elem > 1) {
            radius_derivative[0] = (radius[1] - radius[0]) / (theta[1] - theta[0]);
            radius_derivative[theta.n_elem - 1] =
                (radius[theta.n_elem - 1] - radius[theta.n_elem - 2]) /
                (theta[theta.n_elem - 1] - theta[theta.n_elem - 2]);
            for (arma::uword i = 1; i + 1 < theta.n_elem; ++i) {
                radius_derivative[i] =
                    (radius[i + 1] - radius[i - 1]) /
                    (theta[i + 1] - theta[i - 1]);
            }
        }
        return;
    }

    stop("Unsupported TMM shape geometry.");
}

// Shared collocation-node rule for the spherical-coordinate branch.
static int tmm_collocation_nodes(const std::string& shape,
                                 const std::string& boundary,
                                 int n_terms) {
    if (shape == "Cylinder" && boundary == "fixed_rigid") {
        return std::max(128, 10 * n_terms);
    }
    return std::max(64, 4 * n_terms);
}

// Evaluate one spherical radial-function family for all retained degrees and
// quadrature nodes.
static arma::cx_mat tmm_radial_matrix(const arma::ivec& n_seq,
                                      const arma::vec& argument,
                                      const std::string& fun) {
    arma::cx_mat out(argument.n_elem, n_seq.n_elem);

    for (arma::uword j = 0; j < n_seq.n_elem; ++j) {
        int n = n_seq[j];
        for (arma::uword i = 0; i < argument.n_elem; ++i) {
            double z = argument[i];
            if (fun == "js") {
                out(i, j) = std::complex<double>(js_single_impl(n, z), 0.0);
            } else if (fun == "jsd") {
                out(i, j) = std::complex<double>(js_deriv_single_impl(n, z, 1), 0.0);
            } else if (fun == "hs") {
                out(i, j) = hs_single_impl(n, z);
            } else if (fun == "hsd") {
                out(i, j) = hs_deriv_single_impl(n, z, 1);
            } else {
                stop("Unsupported radial-function request.");
            }
        }
    }

    return out;
}

// Apply the geometric normal derivative on the target surface.
static arma::cx_mat tmm_normal_derivative_matrix(const arma::cx_mat& radial,
                                                 const arma::cx_mat& radial_deriv,
                                                 const arma::mat& angular,
                                                 const arma::mat& angular_theta_deriv,
                                                 double k,
                                                 const arma::vec& radius,
                                                 const arma::vec& radius_derivative) {
    arma::cx_mat out(radial.n_rows, radial.n_cols);

    for (arma::uword i = 0; i < radial.n_rows; ++i) {
        double scale = radius_derivative[i] / (radius[i] * radius[i]);
        for (arma::uword j = 0; j < radial.n_cols; ++j) {
            out(i, j) = k * radial_deriv(i, j) * angular(i, j) -
                radial(i, j) * angular_theta_deriv(i, j) * scale;
        }
    }

    return out;
}

// Incident plane-wave coefficients for one azimuthal block.
static arma::cx_vec tmm_incident_plane_wave_coefficients(int m,
                                                         const arma::ivec& n_seq,
                                                         const arma::vec& p_inc) {
    arma::cx_vec out(n_seq.n_elem);

    for (arma::uword j = 0; j < n_seq.n_elem; ++j) {
        int n = n_seq[j];
        double beta = 1.0;
        if (m > 0) {
            beta = 2.0 * std::exp(
                std::lgamma(static_cast<double>(n - m + 1)) -
                    std::lgamma(static_cast<double>(n + m + 1))
            );
        }
        out[j] = tmm_pow_i(n) * static_cast<double>(2 * n + 1) * beta * p_inc[j];
    }

    return out;
}

// Solve one frequency of the spherical-coordinate TMM branch and return only
// the monostatic backscatter amplitude.
static std::complex<double> tmm_single_frequency_cpp(double frequency,
                                                     double theta_body,
                                                     const std::string& shape,
                                                     const NumericVector& shape_values,
                                                     const std::string& boundary,
                                                     double sound_speed_sw,
                                                     double density_sw,
                                                     double density_body,
                                                     double sound_speed_body,
                                                     int n_max) {
    double mu0 = std::cos(theta_body);
    double k_sw = 2.0 * M_PI * frequency / sound_speed_sw;
    double k_body = std::numeric_limits<double>::quiet_NaN();
    bool penetrable = boundary == "liquid_filled" || boundary == "gas_filled";

    if (penetrable) {
        k_body = 2.0 * M_PI * frequency / sound_speed_body;
    }

    std::complex<double> f_bs(0.0, 0.0);

    for (int m = 0; m <= n_max; ++m) {
        // Axisymmetry decouples the azimuthal orders, so solve one block at a
        // time.
        int n_terms = n_max - m + 1;
        int n_nodes = tmm_collocation_nodes(shape, boundary, n_terms);

        std::vector<double> nodes_vec;
        std::vector<double> weights_vec;
        gauss_legendre(n_nodes, nodes_vec, weights_vec);

        arma::vec mu(nodes_vec);
        arma::vec weights(weights_vec);
        arma::vec theta = arma::acos(mu);

        // Build the collocation geometry and the associated Legendre test
        // basis for this azimuthal block.
        arma::vec radius;
        arma::vec radius_derivative;
        tmm_surface_radius(shape, theta, shape_values, radius, radius_derivative);

        arma::ivec n_seq = arma::regspace<arma::ivec>(m, n_max);
        arma::mat p_mat = tmm_assoc_legendre_table(m, n_max, mu);
        arma::mat dp_dtheta = tmm_assoc_legendre_theta_derivative(m, n_seq, mu, p_mat);

        arma::vec kr_sw = k_sw * radius;
        arma::cx_mat j_sw = tmm_radial_matrix(n_seq, kr_sw, "js");
        arma::cx_mat dj_sw = tmm_radial_matrix(n_seq, kr_sw, "jsd");
        arma::cx_mat h_sw = tmm_radial_matrix(n_seq, kr_sw, "hs");
        arma::cx_mat dh_sw = tmm_radial_matrix(n_seq, kr_sw, "hsd");

        arma::cx_mat reg_normal = tmm_normal_derivative_matrix(
            j_sw, dj_sw, p_mat, dp_dtheta, k_sw, radius, radius_derivative
        );
        arma::cx_mat out_normal = tmm_normal_derivative_matrix(
            h_sw, dh_sw, p_mat, dp_dtheta, k_sw, radius, radius_derivative
        );

        // Assemble the boundary operator in collocation form.
        arma::cx_mat lhs;
        arma::cx_mat rhs;

        if (boundary == "pressure_release") {
            // Pressure-release targets enforce zero total pressure.
            lhs = h_sw % arma::conv_to<arma::cx_mat>::from(p_mat);
            rhs = -(j_sw % arma::conv_to<arma::cx_mat>::from(p_mat));
        } else if (boundary == "fixed_rigid") {
            // Rigid targets enforce zero normal velocity.
            lhs = out_normal;
            rhs = -reg_normal;
        } else {
            // Penetrable targets enforce continuity of pressure and normal
            // velocity across the interface.
            arma::vec kr_body = k_body * radius;
            arma::cx_mat j_body = tmm_radial_matrix(n_seq, kr_body, "js");
            arma::cx_mat dj_body = tmm_radial_matrix(n_seq, kr_body, "jsd");
            arma::cx_mat in_normal = tmm_normal_derivative_matrix(
                j_body, dj_body, p_mat, dp_dtheta, k_body, radius, radius_derivative
            );

            lhs = arma::cx_mat(2 * n_nodes, 2 * n_terms, arma::fill::zeros);
            rhs = arma::cx_mat(2 * n_nodes, n_terms, arma::fill::zeros);

            arma::cx_mat p_cx = arma::conv_to<arma::cx_mat>::from(p_mat);
            lhs.submat(0, 0, n_nodes - 1, n_terms - 1) = h_sw % p_cx;
            lhs.submat(0, n_terms, n_nodes - 1, 2 * n_terms - 1) = -(j_body % p_cx);
            lhs.submat(n_nodes, 0, 2 * n_nodes - 1, n_terms - 1) = out_normal / density_sw;
            lhs.submat(n_nodes, n_terms, 2 * n_nodes - 1, 2 * n_terms - 1) = -in_normal / density_body;

            rhs.submat(0, 0, n_nodes - 1, n_terms - 1) = -(j_sw % p_cx);
            rhs.submat(n_nodes, 0, 2 * n_nodes - 1, n_terms - 1) = -(reg_normal / density_sw);
        }

        // Project the collocation system back onto the retained modal basis.
        arma::cx_mat lhs_proj;
        arma::cx_mat rhs_proj;
        arma::vec surface_weight = radius % arma::sqrt(arma::square(radius) + arma::square(radius_derivative));
        arma::vec projector_weight = weights % surface_weight;
        arma::mat weighted_test = p_mat;
        weighted_test.each_col() %= projector_weight;
        arma::cx_mat wt_h = arma::conv_to<arma::cx_mat>::from(weighted_test.t());

        if (!penetrable) {
            // Rigid and pressure-release cases project to one square block.
            lhs_proj = wt_h * lhs;
            rhs_proj = wt_h * rhs;
        } else {
            // Penetrable cases keep coupled exterior/interior coefficient
            // blocks after projection.
            lhs_proj = arma::cx_mat(2 * n_terms, 2 * n_terms, arma::fill::zeros);
            rhs_proj = arma::cx_mat(2 * n_terms, n_terms, arma::fill::zeros);

            lhs_proj.submat(0, 0, n_terms - 1, n_terms - 1) =
                wt_h * lhs.submat(0, 0, n_nodes - 1, n_terms - 1);
            lhs_proj.submat(0, n_terms, n_terms - 1, 2 * n_terms - 1) =
                wt_h * lhs.submat(0, n_terms, n_nodes - 1, 2 * n_terms - 1);
            lhs_proj.submat(n_terms, 0, 2 * n_terms - 1, n_terms - 1) =
                wt_h * lhs.submat(n_nodes, 0, 2 * n_nodes - 1, n_terms - 1);
            lhs_proj.submat(n_terms, n_terms, 2 * n_terms - 1, 2 * n_terms - 1) =
                wt_h * lhs.submat(n_nodes, n_terms, 2 * n_nodes - 1, 2 * n_terms - 1);

            rhs_proj.submat(0, 0, n_terms - 1, n_terms - 1) =
                wt_h * rhs.submat(0, 0, n_nodes - 1, n_terms - 1);
            rhs_proj.submat(n_terms, 0, 2 * n_terms - 1, n_terms - 1) =
                wt_h * rhs.submat(n_nodes, 0, 2 * n_nodes - 1, n_terms - 1);
        }

        // Solve the projected block, falling back to progressively more
        // forgiving solvers if the fast path is not stable enough.
        arma::cx_mat solution;
        bool ok = arma::solve(solution, lhs_proj, rhs_proj, arma::solve_opts::fast);
        if (!ok || !solution.is_finite()) {
            ok = arma::solve(solution, lhs_proj, rhs_proj);
        }
        if (!ok || !solution.is_finite()) {
            solution = arma::pinv(lhs_proj) * rhs_proj;
            ok = solution.is_finite();
        }
        if (!solution.is_finite()) {
            stop("TMM was unable to solve the projected boundary system.");
        }

        arma::cx_mat t_block = penetrable ? solution.rows(0, n_terms - 1) : solution;

        // Convert the solved block into the outgoing monostatic amplitude.
        arma::vec mu0_vec(1);
        mu0_vec[0] = mu0;
        arma::mat p_inc_mat = tmm_assoc_legendre_table(m, n_max, mu0_vec);
        arma::vec p_inc = p_inc_mat.row(0).t();
        arma::cx_vec a_inc = tmm_incident_plane_wave_coefficients(m, n_seq, p_inc);
        arma::cx_vec coeffs = t_block * a_inc;

        // Evaluate the outgoing expansion in the monostatic receive direction.
        for (arma::uword j = 0; j < coeffs.n_elem; ++j) {
            int n = n_seq[j];
            f_bs += coeffs[j] * (std::pow(-1.0, n) * tmm_pow_i(-(n + 1)) * p_inc[j] / k_sw);
        }
    }

    return f_bs;
}

// Vectorized frequency wrapper used by the R-side spherical TMM path.
// [[Rcpp::export]]
Rcpp::ComplexVector tmm_backscatter_cpp(Rcpp::NumericVector frequency,
                                        double theta_body,
                                        std::string shape,
                                        Rcpp::NumericVector shape_values,
                                        std::string boundary,
                                        double sound_speed_sw,
                                        double density_sw,
                                        double density_body,
                                        double sound_speed_body,
                                        Rcpp::IntegerVector n_max) {
    int n_freq = frequency.size();
    if (n_max.size() != n_freq) {
        stop("'n_max' must match the frequency vector length.");
    }

    Rcpp::ComplexVector out(n_freq);
    for (int i = 0; i < n_freq; ++i) {
        std::complex<double> f_bs = tmm_single_frequency_cpp(
            frequency[i],
            theta_body,
            shape,
            shape_values,
            boundary,
            sound_speed_sw,
            density_sw,
            density_body,
            sound_speed_body,
            n_max[i]
        );
        out[i] = to_Rcomplex(f_bs);
    }

    return out;
}
