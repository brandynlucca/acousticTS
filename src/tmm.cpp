// [[Rcpp::depends(RcppArmadillo, BH)]]
#include <RcppArmadillo.h>
#include <boost/math/special_functions/legendre.hpp>
#include <vector>
#include <complex>
#include <cmath>
#include "bessel_helpers.h"

using namespace Rcpp;

void gauss_legendre(int n, std::vector<double>& nodes, std::vector<double>& weights);

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
    int max_n = n_seq.max();

    for (arma::uword i = 0; i < argument.n_elem; ++i) {
        double z = argument[i];
        if (fun == "js") {
            std::vector<double> seq = js_sequence_miller_impl(max_n, z);
            for (arma::uword j = 0; j < n_seq.n_elem; ++j) {
                out(i, j) = std::complex<double>(seq[n_seq[j]], 0.0);
            }
        } else if (fun == "jsd") {
            std::vector<double> seq = js_deriv_sequence_impl(max_n, z);
            for (arma::uword j = 0; j < n_seq.n_elem; ++j) {
                out(i, j) = std::complex<double>(seq[n_seq[j]], 0.0);
            }
        } else if (fun == "hs") {
            std::vector<std::complex<double>> seq = hs_sequence_impl(max_n, z);
            for (arma::uword j = 0; j < n_seq.n_elem; ++j) {
                out(i, j) = seq[n_seq[j]];
            }
        } else if (fun == "hsd") {
            std::vector<std::complex<double>> seq = hs_deriv_sequence_impl(max_n, z);
            for (arma::uword j = 0; j < n_seq.n_elem; ++j) {
                out(i, j) = seq[n_seq[j]];
            }
        } else {
            stop("Unsupported radial-function request.");
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

struct TmmSphericalBlockSetup {
    int m;
    int n_max;
    int n_terms;
    int n_nodes;
    arma::ivec n_seq;
    arma::mat p_mat;
    arma::mat dp_dtheta;
    arma::cx_mat p_cx;
    arma::vec radius;
    arma::vec radius_derivative;
    arma::cx_mat projector;
    arma::vec p_inc;
    arma::cx_vec a_inc;
};

// Precompute the frequency-independent collocation and projection state for
// one spherical-coordinate azimuthal block.
static TmmSphericalBlockSetup tmm_spherical_block_setup_cpp(int m,
                                                            int n_max,
                                                            double mu0,
                                                            const std::string& shape,
                                                            const NumericVector& shape_values,
                                                            const std::string& boundary,
                                                            bool penetrable) {
    TmmSphericalBlockSetup setup;
    setup.m = m;
    setup.n_max = n_max;
    setup.n_terms = n_max - m + 1;
    setup.n_nodes = tmm_collocation_nodes(shape, boundary, setup.n_terms);

    std::vector<double> nodes_vec;
    std::vector<double> weights_vec;
    gauss_legendre(setup.n_nodes, nodes_vec, weights_vec);

    arma::vec mu(nodes_vec);
    arma::vec weights(weights_vec);
    arma::vec theta = arma::acos(mu);

    tmm_surface_radius(
        shape,
        theta,
        shape_values,
        setup.radius,
        setup.radius_derivative
    );

    setup.n_seq = arma::regspace<arma::ivec>(m, n_max);
    setup.p_mat = tmm_assoc_legendre_table(m, n_max, mu);
    setup.dp_dtheta = tmm_assoc_legendre_theta_derivative(
        m,
        setup.n_seq,
        mu,
        setup.p_mat
    );
    setup.p_cx = arma::conv_to<arma::cx_mat>::from(setup.p_mat);

    arma::vec surface_weight = setup.radius %
        arma::sqrt(arma::square(setup.radius) + arma::square(setup.radius_derivative));
    arma::vec projector_weight = weights % surface_weight;
    arma::mat weighted_test = setup.p_mat;
    weighted_test.each_col() %= projector_weight;
    arma::cx_mat wt_h = arma::conv_to<arma::cx_mat>::from(weighted_test.t());

    if (!penetrable) {
        setup.projector = wt_h;
    } else {
        setup.projector = arma::cx_mat(
            2 * setup.n_terms,
            2 * setup.n_nodes,
            arma::fill::zeros
        );
        setup.projector.submat(
            0,
            0,
            setup.n_terms - 1,
            setup.n_nodes - 1
        ) = wt_h;
        setup.projector.submat(
            setup.n_terms,
            setup.n_nodes,
            2 * setup.n_terms - 1,
            2 * setup.n_nodes - 1
        ) = wt_h;
    }

    arma::vec mu0_vec(1);
    mu0_vec[0] = mu0;
    arma::mat p_inc_mat = tmm_assoc_legendre_table(m, n_max, mu0_vec);
    setup.p_inc = p_inc_mat.row(0).t();
    setup.a_inc = tmm_incident_plane_wave_coefficients(
        m,
        setup.n_seq,
        setup.p_inc
    );

    return setup;
}

// Solve one frequency-dependent radial block using a cached spherical setup.
static arma::cx_mat tmm_solve_spherical_block_cpp(const TmmSphericalBlockSetup& setup,
                                                  const std::string& boundary,
                                                  bool penetrable,
                                                  double k_sw,
                                                  double k_body,
                                                  double density_sw,
                                                  double density_body) {
    arma::vec kr_sw = k_sw * setup.radius;
    arma::cx_mat j_sw = tmm_radial_matrix(setup.n_seq, kr_sw, "js");
    arma::cx_mat dj_sw = tmm_radial_matrix(setup.n_seq, kr_sw, "jsd");
    arma::cx_mat h_sw = tmm_radial_matrix(setup.n_seq, kr_sw, "hs");
    arma::cx_mat dh_sw = tmm_radial_matrix(setup.n_seq, kr_sw, "hsd");

    arma::cx_mat reg_normal = tmm_normal_derivative_matrix(
        j_sw,
        dj_sw,
        setup.p_mat,
        setup.dp_dtheta,
        k_sw,
        setup.radius,
        setup.radius_derivative
    );
    arma::cx_mat out_normal = tmm_normal_derivative_matrix(
        h_sw,
        dh_sw,
        setup.p_mat,
        setup.dp_dtheta,
        k_sw,
        setup.radius,
        setup.radius_derivative
    );

    arma::cx_mat lhs;
    arma::cx_mat rhs;

    if (boundary == "pressure_release") {
        lhs = h_sw % setup.p_cx;
        rhs = -(j_sw % setup.p_cx);
    } else if (boundary == "fixed_rigid") {
        lhs = out_normal;
        rhs = -reg_normal;
    } else {
        arma::vec kr_body = k_body * setup.radius;
        arma::cx_mat j_body = tmm_radial_matrix(setup.n_seq, kr_body, "js");
        arma::cx_mat dj_body = tmm_radial_matrix(setup.n_seq, kr_body, "jsd");
        arma::cx_mat in_normal = tmm_normal_derivative_matrix(
            j_body,
            dj_body,
            setup.p_mat,
            setup.dp_dtheta,
            k_body,
            setup.radius,
            setup.radius_derivative
        );

        lhs = arma::cx_mat(
            2 * setup.n_nodes,
            2 * setup.n_terms,
            arma::fill::zeros
        );
        rhs = arma::cx_mat(
            2 * setup.n_nodes,
            setup.n_terms,
            arma::fill::zeros
        );

        lhs.submat(0, 0, setup.n_nodes - 1, setup.n_terms - 1) =
            h_sw % setup.p_cx;
        lhs.submat(0, setup.n_terms, setup.n_nodes - 1, 2 * setup.n_terms - 1) =
            -(j_body % setup.p_cx);
        lhs.submat(setup.n_nodes, 0, 2 * setup.n_nodes - 1, setup.n_terms - 1) =
            out_normal / density_sw;
        lhs.submat(
            setup.n_nodes,
            setup.n_terms,
            2 * setup.n_nodes - 1,
            2 * setup.n_terms - 1
        ) = -in_normal / density_body;

        rhs.submat(0, 0, setup.n_nodes - 1, setup.n_terms - 1) =
            -(j_sw % setup.p_cx);
        rhs.submat(setup.n_nodes, 0, 2 * setup.n_nodes - 1, setup.n_terms - 1) =
            -(reg_normal / density_sw);
    }

    arma::cx_mat lhs_proj = setup.projector * lhs;
    arma::cx_mat rhs_proj = setup.projector * rhs;

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

    return penetrable ? solution.rows(0, setup.n_terms - 1) : solution;
}

// Convert one solved block into its monostatic backscatter contribution.
static std::complex<double> tmm_backscatter_from_solved_block_cpp(
    const TmmSphericalBlockSetup& setup,
    const arma::cx_mat& t_block,
    double k_sw) {
    arma::cx_vec coeffs = t_block * setup.a_inc;
    std::complex<double> f_bs(0.0, 0.0);

    for (arma::uword j = 0; j < coeffs.n_elem; ++j) {
        int n = setup.n_seq[j];
        f_bs += coeffs[j] *
            (std::pow(-1.0, n) * tmm_pow_i(-(n + 1)) * setup.p_inc[j] / k_sw);
    }

    return f_bs;
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
        TmmSphericalBlockSetup setup = tmm_spherical_block_setup_cpp(
            m,
            n_max,
            mu0,
            shape,
            shape_values,
            boundary,
            penetrable
        );
        arma::cx_mat t_block = tmm_solve_spherical_block_cpp(
            setup,
            boundary,
            penetrable,
            k_sw,
            k_body,
            density_sw,
            density_body
        );
        f_bs += tmm_backscatter_from_solved_block_cpp(setup, t_block, k_sw);
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
    std::vector<int> processed(n_freq, 0);
    double mu0 = std::cos(theta_body);
    bool penetrable = boundary == "liquid_filled" || boundary == "gas_filled";

    for (int start = 0; start < n_freq; ++start) {
        if (processed[start]) {
            continue;
        }

        int bucket_n_max = n_max[start];
        std::vector<int> bucket_indices;
        for (int i = start; i < n_freq; ++i) {
            if (!processed[i] && n_max[i] == bucket_n_max) {
                processed[i] = 1;
                bucket_indices.push_back(i);
            }
        }

        std::vector<std::complex<double>> f_bucket(
            bucket_indices.size(),
            std::complex<double>(0.0, 0.0)
        );

        for (int m = 0; m <= bucket_n_max; ++m) {
            TmmSphericalBlockSetup setup = tmm_spherical_block_setup_cpp(
                m,
                bucket_n_max,
                mu0,
                shape,
                shape_values,
                boundary,
                penetrable
            );

            for (std::size_t j = 0; j < bucket_indices.size(); ++j) {
                int frequency_idx = bucket_indices[j];
                double k_sw = 2.0 * M_PI * frequency[frequency_idx] / sound_speed_sw;
                double k_body = std::numeric_limits<double>::quiet_NaN();
                if (penetrable) {
                    k_body = 2.0 * M_PI * frequency[frequency_idx] / sound_speed_body;
                }

                arma::cx_mat t_block = tmm_solve_spherical_block_cpp(
                    setup,
                    boundary,
                    penetrable,
                    k_sw,
                    k_body,
                    density_sw,
                    density_body
                );
                f_bucket[j] += tmm_backscatter_from_solved_block_cpp(
                    setup,
                    t_block,
                    k_sw
                );
            }
        }

        for (std::size_t j = 0; j < bucket_indices.size(); ++j) {
            out[bucket_indices[j]] = to_Rcomplex(f_bucket[j]);
        }
    }

    return out;
}
