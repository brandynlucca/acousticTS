// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <cmath>
#include <complex>
#include "bessel_helpers.h"

using namespace Rcpp;

static const double vesms_pi = M_PI;
static const double vesms_tol = 1e-12;
static const std::complex<double> vesms_i(0.0, 1.0);

// Solve one dense modal boundary system. Try the fast direct solver first,
// then the full Armadillo solve, and finally fall back to a pseudoinverse if
// the matrix is close to singular.
static arma::cx_vec vesms_solve_system_impl(const arma::cx_mat& M, const arma::cx_vec& F) {
    arma::cx_vec sol;

    bool ok = arma::solve(sol, M, F, arma::solve_opts::fast);
    if (ok && sol.is_finite()) {
        return sol;
    }

    ok = arma::solve(sol, M, F);
    if (ok && sol.is_finite()) {
        return sol;
    }

    arma::cx_mat pinv_M = arma::pinv(M);
    sol = pinv_M * F;
    if (sol.is_finite()) {
        return sol;
    }

    Rcpp::stop(
        "VESMS was unable to solve the modal boundary system. This usually "
        "indicates a numerically singular mode-frequency combination."
    );
}

// Assemble and solve the coupled viscous-layer / elastic-shell / gas-core
// boundary system for one spherical harmonic order. The returned value is the
// exterior scattered-field coefficient for that mode.
static std::complex<double> vesms_single_mode_impl(
    int order,
    double omega,
    double density_sw,
    double density_viscous,
    double density_shell,
    double density_gas,
    double sound_speed_sw,
    double sound_speed_viscous,
    double sound_speed_gas,
    double radius_viscous,
    double radius_shell,
    double radius_gas,
    double shear_viscosity_viscous,
    double bulk_viscosity_viscous,
    double shear_modulus_shell,
    double lambda_shell
) {
    // Convert each physical layer into the compressional and shear
    // wavenumbers used in the modal boundary conditions.
    double k1 = omega / sound_speed_sw;
    double phi2 = bulk_viscosity_viscous + 4.0 * shear_viscosity_viscous / 3.0;
    std::complex<double> kc2 = (omega / sound_speed_viscous) *
        std::pow(
            std::complex<double>(
                1.0,
                -omega * phi2 / (density_viscous * sound_speed_viscous * sound_speed_viscous)
            ),
            -0.5
        );
    std::complex<double> ks2 =
        (1.0 + vesms_i) * std::sqrt(omega * density_viscous / (2.0 * shear_viscosity_viscous));
    double kc3 = omega * std::sqrt(density_shell / (lambda_shell + 2.0 * shear_modulus_shell));
    double ks3 = omega * std::sqrt(density_shell / shear_modulus_shell);
    double k4 = omega / sound_speed_gas;

    // Precompute the spherical Bessel/Hankel values and derivatives at every
    // layer boundary so the linear-system assembly below mirrors the analytic
    // boundary equations directly.
    std::complex<double> j_k1_R2 = js_single_complex_impl(order, k1 * radius_viscous);
    std::complex<double> jd_k1_R2 = js_deriv_single_complex_impl(order, k1 * radius_viscous, 1);
    std::complex<double> h_k1_R2 = hs_single_complex_impl(order, k1 * radius_viscous);
    std::complex<double> hd_k1_R2 = hs_deriv_single_complex_impl(order, k1 * radius_viscous, 1);

    std::complex<double> j_kc2_R2 = js_single_complex_impl(order, kc2 * radius_viscous);
    std::complex<double> jd_kc2_R2 = js_deriv_single_complex_impl(order, kc2 * radius_viscous, 1);
    std::complex<double> jdd_kc2_R2 = js_deriv_single_complex_impl(order, kc2 * radius_viscous, 2);
    std::complex<double> h_kc2_R2 = hs_single_complex_impl(order, kc2 * radius_viscous);
    std::complex<double> hd_kc2_R2 = hs_deriv_single_complex_impl(order, kc2 * radius_viscous, 1);
    std::complex<double> hdd_kc2_R2 = hs_deriv_single_complex_impl(order, kc2 * radius_viscous, 2);

    std::complex<double> j_ks2_R2 = js_single_complex_impl(order, ks2 * radius_viscous);
    std::complex<double> jd_ks2_R2 = js_deriv_single_complex_impl(order, ks2 * radius_viscous, 1);
    std::complex<double> jdd_ks2_R2 = js_deriv_single_complex_impl(order, ks2 * radius_viscous, 2);
    std::complex<double> h_ks2_R2 = hs_single_complex_impl(order, ks2 * radius_viscous);
    std::complex<double> hd_ks2_R2 = hs_deriv_single_complex_impl(order, ks2 * radius_viscous, 1);
    std::complex<double> hdd_ks2_R2 = hs_deriv_single_complex_impl(order, ks2 * radius_viscous, 2);

    std::complex<double> j_kc2_R3 = js_single_complex_impl(order, kc2 * radius_shell);
    std::complex<double> jd_kc2_R3 = js_deriv_single_complex_impl(order, kc2 * radius_shell, 1);
    std::complex<double> jdd_kc2_R3 = js_deriv_single_complex_impl(order, kc2 * radius_shell, 2);
    std::complex<double> h_kc2_R3 = hs_single_complex_impl(order, kc2 * radius_shell);
    std::complex<double> hd_kc2_R3 = hs_deriv_single_complex_impl(order, kc2 * radius_shell, 1);
    std::complex<double> hdd_kc2_R3 = hs_deriv_single_complex_impl(order, kc2 * radius_shell, 2);

    std::complex<double> j_ks2_R3 = js_single_complex_impl(order, ks2 * radius_shell);
    std::complex<double> jd_ks2_R3 = js_deriv_single_complex_impl(order, ks2 * radius_shell, 1);
    std::complex<double> jdd_ks2_R3 = js_deriv_single_complex_impl(order, ks2 * radius_shell, 2);
    std::complex<double> h_ks2_R3 = hs_single_complex_impl(order, ks2 * radius_shell);
    std::complex<double> hd_ks2_R3 = hs_deriv_single_complex_impl(order, ks2 * radius_shell, 1);
    std::complex<double> hdd_ks2_R3 = hs_deriv_single_complex_impl(order, ks2 * radius_shell, 2);

    std::complex<double> j_kc3_R3 = js_single_complex_impl(order, kc3 * radius_shell);
    std::complex<double> jd_kc3_R3 = js_deriv_single_complex_impl(order, kc3 * radius_shell, 1);
    std::complex<double> jdd_kc3_R3 = js_deriv_single_complex_impl(order, kc3 * radius_shell, 2);
    std::complex<double> h_kc3_R3 = hs_single_complex_impl(order, kc3 * radius_shell);
    std::complex<double> hd_kc3_R3 = hs_deriv_single_complex_impl(order, kc3 * radius_shell, 1);
    std::complex<double> hdd_kc3_R3 = hs_deriv_single_complex_impl(order, kc3 * radius_shell, 2);

    std::complex<double> j_ks3_R3 = js_single_complex_impl(order, ks3 * radius_shell);
    std::complex<double> jd_ks3_R3 = js_deriv_single_complex_impl(order, ks3 * radius_shell, 1);
    std::complex<double> jdd_ks3_R3 = js_deriv_single_complex_impl(order, ks3 * radius_shell, 2);
    std::complex<double> h_ks3_R3 = hs_single_complex_impl(order, ks3 * radius_shell);
    std::complex<double> hd_ks3_R3 = hs_deriv_single_complex_impl(order, ks3 * radius_shell, 1);
    std::complex<double> hdd_ks3_R3 = hs_deriv_single_complex_impl(order, ks3 * radius_shell, 2);

    std::complex<double> j_kc3_R4 = js_single_complex_impl(order, kc3 * radius_gas);
    std::complex<double> jd_kc3_R4 = js_deriv_single_complex_impl(order, kc3 * radius_gas, 1);
    std::complex<double> jdd_kc3_R4 = js_deriv_single_complex_impl(order, kc3 * radius_gas, 2);
    std::complex<double> h_kc3_R4 = hs_single_complex_impl(order, kc3 * radius_gas);
    std::complex<double> hd_kc3_R4 = hs_deriv_single_complex_impl(order, kc3 * radius_gas, 1);
    std::complex<double> hdd_kc3_R4 = hs_deriv_single_complex_impl(order, kc3 * radius_gas, 2);

    std::complex<double> j_ks3_R4 = js_single_complex_impl(order, ks3 * radius_gas);
    std::complex<double> jd_ks3_R4 = js_deriv_single_complex_impl(order, ks3 * radius_gas, 1);
    std::complex<double> jdd_ks3_R4 = js_deriv_single_complex_impl(order, ks3 * radius_gas, 2);
    std::complex<double> h_ks3_R4 = hs_single_complex_impl(order, ks3 * radius_gas);
    std::complex<double> hd_ks3_R4 = hs_deriv_single_complex_impl(order, ks3 * radius_gas, 1);
    std::complex<double> hdd_ks3_R4 = hs_deriv_single_complex_impl(order, ks3 * radius_gas, 2);

    std::complex<double> j_k4_R4 = js_single_complex_impl(order, k4 * radius_gas);
    std::complex<double> jd_k4_R4 = js_deriv_single_complex_impl(order, k4 * radius_gas, 1);

    if (order == 0) {
        // The monopole case decouples from the tangential/shear terms, so a
        // reduced 6 x 6 system is sufficient.
        arma::cx_mat M(6, 6, arma::fill::zeros);
        arma::cx_vec F(6, arma::fill::zeros);

        M(0, 0) = k1 * hd_k1_R2;
        M(0, 1) = -kc2 * jd_kc2_R2;
        M(0, 2) = -kc2 * hd_kc2_R2;
        F(0) = -k1 * jd_k1_R2;

        M(1, 0) = -vesms_i * omega * density_sw * h_k1_R2;
        M(1, 1) = -shear_viscosity_viscous * (2.0 * kc2 * kc2 - ks2 * ks2) * j_kc2_R2 -
            2.0 * shear_viscosity_viscous * kc2 * kc2 * jdd_kc2_R2;
        M(1, 2) = -shear_viscosity_viscous * (2.0 * kc2 * kc2 - ks2 * ks2) * h_kc2_R2 -
            2.0 * shear_viscosity_viscous * kc2 * kc2 * hdd_kc2_R2;
        F(1) = vesms_i * omega * density_sw * j_k1_R2;

        M(2, 1) = kc2 * jd_kc2_R3;
        M(2, 2) = kc2 * hd_kc2_R3;
        M(2, 3) = -kc3 * jd_kc3_R3;
        M(2, 4) = -kc3 * hd_kc3_R3;

        M(3, 1) = -vesms_i * omega * shear_viscosity_viscous * (2.0 * kc2 * kc2 - ks2 * ks2) * j_kc2_R3 -
            2.0 * vesms_i * omega * shear_viscosity_viscous * kc2 * kc2 * jdd_kc2_R3;
        M(3, 2) = -vesms_i * omega * shear_viscosity_viscous * (2.0 * kc2 * kc2 - ks2 * ks2) * h_kc2_R3 -
            2.0 * vesms_i * omega * shear_viscosity_viscous * kc2 * kc2 * hdd_kc2_R3;
        M(3, 3) = -shear_modulus_shell * (2.0 * kc3 * kc3 - ks3 * ks3) * j_kc3_R3 -
            2.0 * shear_modulus_shell * kc3 * kc3 * jdd_kc3_R3;
        M(3, 4) = -shear_modulus_shell * (2.0 * kc3 * kc3 - ks3 * ks3) * h_kc3_R3 -
            2.0 * shear_modulus_shell * kc3 * kc3 * hdd_kc3_R3;

        M(4, 3) = kc3 * jd_kc3_R4;
        M(4, 4) = kc3 * hd_kc3_R4;
        M(4, 5) = -k4 * jd_k4_R4;

        M(5, 3) = shear_modulus_shell * (2.0 * kc3 * kc3 - ks3 * ks3) * j_kc3_R4 +
            2.0 * shear_modulus_shell * kc3 * kc3 * jdd_kc3_R4;
        M(5, 4) = shear_modulus_shell * (2.0 * kc3 * kc3 - ks3 * ks3) * h_kc3_R4 +
            2.0 * shear_modulus_shell * kc3 * kc3 * hdd_kc3_R4;
        M(5, 5) = omega * omega * density_gas * j_k4_R4;

        return vesms_solve_system_impl(M, F)(0);
    }

    // Higher orders retain both compressional and shear contributions in the
    // viscous layer and elastic shell, which yields the full 10 x 10 system.
    arma::cx_mat M(10, 10, arma::fill::zeros);
    arma::cx_vec F(10, arma::fill::zeros);
    double mm = -static_cast<double>(order) * (order + 1.0);
    double pp = static_cast<double>(order) * (order + 1.0);

    M(0, 0) = k1 * radius_viscous * hd_k1_R2;
    M(0, 1) = -kc2 * radius_viscous * jd_kc2_R2;
    M(0, 2) = -kc2 * radius_viscous * hd_kc2_R2;
    M(0, 3) = pp * j_ks2_R2;
    M(0, 4) = pp * h_ks2_R2;
    F(0) = -k1 * radius_viscous * std::pow(vesms_i, order) * (2.0 * order + 1.0) * jd_k1_R2;

    M(1, 0) = -vesms_i * omega * density_sw * radius_viscous * radius_viscous * h_k1_R2;
    M(1, 1) = -shear_viscosity_viscous * radius_viscous * radius_viscous * (2.0 * kc2 * kc2 - ks2 * ks2) * j_kc2_R2 -
        2.0 * shear_viscosity_viscous * radius_viscous * radius_viscous * kc2 * kc2 * jdd_kc2_R2;
    M(1, 2) = -shear_viscosity_viscous * radius_viscous * radius_viscous * (2.0 * kc2 * kc2 - ks2 * ks2) * h_kc2_R2 -
        2.0 * shear_viscosity_viscous * radius_viscous * radius_viscous * kc2 * kc2 * hdd_kc2_R2;
    M(1, 3) = -2.0 * shear_viscosity_viscous * mm * ks2 * radius_viscous * jd_ks2_R2 +
        2.0 * shear_viscosity_viscous * mm * j_ks2_R2;
    M(1, 4) = -2.0 * shear_viscosity_viscous * mm * ks2 * radius_viscous * hd_ks2_R2 +
        2.0 * shear_viscosity_viscous * mm * h_ks2_R2;
    F(1) = vesms_i * omega * density_sw * radius_viscous * radius_viscous *
        std::pow(vesms_i, order) * (2.0 * order + 1.0) * j_k1_R2;

    M(2, 1) = 2.0 * kc2 * radius_viscous * jd_kc2_R2 - 2.0 * j_kc2_R2;
    M(2, 2) = 2.0 * kc2 * radius_viscous * hd_kc2_R2 - 2.0 * h_kc2_R2;
    M(2, 3) = -ks2 * ks2 * radius_viscous * radius_viscous * jdd_ks2_R2 + 2.0 * j_ks2_R2 + mm * j_ks2_R2;
    M(2, 4) = -ks2 * ks2 * radius_viscous * radius_viscous * hdd_ks2_R2 + 2.0 * h_ks2_R2 + mm * h_ks2_R2;

    M(3, 1) = kc2 * radius_shell * jd_kc2_R3;
    M(3, 2) = kc2 * radius_shell * hd_kc2_R3;
    M(3, 3) = mm * j_ks2_R3;
    M(3, 4) = mm * h_ks2_R3;
    M(3, 5) = -kc3 * radius_shell * jd_kc3_R3;
    M(3, 6) = -kc3 * radius_shell * hd_kc3_R3;
    M(3, 7) = pp * j_ks3_R3;
    M(3, 8) = pp * h_ks3_R3;

    M(4, 1) = -vesms_i * omega * shear_viscosity_viscous * (2.0 * kc2 * kc2 - ks2 * ks2) * radius_shell * radius_shell * j_kc2_R3 -
        vesms_i * omega * 2.0 * shear_viscosity_viscous * kc2 * kc2 * radius_shell * radius_shell * jdd_kc2_R3;
    M(4, 2) = -vesms_i * omega * shear_viscosity_viscous * (2.0 * kc2 * kc2 - ks2 * ks2) * radius_shell * radius_shell * h_kc2_R3 -
        vesms_i * omega * 2.0 * shear_viscosity_viscous * kc2 * kc2 * radius_shell * radius_shell * hdd_kc2_R3;
    M(4, 3) = -vesms_i * omega * 2.0 * shear_viscosity_viscous * radius_shell * mm * ks2 * jd_ks2_R3 +
        vesms_i * omega * 2.0 * shear_viscosity_viscous * mm * j_ks2_R3;
    M(4, 4) = -vesms_i * omega * 2.0 * shear_viscosity_viscous * radius_shell * mm * ks2 * hd_ks2_R3 +
        vesms_i * omega * 2.0 * shear_viscosity_viscous * mm * h_ks2_R3;
    M(4, 5) = -shear_modulus_shell * (2.0 * kc3 * kc3 - ks3 * ks3) * radius_shell * radius_shell * j_kc3_R3 -
        2.0 * shear_modulus_shell * kc3 * kc3 * radius_shell * radius_shell * jdd_kc3_R3;
    M(4, 6) = -shear_modulus_shell * (2.0 * kc3 * kc3 - ks3 * ks3) * radius_shell * radius_shell * h_kc3_R3 -
        2.0 * shear_modulus_shell * kc3 * kc3 * radius_shell * radius_shell * hdd_kc3_R3;
    M(4, 7) = -2.0 * shear_modulus_shell * mm * ks3 * radius_shell * jd_ks3_R3 +
        2.0 * shear_modulus_shell * mm * j_ks3_R3;
    M(4, 8) = -2.0 * shear_modulus_shell * mm * ks3 * radius_shell * hd_ks3_R3 +
        2.0 * shear_modulus_shell * mm * h_ks3_R3;

    M(5, 1) = -vesms_i * omega * 2.0 * shear_viscosity_viscous * kc2 * radius_shell * jd_kc2_R3 +
        vesms_i * omega * 2.0 * shear_viscosity_viscous * j_kc2_R3;
    M(5, 2) = -vesms_i * omega * 2.0 * shear_viscosity_viscous * kc2 * radius_shell * hd_kc2_R3 +
        vesms_i * omega * 2.0 * shear_viscosity_viscous * h_kc2_R3;
    M(5, 3) = vesms_i * omega * shear_viscosity_viscous * ks2 * ks2 * radius_shell * radius_shell * jdd_ks2_R3 -
        vesms_i * omega * 2.0 * shear_viscosity_viscous * j_ks2_R3 -
        vesms_i * omega * shear_viscosity_viscous * mm * j_ks2_R3;
    M(5, 4) = vesms_i * omega * shear_viscosity_viscous * ks2 * ks2 * radius_shell * radius_shell * hdd_ks2_R3 -
        vesms_i * omega * 2.0 * shear_viscosity_viscous * h_ks2_R3 -
        vesms_i * omega * shear_viscosity_viscous * mm * h_ks2_R3;
    M(5, 5) = -2.0 * shear_modulus_shell * kc3 * radius_shell * jd_kc3_R3 +
        2.0 * shear_modulus_shell * j_kc3_R3;
    M(5, 6) = -2.0 * shear_modulus_shell * kc3 * radius_shell * hd_kc3_R3 +
        2.0 * shear_modulus_shell * h_kc3_R3;
    M(5, 7) = shear_modulus_shell * ks3 * ks3 * radius_shell * radius_shell * jdd_ks3_R3 -
        2.0 * shear_modulus_shell * j_ks3_R3 - shear_modulus_shell * mm * j_ks3_R3;
    M(5, 8) = shear_modulus_shell * ks3 * ks3 * radius_shell * radius_shell * hdd_ks3_R3 -
        2.0 * shear_modulus_shell * h_ks3_R3 - shear_modulus_shell * mm * h_ks3_R3;

    M(6, 1) = j_kc2_R3;
    M(6, 2) = h_kc2_R3;
    M(6, 3) = -j_ks2_R3 - ks2 * radius_shell * jd_ks2_R3;
    M(6, 4) = -h_ks2_R3 - ks2 * radius_shell * hd_ks2_R3;
    M(6, 5) = -j_kc3_R3;
    M(6, 6) = -h_kc3_R3;
    M(6, 7) = j_ks3_R3 + ks3 * radius_shell * jd_ks3_R3;
    M(6, 8) = h_ks3_R3 + ks3 * radius_shell * hd_ks3_R3;

    M(7, 5) = kc3 * radius_gas * jd_kc3_R4;
    M(7, 6) = kc3 * radius_gas * hd_kc3_R4;
    M(7, 7) = mm * j_ks3_R4;
    M(7, 8) = mm * h_ks3_R4;
    M(7, 9) = -k4 * radius_gas * jd_k4_R4;

    M(8, 5) = shear_modulus_shell * (2.0 * kc3 * kc3 - ks3 * ks3) * radius_gas * radius_gas * j_kc3_R4 +
        2.0 * shear_modulus_shell * kc3 * kc3 * radius_gas * radius_gas * jdd_kc3_R4;
    M(8, 6) = shear_modulus_shell * (2.0 * kc3 * kc3 - ks3 * ks3) * radius_gas * radius_gas * h_kc3_R4 +
        2.0 * shear_modulus_shell * kc3 * kc3 * radius_gas * radius_gas * hdd_kc3_R4;
    M(8, 7) = 2.0 * shear_modulus_shell * mm * ks3 * radius_gas * jd_ks3_R4 -
        2.0 * shear_modulus_shell * mm * j_ks3_R4;
    M(8, 8) = 2.0 * shear_modulus_shell * mm * ks3 * radius_gas * hd_ks3_R4 -
        2.0 * shear_modulus_shell * mm * h_ks3_R4;
    M(8, 9) = omega * omega * density_gas * radius_gas * radius_gas * j_k4_R4;

    M(9, 5) = 2.0 * shear_modulus_shell * kc3 * radius_gas * jd_kc3_R4 -
        2.0 * shear_modulus_shell * j_kc3_R4;
    M(9, 6) = 2.0 * shear_modulus_shell * kc3 * radius_gas * hd_kc3_R4 -
        2.0 * shear_modulus_shell * h_kc3_R4;
    M(9, 7) = -shear_modulus_shell * ks3 * ks3 * radius_gas * radius_gas * jdd_ks3_R4 +
        2.0 * shear_modulus_shell * j_ks3_R4 + shear_modulus_shell * mm * j_ks3_R4;
    M(9, 8) = -shear_modulus_shell * ks3 * ks3 * radius_gas * radius_gas * hdd_ks3_R4 +
        2.0 * shear_modulus_shell * h_ks3_R4 + shear_modulus_shell * mm * h_ks3_R4;

    return vesms_solve_system_impl(M, F)(0);
}

// [[Rcpp::export]]
ComplexVector vesms_backscatter_cpp(
    NumericVector frequency,
    IntegerVector m_limit,
    double density_sw,
    double density_viscous,
    double density_shell,
    double density_gas,
    double sound_speed_sw,
    double sound_speed_viscous,
    double sound_speed_gas,
    double radius_viscous,
    double radius_shell,
    double radius_gas,
    double shear_viscosity_viscous,
    double bulk_viscosity_viscous,
    double shear_modulus_shell,
    double lambda_shell
) {
    int n_freq = frequency.size();
    if (m_limit.size() != n_freq) {
        Rcpp::stop("VESMS requires one modal limit per frequency.");
    }

    ComplexVector out(n_freq);
    for (int i = 0; i < n_freq; ++i) {
        double omega = 2.0 * vesms_pi * frequency[i];
        double k_sw = 2.0 * vesms_pi * frequency[i] / sound_speed_sw;
        std::complex<double> f_bs(0.0, 0.0);

        // Sum the retained modal coefficients into the far-field backscatter
        // amplitude for this frequency.
        for (int m = 0; m <= m_limit[i]; ++m) {
            std::complex<double> A1 = vesms_single_mode_impl(
                m,
                omega,
                density_sw,
                density_viscous,
                density_shell,
                density_gas,
                sound_speed_sw,
                sound_speed_viscous,
                sound_speed_gas,
                radius_viscous,
                radius_shell,
                radius_gas,
                shear_viscosity_viscous,
                bulk_viscosity_viscous,
                shear_modulus_shell,
                lambda_shell
            );

            double sign = (m % 2 == 0) ? 1.0 : -1.0;
            f_bs += sign * A1 * hs_single_complex_impl(m, std::complex<double>(k_sw, 0.0));
        }

        out[i] = to_Rcomplex(f_bs);
    }

    return out;
}
