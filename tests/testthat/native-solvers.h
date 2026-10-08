#include <testthat.h>
#include "../../src/svd_solve.h"

// Internal functions are linked directly so zero-order and empty-sequence
// behavior is checked below the R wrappers that normally short-circuit it.
std::vector<double> js_sequence_miller_impl(int, double);
std::vector<double> ys_sequence_upward_impl(int, double);
std::vector<double> js_deriv_sequence_impl(int, double);
std::vector<std::complex<double>> hs_deriv_sequence_impl(int, double);
double js_single_impl(int, double);
double ys_single_impl(int, double);
std::complex<double> hs_single_impl(int, double);
std::complex<double> js_single_complex_impl(int, std::complex<double>);
std::complex<double> ys_single_complex_impl(int, std::complex<double>);
std::complex<double> hs_single_complex_impl(int, std::complex<double>);
double js_deriv_single_impl(int, double, int);
double ys_deriv_single_impl(int, double, int);
std::complex<double> hs_deriv_single_impl(int, double, int);
std::complex<double> js_deriv_single_complex_impl(int, std::complex<double>, int);
std::complex<double> ys_deriv_single_complex_impl(int, std::complex<double>, int);
std::complex<double> hs_deriv_single_complex_impl(int, std::complex<double>, int);
std::complex<double> det6x6_scaled(std::complex<double>[6][6]);

context("Native Bessel and determinant edge cases") {
    test_that("zero-order derivatives retain values and empty orders stay empty") {
        using Complex = std::complex<double>;
        const Complex z(0.7, 0.3);
        for (int n : {0, 1, 4}) {
            expect_true(js_deriv_single_impl(n, 0.7, 0) == js_single_impl(n, 0.7));
            expect_true(ys_deriv_single_impl(n, 0.7, 0) == ys_single_impl(n, 0.7));
            expect_true(hs_deriv_single_impl(n, 0.7, 0) == hs_single_impl(n, 0.7));
            expect_true(js_deriv_single_complex_impl(n, z, 0) == js_single_complex_impl(n, z));
            expect_true(ys_deriv_single_complex_impl(n, z, 0) == ys_single_complex_impl(n, z));
            expect_true(hs_deriv_single_complex_impl(n, z, 0) == hs_single_complex_impl(n, z));
        }
        expect_true(js_sequence_miller_impl(-1, 0.7).empty());
        expect_true(ys_sequence_upward_impl(-1, 0.7).empty());
        expect_true(js_deriv_sequence_impl(-1, 0.7).empty());
        expect_true(hs_deriv_sequence_impl(-1, 0.7).empty());
        auto singular = ys_single_complex_impl(0, Complex(0, 0));
        expect_true(std::isnan(singular.real()));
        expect_true(std::isnan(singular.imag()));
        auto tiny = js_sequence_miller_impl(1, 1e-290);
        expect_true(std::abs(tiny[0] - 1) < 1e-12);
        expect_true(std::abs(tiny[1]) < 1e-280);
        auto converted = to_Rcomplex(Complex(2, -3));
        expect_true(converted.r == 2);
        expect_true(converted.i == -3);
    }

    test_that("scaled determinants preserve products, permutation signs and rank") {
        using Complex = std::complex<double>;
        Complex matrix[6][6] = {};
        Complex expected(1, 0);
        for (int i = 0; i < 6; ++i) {
            matrix[i][i] = Complex(i + 1, 1) * (i % 2 ? 1e30 : 1e-30);
            expected *= matrix[i][i];
        }
        expect_true(std::abs(det6x6_scaled(matrix) / expected - Complex(1, 0)) < 1e-12);
        for (int j = 0; j < 6; ++j) std::swap(matrix[0][j], matrix[1][j]);
        expect_true(std::abs(det6x6_scaled(matrix) / expected + Complex(1, 0)) < 1e-12);
        for (int j = 0; j < 6; ++j) matrix[0][j] = matrix[1][j];
        expect_true(det6x6_scaled(matrix) == Complex(0, 0));
    }
}

context("Minimum-norm SVD solutions") {
    test_that("rectangular complex systems recover known solutions") {
        arma::cx_mat matrix(3, 2);
        matrix.col(0) = arma::cx_vec({{1, 1}, {2, -1}, {0, 1}});
        matrix.col(1) = arma::cx_vec({{0, 2}, {1, 0}, {3, -1}});
        arma::cx_mat expected(2, 2);
        expected.col(0) = arma::cx_vec({{2, 1}, {-1, 2}});
        expected.col(1) = arma::cx_vec({{0, -1}, {3, 1}});
        arma::cx_mat rhs = matrix * expected;
        arma::cx_mat actual;
        expect_true(acousticts_svd_solve(actual, matrix, rhs));
        expect_true(arma::norm(actual - expected, "fro") < 1e-12);

        arma::cx_mat wide = matrix.t();
        arma::cx_vec minimum = matrix * expected.col(0);
        arma::cx_vec wide_rhs = wide * minimum;
        arma::cx_vec wide_result;
        expect_true(acousticts_svd_solve(wide_result, wide, wide_rhs));
        expect_true(arma::norm(wide_result - minimum) < 1e-12);
    }

    test_that("rank deficiency selects the minimum-norm solution") {
        arma::cx_mat matrix(2, 2, arma::fill::ones);
        arma::cx_vec rhs = {{2, 2}, {2, 2}};
        arma::cx_vec actual;
        expect_true(acousticts_svd_solve(actual, matrix, rhs));
        arma::cx_vec expected = {{1, 1}, {1, 1}};
        expect_true(arma::norm(actual - expected) < 1e-12);

        matrix.zeros(2, 3);
        expect_true(acousticts_svd_solve(actual, matrix, rhs));
        expect_true(actual.n_elem == 3);
        expect_true(actual.is_finite());
        expect_true(arma::norm(actual) == 0);
    }

    test_that("empty, mismatched and nonfinite systems are rejected") {
        arma::cx_mat result;
        arma::cx_mat empty;
        arma::cx_mat matrix(2, 2, arma::fill::eye);
        arma::cx_mat rhs(3, 1, arma::fill::ones);
        expect_false(acousticts_svd_solve(result, empty, rhs));
        expect_false(acousticts_svd_solve(result, matrix, rhs));
        rhs.set_size(2, 1);
        matrix(0, 0) = std::numeric_limits<double>::quiet_NaN();
        expect_false(acousticts_svd_solve(result, matrix, rhs));
    }
}

template<typename T>
void check_pivoted_solver() {
    using Complex = std::complex<T>;
    const std::vector<Complex> matrix = {{0, 0}, {4, 1}, {1, -1}, {2, 0}};
    const std::vector<Complex> expected = {{2, 1}, {-1, 2}};
    const auto rhs = matvec_product(matrix, expected, 2);
    auto factors = lup_decompose(matrix, 2);
    expect_false(factors.singular);
    auto result = solve_linear_system_lup_refined(matrix, rhs, 2, 3);
    for (int i = 0; i < 2; ++i) {
        expect_true(static_cast<double>(complex_abs_sq(result[i] - expected[i])) < 1e-24);
    }
    const std::vector<Complex> singular(4, Complex(1, 0));
    expect_true(lup_decompose(singular, 2).singular);
    expect_error(solve_linear_system_lup_refined(singular, rhs, 2));
    expect_error(solve_fluid_t_blocks_native<T>({matrix}, {singular}));
}

context("Pivoted spheroidal kernel solves") {
    test_that("residual refinement retains accuracy for large complex solutions") {
        using Complex = std::complex<double>;
        const std::vector<Complex> matrix = {
            {0.1, 0.3}, {0.2, -0.1}, {0.3, 0.7},
            {0.7, -0.2}, {0.11, 0.4}, {0.13, -0.9},
            {0.17, 0.6}, {0.19, -0.2}, {0.23, 0.5}
        };
        const std::vector<Complex> expected = {
            {1e14, 2e14}, {-3e14, 1e14}, {2e14, -1e14}
        };
        auto rhs = matvec_product(matrix, expected, 3);
        auto actual = solve_linear_system_lup_refined(matrix, rhs, 3, 3);
        auto recovered = matvec_product(matrix, actual, 3);
        for (int i = 0; i < 3; ++i) {
            expect_true(std::abs(actual[i] / expected[i] - Complex(1, 0)) < 1e-13);
            expect_true(std::abs(recovered[i] / rhs[i] - Complex(1, 0)) < 1e-13);
        }
    }

    test_that("complete pivoting preserves complex systems and rejects singularity") {
        // The double T-block specialization permits rank-deficient systems.
        // Exercise its LUP counterpart with the precision-independent kernel.
        check_pivoted_solver<long double>();
#if ACOUSTICTS_HAVE_QUADMATH
        check_pivoted_solver<acousticts_quad_t>();
#endif
    }

    test_that("double precision SVD handles dependent modal rows") {
        using Complex = std::complex<double>;
        const std::vector<std::vector<Complex>> matrix = {
            {{1, 0}, {1, 0}, {1, 0}, {1, 0}}
        };
        const std::vector<std::vector<Complex>> rhs = {{{2, 2}, {2, 2}}};
        auto actual = solve_fluid_Amn_divide_and_conquer(rhs, matrix);
        expect_true(std::abs(actual[0][0] - Complex(1, 1)) < 1e-12);
        expect_true(std::abs(actual[0][1] - Complex(1, 1)) < 1e-12);
    }
}

context("Packed spheroidal radial values") {
    test_that("mantissas, exponents and wave conventions are preserved") {
        ProfcnResult<double> packed;
        packed.r1c = {2};
        packed.r1dc = {-4};
        packed.ir1e = {2};
        packed.ir1de = {0};
        packed.r2c = {6};
        packed.r2dc = {8};
        packed.ir2e = {-1};
        packed.ir2de = {-2};
        ProfcnBatchResult<double> block;
        block.lnum = 1;
        block.r1c = {200};
        block.r1dc = {-4};
        block.r2c = {0.6};
        block.r2dc = {0.08};
        for (int kind = 1; kind <= 4; ++kind) {
            auto single = extract_radial_from_batch(
                packed, packed, 0, kind, 0, 0, 1.0, 1.5
            );
            auto batched = extract_radial_from_mblock(block, 0, 0, kind);
            expect_true(std::abs(single.val_real - (kind == 2 ? 0.6 : 200)) < 1e-12);
            expect_true(std::abs(single.der_real - (kind == 2 ? 0.08 : -4)) < 1e-12);
            expect_true(std::abs(single.val_imag - (kind < 3 ? 0 : kind == 3 ? 0.6 : -0.6)) < 1e-12);
            expect_true(std::abs(single.der_imag - (kind < 3 ? 0 : kind == 3 ? 0.08 : -0.08)) < 1e-12);
            expect_true(std::abs(batched.val_real - single.val_real) < 1e-12);
            expect_true(std::abs(batched.der_real - single.der_real) < 1e-12);
            expect_true(std::abs(batched.val_imag - single.val_imag) < 1e-12);
            expect_true(std::abs(batched.der_imag - single.der_imag) < 1e-12);
        }
    }

    test_that("missing and underflowed radial values remain invalid") {
        ProfcnResult<double> packed;
        ProfcnBatchResult<double> block;
        block.lnum = 1;
        for (int kind = 1; kind <= 4; ++kind) {
            auto single = extract_radial_from_batch(
                packed, packed, 0, kind, 0, 0, 1.0, 1.5
            );
            auto batched = extract_radial_from_mblock(block, 0, 0, kind);
            expect_true(std::isnan(single.val_real));
            expect_true(std::isnan(single.der_real));
            expect_true(std::isnan(batched.val_real));
            expect_true(std::isnan(batched.der_real));
        }
        packed.r2c = {0};
        packed.r2dc = {0};
        block.r2c = {0};
        block.r2dc = {0};
        for (int kind = 2; kind <= 4; ++kind) {
            auto single = extract_radial_from_batch(
                packed, packed, 0, kind, 0, 0, 1.0, 1.5
            );
            auto batched = extract_radial_from_mblock(block, 0, 0, kind);
            expect_true(std::isnan(kind == 2 ? single.val_real : single.val_imag));
            expect_true(std::isnan(kind == 2 ? batched.der_real : batched.der_imag));
        }
    }

    test_that("native radial layouts retain invalid lower-triangle modes") {
        for (int kind = 1; kind <= 4; ++kind) {
            expect_error(Rmn_matrix<double>({1}, {0, 1, 3}, 2, 1.5, kind));
            auto outer = Rmn_matrix<double>({0, 2}, {0, 1, 3}, 2, 1.5, kind);
            expect_true(std::isnan(outer.value[1][0].real()));
            expect_true(std::isnan(outer.derivative[1][1].real()));
            auto paired = Rmn_matrix<double>({0, 1}, {2, 3}, 2, 1.5, kind);
            auto expected = Rmn_scalar<double>(1, 3, 2, 1.5, kind);
            expect_true(std::abs(paired.value[1][1] - expected.first) < 1e-10);
        }
    }
}

context("Modal summation invariants") {
    test_that("adaptive and complete sums recover a finite geometric series") {
        const int order = 24;
        const double ratio = 0.1;
        using Complex = std::complex<double>;
        const Complex amplitude(2, -1);
        std::vector<std::vector<double>> angular(order + 1, std::vector<double>(order + 1, 1));
        std::vector<std::vector<Complex>> coefficients(order + 1, std::vector<Complex>(order + 1));
        std::vector<std::vector<Complex>> triangular(order + 1);
        for (int m = 0; m <= order; ++m) {
            for (int n = m; n <= order; ++n) {
                coefficients[m][n] = amplitude * std::pow(ratio, m + n);
                triangular[m].push_back(coefficients[m][n]);
            }
        }
        auto geometric = [&](double q) {
            return 1 + 2 * q * (1 - std::pow(q, order)) / (1 - q);
        };
        Complex expected = amplitude * (
            geometric(-ratio * ratio) - std::pow(-ratio, order + 1) * geometric(ratio)
        ) / (1 + ratio);
        for (bool adaptive : {false, true}) {
            expect_true(std::abs(compute_fbs_backscatter(
                order, order, angular, coefficients, adaptive
            ) - expected) < 1e-10);
            expect_true(std::abs(compute_fbs_backscatter_triangular(
                order, angular, triangular, adaptive
            ) - expected) < 1e-10);
        }
        auto azimuth = compute_azimuth(order, 0.0, std::acos(-1.0));
        auto reflected = reflect_smn_matrix(angular);
        expect_true(std::abs(compute_fbs(
            order, order, azimuth, angular, reflected, coefficients
        ) - expected) < 1e-12);
        expect_error(compute_fbs(
            order, order, std::vector<double>(), angular, reflected, coefficients
        ));
    }

    test_that("invalid angular modes do not contaminate the retained sum") {
        using Complex = std::complex<double>;
        const double nan = std::numeric_limits<double>::quiet_NaN();
        const std::vector<std::vector<double>> angular = {{1, nan}, {nan, 1}};
        const std::vector<std::vector<Complex>> coefficients = {{{1, 0}, {2, 0}}, {{nan, 0}, {2, 0}}};
        expect_true(compute_fbs_backscatter(1, 1, angular, coefficients) == Complex(-3, 0));
        expect_true(compute_fbs(1, 1, {1.0, -1.0}, angular, angular, coefficients) == Complex(-3, 0));
        auto expanded = expand_Amn_triangular<double>(1, 1, {{{1, 0}, {nan, 0}, {3, 0}}, {{2, 0}}});
        expect_true(std::isnan(expanded[0][1].real()));
        expect_true(expanded[1][1] == Complex(2, 0));
        expect_true(std::isnan(reflect_smn_matrix(angular)[0][1]));
    }
}

context("Packed angular data and native validation") {
    test_that("missing exponents and modes preserve missingness") {
        const double nan = std::numeric_limits<double>::quiet_NaN();
        ProfcnResult<double> packed;
        packed.s1c = {2, nan};
        packed.is1e = {2, NA_INTEGER};
        expect_true(extract_angular_value_from_batch(packed, 0) == 200);
        expect_true(std::isnan(extract_angular_value_from_batch(packed, 1)));
        expect_true(std::isnan(extract_angular_value_from_batch(packed, 3)));
        ProfcnBatchResult<double> block;
        block.lnum = 2;
        block.narg = 1;
        block.s1c = {200, nan};
        expect_true(extract_angular_value_from_mblock(block, 0, 0) == 200);
        expect_true(std::isnan(extract_angular_value_from_mblock(block, 0, 1)));
        expect_true(std::isnan(extract_angular_value_from_mblock(block, 1, 0)));
        std::vector<double> values = {2, 3, nan};
        std::vector<int> exponents = {2, NA_INTEGER, 4};
        scale_profcn_component(values, exponents);
        expect_true(values[0] == 200);
        expect_true(values[1] == 3);
        expect_true(std::isnan(values[2]));
        expect_true(exponents[1] == NA_INTEGER);
        expect_error(Smn_matrix<double>({}, {0}, 1, {0}));
        expect_error(Smn_matrix<double>({0}, {0}, 1, {}));
    }
}
