#include <testthat.h>
#include "../../src/svd_solve.h"

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
    const std::vector<Complex> matrix = {{0, 0}, {2, 1}, {1, -1}, {3, 0}};
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
