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

context("Portable precision sentinels") {
    test_that("broadside symmetry preserves adjacent public angle values") {
        const double half_pi = std::acos(-1.0) / 2;
        auto check = [half_pi](auto scalar) {
            using T = decltype(scalar);
            expect_true((psms_direction_cosine(T(half_pi)) == T(0)));
            for (double angle : {std::nextafter(half_pi, 0.0),
                    std::nextafter(half_pi, 2.0)}) {
                expect_true((psms_direction_cosine(T(angle)) == preccos(T(angle))));
                expect_true((psms_direction_cosine(T(angle)) != T(0)));
            }
        };
        check(double(0));
#if ACOUSTICTS_HAVE_QUADMATH
        check(acousticts_quad_t(0));
#endif
    }

    test_that("coupled gas cache respects physical inputs and angular parity") {
        auto check = [](auto scalar) {
            using T = decltype(scalar);
            auto evaluate = [](int m, int n, T c, T ci, T xi, T density, T eta,
                    T receive = T(0.4), T phi = T(1.3)) {
                return psms_full_coupled_field<T>(m, n, c, ci, xi,
                    eta, receive, phi, density, {}, {}, false, true);
            };
            auto invalidate = [&]() {
                evaluate(0, 2, T(0.2), T(0.7), T(1.1), T(0.01), T(0.4));
            };
            const T tolerance = is_quad_precision_v<T> ? T(1e-27) : T(1e-11);
            for (int changed = 0; changed < 10; ++changed) {
                const int m = changed == 0 ? 3 : 4, n = changed == 1 ? 19 : 18;
                const T c = changed == 2 ? T(4.1) : T(4);
                const T ci = changed == 3 ? T(17.1) : T(17);
                const T xi = changed == 4 ? T(1.02) : T(1.01);
                const T density = changed == 5 ? T(0.0013) : T(0.0012);
                const T eta = changed == 6 ? T(0) : T(0.6);
                const T receive = changed == 8 ? T(0.7) : T(0.4);
                const T phi = changed == 9 ? T(2.3) : T(1.3);
                evaluate(4, 18, T(4), T(17), T(1.01), T(0.0012), T(0.3));
                const auto reused = evaluate(m, n, c, ci, xi, density, eta, receive, phi);
                invalidate();
                const auto rebuilt = evaluate(m, n, c, ci, xi, density, eta, receive, phi);
                expect_true((complex_abs_sq(reused - rebuilt) <=
                    tolerance * tolerance * complex_abs_sq(rebuilt)));
            }
        };
        check(double(0));
#if ACOUSTICTS_HAVE_QUADMATH
        check(acousticts_quad_t(0));
#endif
    }

    test_that("zero angular modes preserve the next degree's convergence history") {
        auto check = [](auto scalar) {
            using T = decltype(scalar);
            const T tolerance = is_quad_precision_v<T> ? T(1e-25) : T(1e-10);
            for (int m : {0, 38}) {
                auto result = cprofcn_mblock_coefficients<T>(
                    T(58.93344), m, 1, 88, {T(0)}, 0, 1, T(0.01)
                );
                // Evaluate the normalized Ferrers basis at zero independently.
                T seed = T(1) / precsqrt(T(2));
                for (int order = 1; order <= m; ++order) {
                    seed *= precsqrt(T(2 * order + 1) / T(2 * order));
                }
                for (int index : {0, 1, 2, 3, 10, 11, 40, 41, 86, 87}) {
                    T basis = seed, expected = T(0);
                    for (int offset = 0; offset < result.lengths[index]; offset += 2) {
                        if (offset != 0) {
                            const T degree = T(m + offset);
                            basis *= -precsqrt(
                                ((degree + T(m) - T(1)) * (degree - T(m) - T(1)) *
                                    (T(2) * degree + T(1))) /
                                ((degree - T(m)) * (degree + T(m)) *
                                    (T(2) * degree - T(3)))
                            );
                        }
                        expected += result.coefficients[result.capacity * index + offset] * basis;
                    }
                    const T actual = extract_angular_value_from_mblock<T>(
                        result.functions, 0, index, 0
                    );
                    expect_true((precabs(actual - expected) < tolerance));
                    if (index % 2 != 0) expect_true((actual == T(0)));
                }
            }
        };
        check(double(0));
#if ACOUSTICTS_HAVE_QUADMATH
        check(acousticts_quad_t(0));
#endif
    }

    test_that("one-degree radial batches preserve the joining-factor recurrence") {
        auto check = [](auto scalar) {
            using T = decltype(scalar);
            const T c = T(11.2), x1 = T(0.01);
            const T tolerance = is_quad_precision_v<T> ? T(1e-25) : T(1e-10);
            auto single = cprofcn_mblock<T>(c, 48, 2, 1, {T(0.2)}, 2, 1, 1, x1);
            auto padded = cprofcn_mblock<T>(c, 48, 2, 2, {T(0.2)}, 2, 1, 1, x1);
            auto coefficients = cprofcn_mblock_coefficients<T>(
                c, 48, 2, 1, {T(0.2)}, 2, 1, x1
            );
            expect_true((single.lnum == 1 && single.r2c.size() == 2));
            expect_true((coefficients.functions.r2c == single.r2c));
            expect_true((coefficients.functions.r2dc == single.r2dc));
            const T wronskian = T(1) / (c * x1 * (x1 + T(2)));
            for (int m = 0; m < 2; ++m) {
                // The first degree must agree with a two-degree native request.
                expect_true((precabs(single.r1c[m] / padded.r1c[2 * m] - T(1)) < tolerance));
                expect_true((precabs(single.r1dc[m] / padded.r1dc[2 * m] - T(1)) < tolerance));
                expect_true((precabs(single.r2c[m] / padded.r2c[2 * m] - T(1)) < tolerance));
                expect_true((precabs(single.r2dc[m] / padded.r2dc[2 * m] - T(1)) < tolerance));
                const T first_product = single.r1c[m] * single.r2dc[m] *
                    pow10_typed<T>(single.ir1e[m] + single.ir2de[m]);
                const T second_product = single.r1dc[m] * single.r2c[m] *
                    pow10_typed<T>(single.ir1de[m] + single.ir2e[m]);
                // Check the Wronskian against the backend's reported decimal accuracy.
                expect_true((single.naccr[m] > 0));
                const T reported_radial_error = pow10_typed<T>(-single.naccr[m]);
                expect_true((precabs((first_product - second_product) / wronskian - T(1)) < reported_radial_error));
            }
        };
        check(double(0));
#if ACOUSTICTS_HAVE_QUADMATH
        check(acousticts_quad_t(0));
#endif
    }

    test_that("spheroidal expansion coefficients preserve the orthonormal basis") {
        auto check = [](auto scalar) {
            using T = decltype(scalar);
            for (T c : {T(0.2), T(8), T(80)}) {
                auto coefficients = cprofcn_mblock_coefficients<T>(
                    c, 0, 3, 12, {T(-0.6), T(0), T(0.6)}, 0, 1, T(0)
                );
                auto legacy = cprofcn_mblock<T>(
                    c, 0, 3, 12, {T(-0.6), T(0), T(0.6)}, 0, 1, 1, T(0)
                );
                expect_true((coefficients.functions.s1c == legacy.s1c));
                expect_true((coefficients.functions.naccs == legacy.naccs));
                // Unit-normalized angular eigenfunctions have an identity Gram matrix.
                T largest_error = 0;
                for (int m = 0; m < 3; ++m) {
                    for (int row = 0; row < 12; ++row) {
                        for (int col = 0; col < 12; ++col) {
                            const T value = spheroidal_coefficient_overlap<T>(
                                coefficients, coefficients, m, row, col
                            );
                            largest_error = std::max(largest_error,
                                precabs(value - (row == col ? T(1) : T(0))));
                        }
                    }
                }
                const T tolerance = is_quad_precision_v<T> ? T(1e-27) : T(1e-12);
                expect_true((largest_error < tolerance));
            }
        };
        check(double(0));
#if ACOUSTICTS_HAVE_QUADMATH
        check(acousticts_quad_t(0));
#endif
    }

    test_that("radial quadrature reuse is independent of scattering angle") {
        auto check = [](auto scalar) {
            using T = decltype(scalar);
            auto first = cprofcn_mblock<T>(T(40), 0, 3, 10, {T(0.2)}, 2, 1, 1, T(0.01));
            auto second = cprofcn_mblock<T>(T(40), 0, 3, 10, {T(0.4)}, 2, 1, 1, T(0.01));
            expect_true((first.r1c == second.r1c));
            expect_true((first.r2c == second.r2c));
            expect_true((first.r1dc == second.r1dc));
            expect_true((first.r2dc == second.r2dc));
            expect_true((first.ir2e == second.ir2e));
            expect_true((first.s1c != second.s1c));
            // Check the radial Wronskian before any conversion back to R doubles.
            const T wronskian = T(1) / (T(40) * T(0.01) * (T(0.01) + T(2)));
            const T tolerance = is_quad_precision_v<T> ? T(1e-25) : T(1e-10);
            for (size_t i = 0; i < first.r1c.size(); ++i) {
                const T first_product = first.r1c[i] * first.r2dc[i] *
                    pow10_typed<T>(first.ir1e[i] + first.ir2de[i]);
                const T second_product = first.r1dc[i] * first.r2c[i] *
                    pow10_typed<T>(first.ir1de[i] + first.ir2e[i]);
                expect_true((precabs((first_product - second_product) / wronskian - T(1)) < tolerance));
            }
        };
        check(double(0));
#if ACOUSTICTS_HAVE_QUADMATH
        check(acousticts_quad_t(0));
#endif
    }

    test_that("values-only spheroidal batches preserve angular normalization") {
        auto check = [](auto scalar) {
            using T = decltype(scalar);
            const std::vector<T> eta = {T(-0.7), T(0), T(0.3), T(1)};
            auto reference = cprofcn_mblock<T>(T(8), 0, 3, 8, eta, 2, 1, 2, T(0.15));
            for (int repeat = 0; repeat < 2; ++repeat) {
                auto actual = cprofcn_mblock<T>(T(8), 0, 3, 8, eta, 2, 1, 1, T(0.15));
                for (size_t i = 0; i < actual.s1c.size(); ++i) {
                    expect_true((precabs(actual.s1c[i] - reference.s1c[i]) < T(1e-12)));
                    expect_true((actual.s1dc[i] == T(0)));
                }
            }
        };
        check(double(0));
#if ACOUSTICTS_HAVE_QUADMATH
        check(acousticts_quad_t(0));
#endif
    }

    test_that("missing values remain NaN across native precision conversions") {
        expect_true(is_na_real(precnan<double>()));
        expect_true(is_na_real(precnan<acousticts_quad_t>()));
        expect_true(std::isnan(quad_to_double(precnan<acousticts_quad_t>())));
        expect_true(is_na_real(double_to_quad(precnan<double>())));
    }
}

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
    test_that("column scaling preserves small coefficients and deficient norms") {
        arma::cx_mat matrix(3, 3, arma::fill::eye);
        matrix(0, 1) = {0.01, 0.02};
        matrix(2, 0) = {-0.02, 0.01};
        // The unscaled SVD drops a physical mode at this radial scale ratio.
        matrix.col(1) *= 1e16;
        matrix.col(2) *= 1e8;
        arma::cx_mat expected(3, 2);
        expected.col(0) = arma::cx_vec({{2, 1}, {1e-16, -2e-16}, {3e-8, 1e-8}});
        expected.col(1) = arma::cx_vec({{1, -2}, {3e-16, 1e-16}, {-2e-8, 1e-8}});
        arma::cx_mat rhs = matrix * expected;
        arma::cx_mat actual = solve_fluid_system_svd(matrix, rhs);
        for (arma::uword i = 0; i < actual.n_elem; ++i) {
            expect_true(std::abs(actual[i] / expected[i] - 1.0) < 1e-12);
        }

        matrix.ones(2, 2);
        matrix.col(1) *= 2.0;
        rhs.ones(2, 1);
        actual = solve_fluid_system_svd(matrix, rhs);
        expect_true(std::abs(actual[0] - 0.2) < 1e-12);
        expect_true(std::abs(actual[1] - 0.4) < 1e-12);
        matrix.zeros();
        actual = solve_fluid_system_svd(matrix, rhs);
        expect_true(arma::norm(actual) == 0.0);
    }

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
