#include <testthat.h>
#include <cstdio>
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
void gauss_legendre(int, std::vector<double>&, std::vector<double>&);
extern "C" void profcn_auxiliary_cache_audit(int*);
#if ACOUSTICTS_HAVE_QUADMATH
extern "C" void profcn_auxiliary_cache_audit_quad(int*);
#endif


context("Checked gas recovery linear algebra") {
    test_that("factorized Robin solves match an independent manufactured dense operator") {
        auto check = [](auto scalar) {
            using T = decltype(scalar);
            using Complex = std::complex<T>;
            const int rows = 5, columns = 3;
            const std::vector<T> overlap{T(1), T(1)/4, T(-1)/8, T(1)/2, T(1), T(1)/4,
                T(-1)/4, T(1)/8, T(1), T(1)/8, T(-1)/2, T(1)/4, T(1)/4, T(1)/4, T(-1)/2};
            const std::vector<Complex> response{{T(1)/8, T(1)/16}, {T(-1)/16, T(1)/32},
                {T(1)/16, T(-1)/8}, {T(1)/32, T(1)/16}, {T(-1)/8, T(-1)/16}};
            const std::vector<Complex> diagonal{{T(2), T(1)/4}, {T(3), T(-1)/2}, {T(4), T(1)/8}};
            const std::vector<Complex> expected{{T(1)/2, T(1)/4}, {T(-1), T(1)/8}, {T(2), T(-1)/2}};
            std::vector<Complex> matrix(columns * columns), rhs(columns);
            T matrix_norm = 0;
            // Dyadic inputs make this independently assembled operator exactly representable.
            for (int row = 0; row < columns; ++row) {
                T row_norm = 0;
                for (int col = 0; col < columns; ++col) {
                    Complex entry = row == col ? diagonal[row] : Complex(0, 0);
                    for (int inner = 0; inner < rows; ++inner)
                        entry -= overlap[size_t(inner) * columns + row] * response[inner] *
                            overlap[size_t(inner) * columns + col];
                    matrix[size_t(row) * columns + col] = entry;
                    rhs[row] += entry * expected[col];
                    row_norm += precabs(entry.real()) + precabs(entry.imag());
                }
                matrix_norm = std::max(matrix_norm, row_norm);
            }
            const auto dense = solve_linear_system_lup_refined<T>(matrix, rhs, columns, 2);
            const auto actual = psms_robin_solver::solve(overlap, response, diagonal, rhs, rows, columns);
            expect_true(bool(actual.success && actual.status == "converged_true_residual"));
            expect_true(bool(actual.operator_norm_is_bound && actual.operator_norm >= matrix_norm));
            const T allowance = T(2048) * preceps<T>();
            const auto mapped = psms_robin_solver::apply_map(overlap, response, expected, rows, columns);
            if (actual.solution.size() == expected.size()) {
                for (int row = 0; row < columns; ++row) {
                    const auto known_error = actual.solution[row] - expected[row];
                    const auto dense_error = actual.solution[row] - dense[row];
                    const auto map_error = diagonal[row] * expected[row] - mapped[row] - rhs[row];
                    expect_true(bool(precabs(known_error.real()) + precabs(known_error.imag()) <= allowance));
                    expect_true(bool(precabs(dense_error.real()) + precabs(dense_error.imag()) <= allowance));
                    expect_true(bool(precabs(map_error.real()) + precabs(map_error.imag()) <= allowance));
                }
            } else expect_true(false);
            expect_true(bool(actual.rhs_relative_residual <= T(64) * preceps<T>()));
            // Iteration exhaustion and an unusable preconditioner remain explicit failures.
            psms_robin_solver::Options<T> exhausted;
            exhausted.maximum_iterations = 0;
            const auto limited = psms_robin_solver::solve(overlap, response, diagonal, rhs, rows, columns, exhausted);
            expect_true(bool(!limited.success && limited.status == "iteration_budget_exhausted"));
            const auto unavailable = psms_robin_solver::solve(std::vector<T>{T(1)},
                std::vector<Complex>{{T(1), T(0)}}, std::vector<Complex>{{T(1), T(0)}},
                std::vector<Complex>{{T(1), T(0)}}, 1, 1);
            expect_true(bool(!unavailable.success && unavailable.status == "diagonal_preconditioner_unavailable"));
            const auto zero = psms_robin_solver::solve(overlap, response, diagonal,
                std::vector<Complex>(columns), rows, columns);
            expect_true(bool(zero.success && zero.iterations == 0 && zero.residual_norm == T(0)));
            const auto mismatch = psms_robin_solver::solve(overlap, response, diagonal, rhs, rows + 1, columns);
            expect_true(bool(!mismatch.success && mismatch.status == "dimension_mismatch"));
        };
        check(double(0));
#if ACOUSTICTS_HAVE_QUADMATH
        check(acousticts_quad_t(0));
#endif
    }

    test_that("componentwise refinement preserves accepted values and corrects a tiny forcing") {
        auto check = [](auto scalar) {
            using T = decltype(scalar);
            using Complex = std::complex<T>;
            const std::vector<Complex> matrix{{T(0), T(0)}, {T(2), T(0)}, {T(1), T(-1)},
                {T(3), T(1)}, {T(1), T(0)}, {T(0), T(0)},
                {T(1), T(0)}, {T(0), T(0)}, {T(4), T(-1)}};
            const T scale = preceps<T>() * preceps<T>();
            const std::vector<Complex> expected{Complex(T(1)/2, T(1)/4) * scale,
                Complex(T(-1), T(1)/8) * scale, Complex(T(2), T(-1)/2) * scale};
            std::vector<Complex> rhs(3);
            for (int row = 0; row < 3; ++row)
                for (int col = 0; col < 3; ++col) rhs[row] += matrix[size_t(row) * 3 + col] * expected[col];
            const auto factors = lup_decompose<T>(matrix, 3);
            expect_true(bool(!factors.singular));
            auto accepted = expected;
            const auto initial = psms_direct_refinement::refine(matrix, rhs, accepted, factors, 3);
            expect_true(bool(initial.accepted && initial.correction_count == 0 && accepted == expected));
            // The absolute defect is tiny, but its componentwise relative error requires correction.
            auto perturbed = expected;
            perturbed[0] += Complex(T(1)/8, T(-1)/16) * scale;
            perturbed[2] -= Complex(T(1)/16, T(1)/32) * scale;
            const auto corrected = psms_direct_refinement::refine(matrix, rhs, perturbed, factors, 3);
            expect_true(bool(corrected.accepted && corrected.correction_count > 0 && corrected.correction_count <= 5));
            expect_true(bool(corrected.initial_backward_error > T(64) * preceps<T>() &&
                corrected.backward_error <= T(64) * preceps<T>()));
            for (int index = 0; index < 3; ++index) {
                const auto error = (perturbed[index] - expected[index]) / scale;
                expect_true(bool(precabs(error.real()) + precabs(error.imag()) <= T(2048) * preceps<T>()));
            }
        };
        check(double(0));
#if ACOUSTICTS_HAVE_QUADMATH
        check(acousticts_quad_t(0));
#endif
    }

    test_that("native forward error detects sensitivity after backward convergence") {
        auto check = [](auto scalar) {
            using T = decltype(scalar);
            using Complex = std::complex<T>;
            const T epsilon = preceps<T>(), allowance = T(4096) * epsilon;
            const std::vector<Complex> identity{{T(1), T(0)}, {T(0), T(0)},
                {T(0), T(0)}, {T(1), T(0)}};
            const std::vector<Complex> expected{{T(1), T(1)}, {T(-1)/2, T(1)/4}};
            const auto identity_factors = lup_decompose<T>(identity, 2);
            auto solution = expected;
            const auto ordinary = psms_direct_refinement::refine(
                identity, expected, solution, identity_factors, 2);
            expect_true(bool(ordinary.accepted && !ordinary.forward_checked && solution == expected));
            const auto checked = psms_direct_refinement::refine(
                identity, expected, solution, identity_factors, 2, T(64) * epsilon);
            // For identity, the largest ZGERFS row weight divided by ||x|| is six eps.
            expect_true(bool(checked.accepted && checked.forward_checked && solution == expected));
            expect_true(bool(precabs(checked.forward_error / epsilon - T(6)) <= allowance));

            const T separation = T(1024) * epsilon;
            const std::vector<Complex> sensitive{{T(1), T(0)}, {T(1), T(0)},
                {T(1), T(0)}, {T(1) + separation, T(0)}};
            const std::vector<Complex> forcing{{T(0), T(0)}, {-separation, T(0)}};
            const auto sensitive_factors = lup_decompose<T>(sensitive, 2);
            expect_true(bool(!sensitive_factors.singular));
            solution = {{T(1), T(0)}, {T(-1), T(0)}};
            const auto unguarded = psms_direct_refinement::refine(
                sensitive, forcing, solution, sensitive_factors, 2);
            expect_true(bool(unguarded.accepted && unguarded.backward_error == T(0)));
            const auto guarded = psms_direct_refinement::refine(
                sensitive, forcing, solution, sensitive_factors, 2, T(64) * epsilon);
            expect_true(bool(!guarded.accepted && guarded.status == "forward_error_budget" &&
                guarded.forward_error > T(64) * epsilon && guarded.backward_error == T(0)));

            // A nonsymmetric complex example exercises both complete-pivot permutations.
            const std::vector<Complex> matrix{{T(2), T(1)}, {T(-1), T(2)}, {T(1)/2, T(-3)/4},
                {T(-2), T(1)/4}, {T(3), T(-1)}, {T(1), T(3)/2},
                {T(1)/4, T(1)/2}, {T(8), T(-3)}, {T(4), T(2)}};
            const std::vector<Complex> adjoint_expected{{T(1)/2, T(-3)/4},
                {T(-5)/4, T(1)/4}, {T(2), T(3)/2}};
            const auto factors = lup_decompose<T>(matrix, 3);
            expect_true(bool(!factors.singular && factors.row_piv[0] != 0 && factors.col_piv[0] != 0));
            std::vector<Complex> adjoint_rhs(3);
            for (int row = 0; row < 3; ++row)
                for (int column = 0; column < 3; ++column)
                    adjoint_rhs[row] += std::conj(matrix[size_t(column) * 3 + row]) * adjoint_expected[column];
            const auto adjoint = psms_direct_refinement::adjoint_solve(factors, adjoint_rhs);
            for (int index = 0; index < 3; ++index)
                expect_true(bool(psms_direct_refinement::abs1(adjoint[index] - adjoint_expected[index]) <= allowance));

            bool nonfinite_rejected = false, dimensions_rejected = false;
            try {
                psms_direct_refinement::forward_error_estimate(identity_factors,
                    std::vector<T>{precnan<T>(), T(1)}, T(1));
            } catch (const psms_direct_refinement::ForwardArithmeticError&) { nonfinite_rejected = true; }
            try {
                psms_direct_refinement::forward_error_estimate(identity_factors,
                    std::vector<T>{T(1)}, T(1));
            } catch (const std::invalid_argument&) { dimensions_rejected = true; }
            expect_true(bool(nonfinite_rejected && dimensions_rejected));
        };
        check(double(0));
#if ACOUSTICTS_HAVE_QUADMATH
        check(acousticts_quad_t(0));
#endif
    }

    test_that("matrix-free gas coupling retains general receiver geometry and boundary traces") {
        auto check = [](auto scalar) {
            using T = decltype(scalar);
            using Complex = std::complex<T>;
            const int m = 1, degree = 32, interior_degree = 40;
            const T chi = T(4), interior_chi = T(17), xi = T(3)/T(2), density = T(1)/T(800);
            const auto incidence = psms_direction_from_angle(T(0.61));
            const auto reception = psms_direction_from_angle(T(1.04));
            psms_automatic::LegendreWorkspace<T> workspace;
            auto exterior = cprofcn_mblock_coefficients<T>(chi, m, 1, degree - m + 1,
                {incidence.cosine, reception.cosine}, 2, 0, xi - T(1), true);
            psms_automatic::evaluate_coefficient_angular_values(exterior, m,
                std::vector<SpheroidalDirection<T>>{incidence, reception}, workspace);
            const auto dense = psms_automatic::solve_order<T>(m, degree, interior_degree,
                interior_chi, xi, incidence.cosine, T(0.7), density, false, true,
                exterior, nullptr, nullptr, true, true);
            const auto actual = psms_automatic::solve_order<T>(m, degree, interior_degree,
                interior_chi, xi, incidence.cosine, T(0.7), density, false, true,
                exterior, nullptr, nullptr, true, true, false, true);
            expect_true(bool(!dense.singular && !actual.singular));
            if (dense.singular || actual.singular) return;
            auto compare = [](const std::vector<Complex>& observed, const std::vector<Complex>& expected) {
                expect_true(bool(observed.size() == expected.size()));
                if (observed.size() != expected.size()) return;
                T scale = T(1), error = 0;
                bool finite = true;
                for (size_t index = 0; index < observed.size(); ++index) {
                    finite = finite && is_finite_real(observed[index].real()) &&
                        is_finite_real(observed[index].imag()) && is_finite_real(expected[index].real()) &&
                        is_finite_real(expected[index].imag());
                    scale = std::max(scale, precabs(expected[index].real()) + precabs(expected[index].imag()));
                    const auto difference = observed[index] - expected[index];
                    error = std::max(error, precabs(difference.real()) + precabs(difference.imag()));
                }
                expect_true(bool(finite && error <= T(8192) * preceps<T>() * scale));
            };
            compare({actual.field}, {dense.field});
            compare(actual.coefficients, dense.coefficients);
            compare(actual.outgoing_terms, dense.outgoing_terms);
            expect_true(bool(actual.radial_exponents == dense.radial_exponents));
            expect_true(bool(actual.outgoing_terms.size() == size_t(degree - m + 1)));
            for (size_t trace = 0; trace < actual.traces.size(); ++trace)
                compare(actual.traces[trace], dense.traces[trace]);
        };
        check(double(0));
#if ACOUSTICTS_HAVE_QUADMATH
        check(acousticts_quad_t(0));
#endif
    }
}

context("Automatic spheroidal boundary convergence") {
    test_that("native spheroidal sizing rejects overflowing expressions before allocation") {
        auto check = [](auto scalar) {
            using T = decltype(scalar);
            const bool quad = is_quad_precision_v<T>;
            const auto angular = preflight_profcn_dimensions(T(4), 2, 3, 8, 5, 0, 2, T(0));
            expect_true(bool(angular.maxm == 4 && angular.maxint == (quad ? 116 : 62)));
            expect_true(bool(angular.maxp == (quad ? 137 : 65) && angular.maxq == (quad ? 124 : 70)));
            expect_true(bool(angular.radial_size == 24 && angular.angular_size == 120));
            const auto radial = preflight_profcn_dimensions(T(20), 2, 1, 12, 1, 2, 0, T(2));
            expect_true(bool(radial.maxn == (quad ? 1457 : 1277)));
            expect_true(bool(radial.maxp == radial.maxn && radial.ngau == 200));
            // PLEG's integer 2*m*m+1 changes from representable to overflowing at this adjacent pair.
            expect_true(bool(preflight_profcn_dimensions(T(4), 32767, 1, 2, 1, 2, 0, T(0.1)).maxm == 32767));
            expect_error(preflight_profcn_dimensions(T(4), 32768, 1, 2, 1, 2, 0, T(0.1)));
            expect_error(preflight_profcn_dimensions(T(4), INT_MAX, 1, 1, 1, 0, 0, T(0)));
            expect_error(preflight_profcn_dimensions(T(4), 0, 1, INT_MAX, 1, 0, 0, T(0)));
            expect_error(preflight_profcn_dimensions(T(4), 0, 1, 2, size_t(INT_MAX) + 1, 0, 0, T(0)));
            expect_error(preflight_profcn_dimensions(T(1100000000), 0, 1, 2, 1, 2, 0, T(0.01)));
            expect_error(preflight_profcn_dimensions(T(100000000), 0, 1, 2, 1, 2, 0, T(0.1)));
            expect_error(preflight_profcn_dimensions(T(1000000000), 0, 1, 2, 1, 1, 0, T(2)));
            // All production wrapper families must reject before allocating their output rectangles.
            expect_error(cprofcn<T>(T(4), INT_MAX, 1, {T(0)}, 0, 1, 2, T(0)));
            expect_error(cprofcn_mblock<T>(T(4), INT_MAX, 1, 1, {T(0)}, 0, 1, 2, T(0)));
            expect_error(cprofcn_mblock_coefficients<T>(T(4), INT_MAX, 1, 1, {T(0)}, 0, 0, T(0)));
        };
        check(double(0));
#if ACOUSTICTS_HAVE_QUADMATH
        check(acousticts_quad_t(0));
#endif
        expect_error(checked_profcn_product(std::numeric_limits<size_t>::max(), 2,
            std::numeric_limits<size_t>::max(), "test product"));
    }

    test_that("inclusive spheroidal degree counts reject before signed overflow or allocation") {
        const int largest_degree = (INT_MAX - 1) / 2;
        expect_true(bool(checked_profcn_integer(2LL * largest_degree + 1,
            "native maximum degree") == INT_MAX));
        expect_error(checked_profcn_integer(2LL * (static_cast<long long>(largest_degree) + 1) + 1,
            "native maximum degree"));
        expect_true(bool(psms_automatic::grow_degree(700, 37) == 866));
        expect_error(psms_automatic::grow_degree(largest_degree, 0));
        expect_error(psms_automatic::grow_degree(INT_MAX - 2, INT_MAX - 2));
        auto check = [](auto scalar) {
            using T = decltype(scalar);
            // Small request vectors exercise every radial layout without constructing a large result.
            expect_error(Rmn_scalar<T>(0, INT_MAX, T(4), T(1.5), 1));
            expect_error(Rmn_matrix<T>({0}, {0, INT_MAX}, T(4), T(1.5), 1));
            expect_error(Rmn_matrix<T>({0, 0}, {0, INT_MAX}, T(4), T(1.5), 1));
            expect_error(Rmn_matrix<T>({0, 1}, {0, 1, INT_MAX}, T(4), T(1.5), 1));
            expect_error(compute_smn_direction_matrices<T>(0, INT_MAX, T(4), {}));
            for (int degree : {INT_MAX - 1, INT_MAX}) {
                const auto fluid = psms_automatic_full_fluid_field<T>(0, degree,
                    T(4), T(5), T(1.5), T(0.3), T(0.4), T(1), T(1.1), false);
                expect_true(bool(!fluid.converged && fluid.orders.empty()));
                expect_true(bool(fluid.detail.find("integer range exceeded") != std::string::npos));
                for (bool rigid : {false, true}) {
                    const auto boundary_result = psms_automatic_rigid_release_field<T>(0, degree,
                        T(4), T(1.5), T(0.3), T(0.4), T(1), false, rigid);
                    expect_true(bool(!boundary_result.converged && boundary_result.orders.empty()));
                    expect_true(bool(boundary_result.detail.find("integer range exceeded") != std::string::npos));
                }
            }
        };
        check(double(0));
#if ACOUSTICTS_HAVE_QUADMATH
        check(acousticts_quad_t(0));
#endif
    }

    test_that("bounded native auxiliary caches preserve cold warm and uncached results") {
        const std::array<const char*, 15> checks = {"checked dimension arithmetic",
            "bounded nonempty retention", "oversize admission bypass", "budget eviction",
            "cold versus warm native values", "cold versus uncached native values",
            "zero-budget bypass", "adjacent signed native arguments remain distinct",
            "changed PLEG arguments miss the cache", "changed PLEG cached versus disabled values",
            "changed QLEG arguments miss the cache", "changed QLEG cached versus disabled values",
            "PINT even-to-odd prefix versus cold values", "PINT odd-to-even prefix versus cold values",
            "PINT after clearing cache versus cold values"};
        auto verify = [&](void (*audit)(int*), const char* precision) {
            std::array<int, 15> status{};
            audit(status.data());
            for (size_t index = 0; index < status.size(); ++index) {
                CATCH_INFO(precision << ": " << checks[index]);
                expect_true(bool(status[index] == 1));
            }
        };
        verify(profcn_auxiliary_cache_audit, "double");
#if ACOUSTICTS_HAVE_QUADMATH
        verify(profcn_auxiliary_cache_audit_quad, "quad");
#endif
    }

    test_that("public backscatter detection preserves adjacent receiver directions") {
        const double pi = std::acos(-1.0);
        for (double incident : {0.0, 0.61, pi / 2, pi}) {
            const double receiver = pi - incident;
            for (double azimuth : {0.0, 0.23, pi, -pi}) {
                for (double opposite : {azimuth + pi, azimuth - pi}) {
                    expect_true(bool(psms_public_backscatter_angles(incident, receiver, azimuth, opposite)));
                    for (double toward : {-std::numeric_limits<double>::infinity(),
                            std::numeric_limits<double>::infinity()}) {
                        expect_true(bool(!psms_public_backscatter_angles(incident,
                            std::nextafter(receiver, toward), azimuth, opposite)));
                        expect_true(bool(!psms_public_backscatter_angles(incident,
                            receiver, azimuth, std::nextafter(opposite, toward))));
                    }
                }
            }
        }
    }

    test_that("public near-monostatic fields use the specified general receiver direction") {
        auto check = [](auto scalar) {
            using T = decltype(scalar);
            const double pi = std::acos(-1.0), incident = 0.61, receiver = pi - incident + 1e-14;
            const std::string precision = is_quad_precision_v<T> ? "quad" : "double";
            const auto reference = psms_fixed_support_field<T>(2, 12, T(3), T(2.5), T(1.5),
                psms_direction_from_angle(T(incident)), psms_direction_from_angle(T(receiver)),
                T(pi), T(1.1), false, true, false, false);
            const auto acoustics = Rcpp::DataFrame::create(Rcpp::_["chi_sw"] = Rcpp::NumericVector::create(3),
                Rcpp::_["chi_body"] = Rcpp::NumericVector::create(2.5),
                Rcpp::_["m_max"] = Rcpp::IntegerVector::create(2),
                Rcpp::_["n_max"] = Rcpp::IntegerVector::create(12));
            const auto body = Rcpp::DataFrame::create(Rcpp::_["xi"] = Rcpp::NumericVector::create(1.5),
                Rcpp::_["theta_body"] = Rcpp::NumericVector::create(incident),
                Rcpp::_["theta_scatter"] = Rcpp::NumericVector::create(receiver),
                Rcpp::_["phi_body"] = Rcpp::NumericVector::create(0),
                Rcpp::_["phi_scatter"] = Rcpp::NumericVector::create(pi),
                Rcpp::_["density"] = Rcpp::NumericVector::create(1100));
            const auto medium = Rcpp::DataFrame::create(Rcpp::_["density"] = Rcpp::NumericVector::create(1000));
            const auto integration = Rcpp::List::create(Rcpp::_["nodes"] = Rcpp::NumericVector(),
                Rcpp::_["weights"] = Rcpp::NumericVector());
            const auto actual = prolate_spheroid_fbs(acoustics, body, medium, integration,
                precision, "Amn_pressure_release", false, false);
            const Rcomplex expected = to_Rcomplex(reference), observed = actual[0];
            expect_true(bool(observed.r == expected.r && observed.i == expected.i));
        };
        check(double(0));
#if ACOUSTICTS_HAVE_QUADMATH
        check(acousticts_quad_t(0));
#endif
    }

    test_that("native directions preserve near-axis and adjacent broadside Ferrers values") {
        auto check = [](auto scalar) {
            using T = decltype(scalar);
            const double public_pi = std::acos(-1.0), public_broadside = public_pi / 2;
            const T small = pow10_typed<T>(is_quad_precision_v<T> ? -20 : -8);
            const auto near_axis = psms_direction_from_angle(small);
            expect_true(bool(near_axis.cosine == T(1)));
            expect_true(bool(near_axis.sine > T(0) && !near_axis.exact_axis));
            expect_true(bool(psms_direction_from_angle(T(0)).exact_axis));
            expect_true(bool(psms_direction_from_angle(T(public_pi)).exact_axis));
            psms_automatic::LegendreWorkspace<T> workspace;
            const std::array<T, 5> angles = {small, T(public_broadside),
                T(std::nextafter(public_broadside, 0.0)), T(std::nextafter(public_broadside, public_pi)),
                T(std::nextafter(public_pi, 0.0))};
            for (T angle : angles) {
                const auto direction = psms_direction_from_angle(angle);
                const auto first = psms_automatic::normalized_ferrers_values(1, direction, 2, workspace);
                const auto second = psms_automatic::normalized_ferrers_values(2, direction, 1, workspace);
                // Closed polynomials retain their sine factor even when cosine rounds to an endpoint.
                const T p11 = precsqrt(T(3) / T(4)) * direction.sine;
                const T p21 = precsqrt(T(15) / T(4)) * direction.cosine * direction.sine;
                const T p22 = precsqrt(T(15) / T(16)) * direction.sine * direction.sine;
                const T tolerance = T(128) * preceps<T>();
                expect_true(bool(precabs(first[0] - p11) <= tolerance * precabs(p11)));
                expect_true(bool(precabs(first[1] - p21) <= tolerance * precabs(p21)));
                expect_true(bool(precabs(second[0] - p22) <= tolerance * precabs(p22)));
                expect_true(bool(first[0] > T(0)));
                if (angle == T(public_broadside)) expect_true(bool(first[1] == T(0)));
                else expect_true(bool(first[1] != T(0)));
            }
            const T xi = T(7) / precsqrt(T(48));
            const auto omitted = psms_automatic::omitted_incident_bound(0, T(120), xi, near_axis);
            expect_true(bool(omitted.first > T(0) && omitted.second > T(0)));
            const auto axial = psms_automatic::omitted_incident_bound(0, T(120), xi,
                psms_direction_from_angle(T(0)));
            expect_true(bool(axial.first == T(0) && axial.second == T(0)));
        };
        check(double(0));
#if ACOUSTICTS_HAVE_QUADMATH
        check(acousticts_quad_t(0));
#endif
    }

    test_that("coefficient angular matrices preserve ordinary values and neighboring broadside modes") {
        auto check = [](auto scalar) {
            using T = decltype(scalar);
            const double half_pi = std::acos(-1.0) / 2;
            const std::vector<SpheroidalDirection<T>> directions = {
                psms_direction_from_angle(T(0.61)), psms_direction_from_angle(T(half_pi)),
                psms_direction_from_angle(T(std::nextafter(half_pi, 0.0))),
                psms_direction_from_angle(T(std::nextafter(half_pi, 2.0)))};
            const auto matrices = compute_smn_direction_matrices<T>(2, 6, T(3), directions);
            const auto native = compute_smn_matrix<T>(2, 6, T(3), directions[0].cosine, true);
            const T tolerance = pow10_typed<T>(is_quad_precision_v<T> ? -27 : -12);
            for (int m = 0; m <= 2; ++m) {
                for (int degree = m; degree <= 6; ++degree)
                    expect_true(bool(precabs(matrices[0][m][degree] - native[m][degree]) <=
                        tolerance * std::max(T(1), precabs(native[m][degree]))));
                expect_true(bool(matrices[1][m][m + 1] == T(0)));
                expect_true(bool(matrices[2][m][m + 1] != T(0)));
                expect_true(bool(matrices[3][m][m + 1] != T(0)));
                // An odd angular mode reverses across broadside without a finite zeroing window.
                expect_true(bool((matrices[2][m][m + 1] > T(0)) != (matrices[3][m][m + 1] > T(0))));
            }
        };
        check(double(0));
#if ACOUSTICTS_HAVE_QUADMATH
        check(acousticts_quad_t(0));
#endif
    }

    test_that("near-axis fluid, rigid, and release fields preserve bistatic reciprocity and nonzero odd orders") {
        auto check = [](auto scalar) {
            using T = decltype(scalar);
            const auto incidence = psms_direction_from_angle(T(1e-8));
            const auto reception = psms_direction_from_angle(T(37) * precpi<T>() / T(180));
            const SpheroidalDirection<T> reverse_incidence = {-reception.cosine, reception.sine, false};
            const SpheroidalDirection<T> reverse_reception = {-incidence.cosine, incidence.sine, false};
            const T xi = T(3) / precsqrt(T(8)), phi = T(0.37);
            for (bool fluid : {true, false}) {
                auto solve = [&](const SpheroidalDirection<T>& incoming,
                        const SpheroidalDirection<T>& outgoing, T azimuth) {
                    return fluid ? psms_automatic_full_fluid_field<T>(4, 24, T(2), T(8.5), xi,
                        incoming, outgoing, azimuth, T(1) / T(800), false) :
                        psms_automatic_rigid_release_field<T>(4, 24, T(2), xi,
                            incoming, outgoing, azimuth, false, false);
                };
                const auto forward = solve(incidence, reception, phi);
                const auto reciprocal = solve(reverse_incidence, reverse_reception, -phi);
                const auto opposite = solve(incidence, reception, phi + precpi<T>());
                CATCH_INFO("precision=" << (is_quad_precision_v<T> ? "quad" : "double") <<
                    " fluid=" << fluid << " forward=" << forward.detail << " reciprocal=" << reciprocal.detail <<
                    " opposite=" << opposite.detail);
                expect_true(bool(forward.converged && reciprocal.converged && opposite.converged));
                if (!forward.converged || !reciprocal.converged || !opposite.converged) continue;
                const T scale = precsqrt(std::max(complex_abs_sq(forward.field), complex_abs_sq(reciprocal.field)));
                expect_true(bool(precsqrt(complex_abs_sq(forward.field - reciprocal.field)) <=
                    T(8) * psms_automatic::numerical_budget(scale)));
                expect_true(bool(forward.orders.size() > 1));
                if (forward.orders.size() > 1)
                    expect_true(bool(psms_automatic::norm(forward.orders[1].coefficients) > T(0)));
                // Oblique reception exposes the O(sin(theta)) response hidden by monostatic symmetry.
                expect_true(bool(precsqrt(complex_abs_sq(forward.field - opposite.field)) >
                    T(100) * preceps<T>() * scale));
            }
        };
        check(double(0));
#if ACOUSTICTS_HAVE_QUADMATH
        check(acousticts_quad_t(0));
#endif
    }

    test_that("retained T-matrix single, grid and point evaluators preserve near-axis reception") {
        auto check = [](auto scalar) {
            using T = decltype(scalar);
            using Complex = std::complex<T>;
            // A single nonzero m=1 block gives an independent exact-axis zero and near-axis response.
            const std::vector<std::vector<Complex>> blocks = {
                std::vector<Complex>(4, Complex(0, 0)), {Complex(1, 0)}};
            const auto retained = prolate_tmatrix_blocks_to_rcpp<T>(1, 1, blocks);
            const T tiny = pow10_typed<T>(is_quad_precision_v<T> ? -20 : -8), incident = T(0.61), phi = T(0.23);
            const Complex single = prolate_scattering_from_retained_blocks_single<T>(
                retained, T(3), T(1), incident, T(0), tiny, phi);
            const Complex axial = prolate_scattering_from_retained_blocks_single<T>(
                retained, T(3), T(1), incident, T(0), T(0), phi);
            const auto grid = prolate_scattering_grid_from_retained_blocks_single<T>(
                retained, T(3), T(1), incident, T(0), {T(0), tiny}, {phi});
            const auto points = prolate_scattering_points_from_retained_blocks_single<T>(
                retained, T(3), T(1), incident, T(0), {T(0), tiny}, {phi, phi});
            expect_true(bool(single != Complex(0, 0)));
            expect_true(bool(axial == Complex(0, 0)));
            const Rcomplex grid_zero = grid(0, 0), point_zero = points[0];
            const Rcomplex grid_value = grid(1, 0), point_value = points[1];
            const std::complex<double> reference(static_cast<double>(single.real()), static_cast<double>(single.imag()));
            const double tolerance = 128 * std::numeric_limits<double>::epsilon() * std::abs(reference);
            expect_true(bool(grid_zero.r == 0 && grid_zero.i == 0 && point_zero.r == 0 && point_zero.i == 0));
            expect_true(bool(std::abs(std::complex<double>(grid_value.r, grid_value.i) - reference) <= tolerance));
            expect_true(bool(std::abs(std::complex<double>(point_value.r, point_value.i) - reference) <= tolerance));
        };
        check(double(0));
#if ACOUSTICTS_HAVE_QUADMATH
        check(acousticts_quad_t(0));
#endif
    }

    test_that("interior support brackets are reused without losing an independent comparison") {
        const auto unchanged = psms_automatic::reuse_interior_bracket(300, 0, 400, 500);
        const auto clipped = psms_automatic::reuse_interior_bracket(450, 0, 400, 500);
        const auto merged = psms_automatic::reuse_interior_bracket(500, 0, 400, 500);
        expect_true(bool(unchanged.first == 400 && unchanged.second == 500));
        expect_true(bool(clipped.first == 450 && clipped.second == 500));
        expect_true(bool(merged.first == 500 && merged.second > merged.first));
        expect_error(psms_automatic::reuse_interior_bracket(300, 0, 500, 500));
    }

    test_that("aggregate overflow cannot satisfy automatic numerical tolerances") {
        auto check = [](auto scalar) {
            using T = decltype(scalar);
            using Complex = std::complex<T>;
            const T infinity = T(std::numeric_limits<double>::infinity());
            T largest = std::numeric_limits<T>::max();
#if ACOUSTICTS_HAVE_QUADMATH
            if constexpr (std::is_same_v<T, acousticts_quad_t>) largest = nextafterq(infinity, T(0));
#endif
            const std::vector<Complex> overflowing = {Complex(largest, 0), Complex(largest, 0)};
            expect_error(psms_automatic::norm(overflowing));
            expect_error(psms_automatic::combine_norms(largest, largest));
            expect_error(psms_automatic::coefficient_budget(1, T(0), largest));
            psms_automatic::TailNorm<T> tail;
            for (const auto& value : overflowing) tail.append(value);
            expect_error(tail.value());
            Complex field(largest, 0), correction(0, 0);
            expect_error(psms_automatic::add_field(field, correction, Complex(largest, 0)));
            psms_automatic::BoundaryCheck<T> boundary;
            boundary.pressure_error = infinity;
            const std::pair<T, T> infinite_budgets = {infinity, infinity};
            expect_error(psms_automatic::boundary_passes(boundary, infinite_budgets));
            boundary.pressure_reference = infinity;
            expect_error(psms_automatic::boundary_budgets(boundary, 0, T(1), T(1)));
            // A temporarily unavailable tail bound is distinct from an overflowing computed norm.
            expect_true(bool(psms_automatic::combine_omitted_bound(T(1), infinity) == infinity));
            expect_error(psms_automatic::combine_omitted_bound(T(1), precnan<T>()));
            const T finite_combination = psms_automatic::combine_norms(largest / T(4), largest / T(4));
            expect_true(bool(is_finite_real(finite_combination) && finite_combination > largest / T(4)));
            expect_true(bool(is_finite_real(psms_automatic::numerical_budget(largest))));
        };
        check(double(0));
#if ACOUSTICTS_HAVE_QUADMATH
        check(acousticts_quad_t(0));
#endif
    }

    test_that("gas refinement recovers off-axis and overflow-range complex fields") {
        const double focal = std::sqrt(0.07 * 0.07 - 0.01 * 0.01);
        const double xi = 0.07 / focal, density = 1.24 / 1026.8;
        for (int high = 0; high < 2; ++high) {
            const double chi = high ? 120.0 : 20.0;
            const double eta = high ? 1.0 : std::cos(37.0 * std::acos(-1.0) / 180.0);
            const auto result = psms_automatic_full_fluid_field<double>(high ? 35 : 11,
                high ? 398 : 73, chi, chi * 1477.3 / 345.0, xi, eta, -eta, 0.0, density, true);
            expect_true(bool(result.converged));
            if (!result.converged) continue;
            // Stored quad references agree with enlarged original and Robin boundary equations.
            const std::complex<double> reference = high ?
                std::complex<double>(0.000758862317322051, -0.000210950066298305) :
                std::complex<double>(-0.00174320476740184, 0.000722475327611388);
            const auto amplitude = std::complex<double>(0, -2 * focal / chi) * result.field;
            expect_true(bool(std::abs(amplitude - reference) <= 1e-8 * std::abs(reference)));
            expect_true(bool(result.pressure_relative < 1.01e-8));
            expect_true(bool(result.derivative_relative < 1.01e-8));
            expect_true(bool(!result.orders.empty()));
            if (high) {
                expect_true(bool(result.orders.size() == 1));
                expect_true(bool(result.orders.front().robin));
                expect_true(bool(result.orders.front().interior_comparison_lower >= result.orders.front().degree));
                expect_true(bool(result.orders.front().interior_comparison_upper == result.orders.front().interior_degree));
                expect_true(bool(result.orders.front().interior_comparison_lower <
                    result.orders.front().interior_comparison_upper));
                expect_true(bool(result.omitted_pressure_bound == 0));
                expect_true(bool(result.omitted_derivative_bound == 0));
                bool beyond_double_radial_range = false;
                for (long long exponent : result.orders.front().radial_exponents)
                    beyond_double_radial_range |= exponent > std::numeric_limits<double>::max_exponent10;
                expect_true(bool(beyond_double_radial_range));
            } else {
                expect_true(bool(result.last_order > 11));
                expect_true(bool(result.last_degree > 73));
            }
        }
    }

    test_that("fixed-rigid and pressure-release boundaries agree with independently enlarged modal sums off-axis") {
        const double xi = 0.07 / std::sqrt(0.07 * 0.07 - 0.01 * 0.01);
        const double eta = std::cos(37.0 * std::acos(-1.0) / 180.0);
        for (bool rigid : {false, true}) {
            const auto result = psms_automatic_rigid_release_field<double>(8, 28,
                20.0, xi, eta, -eta, 0.0, true, rigid);
            expect_true(bool(result.converged));
            if (!result.converged) continue;
            auto fixed = [&](int m, int n) {
                return rigid ? psms_backscatter_fixed_rigid<double>(m, n, 20.0, xi, eta, false) :
                    psms_backscatter_pressure_release<double>(m, n, 20.0, xi, eta, false);
            };
            // Two enlargements establish the legacy reference independently of automatic support.
            const auto reference = fixed(24, 72), enlarged = fixed(32, 96);
            expect_true(bool(std::abs(reference - enlarged) <= 1e-10 * std::abs(enlarged)));
            expect_true(bool(std::abs(result.field - enlarged) <= 1e-8 * std::abs(enlarged)));
            if (rigid) {
                expect_true(bool(result.derivative_relative < 1.01e-8));
                expect_true(bool(std::isnan(result.pressure_relative)));
            } else {
                expect_true(bool(result.pressure_relative < 1.01e-8));
                expect_true(bool(std::isnan(result.derivative_relative)));
            }
        }
    }

    test_that("matched media return analytic zero in each native precision") {
        auto check = [](auto scalar) {
            using T = decltype(scalar);
            const auto result = psms_automatic_full_fluid_field<T>(4, 10, T(20), T(20),
                T(1.1), T(0.3), T(-0.6), T(0.7), T(1), false);
            expect_true(bool(result.converged && result.matched_medium));
            expect_true(bool(result.field == std::complex<T>(0, 0)));
            expect_true(bool(result.orders.empty()));
            expect_true(bool(result.peak_working_bytes == 0));
        };
        check(double(0));
#if ACOUSTICTS_HAVE_QUADMATH
        check(acousticts_quad_t(0));
#endif
    }

    test_that("discarded Legendre tails remain included in the complete residual") {
        auto check = [](auto scalar) {
            using T = decltype(scalar);
            psms_automatic::OrderSolution<T> solution;
            for (auto& trace : solution.traces) trace.resize(13);
            const T epsilon = psms_modal_abs_tol<T>() / T(16);
            solution.traces[4][0] = T(1) / precsqrt(T(2));
            solution.traces[0][12] = epsilon;
            psms_automatic::Quadrature<T> quadrature;
            // At chi=0 the incident trace is exactly 1/2 and the residual is epsilon times normalized P12.
            psms_automatic::check_boundary_quadrature(solution, T(0), T(1.1), T(0.3),
                T(1) / precsqrt(T(2)), T(0), quadrature, true, false);
            expect_true(bool(solution.boundary.quadrature_converged));
            expect_true(bool(solution.boundary.trace_prefix == 1));
            expect_true(bool(solution.boundary.pressure_tail_bound == epsilon));
            expect_true(bool(solution.boundary.pressure_error >= epsilon));
            expect_true(bool(solution.boundary.incident_pressure_tail_bound == T(0)));
            expect_true(bool(solution.traces[0].size() == 13));
            expect_true(bool(solution.traces[0][12].real() == epsilon));
        };
        check(double(0));
#if ACOUSTICTS_HAVE_QUADMATH
        check(acousticts_quad_t(0));
#endif
    }

    test_that("native incident Miller sequences preserve anchors, recurrences and closure") {
        auto check = [](auto scalar) {
            using T = decltype(scalar);
            const std::array<T, 4> arguments = {T(1), precpi<T>(), T(20), T(120)};
            const std::array<int, 4> degrees = {32, 64, 128, 256};
            for (size_t case_index = 0; case_index < arguments.size(); ++case_index) {
                const T argument = arguments[case_index];
                const int degree = degrees[case_index];
                psms_automatic::IncidentSequence<T> sequence;
                sequence.reset(argument);
                const T budget = psms_modal_abs_tol<T>() / T(32);
                psms_automatic::prepare_incident_sequence(sequence, degree,
                    T(1), T(0), budget, budget, true, true);
                const T tolerance = T(1024 * (degree + 1)) * preceps<T>();
                const T first = psms_automatic::sine(argument) / argument;
                const T second = (first - preccos(argument)) / argument;
                expect_true(bool(precabs(sequence.values[0] - first) <= tolerance));
                expect_true(bool(precabs(sequence.values[1] - second) <= tolerance));
                expect_true(bool(sequence.derivatives[0] == -sequence.values[1]));
                expect_true(bool(sequence.pressure_change <= budget));
                expect_true(bool(sequence.radial_change <= budget));
                T value_norm = 0, value_correction = 0, derivative_norm = 0, derivative_correction = 0;
                for (int index = 0; index <= degree; ++index) {
                    const T weight = T(2) * index + T(1);
                    psms_automatic::add(value_norm, value_correction,
                        weight * sequence.values[index] * sequence.values[index]);
                    psms_automatic::add(derivative_norm, derivative_correction,
                        weight * sequence.derivatives[index] * sequence.derivatives[index]);
                    if (index == 0 || index == degree) continue;
                    const T recurrence = sequence.values[index - 1] + sequence.values[index + 1] -
                        weight * sequence.values[index] / argument;
                    const T derivative = sequence.values[index - 1] -
                        T(index + 1) * sequence.values[index] / argument;
                    expect_true(bool(precabs(recurrence) <= tolerance));
                    expect_true(bool(precabs(sequence.derivatives[index] - derivative) <= tolerance));
                }
                // DLMF 10.60.12 and 10.60.14 supply independent infinite-series closure identities.
                const T value_tail = psms_automatic::incident_degree_tail(sequence, degree + 1);
                const T derivative_tail = psms_automatic::incident_degree_tail(sequence, degree);
                expect_true(bool(precabs(value_norm - T(1)) <= tolerance + T(2) * value_tail * value_tail));
                expect_true(bool(precabs(derivative_norm - T(1) / T(3)) <=
                    tolerance + T(2) * derivative_tail * derivative_tail));
            }
        };
        check(double(0));
#if ACOUSTICTS_HAVE_QUADMATH
        check(acousticts_quad_t(0));
#endif
    }

    test_that("analytic incident boundary norms agree with independent refined quadrature") {
        auto check = [](auto scalar) {
            using T = decltype(scalar);
            using Complex = std::complex<T>;
            const T xi = T(7) / precsqrt(T(48));
            const T oblique = preccos(T(37) * precpi<T>() / T(180));
            const std::array<T, 6> parameters = {T(20), T(120), T(120), T(120), T(20), T(120)};
            const std::array<SpheroidalDirection<T>, 6> directions = {
                psms_direction_from_cosine(oblique), psms_direction_from_cosine(T(1)),
                psms_direction_from_cosine(T(0)), psms_direction_from_cosine(oblique),
                psms_direction_from_cosine(-oblique), psms_direction_from_angle(T(1e-8))};
            const std::array<int, 6> orders = {1, 0, 2, 1, 1, 1};
            for (size_t case_index = 0; case_index < parameters.size(); ++case_index) {
                const T chi = parameters[case_index];
                const auto direction = directions[case_index];
                const int m = orders[case_index];
                const T sine_squared = direction.sine * direction.sine;
                const T pressure_reference = T(1) / precsqrt(T(2));
                const T derivative_reference = chi * precsqrt((direction.cosine * direction.cosine +
                    xi * xi * sine_squared / ((xi - T(1)) * (xi + T(1)))) / T(6));
                const auto budgets = psms_automatic::boundary_budgets(psms_automatic::BoundaryCheck<T>(),
                    m, pressure_reference, derivative_reference);
                psms_automatic::Quadrature<T> coefficient_workspace, quadrature_workspace;
                const auto incident = psms_automatic::incident_coefficients(m, chi, xi, direction,
                    budgets.first, budgets.second, coefficient_workspace);
                psms_automatic::OrderSolution<T> analytic;
                analytic.m = m;
                for (size_t index = 0; index < analytic.traces.size(); ++index)
                    analytic.traces[index] = index % 2 ? incident.derivative : incident.pressure;
                // A known orthonormal perturbation exercises nonzero pressure and derivative residuals.
                const T pressure_error = budgets.first / T(8), derivative_error = budgets.second / T(8);
                analytic.traces[0][0] += Complex(T(3) * pressure_error / T(5), T(4) * pressure_error / T(5));
                analytic.traces[1][0] += Complex(T(4) * derivative_error / T(5), T(3) * derivative_error / T(5));
                auto quadrature = analytic;
                const auto original_traces = analytic.traces;
                psms_automatic::check_boundary_coefficients(analytic, chi, xi, direction,
                    pressure_reference, derivative_reference, coefficient_workspace);
                psms_automatic::check_boundary_quadrature(quadrature, chi, xi, direction,
                    pressure_reference, derivative_reference, quadrature_workspace);
                expect_true(bool(analytic.boundary.quadrature_converged));
                expect_true(bool(quadrature.boundary.quadrature_converged));
                expect_true(bool(analytic.boundary.residual_method == "analytic_coefficients"));
                expect_true(bool(quadrature.boundary.residual_method == "quadrature"));
                expect_true(bool(analytic.boundary.quadrature_order == 0));
                expect_true(bool(analytic.boundary.incident_degree >= m));
                expect_true(bool(analytic.boundary.incident_recurrence_start > analytic.boundary.incident_degree));
                expect_true(bool(analytic.traces == original_traces));
                expect_true(bool(quadrature.traces == original_traces));
                const T rounding = T(10000) * preceps<T>();
                const T pressure_allowance = analytic.boundary.incident_pressure_tail_bound +
                    analytic.boundary.incident_pressure_recurrence_error +
                    quadrature.boundary.pressure_tail_bound + quadrature.boundary.quadrature_pressure_change + rounding;
                const T derivative_allowance = analytic.boundary.incident_derivative_tail_bound +
                    analytic.boundary.incident_derivative_recurrence_error +
                    quadrature.boundary.derivative_tail_bound + quadrature.boundary.quadrature_derivative_change +
                    rounding * std::max(T(1), derivative_reference);
                const T analytic_pressure_center = std::max(T(0), analytic.boundary.pressure_error -
                    analytic.boundary.incident_pressure_tail_bound - analytic.boundary.incident_pressure_recurrence_error);
                const T analytic_derivative_center = std::max(T(0), analytic.boundary.derivative_error -
                    analytic.boundary.incident_derivative_tail_bound - analytic.boundary.incident_derivative_recurrence_error);
                // Compare norm centers so each already recorded tail allowance is counted once.
                auto quadrature_center = [](T upper, T tail, T change) {
                    const T unpadded = std::max(T(0), upper - tail);
                    return precsqrt(std::max(T(0), unpadded * unpadded - change * change));
                };
                const T quadrature_pressure_center = quadrature_center(quadrature.boundary.pressure_error,
                    quadrature.boundary.pressure_tail_bound, quadrature.boundary.quadrature_pressure_change);
                const T quadrature_derivative_center = quadrature_center(quadrature.boundary.derivative_error,
                    quadrature.boundary.derivative_tail_bound, quadrature.boundary.quadrature_derivative_change);
                auto native_text = [](T value) {
                    char buffer[128];
#if ACOUSTICTS_HAVE_QUADMATH
                    if constexpr (std::is_same_v<T, acousticts_quad_t>)
                        quadmath_snprintf(buffer, sizeof(buffer), "%.36Qg", value);
                    else
#endif
                        std::snprintf(buffer, sizeof(buffer), "%.17g", static_cast<double>(value));
                    return std::string(buffer);
                };
                CATCH_INFO("precision=" << (is_quad_precision_v<T> ? "quad" : "double") <<
                    " chi=" << native_text(chi) << " direction=" << native_text(direction.cosine) <<
                    " sine=" << native_text(direction.sine) << " m=" << m <<
                    " analytic_derivative=" << native_text(analytic.boundary.derivative_error) <<
                    " quadrature_derivative=" << native_text(quadrature.boundary.derivative_error) <<
                    " analytic_center=" << native_text(analytic_derivative_center) <<
                    " quadrature_center=" << native_text(quadrature_derivative_center) <<
                    " expected_derivative=" << native_text(derivative_error) <<
                    " allowance=" << native_text(derivative_allowance) <<
                    " incident_tail=" << native_text(analytic.boundary.incident_derivative_tail_bound) <<
                    " recurrence=" << native_text(analytic.boundary.incident_derivative_recurrence_error) <<
                    " trace_tail=" << native_text(quadrature.boundary.derivative_tail_bound) <<
                    " quadrature_change=" << native_text(quadrature.boundary.quadrature_derivative_change) <<
                    " quadrature_order=" << quadrature.boundary.quadrature_order <<
                    " trace_prefix=" << quadrature.boundary.trace_prefix <<
                    " incident_degree=" << analytic.boundary.incident_degree);
                expect_true(bool(precabs(analytic_pressure_center - quadrature_pressure_center) <=
                    pressure_allowance));
                expect_true(bool(precabs(analytic_derivative_center - quadrature_derivative_center) <=
                    derivative_allowance));
                expect_true(bool(precabs(analytic.boundary.pressure_error - pressure_error) <= pressure_allowance));
                expect_true(bool(precabs(analytic.boundary.derivative_error - derivative_error) <= derivative_allowance));
            }
        };
        check(double(0));
#if ACOUSTICTS_HAVE_QUADMATH
        check(acousticts_quad_t(0));
#endif
    }

    test_that("an unavailable analytic incident expansion reports the quadrature fallback") {
        auto check = [](auto scalar) {
            using T = decltype(scalar);
            psms_automatic::OrderSolution<T> solution;
            for (auto& trace : solution.traces) trace.resize(1);
            solution.traces[0][0] = solution.traces[2][0] = solution.traces[4][0] = T(1) / precsqrt(T(2));
            psms_automatic::Quadrature<T> workspace;
            psms_automatic::check_boundary(solution, T(0), T(1.1), T(0.3),
                T(1) / precsqrt(T(2)), T(0), workspace, true, false);
            expect_true(bool(solution.boundary.quadrature_converged));
            expect_true(bool(solution.boundary.residual_method == "quadrature"));
            expect_true(bool(!solution.boundary.analytic_failure.empty()));
            expect_true(bool(solution.boundary.incident_degree == -1));
        };
        check(double(0));
#if ACOUSTICTS_HAVE_QUADMATH
        check(acousticts_quad_t(0));
#endif
    }

    test_that("high-order native quadrature preserves analytic moments and distinct nodes") {
        auto check = [](auto scalar) {
            using T = decltype(scalar);
            const int count = 257;
            psms_automatic::Quadrature<T> quadrature;
            quadrature.prepare(count);
            const T tolerance = T(4 * count) * preceps<T>();
            for (int index = 0; index < count; ++index) {
                expect_true(bool(quadrature.nodes[index] > T(-1) && quadrature.nodes[index] < T(1)));
                expect_true(bool(quadrature.weights[index] > T(0)));
                expect_true(bool(is_finite_real(quadrature.weights[index])));
                if (index > 0) expect_true(bool(quadrature.nodes[index] > quadrature.nodes[index - 1]));
            }
            expect_true(bool(quadrature.nodes[count / 2] == T(0)));
            // Analytic integrals check final native roots and weights independently of their initial guesses.
            for (int degree : {0, 2, 4, 8}) {
                T integral = 0, correction = 0;
                for (int index = 0; index < count; ++index) {
                    T polynomial = 1;
                    for (int power = 0; power < degree; ++power) polynomial *= quadrature.nodes[index];
                    psms_automatic::add(integral, correction, quadrature.weights[index] * polynomial);
                }
                expect_true(bool(precabs(integral - T(2) / T(degree + 1)) <= tolerance));
            }
        };
        check(double(0));
#if ACOUSTICTS_HAVE_QUADMATH
        check(acousticts_quad_t(0));
#endif
    }

    test_that("adaptive Gauss exactness agrees with complete finite Ferrers coefficient overlaps") {
        auto check = [](auto scalar) {
            using T = decltype(scalar);
            auto make_batch = [](int capacity, std::vector<int> lengths, int shift) {
                ProfcnCoefficientBatchResult<T> batch;
                batch.capacity = capacity;
                batch.functions.lnum = int(lengths.size());
                batch.lengths = std::move(lengths);
                batch.coefficients.resize(size_t(capacity) * batch.functions.lnum);
                for (int mode = 0; mode < batch.functions.lnum; ++mode)
                    for (int index = mode % 2; index < batch.lengths[mode]; index += 2) {
                        const T sign = (index / 2 + mode + shift) % 2 ? T(-1) : T(1);
                        batch.coefficients[size_t(mode) * capacity + index] =
                            sign * T((index + 2 * mode + shift) % 7 + 1) / T(8);
                    }
                return batch;
            };
            const auto exterior = make_batch(9, {9, 8, 7, 6}, 1);
            const auto interior = make_batch(12, {11, 12, 9, 10}, 2);
            const std::vector<T> empty;
            psms_automatic::OverlapQuadrature<T> adaptive;
            adaptive.supplied_nodes = &empty;
            adaptive.supplied_weights = &empty;
            adaptive.adaptive = true;
            adaptive.validate();
            for (int m : {0, 1, 2, 3}) {
                adaptive.prepare(m, 9, 12);
                const int product_degree = 2 * m + 9 + 12 - 2;
                expect_true(bool(2 * adaptive.exact_order - 1 >= product_degree));
                expect_true(bool(2 * (adaptive.exact_order - 1) - 1 < product_degree));
                expect_true(bool(adaptive.bytes() == 0 && adaptive.nodes().empty() && adaptive.weights().empty()));
                // Evaluate actual native Gauss nodes independently of the coefficient identity.
                psms_automatic::Quadrature<T> gauss;
                gauss.prepare(adaptive.exact_order);
                psms_automatic::OverlapQuadrature<T> prescribed;
                prescribed.supplied_nodes = &gauss.nodes;
                prescribed.supplied_weights = &gauss.weights;
                prescribed.validate();
                const auto sampled = psms_automatic::sample_overlaps(m, exterior, interior, prescribed, false);
                for (int row = 0; row < 4; ++row) {
                    for (int col = 0; col < 4; ++col) {
                        T expected = 0, product_norm = 0;
                        for (int index = 0; index < std::min(exterior.lengths[row], interior.lengths[col]); ++index) {
                            const T term = exterior.coefficients[size_t(row) * 9 + index] *
                                interior.coefficients[size_t(col) * 12 + index];
                            expected += term;
                            product_norm += precabs(term);
                        }
                        const T identity = spheroidal_coefficient_overlap(exterior, interior, 0, row, col);
                        const T tolerance = T(256) * preceps<T>() * std::max(T(1), product_norm);
                        expect_true(bool(identity == expected));
                        expect_true(bool(precabs(sampled.integral(row, col) - identity) <= tolerance));
                    }
                }
            }
            // The highest retained basis polynomial aliases to zero with one too few Gauss nodes.
            ProfcnCoefficientBatchResult<T> highest;
            highest.capacity = 5;
            highest.functions.lnum = 1;
            highest.lengths = {5};
            highest.coefficients = {T(0), T(0), T(0), T(0), T(1)};
            adaptive.prepare(0, 5, 5);
            expect_true(bool(adaptive.exact_order == 5));
            for (int points : {4, 5}) {
                psms_automatic::Quadrature<T> gauss;
                gauss.prepare(points);
                psms_automatic::OverlapQuadrature<T> prescribed;
                prescribed.supplied_nodes = &gauss.nodes;
                prescribed.supplied_weights = &gauss.weights;
                const auto sampled = psms_automatic::sample_overlaps(0, highest, highest, prescribed, false);
                expect_true(bool(precabs(sampled.integral(0, 0) - T(points == 5 ? 1 : 0)) <=
                    T(256) * preceps<T>()));
            }
            // Formal exactness orders are checked in wide integers without allocating a rule.
            adaptive.prepare(0, INT_MAX - 1, INT_MAX - 1);
            expect_true(bool(adaptive.exact_order == INT_MAX - 1 && adaptive.bytes() == 0));
            expect_error(adaptive.prepare(0, INT_MAX, INT_MAX));
            expect_error(adaptive.prepare(INT_MAX, 1, 1));
            expect_error(adaptive.prepare(-1, 1, 1));
            expect_error(adaptive.prepare(0, 0, 1));
        };
        check(double(0));
#if ACOUSTICTS_HAVE_QUADMATH
        check(acousticts_quad_t(0));
#endif
    }

    test_that("overlap reuse verifies complete native columns and invalidates changed support") {
        auto check = [](auto scalar) {
            using T = decltype(scalar);
            auto batch = [](int modes, int capacity, int shift) {
                ProfcnCoefficientBatchResult<T> result;
                result.capacity = capacity;
                result.functions.lnum = modes;
                result.functions.m_count = 1;
                result.lengths.resize(modes, 7);
                result.coefficients.resize(size_t(modes) * capacity);
                for (int mode = 0; mode < modes; ++mode)
                    for (int index = mode % 2; index < 7; index += 2)
                        result.coefficients[size_t(mode) * capacity + index] =
                            T((mode + index + shift) % 9 - 4) / T(8);
                return result;
            };
            auto exterior = batch(4, 9, 1), interior = batch(4, 10, 2);
            SpheroidalCoefficientOverlapCache<T> cache;
            auto compare = [&](int m, int parity, const auto& external, const auto& internal,
                    size_t reused, size_t computed) {
                const int columns = (external.functions.lnum + 1 - parity) / 2;
                const int rows = (internal.functions.lnum + 1 - parity) / 2;
                std::vector<T> values(size_t(rows) * columns);
                cache.fill(m, parity, external, internal, values);
                expect_true(bool(cache.reused_entries == reused && cache.computed_entries == computed));
                expect_true(bool(cache.peak_bytes >= cache.bytes()));
                expect_true(bool(cache.bytes() <= cache.budget));
                for (int row = 0; row < rows; ++row) {
                    for (int col = 0; col < columns; ++col) {
                        T expected = 0;
                        const int ei = parity + 2 * col, ii = parity + 2 * row;
                        for (int index = parity; index < std::min(external.lengths[ei], internal.lengths[ii]); index += 2)
                            expected += external.coefficients[size_t(ei) * external.capacity + index] *
                                internal.coefficients[size_t(ii) * internal.capacity + index];
                        expect_true(bool(precabs(values[size_t(row) * columns + col] - expected) <=
                            T(8) * preceps<T>() * std::max(T(1), precabs(expected))));
                        const T uncached = spheroidal_coefficient_overlap(external, internal, 0, ei, ii);
                        expect_true(bool(std::memcmp(&values[size_t(row) * columns + col], &uncached, sizeof(T)) == 0));
                    }
                }
            };
            compare(2, 0, exterior, interior, 0, 4);
            expect_true(bool(cache.blocks[1].overlaps.empty()));
            compare(2, 0, exterior, interior, 4, 0);
            // Repacking capacity preserves exact columns, while new degrees add only new intersections.
            exterior = batch(6, 12, 1);
            interior = batch(6, 13, 2);
            compare(2, 0, exterior, interior, 4, 5);
            compare(2, 1, exterior, interior, 0, 9);
            compare(2, 1, exterior, interior, 9, 0);
            exterior.coefficients[size_t(2) * exterior.capacity + 6] += T(1) / T(8);
            compare(2, 0, exterior, interior, 6, 3);
            interior.coefficients[size_t(4) * interior.capacity + 6] -= T(1) / T(8);
            compare(2, 0, exterior, interior, 6, 3);
            // Even an exact-zero extension changes the exported support and invalidates that column.
            exterior.lengths[0] = 9;
            compare(2, 0, exterior, interior, 6, 3);
            compare(2, 0, exterior, interior, 9, 0);
            compare(3, 0, exterior, interior, 0, 9);
            expect_true(bool(cache.blocks[1].overlaps.empty()));
            cache.clear();
            expect_true(bool(cache.order == -1 && cache.bytes() == sizeof(cache)));
            compare(3, 0, exterior, interior, 0, 9);
            auto changed = exterior;
            changed.coefficients[0] += preceps<T>();
            compare(3, 0, changed, interior, 6, 3);
            // A zero retention budget computes the same complete products without retaining arrays.
            SpheroidalCoefficientOverlapCache<T> disabled(0);
            std::vector<T> values(9), repeated(9);
            disabled.fill(3, 0, exterior, interior, values);
            disabled.fill(3, 0, exterior, interior, repeated);
            expect_true(bool(values == repeated && disabled.reused_entries == 0 && disabled.computed_entries == 9));
            expect_true(bool(disabled.bytes() == sizeof(disabled)));
            expect_error(disabled.fill(-1, 0, exterior, interior, values));
            expect_error(disabled.fill(3, 2, exterior, interior, values));
            changed.lengths[0] = changed.capacity + 1;
            expect_error(disabled.fill(3, 0, changed, interior, values));
            size_t bytes = 0;
            expect_true(bool(!SpheroidalCoefficientOverlapCache<T>::reserve_bytes(
                std::numeric_limits<size_t>::max(), sizeof(T), 1024, bytes)));
        };
        check(double(0));
#if ACOUSTICTS_HAVE_QUADMATH
        check(acousticts_quad_t(0));
#endif
    }

    test_that("coupled orders preserve fields with overlap reuse and bypass prescribed quadrature") {
        auto check = [](auto scalar) {
            using T = decltype(scalar);
            const T chi = T(4), interior_chi = T(7) / T(2), xi = T(3) / T(2);
            psms_automatic::Quadrature<T> workspace;
            auto exterior = cprofcn_mblock_coefficients<T>(chi, 1, 1, 16, {T(0)}, 2, 0, xi - T(1), true);
            psms_automatic::evaluate_coefficient_angular_values(exterior, 1,
                std::vector<SpheroidalDirection<T>>{psms_direction_from_cosine(T(0))}, workspace.legendre);
            SpheroidalCoefficientOverlapCache<T> cache;
            const auto cold = psms_automatic::solve_order<T>(1, 16, 16, interior_chi, xi,
                T(0), T(0), T(11) / T(10), true, false, exterior, nullptr, &cache);
            const auto warm = psms_automatic::solve_order<T>(1, 16, 16, interior_chi, xi,
                T(0), T(0), T(11) / T(10), true, false, exterior, nullptr, &cache);
            const auto uncached = psms_automatic::solve_order<T>(1, 16, 16, interior_chi, xi,
                T(0), T(0), T(11) / T(10), true, false, exterior);
            expect_true(bool(!cold.singular && !warm.singular && !uncached.singular));
            expect_true(bool(cold.field == uncached.field && warm.field == uncached.field));
            expect_true(bool(warm.coefficients == uncached.coefficients && warm.traces == uncached.traces));
            expect_true(bool(cold.overlap_entries_computed == 64 && warm.overlap_entries_reused == 64));
            expect_true(bool(warm.overlap_entries_computed == 0 && cache.blocks[1].overlaps.empty()));
            expect_true(bool(warm.working_bytes >= uncached.working_bytes + cache.bytes()));
            // Omitting diagnostic traces preserves every solved coefficient and the complete field.
            const auto field_only = psms_automatic::solve_order<T>(1, 16, 16, interior_chi, xi,
                T(0), T(0), T(11) / T(10), true, false, exterior, nullptr, nullptr, false, true);
            expect_true(bool(!field_only.singular && field_only.field == uncached.field));
            expect_true(bool(field_only.coefficients == uncached.coefficients &&
                field_only.scaled_coefficients == uncached.scaled_coefficients &&
                field_only.radial_exponents == uncached.radial_exponents));
            expect_true(bool(field_only.outgoing_terms.size() == 16 &&
                field_only.working_bytes < uncached.working_bytes));
            for (const auto& trace : field_only.traces) expect_true(bool(trace.empty()));
            // Direction changes preserve overlap keys while rebuilding forcing and both parity solves.
            psms_automatic::evaluate_coefficient_angular_values(exterior, 1,
                std::vector<SpheroidalDirection<T>>{psms_direction_from_cosine(T(0.3))}, workspace.legendre);
            const auto oblique = psms_automatic::solve_order<T>(1, 16, 16, interior_chi, xi,
                T(0.3), T(0), T(11) / T(10), true, false, exterior, nullptr, &cache);
            const auto oblique_reference = psms_automatic::solve_order<T>(1, 16, 16, interior_chi, xi,
                T(0.3), T(0), T(11) / T(10), true, false, exterior);
            expect_true(bool(!oblique.singular && !oblique_reference.singular));
            expect_true(bool(oblique.field == oblique_reference.field &&
                oblique.coefficients == oblique_reference.coefficients && oblique.traces == oblique_reference.traces));
            expect_true(bool(oblique.overlap_entries_reused == 64 && oblique.overlap_entries_computed == 64));
            const auto robin = psms_automatic::solve_order<T>(1, 16, 20, interior_chi, xi,
                T(0.3), T(0), T(11) / T(10), true, true, exterior, nullptr, &cache);
            const auto robin_reference = psms_automatic::solve_order<T>(1, 16, 20, interior_chi, xi,
                T(0.3), T(0), T(11) / T(10), true, true, exterior);
            expect_true(bool(!robin.singular && !robin_reference.singular));
            expect_true(bool(robin.field == robin_reference.field &&
                robin.coefficients == robin_reference.coefficients && robin.traces == robin_reference.traces));
            // Enlarged native blocks may change old columns, so compare results without assuming reuse.
            expect_true(bool(robin.overlap_entries_reused + robin.overlap_entries_computed == 160));
            workspace.prepare(24);
            psms_automatic::OverlapQuadrature<T> prescribed;
            prescribed.supplied_nodes = &workspace.nodes;
            prescribed.supplied_weights = &workspace.weights;
            cache.clear();
            const auto sampled = psms_automatic::solve_order<T>(1, 16, 16, interior_chi, xi,
                T(0.3), T(0), T(11) / T(10), true, false, exterior, &prescribed, &cache);
            const auto reference = psms_automatic::solve_order<T>(1, 16, 16, interior_chi, xi,
                T(0.3), T(0), T(11) / T(10), true, false, exterior, &prescribed);
            expect_true(bool(!sampled.singular && !reference.singular && sampled.field == reference.field));
            expect_true(bool(sampled.coefficients == reference.coefficients && sampled.traces == reference.traces));
            expect_true(bool(sampled.overlap_entries_reused == 0 && sampled.overlap_entries_computed == 128));
            expect_true(bool(cache.order == -1 && cache.blocks[0].overlaps.empty()));
            const auto sampled_field_only = psms_automatic::solve_order<T>(1, 16, 16, interior_chi, xi,
                T(0.3), T(0), T(11) / T(10), true, false, exterior, &prescribed, nullptr, false);
            expect_true(bool(!sampled_field_only.singular && sampled_field_only.field == reference.field));
            expect_true(bool(sampled_field_only.coefficients == reference.coefficients));
            for (const auto& trace : sampled_field_only.traces) expect_true(bool(trace.empty()));
            // Both boundary types retain packed coefficients without trace exports.
            for (bool rigid : {false, true}) {
                const auto complete = psms_automatic::solve_rigid_release_order(1, 16, T(0.3), T(0),
                    true, rigid, exterior);
                const auto field_only_boundary = psms_automatic::solve_rigid_release_order(1, 16, T(0.3), T(0),
                    true, rigid, exterior, false, true);
                expect_true(bool(field_only_boundary.field == complete.field &&
                    field_only_boundary.coefficients == complete.coefficients &&
                    field_only_boundary.scaled_coefficients == complete.scaled_coefficients &&
                    field_only_boundary.radial_exponents == complete.radial_exponents));
                expect_true(bool(field_only_boundary.outgoing_terms.size() == 16 &&
                    field_only_boundary.working_bytes < complete.working_bytes));
                for (const auto& trace : field_only_boundary.traces) expect_true(bool(trace.empty()));
            }
        };
        check(double(0));
#if ACOUSTICTS_HAVE_QUADMATH
        check(acousticts_quad_t(0));
#endif
    }

    test_that("adaptive degree tails retain late contributions and cancellation") {
        auto check = [](auto scalar) {
            using T = decltype(scalar);
            using Complex = std::complex<T>;
            using namespace psms_automatic;
            const T epsilon = psms_modal_abs_tol<T>() / T(1024);
            auto complete = [](OrderSolution<T>& order) {
                Complex correction(0, 0);
                order.field = Complex(0, 0);
                for (const auto& term : order.outgoing_terms) add_field(order.field, correction, term);
            };
            OrderSolution<T> order;
            order.outgoing_terms.assign(12, Complex(epsilon, epsilon));
            // A tiny leading prefix must not conceal a later real and imaginary lobe.
            order.outgoing_terms[10] = Complex(T(1), T(-2));
            complete(order);
            expect_true(bool(adaptive_outgoing_band(order) == order.field));
            // Opposing late terms require a budget based on the cancelled complete field.
            order.outgoing_terms[0] = Complex(T(-1), T(2));
            complete(order);
            expect_true(bool(adaptive_outgoing_band(order) == order.field));
            // A small omitted suffix may still be accepted within the unchanged native budget.
            order.outgoing_terms.assign(12, Complex(epsilon, T(0)));
            order.outgoing_terms[0] = Complex(T(1), T(0));
            order.outgoing_terms[10] = Complex(psms_modal_rel_tol<T>() / T(4), T(0));
            complete(order);
            const auto truncated = adaptive_outgoing_band(order);
            expect_true(bool(truncated != order.field));
            const auto difference = truncated - order.field;
            expect_true(bool(combine_norms(precabs(difference.real()), precabs(difference.imag())) <=
                numerical_budget(combine_norms(precabs(order.field.real()), precabs(order.field.imag())))));
            order.outgoing_terms.clear();
            expect_true(bool(adaptive_outgoing_band(order) == order.field));
        };
        check(double(0));
#if ACOUSTICTS_HAVE_QUADMATH
        check(acousticts_quad_t(0));
#endif
    }

    test_that("adaptive azimuthal sums guard later bands and accumulated degree omissions") {
        auto check = [](auto scalar) {
            using T = decltype(scalar);
            using Complex = std::complex<T>;
            using namespace psms_automatic;
            const T epsilon = psms_modal_abs_tol<T>() / T(1024);
            std::vector<OrderSolution<T>> orders(10);
            auto complete = [&]() {
                Complex value(0, 0), correction(0, 0);
                for (const auto& order : orders) add_field(value, correction, order.field);
                return value;
            };
            for (size_t index = 0; index < orders.size(); ++index) {
                orders[index].m = int(index);
                orders[index].field = Complex(epsilon, epsilon);
            }
            // Small early orders do not bound an unexamined later azimuthal contribution.
            orders.back().field = Complex(T(1), T(-2));
            expect_true(bool(adaptive_outgoing_sum(orders) == complete()));
            orders.front().field = -orders.back().field;
            expect_true(bool(adaptive_outgoing_sum(orders) == complete()));
            // Individually acceptable degree omissions can exceed the budget after band cancellation.
            for (size_t index = 0; index < orders.size(); ++index) {
                auto& order = orders[index];
                order.outgoing_terms.assign(12, Complex(epsilon, T(0)));
                order.outgoing_terms[0] = Complex(index % 2 ? T(-1) : T(1), T(0));
                order.outgoing_terms[10] = Complex(psms_modal_rel_tol<T>() / T(4), T(0));
                order.field = Complex(0, 0);
                Complex correction(0, 0);
                for (const auto& term : order.outgoing_terms) add_field(order.field, correction, term);
                expect_true(bool(adaptive_outgoing_band(order) != order.field));
            }
            expect_true(bool(adaptive_outgoing_sum(orders) == complete()));
            orders.clear();
            expect_true(bool(adaptive_outgoing_sum(orders) == Complex(0, 0)));
        };
        check(double(0));
#if ACOUSTICTS_HAVE_QUADMATH
        check(acousticts_quad_t(0));
#endif
    }

    test_that("axial h0=60 fixed-rigid and pressure-release boundaries retain fixed double and converged quad fields") {
        auto check = [](auto scalar) {
            using T = decltype(scalar);
            using namespace psms_automatic;
            const T chi = T(60), xi = T(7) / precsqrt(T(48));
            const std::vector<T> empty;
            // Quad support is checked against two independently enlarged fixed series.
            const int reference_degree = is_quad_precision_v<T> ? 160 : 49;
            for (const std::string method : {"Amn_pressure_release", "Amn_fixed_rigid"}) {
                const auto complete = psms_public_field<T>(18, reference_degree, chi, chi, xi,
                    T(0), T(0), T(0), T(0), T(1), T(1), empty, empty,
                    method, true, false, false);
                if constexpr (is_quad_precision_v<T>) {
                    const auto enlarged = psms_public_field<T>(18, 200, chi, chi, xi,
                        T(0), T(0), T(0), T(0), T(1), T(1), empty, empty,
                        method, true, false, false);
                    const T scale = combine_norms(precabs(enlarged.real()), precabs(enlarged.imag()));
                    const auto difference = complete - enlarged;
                    expect_true(bool(combine_norms(precabs(difference.real()), precabs(difference.imag())) <=
                        numerical_budget(scale)));
                }
                const auto adaptive = psms_public_field<T>(18, 49, chi, chi, xi,
                    T(0), T(0), T(0), T(0), T(1), T(1), empty, empty,
                    method, true, true, false);
                const T reference = combine_norms(precabs(complete.real()), precabs(complete.imag()));
                const auto difference = adaptive - complete;
                expect_true(bool(reference > numerical_budget(reference)));
                expect_true(bool(combine_norms(precabs(difference.real()), precabs(difference.imag())) <=
                    numerical_budget(reference)));
            }
        };
        check(double(0));
#if ACOUSTICTS_HAVE_QUADMATH
        check(acousticts_quad_t(0));
#endif
    }

    test_that("prescribed PSMS ceilings remain fixed outside adaptive quad backscatter") {
        auto check = [](auto scalar) {
            using T = decltype(scalar);
            const T chi = T(0.4), internal_chi = T(0.7), xi = T(1.5), density = T(1.1);
            psms_automatic::Quadrature<T> quadrature;
            quadrature.prepare(16);
            const double pi = std::acos(-1.0), receive_angle = 0.93;
            const auto reception = psms_direction_from_angle(T(receive_angle));
            for (double angle : {0.0, 0.61, pi / 2, pi}) {
                const auto incidence = psms_direction_from_angle(T(angle));
                const T cosine = incidence.cosine;
                auto exterior = cprofcn_mblock_coefficients<T>(chi, 0, 1, 1,
                    {cosine, reception.cosine}, 2, 0, xi - T(1), true);
                psms_automatic::evaluate_coefficient_angular_values(exterior, 0,
                    std::vector<SpheroidalDirection<T>>{incidence, reception}, quadrature.legendre);
                for (bool backscatter : {false, true}) {
                    // The single retained degree gives an independent support-sensitive reference.
                    const auto fluid = psms_automatic::solve_order<T>(0, 0, 0, internal_chi, xi,
                        cosine, T(0.7), density, backscatter, false, exterior);
                    expect_true(bool(!fluid.singular));
                    for (bool adaptive : {false, true}) {
                        // Adaptive quad backscatter selects support beyond the supplied seeds.
                        if (adaptive && is_quad_precision_v<T> && backscatter) continue;
                        const auto actual = psms_public_field<T>(0, 0, chi, internal_chi, xi,
                            T(angle), T(receive_angle), T(0), T(0.7), density, T(1),
                            quadrature.nodes, quadrature.weights, "Amn_fluid_gas", backscatter, adaptive, false);
                        expect_true(bool(actual == fluid.field));
                        if (adaptive && (!is_quad_precision_v<T> || !backscatter)) {
                            const auto liquid = psms_public_field<T>(0, 0, chi, internal_chi, xi,
                                T(angle), T(receive_angle), T(0), T(0.7), density, T(1),
                                quadrature.nodes, quadrature.weights, "Amn_fluid", backscatter, true, false);
                            expect_true(bool(liquid == fluid.field));
                        }
                        for (bool rigid : {false, true}) {
                            const auto boundary_result = psms_automatic::solve_rigid_release_order(0, 0, cosine,
                                T(0.7), backscatter, rigid, exterior);
                            const auto value = psms_public_field<T>(0, 0, chi, internal_chi, xi,
                                T(angle), T(receive_angle), T(0), T(0.7), density, T(1),
                                quadrature.nodes, quadrature.weights,
                                rigid ? "Amn_fixed_rigid" : "Amn_pressure_release", backscatter, adaptive, false);
                            expect_true(bool(value == boundary_result.field));
                        }
                    }
                    psms_automatic::OverlapQuadrature<T> prescribed;
                    prescribed.supplied_nodes = &quadrature.nodes;
                    prescribed.supplied_weights = &quadrature.weights;
                    const auto liquid = psms_automatic::solve_order<T>(0, 0, 0, internal_chi, xi,
                        cosine, T(0.7), density, backscatter, false, exterior, &prescribed);
                    const auto value = psms_public_field<T>(0, 0, chi, internal_chi, xi,
                        T(angle), T(receive_angle), T(0), T(0.7), density, T(1),
                        quadrature.nodes, quadrature.weights, "Amn_fluid", backscatter, false, false);
                    expect_true(bool(!liquid.singular && value == liquid.field));
                }
            }
            const auto incidence = psms_direction_from_cosine(T(0.3));
            const std::vector<T> empty;
            expect_error(psms_fixed_support_field<T>(0, 0, chi, internal_chi, xi,
                incidence, reception, T(0), density, true, false, false, false, &empty, &empty));
            const auto automatic_quadrature = psms_fixed_support_field<T>(0, 0, chi, internal_chi, xi,
                incidence, reception, T(0), density, true, false, false, true, &empty, &empty);
            const auto coefficient_overlap = psms_fixed_support_field<T>(0, 0, chi, internal_chi, xi,
                incidence, reception, T(0), density, true, false, false, true);
            expect_true(bool(automatic_quadrature == coefficient_overlap));
        };
        check(double(0));
#if ACOUSTICTS_HAVE_QUADMATH
        check(acousticts_quad_t(0));
#endif
    }

    test_that("prescribed overlap rules retain arbitrary signed weights and node order") {
        auto check = [](auto scalar) {
            using T = decltype(scalar);
            const std::vector<T> nodes = {T(1) / T(2), T(-1) / T(2), T(0), T(1),
                T(-1), T(1) / T(4), T(-1) / T(4), T(1) / T(2)};
            const std::vector<T> weights = {T(2), T(-1) / T(2), T(0), T(1) / T(8),
                T(3) / T(8), T(1) / T(4), T(-3) / T(4), T(5) / T(8)};
            ProfcnCoefficientBatchResult<T> batch;
            batch.capacity = 3;
            batch.functions.lnum = 3;
            batch.lengths = {1, 2, 3};
            batch.coefficients = {T(1), T(0), T(0), T(0), T(1), T(0),
                T(-1) / T(2), T(0), T(3) / T(4)};
            psms_automatic::OverlapQuadrature<T> prescribed;
            prescribed.supplied_nodes = &nodes;
            prescribed.supplied_weights = &weights;
            prescribed.validate();
            const auto samples = psms_automatic::sample_overlaps(0, batch, batch, prescribed, false);
            expect_true(bool(samples.count == nodes.size()));
            expect_true(bool(&prescribed.nodes() == &nodes && &prescribed.weights() == &weights));
            expect_true(bool(prescribed.rank_bound(0, 0) == 3 && prescribed.rank_bound(1, 0) == 2));
            std::array<std::vector<T>, 3> expected;
            for (auto& values : expected) values.resize(nodes.size());
            // Closed normalized P_0, P_1 and P_2 give independent values at every supplied point.
            const T tolerance = T(128) * preceps<T>();
            for (size_t node = 0; node < nodes.size(); ++node) {
                const T eta = nodes[node];
                expected[0][node] = T(1) / precsqrt(T(2));
                expected[1][node] = precsqrt(T(3) / T(2)) * eta;
                expected[2][node] = -expected[0][node] / T(2) + T(3) / T(4) *
                    precsqrt(T(5) / T(8)) * (T(3) * eta * eta - T(1));
                for (int mode = 0; mode < 3; ++mode) {
                    expect_true(bool(precabs(samples.exterior[size_t(mode) * nodes.size() + node] -
                        expected[mode][node]) <= tolerance));
                    expect_true(bool(precabs(samples.weighted_interior[size_t(mode) * nodes.size() + node] -
                        weights[node] * expected[mode][node]) <= tolerance));
                }
            }
            for (int row = 0; row < 3; ++row)
                for (int col = 0; col < 3; ++col) {
                    T integral = 0;
                    for (size_t node = 0; node < nodes.size(); ++node)
                        integral += weights[node] * expected[row][node] * expected[col][node];
                    expect_true(bool(precabs(samples.integral(row, col) - integral) <= tolerance));
                }
            expect_true(bool(precabs(samples.integral(0, 1)) > T(1) / T(4)));
        };
        check(double(0));
#if ACOUSTICTS_HAVE_QUADMATH
        check(acousticts_quad_t(0));
#endif
    }

    test_that("boundary reconstruction preserves normalized Ferrers values and orthogonality") {
        auto check = [](auto scalar) {
            using T = decltype(scalar);
            using Complex = std::complex<T>;
            psms_automatic::Quadrature<T> quadrature;
            quadrature.prepare(8);
            const T tolerance = T(100) * preceps<T>();
            for (int m : {0, 1}) {
                std::array<std::vector<Complex>, 6> traces;
                for (auto& trace : traces) trace.resize(3);
                traces[0][0] = T(1);
                traces[1][1] = T(1);
                traces[2][2] = T(1);
                quadrature.legendre.prepare(m, 3);
                const T eta = T(0.3), transverse = precsqrt(T(1) - eta * eta);
                const auto actual = psms_automatic::evaluate_traces(m, eta, traces, quadrature.legendre);
                // Closed low-degree Ferrers polynomials fix normalization and the native phase sign.
                const T first = m == 0 ? T(1) / precsqrt(T(2)) : precsqrt(T(3) / T(4)) * transverse;
                const T second = m == 0 ? precsqrt(T(3) / T(2)) * eta : precsqrt(T(15) / T(4)) * eta * transverse;
                expect_true(bool(precabs(actual[0].real() - first) <= tolerance));
                expect_true(bool(precabs(actual[1].real() - second) <= tolerance));
                T first_norm = 0, second_norm = 0, overlap = 0;
                for (size_t node = 0; node < quadrature.nodes.size(); ++node) {
                    const auto values = psms_automatic::evaluate_traces(m,
                        quadrature.nodes[node], traces, quadrature.legendre);
                    first_norm += quadrature.weights[node] * complex_abs_sq(values[0]);
                    second_norm += quadrature.weights[node] * complex_abs_sq(values[1]);
                    overlap += quadrature.weights[node] * values[0].real() * values[1].real();
                }
                expect_true(bool(precabs(first_norm - T(1)) <= tolerance));
                expect_true(bool(precabs(second_norm - T(1)) <= tolerance));
                expect_true(bool(precabs(overlap) <= tolerance));
            }
        };
        check(double(0));
#if ACOUSTICTS_HAVE_QUADMATH
        check(acousticts_quad_t(0));
#endif
    }
}

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

    test_that("quad liquid cache respects quadrature, formulation and angular inputs") {
#if ACOUSTICTS_HAVE_QUADMATH
        using T = acousticts_quad_t;
        auto quadrature = [](int count) {
            std::vector<double> nodes, weights;
            gauss_legendre(count, nodes, weights);
            return std::make_pair(std::vector<T>(nodes.begin(), nodes.end()),
                std::vector<T>(weights.begin(), weights.end()));
        };
        const auto base = quadrature(16);
        auto evaluate = [](int m, int n, T c, T ci, T xi, T density,
                const std::vector<T>& nodes, const std::vector<T>& weights,
                T eta, T receive, T phi, bool backscatter, bool coefficients) {
            return psms_full_coupled_field<T>(m, n, c, ci, xi, eta,
                receive, phi, density, nodes, weights, backscatter, coefficients);
        };
        auto invalidate = [&]() {
            evaluate(0, 2, T(0.2), T(0.7), T(1.1), T(0.01), {}, {},
                T(0.4), T(0.4), T(1.3), false, true);
        };
        for (int changed = 0; changed < 21; ++changed) {
            auto rule = changed == 6 ? quadrature(18) : base;
            auto& nodes = rule.first;
            auto& weights = rule.second;
            if (changed == 1 || changed == 3) nodes.front() += T(0.001);
            if (changed == 1) nodes.back() -= T(0.001);
            if (changed == 2 || changed == 4) weights.front() *= T(1.001);
            if (changed == 2) weights.back() *= T(1.001);
            if (changed == 5) {
                std::reverse(nodes.begin(), nodes.end());
                std::reverse(weights.begin(), weights.end());
            }
            const int m = changed == 15 ? 1 : 2, n = changed == 14 ? 9 : 8;
            const T c = changed == 10 ? T(4.2) : T(4);
            const T ci = changed == 11 ? T(4.3) : T(4.1);
            const T xi = changed == 12 ? T(1.2) : T(1.1);
            const T density = changed == 13 ? T(1.1) : T(1.01);
            T eta = changed == 7 ? T(0) : T(0.6);
            if (changed == 17 || changed == 18) {
                const double half_pi = std::acos(-1.0) / 2;
                eta = psms_direction_cosine(T(std::nextafter(half_pi,
                    changed == 17 ? 0.0 : 2.0)));
            }
            if (changed == 19) eta = T(1);
            if (changed == 20) eta = preccos(T(1e-8));
            const T receive = changed == 8 ? T(0.7) : T(0.4);
            const T phi = changed == 9 ? T(2.3) : T(1.3);
            const bool backscatter = changed == 16;
            // Prime a different incident direction before comparing with a cold solve.
            evaluate(2, 8, T(4), T(4.1), T(1.1), T(1.01), base.first, base.second,
                T(0.3), T(0.4), T(1.3), true, false);
            const auto reused = evaluate(m, n, c, ci, xi, density, nodes, weights,
                eta, receive, phi, backscatter, false);
            invalidate();
            const auto rebuilt = evaluate(m, n, c, ci, xi, density, nodes, weights,
                eta, receive, phi, backscatter, false);
            expect_true((complex_abs_sq(reused - rebuilt) <=
                T(1e-54) * complex_abs_sq(rebuilt)));
        }
        // Identical physical inputs must not reuse a different overlap formulation.
        for (bool coefficients : {false, true}) {
            evaluate(2, 8, T(4), T(4.1), T(1.1), T(1.01), base.first, base.second,
                T(0.3), T(0.4), T(1.3), false, !coefficients);
            const auto changed = evaluate(2, 8, T(4), T(4.1), T(1.1), T(1.01),
                base.first, base.second, T(0.6), T(0.4), T(1.3), false, coefficients);
            invalidate();
            const auto rebuilt = evaluate(2, 8, T(4), T(4.1), T(1.1), T(1.01),
                base.first, base.second, T(0.6), T(0.4), T(1.3), false, coefficients);
            expect_true((complex_abs_sq(changed - rebuilt) <=
                T(1e-54) * complex_abs_sq(rebuilt)));
        }
#endif
    }

    test_that("real boundary forcing equals the assembled complex matrix component") {
        auto check = [](auto scalar) {
            using T = decltype(scalar);
            const std::complex<T> value(T(0.7), T(-1.3));
            const std::complex<T> derivative(T(-0.9), T(2.1));
            for (T overlap : {T(-1e-100), T(0), T(0.8), T(1e100)}) {
                for (T interior : {T(-0.3), T(0), T(1.4)}) {
                    const T interior_derivative = T(0.6);
                    const auto entry = overlap *
                        (value * interior_derivative - interior * derivative);
                    const T forcing = overlap *
                        (value.real() * interior_derivative - interior * derivative.real());
                    expect_true((entry.real() == forcing));
                }
            }
        };
        check(double(0));
#if ACOUSTICTS_HAVE_QUADMATH
        check(acousticts_quad_t(0));
#endif
    }

    test_that("angular normalization cache collisions preserve exact modal values") {
        auto check = [](auto scalar) {
            using T = decltype(scalar);
            const std::vector<T> eta = {T(-0.4), T(0.2), T(0.7)};
            const auto original = cprofcn_mblock<T>(
                T(2.3), 0, 1, 1, eta, 0, 1, 2, T(0));
            // The triangular index for (m,l)=(94,180) replaces the (0,0) slot.
            cprofcn_mblock<T>(T(2.3), 94, 1, 87, eta, 0, 1, 2, T(0));
            const auto repeated = cprofcn_mblock<T>(
                T(2.3), 0, 1, 1, eta, 0, 1, 2, T(0));
            expect_true((original.s1c == repeated.s1c));
            expect_true((original.is1e == repeated.is1e));
            expect_true((original.s1dc == repeated.s1dc));
            expect_true((original.is1de == repeated.is1de));
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

    test_that("rejected Neumann evaluation restores a valid second-kind radial method") {
        auto check = [](auto scalar) {
            using T = decltype(scalar);
            const T c = T(20), xi = T(1.5), tolerance =
                is_quad_precision_v<T> ? T(1e-25) : T(1e-12);
            const auto first = Rmn_scalar<T>(0, 0, c, xi, 1);
            const auto second = Rmn_scalar<T>(0, 0, c, xi, 2);
            const auto batch = cprofcn_mblock<T>(c, 0, 1, 36, {T(1)}, 2, 0, 0, xi - T(1));
            const auto coefficients = cprofcn_mblock_coefficients<T>(
                c, 0, 1, 36, {T(1)}, 2, 0, xi - T(1));
            const auto packed = cprofcn_mblock_coefficients<T>(
                c, 0, 1, 36, {T(1)}, 2, 0, xi - T(1), true);
            const T value = second.first.real(), derivative = second.second.real();
            expect_true(bool(is_finite_real(value) && is_finite_real(derivative)));
            expect_true(bool(value != T(0) && derivative != T(0)));
            // This case formerly disabled every radial method after joining-factor cancellation.
            const std::array<std::array<T, 2>, 3> alternate = {{
                {batch.r2c[0], batch.r2dc[0]},
                {coefficients.functions.r2c[0], coefficients.functions.r2dc[0]},
                {psms_packed::scale_decimal(packed.functions.r2c[0], packed.functions.ir2e[0]),
                    psms_packed::scale_decimal(packed.functions.r2dc[0], packed.functions.ir2de[0])}
            }};
            for (const auto& radial : alternate) {
                expect_true(bool(precabs(radial[0] / value - T(1)) <= tolerance));
                expect_true(bool(precabs(radial[1] / derivative - T(1)) <= tolerance));
            }
            const T expected_wronskian = T(1) / (c * (xi - T(1)) * (xi + T(1)));
            const T actual_wronskian = first.first.real() * derivative - first.second.real() * value;
            expect_true(bool(precabs(actual_wronskian / expected_wronskian - T(1)) <= tolerance));
            return std::array<T, 2>{value, derivative};
        };
        const auto ordinary = check(double(0));
#if ACOUSTICTS_HAVE_QUADMATH
        using T = acousticts_quad_t;
        const auto extended = check(T(0));
        for (size_t index = 0; index < ordinary.size(); ++index)
            expect_true(bool(precabs(T(ordinary[index]) / extended[index] - T(1)) <= T(1e-12)));
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

    test_that("high parameter angular normalization retains later coefficient lobes") {
        auto check = [](auto scalar) {
            using T = decltype(scalar);
            const T c = T(std::nextafter(700.0, std::numeric_limits<double>::infinity()));
            const T x1 = T(7) / precsqrt(T(48)) - T(1);
            const auto coefficients = cprofcn_mblock_coefficients<T>(
                c, 2, 1, 379, {T(0.3)}, 1, 0, x1);
            // The former single-term stop gave a squared norm of 1.567379 for degree 376.
            for (int degree : {374, 375, 376, 377, 378, 379, 380}) {
                const int index = degree - 2;
                const T squared_norm = spheroidal_coefficient_overlap<T>(
                    coefficients, coefficients, 0, index, index);
                const T tolerance = T(32) * T(coefficients.lengths[index]) * preceps<T>();
                CATCH_INFO((is_quad_precision_v<T> ? "quad" : "double") << " degree=" << degree
                    << " norm=" << static_cast<double>(squared_norm)
                    << " allowance=" << static_cast<double>(tolerance));
                expect_true(bool(is_finite_real(squared_norm)));
                expect_true(bool(precabs(squared_norm - T(1)) <= tolerance));
                expect_true(bool((coefficients.lengths[index] - 1) % 2 == index % 2));
                bool opposite_parity_is_zero = true;
                for (int offset = 1 - index % 2; offset < coefficients.lengths[index]; offset += 2)
                    opposite_parity_is_zero = opposite_parity_is_zero && coefficients.coefficients[
                        size_t(coefficients.capacity) * index + offset] == T(0);
                expect_true(bool(opposite_parity_is_zero));
            }
        };
        check(double(0));
#if ACOUSTICTS_HAVE_QUADMATH
        check(acousticts_quad_t(0));
#endif
    }

    test_that("complete normalization agrees with public angular evaluation in both parities") {
        auto check = [](auto scalar) {
            using T = decltype(scalar);
            const std::vector<T> arguments = {T(-0.6), T(0.2), T(0.7)};
            const std::vector<int> degrees = {2, 3, 4, 5, 12, 13, 14, 15};
            for (T c : {T(3), T(20)}) {
                auto coefficients = cprofcn_mblock_coefficients<T>(
                    c, 2, 1, 14, arguments, 0, 1, T(0));
                const auto angular = Smn_matrix<T>({2}, degrees, c, arguments, true);
                psms_automatic::LegendreWorkspace<T> workspace;
                std::vector<SpheroidalDirection<T>> directions;
                for (T eta : arguments) directions.push_back(psms_direction_from_cosine(eta));
                psms_automatic::evaluate_coefficient_angular_values(coefficients, 2, directions, workspace);
                // The public evaluator sums the native angular expansion independently of C++ Ferrers recurrence.
                for (size_t row = 0; row < degrees.size(); ++row) {
                    const int index = degrees[row] - 2;
                    const T squared_norm = spheroidal_coefficient_overlap<T>(
                        coefficients, coefficients, 0, index, index);
                    const T tolerance = T(64) * T(coefficients.lengths[index]) * preceps<T>();
                    expect_true(bool(precabs(squared_norm - T(1)) <= tolerance));
                    for (size_t angle = 0; angle < arguments.size(); ++angle) {
                        const T reconstructed = extract_angular_value_from_mblock<T>(
                            coefficients.functions, 0, index, static_cast<int>(angle));
                        expect_true(bool(precabs(reconstructed - angular.value[row][angle]) <=
                            tolerance * (T(1) + precabs(angular.value[row][angle]))));
                    }
                }
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
    test_that("retained block equilibration preserves nonsymmetric incident columns and outgoing rows") {
        auto check = [](auto scalar, int exponent) {
            using T = decltype(scalar);
            using Complex = std::complex<T>;
            const std::vector<Complex> unscaled = {
                {4, 1}, {1, -2}, {-1, 1},
                {2, T(0.5)}, {5, -1}, {1, T(0.5)},
                {-1, T(0.5)}, {2, 1}, {6, 2}
            };
            const std::vector<Complex> expected = {
                {1, 2}, {-3, 1}, {2, -1},
                {-2, -1}, {1, 3}, {4, 1},
                {3, -2}, {-1, T(0.5)}, {-2, 4}
            };
            const T high = pow10_typed<T>(exponent), low = pow10_typed<T>(-exponent);
            const std::array<T, 3> rows{low, high, T(1)}, columns{high, T(1), low};
            std::vector<Complex> matrix(9), rhs(9);
            // Form A=Dr*A0*Dc and B=Dr*A0*T0, whose exact solution is Dc^-1*T0.
            for (int row = 0; row < 3; ++row) {
                for (int col = 0; col < 3; ++col) {
                    matrix[row * 3 + col] = (unscaled[row * 3 + col] * rows[row]) * columns[col];
                    Complex sum(0, 0);
                    for (int inner = 0; inner < 3; ++inner)
                        sum += unscaled[row * 3 + inner] * expected[inner * 3 + col];
                    rhs[row * 3 + col] = sum * rows[row];
                }
            }
            const auto solved = solve_matrix_system_lup_equilibrated(matrix, rhs, 3);
            const T tolerance = T(512) * preceps<T>();
            for (int row = 0; row < 3; ++row) {
                for (int col = 0; col < 3; ++col) {
                    const Complex recovered = solved[row * 3 + col] * columns[row];
                    const Complex error = recovered - expected[row * 3 + col];
                    expect_true(bool(complex_abs_sq(error) <= tolerance * tolerance *
                        complex_abs_sq(expected[row * 3 + col])));
                }
            }
            if constexpr (is_quad_precision_v<T>) {
                const auto retained = solve_fluid_t_blocks_native<T>({rhs}, {matrix});
                for (size_t index = 0; index < solved.size(); ++index)
                    expect_true(bool(retained[0][index] == solved[index]));
            }
            // Equilibration must preserve rejection of genuine rank deficiency and nonfinite input.
            const std::vector<Complex> singular(9, Complex(1, 0));
            expect_error(solve_matrix_system_lup_equilibrated(singular, rhs, 3));
            matrix[0] = Complex(precnan<T>(), T(0));
            expect_error(solve_matrix_system_lup_equilibrated(matrix, rhs, 3));
        };
        check(double(0), 100);
#if ACOUSTICTS_HAVE_QUADMATH
        check(acousticts_quad_t(0), 1000);
#endif
    }

    test_that("enlarged quad retained fluid blocks agree with the direct coefficient equations") {
#if ACOUSTICTS_HAVE_QUADMATH
        using T = acousticts_quad_t;
        psms_automatic::Quadrature<T> quadrature;
        quadrature.prepare(128);
        const T chi = T(5), interior_chi = T(35) / T(8), xi = T(3) / T(2);
        const auto blocks = psms_tmatrix_blocks<T>(0, 48, chi, interior_chi, xi,
            T(1100), T(1000), quadrature.nodes, quadrature.weights, "Amn_fluid_gas");
        const auto incident = compute_smn_direction_matrix<T>(0, 48, chi,
            psms_direction_from_cosine(T(1)));
        const auto receiving = compute_smn_direction_matrix<T>(0, 48, chi,
            psms_direction_from_cosine(T(-1)));
        const auto retained = compute_fbs_from_tmatrix_blocks<T>(0, 48,
            std::vector<T>{T(1)}, incident, receiving, blocks);
        // The independent assembly uses normalized coefficient overlaps instead of quadrature.
        const auto reference = psms_full_coupled_field<T>(0, 48, chi, interior_chi, xi,
            T(1), T(-1), T(0), T(11) / T(10), {}, {}, true, true);
        const T tolerance = T(1e-24);
        expect_true(bool(complex_abs_sq(retained - reference) <=
            tolerance * tolerance * complex_abs_sq(reference)));
#endif
    }

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
    test_that("scalar outgoing and incoming kinds preserve the matrix second-kind guard") {
        auto check = [](auto scalar) {
            using T = decltype(scalar);
            const T tiny = T(1e-300), nan = precnan<T>();
            const T infinity = T(std::numeric_limits<double>::infinity());
            const std::array<std::array<T, 2>, 11> pairs = {{{T(0), T(0)},
                {tiny / T(2), -tiny / T(3)}, {nan, T(2)}, {T(2), nan},
                {infinity, T(2)}, {T(2), infinity}, {T(0), T(2)},
                {T(2), T(0)}, {tiny, T(0)}, {T(0), tiny}, {T(2), T(3)}}};
            ProfcnResult<T> single;
            ProfcnBatchResult<T> block;
            block.lnum = block.narg = block.m_count = 1;
            single.r1c = block.r1c = {T(1)};
            single.r1dc = block.r1dc = {T(-2)};
            single.ir1e = single.ir1de = single.ir2e = single.ir2de = {0};
            block.ir1e = block.ir1de = block.ir2e = block.ir2de = {0};
            for (size_t index = 0; index < pairs.size(); ++index) {
                const auto& pair = pairs[index];
                single.r2c = block.r2c = {pair[0]};
                single.r2dc = block.r2dc = {pair[1]};
                const RmnResult<T> ingredients{T(1), T(-2), pair[0], pair[1], true, true};
                for (int kind : {3, 4}) {
                    const auto combined = combine_higher_radial_kind(ingredients, kind);
                    const auto extracted = extract_radial_from_batch(single, single,
                        0, kind, 0, 0, T(2), T(1.5));
                    const auto batched = extract_radial_from_mblock(block, 0, 0, kind);
                    expect_true(bool(combined.first.real() == T(1)));
                    expect_true(bool(combined.second.real() == T(-2)));
                    if (index < 6) {
                        expect_true(bool(is_na_real(combined.first.imag())));
                        expect_true(bool(is_na_real(combined.second.imag())));
                        expect_true(bool(is_na_real(extracted.val_imag) && is_na_real(extracted.der_imag)));
                        expect_true(bool(is_na_real(batched.val_imag) && is_na_real(batched.der_imag)));
                    } else {
                        const T sign = kind == 3 ? T(1) : T(-1);
                        expect_true(bool(combined.first.imag() == sign * pair[0]));
                        expect_true(bool(combined.second.imag() == sign * pair[1]));
                        expect_true(bool(combined.first.imag() == extracted.val_imag));
                        expect_true(bool(combined.second.imag() == extracted.der_imag));
                        expect_true(bool(combined.first.imag() == batched.val_imag));
                        expect_true(bool(combined.second.imag() == batched.der_imag));
                    }
                }
            }
            if constexpr (is_quad_precision_v<T>) {
                // A finite quad pair beyond double range must remain a valid complex pair.
                const T large = pow10_typed<T>(400);
                const RmnResult<T> ingredients{T(1), T(-2), large, -large, true, true};
                for (int kind : {3, 4}) {
                    const auto combined = combine_higher_radial_kind(ingredients, kind);
                    expect_true(bool(is_finite_real(combined.first.imag())));
                    expect_true(bool(is_finite_real(combined.second.imag())));
                    expect_true(bool(combined.first.imag() == (kind == 3 ? large : -large)));
                    expect_true(bool(combined.second.imag() == (kind == 3 ? -large : large)));
                }
            }
        };
        check(double(0));
#if ACOUSTICTS_HAVE_QUADMATH
        check(acousticts_quad_t(0));
#endif
    }

    test_that("radial extraction preserves finite values in the requested precision") {
        auto check = [](auto scalar, int exponent) {
            using T = decltype(scalar);
            const T scale = pow10_typed<T>(exponent);
            const T expected_value = -T(1.25) * scale;
            const T expected_derivative = T(2.5) * scale;
            ProfcnResult<T> packed;
            packed.r1c = {2};
            packed.r1dc = {-4};
            packed.ir1e = {0};
            packed.ir1de = {0};
            packed.r2c = {-T(1.25)};
            packed.r2dc = {T(2.5)};
            packed.ir2e = {exponent};
            packed.ir2de = {exponent};
            ProfcnBatchResult<T> block;
            block.lnum = 1;
            block.r1c = {2};
            block.r1dc = {-4};
            block.r2c = {expected_value};
            block.r2dc = {expected_derivative};
            // Runtime pow and compiler-folded powers can differ by a few native rounding units.
            const auto same_native_value = [](T actual, T expected) {
                return is_finite_real(actual) && (actual == expected ||
                    (expected != T(0) && precabs((actual - expected) / expected) <= T(8) * preceps<T>()));
            };
            // Keep signs and both wave conventions beyond the double exponent range.
            for (int kind = 2; kind <= 4; ++kind) {
                const auto single = extract_radial_from_batch(
                    packed, packed, 0, kind, 0, 0, T(1), T(1.5));
                const auto batched = extract_radial_from_mblock(block, 0, 0, kind);
                const T sign = kind == 4 ? T(-1) : T(1);
                expect_true(same_native_value(kind == 2 ? single.val_real : single.val_imag, sign * expected_value));
                expect_true(same_native_value(kind == 2 ? single.der_real : single.der_imag, sign * expected_derivative));
                expect_true(same_native_value(batched.val_real, single.val_real));
                expect_true(same_native_value(batched.der_real, single.der_real));
                expect_true(same_native_value(batched.val_imag, single.val_imag));
                expect_true(same_native_value(batched.der_imag, single.der_imag));
            }
            // A nonfinite component still invalidates the complete second-kind pair.
            for (const T invalid : {T(std::numeric_limits<double>::infinity()), precnan<T>()}) {
                for (int component = 0; component < 2; ++component) {
                    packed.ir2e = packed.ir2de = {0};
                    packed.r2c = block.r2c = {component == 0 ? invalid : T(2)};
                    packed.r2dc = block.r2dc = {component == 1 ? invalid : T(3)};
                    for (int kind = 2; kind <= 4; ++kind) {
                        const auto single = extract_radial_from_batch(
                            packed, packed, 0, kind, 0, 0, T(1), T(1.5));
                        const auto batched = extract_radial_from_mblock(block, 0, 0, kind);
                        expect_true(is_na_real(kind == 2 ? single.val_real : single.val_imag));
                        expect_true(is_na_real(kind == 2 ? single.der_real : single.der_imag));
                        expect_true(is_na_real(kind == 2 ? batched.val_real : batched.val_imag));
                        expect_true(is_na_real(kind == 2 ? batched.der_real : batched.der_imag));
                    }
                }
            }
        };
        check(double(0), 300);
#if ACOUSTICTS_HAVE_QUADMATH
        check(acousticts_quad_t(0), 400);
#endif
    }

    test_that("large quad radial modes retain their native Wronskian") {
#if ACOUSTICTS_HAVE_QUADMATH
        using T = acousticts_quad_t;
        const T c = T(120);
        const T xi = T(0.07 / std::sqrt(0.07 * 0.07 - 0.01 * 0.01));
        const int degree = 450;
        const auto block = cprofcn_mblock<T>(c, 0, 1, degree + 1,
            {T(0)}, 2, 1, 1, xi - T(1));
        const T second = block.r2c[degree], derivative = block.r2dc[degree];
        expect_true(bool(finiteq(second) && finiteq(derivative)));
        expect_true(bool(precabs(second) > T(std::numeric_limits<double>::max()) ||
            precabs(derivative) > T(std::numeric_limits<double>::max())));
        // Compare the physical pair before its intentional conversion at the R boundary.
        for (int kind = 2; kind <= 4; ++kind) {
            const auto extracted = extract_radial_from_mblock(block, 0, degree, kind);
            const T sign = kind == 4 ? T(-1) : T(1);
            expect_true(bool((kind == 2 ? extracted.val_real : extracted.val_imag) == sign * second));
            expect_true(bool((kind == 2 ? extracted.der_real : extracted.der_imag) == sign * derivative));
        }
        const auto scalar = Rmn_scalar<T>(0, degree, c, xi, 2);
        expect_true(bool(finiteq(scalar.first.real()) && finiteq(scalar.second.real())));
        const T wronskian = T(1) / (c * (xi - T(1)) * (xi + T(1)));
        const T actual = block.r1c[degree] * derivative - block.r1dc[degree] * second;
        expect_true(block.naccr[degree] >= 8);
        // Allow one reported decimal digit for reconstruction and subtraction rounding.
        expect_true((precabs(actual / wronskian - T(1)) <
            T(10) * pow10_typed<T>(-block.naccr[degree])));
#endif
    }

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

#if ACOUSTICTS_HAVE_EXACT_QUAD_DOT
namespace exact_dot_regression {

using Quad = acousticts_quad_t;
using Accumulator = quad_dot_product::Accumulator;

Quad bits(std::uint64_t low, std::uint64_t high) {
    const std::uint64_t words[2]{low, high};
    Quad value;
    std::memcpy(&value, words, sizeof(value));
    return value;
}

bool identical(Quad first, Quad second) {
    return std::memcmp(&first, &second, sizeof(first)) == 0;
}

struct RestoreRounding {
    int original = std::fegetround();
    ~RestoreRounding() { std::fesetround(original); }
};

Quad evaluate(const std::vector<Quad>& first, const std::vector<Quad>& second,
    Accumulator& accumulator) {
    return quad_dot_product::dot(first.data(), second.data(), first.size(), accumulator);
}

} // namespace exact_dot_regression

context("Exact binary128 coefficient dot products") {
    test_that("full products and halfway values are rounded once to even") {
        using namespace exact_dot_regression;
        RestoreRounding rounding;
        expect_true(std::fesetround(FE_TONEAREST) == 0);
        Accumulator accumulator;
        const Quad one = 1;
        const Quad half_ulp = scalbnq(one, -113);
        const Quad minimum = bits(1, 0);
        const Quad one_ulp = bits(1, UINT64_C(0x3fff000000000000));
        const Quad two_ulps = bits(2, UINT64_C(0x3fff000000000000));
        // The final tiny product distinguishes a tie from either adjacent side.
        expect_true(identical(evaluate({one, one}, {one, half_ulp}, accumulator), one));
        expect_true(identical(evaluate({one_ulp, one}, {one, half_ulp}, accumulator), two_ulps));
        expect_true(identical(evaluate({one, one, minimum},
            {one, half_ulp, minimum}, accumulator), one_ulp));
        expect_true(identical(evaluate({one, one, minimum},
            {one, half_ulp, -minimum}, accumulator), one));
        expect_true(identical(evaluate({-one, -one, -minimum},
            {one, half_ulp, minimum}, accumulator), -one_ulp));
        // Products retain their low bits even when the leading terms cancel.
        expect_true(identical(evaluate({one_ulp, -one, -scalbnq(one, -111)},
            {one_ulp, one, one}, accumulator), scalbnq(one, -224)));
    }

    test_that("subnormal boundaries preserve ties and signed underflow zero") {
        using namespace exact_dot_regression;
        RestoreRounding rounding;
        expect_true(std::fesetround(FE_TONEAREST) == 0);
        Accumulator accumulator;
        const Quad minimum = bits(1, 0);
        const Quad largest_subnormal = bits(UINT64_MAX, UINT64_C(0x0000ffffffffffff));
        const Quad smallest_normal = bits(0, UINT64_C(0x0001000000000000));
        expect_true(identical(evaluate({minimum}, {Quad(0.5)}, accumulator), bits(0, 0)));
        expect_true(identical(evaluate({-minimum}, {Quad(0.5)}, accumulator),
            bits(0, UINT64_C(0x8000000000000000))));
        expect_true(identical(evaluate({minimum}, {Quad(1.5)}, accumulator), bits(2, 0)));
        expect_true(identical(evaluate({minimum}, {Quad(2.5)}, accumulator), bits(2, 0)));
        expect_true(identical(evaluate({largest_subnormal, minimum},
            {Quad(1), Quad(0.5)}, accumulator), smallest_normal));
    }

    test_that("overflowing products can cancel before final rounding") {
        using namespace exact_dot_regression;
        RestoreRounding rounding;
        expect_true(std::fesetround(FE_TONEAREST) == 0);
        Accumulator accumulator;
        const Quad maximum = bits(UINT64_MAX, UINT64_C(0x7ffeffffffffffff));
        const Quad minimum = bits(1, 0);
        const Quad half_ulp = scalbnq(Quad(1), 16270);
        const Quad infinity = bits(0, UINT64_C(0x7fff000000000000));
        expect_true(identical(evaluate({maximum, half_ulp}, {Quad(1), Quad(1)}, accumulator), infinity));
        expect_true(identical(evaluate({maximum, half_ulp, minimum},
            {Quad(1), Quad(1), Quad(-1)}, accumulator), maximum));
        expect_true(identical(evaluate({-maximum, -half_ulp},
            {Quad(1), Quad(1)}, accumulator), -infinity));
        expect_true(identical(evaluate({maximum, -maximum, Quad(1)},
            {maximum, maximum, Quad(1)}, accumulator), Quad(1)));
        expect_true(identical(evaluate({maximum, minimum, -maximum},
            {Quad(1), Quad(1), Quad(1)}, accumulator), minimum));
    }

    test_that("reused accumulators and unequal strides retain the intended terms") {
        using namespace exact_dot_regression;
        RestoreRounding rounding;
        expect_true(std::fesetround(FE_TONEAREST) == 0);
        Accumulator accumulator;
        const Quad maximum = bits(UINT64_MAX, UINT64_C(0x7ffeffffffffffff));
        const Quad minimum = bits(1, 0);
        const Quad infinity = bits(0, UINT64_C(0x7fff000000000000));
        const Quad nan = bits(1, UINT64_C(0x7fff800000000000));
        // Moving between exponent extremes exercises touched-range clearing.
        for (Quad value : {-maximum, minimum, -minimum, Quad(3), maximum, Quad(0)}) {
            expect_true(identical(evaluate({value}, {Quad(1)}, accumulator), value));
        }
        const Quad first[4]{Quad(1), infinity, Quad(3), infinity};
        const Quad second[6]{Quad(2), nan, nan, Quad(4), nan, nan};
        expect_true(identical(quad_dot_product::dot(first, second, 2,
            accumulator, 2, 3), Quad(14)));
        expect_true(identical(quad_dot_product::dot(nullptr, nullptr, 0, accumulator), Quad(0)));
    }

    test_that("forced carry passes preserve extreme cancellation and a sticky product") {
        using namespace exact_dot_regression;
        RestoreRounding rounding;
        expect_true(std::fesetround(FE_TONEAREST) == 0);
        Accumulator accumulator(7);
        const Quad maximum = bits(UINT64_MAX, UINT64_C(0x7ffeffffffffffff));
        const Quad minimum = bits(1, 0);
        const Quad half_ulp = scalbnq(Quad(1), -113);
        const Quad successor = bits(1, UINT64_C(0x3fff000000000000));
        for (int sign : {-1, 1}) {
            std::vector<Quad> first(16), second(16, maximum);
            // Two carry passes occur before the overflowing products finish cancelling.
            for (int index = 0; index < 16; ++index) {
                first[index] = index < 8 ? Quad(sign) * maximum : -Quad(sign) * maximum;
            }
            first.push_back(minimum);
            second.push_back(Quad(1));
            expect_true(identical(evaluate(first, second, accumulator), minimum));
            first.resize(16);
            second.resize(16);
            // The minimum squared product must survive as a sticky bit above the tie.
            first.insert(first.end(), {Quad(1), Quad(1), minimum});
            second.insert(second.end(), {Quad(1), half_ulp, minimum});
            expect_true(identical(evaluate(first, second, accumulator), successor));
        }
    }

    test_that("nonfinite values and directed rounding use the established fallback") {
        using namespace exact_dot_regression;
        RestoreRounding rounding;
        Accumulator accumulator;
        const Quad infinity = bits(0, UINT64_C(0x7fff000000000000));
        const Quad nan = bits(1, UINT64_C(0x7fff800000000000));
        expect_true(std::fesetround(FE_TONEAREST) == 0);
        expect_true(identical(evaluate({infinity}, {Quad(1)}, accumulator), infinity));
        expect_true(isnanq(evaluate({infinity}, {Quad(0)}, accumulator)));
        expect_true(isnanq(evaluate({nan}, {Quad(1)}, accumulator)));
        // Runtime operands prevent folding the sum under nearest-rounding assumptions.
        volatile Quad runtime_half_ulp = scalbnq(Quad(1), -113);
        const Quad half_ulp = runtime_half_ulp;
        // Upward rounding must produce the successor instead of the nearest-even result.
        const int upward = std::fesetround(FE_UPWARD);
        expect_true(upward == 0);
        if (upward == 0) {
            expect_true(identical(evaluate({Quad(1), Quad(1)}, {Quad(1), half_ulp}, accumulator),
                bits(1, UINT64_C(0x3fff000000000000))));
        }
        const int downward = std::fesetround(FE_DOWNWARD);
        expect_true(downward == 0);
        if (downward == 0) {
            expect_true(identical(evaluate({Quad(1), Quad(1)},
                {Quad(1), half_ulp}, accumulator), Quad(1)));
        }
    }
}
#endif

#if ACOUSTICTS_HAVE_EXACT_QUAD_DOT
namespace lup_pivot_regression {

using Quad = acousticts_quad_t;
using Complex = std::complex<Quad>;

LUPFactorization<Quad> reference(const std::vector<Complex>& matrix, int n) {
    LUPFactorization<Quad> result;
    result.lu = matrix;
    result.row_piv.resize(n);
    result.col_piv.resize(n);
    for (int i = 0; i < n; ++i) result.row_piv[i] = result.col_piv[i] = i;
    result.singular = false;
    Quad maximum = 0;
    for (const auto& value : matrix) {
        const Quad norm = complex_abs_sq(value);
        if (norm > maximum) maximum = norm;
    }
    result.tolerance = std::max(Quad(1), precsqrt(maximum)) * Quad(64) * preceps<Quad>();
    for (int k = 0; k < n; ++k) {
        int row = k, col = k;
        Quad pivot_norm = complex_abs_sq(result.lu[k * n + k]);
        Quad half = std::max(precabs(result.lu[k * n + k].real()),
            precabs(result.lu[k * n + k].imag())) * Quad(0.5);
        auto consider = [&](int i, int j) {
            const auto& value = result.lu[i * n + j];
            // Retain the previous floating half-bound and strict norm comparison.
            if (precabs(value.real()) <= half && precabs(value.imag()) <= half) return;
            const Quad norm = complex_abs_sq(value);
            if (norm > pivot_norm) {
                row = i;
                col = j;
                pivot_norm = norm;
                half = std::max(precabs(value.real()), precabs(value.imag())) * Quad(0.5);
            }
        };
        for (int i = k + 1; i < n; ++i) {
            for (int j = k; j < n; ++j) consider(i, j);
        }
        for (int j = k; j < n; ++j) consider(k, j);
        if (pivot_norm <= result.tolerance * result.tolerance) {
            result.singular = true;
            return result;
        }
        if (row != k) {
            for (int j = 0; j < n; ++j) std::swap(result.lu[k * n + j], result.lu[row * n + j]);
            std::swap(result.row_piv[k], result.row_piv[row]);
        }
        if (col != k) {
            for (int i = 0; i < n; ++i) std::swap(result.lu[i * n + k], result.lu[i * n + col]);
            std::swap(result.col_piv[k], result.col_piv[col]);
        }
        const Complex diagonal = result.lu[k * n + k];
        for (int i = k + 1; i < n; ++i) {
            result.lu[i * n + k] /= diagonal;
            const Complex factor = result.lu[i * n + k];
            for (int j = k + 1; j < n; ++j) {
                result.lu[i * n + j] -= factor * result.lu[k * n + j];
            }
        }
    }
    return result;
}

bool identical(Quad first, Quad second) {
    return std::memcmp(&first, &second, sizeof(first)) == 0;
}

bool matches_reference(const std::vector<Complex>& matrix, int n) {
    const auto expected = reference(matrix, n), actual = lup_decompose(matrix, n);
    if (actual.singular != expected.singular || actual.row_piv != expected.row_piv ||
        actual.col_piv != expected.col_piv || !identical(actual.tolerance, expected.tolerance)) return false;
    for (size_t i = 0; i < matrix.size(); ++i) {
        if (!identical(actual.lu[i].real(), expected.lu[i].real()) ||
            !identical(actual.lu[i].imag(), expected.lu[i].imag())) return false;
    }
    return true;
}

} // namespace lup_pivot_regression

context("Binary128 complete-pivot screening") {
    test_that("component dominance preserves pivot ties and rounded norm ties") {
        using namespace lup_pivot_regression;
        const std::vector<Complex> ties{{3, 4}, {4, 3}, {-4, 3}, {0, -5}};
        expect_true(matches_reference(ties, 2));
        const auto tied = lup_decompose(ties, 2);
        expect_true(bool(tied.row_piv[0] == 0 && tied.col_piv[0] == 0));
        const Quad tiny = scalbnq(Quad(1), -60);
        const Complex rounded_tie(Quad(1), tiny);
        expect_true(bool(complex_abs_sq(rounded_tie) == Quad(1)));
        expect_true(matches_reference({{1, 0}, rounded_tie, {0, 0}, {1, 0}}, 2));
        expect_true(matches_reference({{1, 0}, {nextafterq(Quad(1), Quad(2)), 0},
            {0, 0}, {1, 0}}, 2));
        // Swapped, signed components have the same squared norm and stay tied.
        LUPPivotBounds<Quad> bounds(Complex(3, 4));
        expect_true(bounds.rejects(Complex(-4, -3)));
        expect_true(!bounds.rejects(Complex(4, 4)));
    }

    test_that("deterministic dense systems retain every factor and permutation") {
        using namespace lup_pivot_regression;
        std::uint64_t state = UINT64_C(0xf372b569c1d428af);
        auto random = [&]() {
            state = state * UINT64_C(6364136223846793005) + UINT64_C(1442695040888963407);
            const int mantissa = static_cast<int>((state >> 32) % 2049) - 1024;
            const int exponent = static_cast<int>((state >> 48) % 17) - 8;
            return scalbnq(Quad(mantissa) / Quad(1024), exponent);
        };
        for (int n : {1, 2, 3, 6, 10}) {
            for (int trial = 0; trial < 8; ++trial) {
                std::vector<Complex> matrix(n * n);
                for (auto& value : matrix) value = Complex(random(), random());
                expect_true(matches_reference(matrix, n));
            }
        }
    }

    test_that("extreme exponents and signed zero preserve singularity decisions") {
        using namespace lup_pivot_regression;
        for (int exponent : {-16494, -8200, -8000, 0, 8000, 8200, 16380}) {
            const Quad scale = scalbnq(Quad(1), exponent);
            expect_true(matches_reference({{scale, 0}, {0, -scale},
                {scale / Quad(2), scale / Quad(2)}, {-scale, scale}}, 2));
        }
        const Quad negative_zero = -Quad(0);
        expect_true(matches_reference({{0, negative_zero}, {negative_zero, 0},
            {0, 0}, {negative_zero, negative_zero}}, 2));
        expect_true(matches_reference({{1, 1}, {2, 2}, {2, 2}, {4, 4}}, 2));
    }

    test_that("nonfinite components retain the original floating comparison path") {
        using namespace lup_pivot_regression;
        const Quad infinity = static_cast<Quad>(std::numeric_limits<double>::infinity());
        const Quad nan = precnan<Quad>();
        const std::vector<Complex> values{{1, 2}, {0, -Quad(0)}, {infinity, 1},
            {1, -infinity}, {nan, 1}, {1, nan}, {nan, infinity}};
        for (const auto& pivot : values) {
            const Quad half = std::max(precabs(pivot.real()), precabs(pivot.imag())) * Quad(0.5);
            LUPPivotBounds<Quad> bounds(pivot);
            for (const auto& candidate : values) {
                if (finiteq(pivot.real()) && finiteq(pivot.imag()) &&
                    finiteq(candidate.real()) && finiteq(candidate.imag())) continue;
                const bool expected = precabs(candidate.real()) <= half &&
                    precabs(candidate.imag()) <= half;
                expect_true(bounds.rejects(candidate) == expected);
            }
            expect_true(matches_reference({pivot, {1, 0}, {0, 1}, {2, 0}}, 2));
        }
    }
}
#endif


context("Embedded fluid degree support") {
    test_that("adaptive quad gas recovers a degree allowance invalidated by cancellation") {
#if ACOUSTICTS_HAVE_QUADMATH
        using T = acousticts_quad_t;
        // Preserve public binary64 inputs for the 188 kHz, 45 degree gas case.
        const T chi = T(0x1.bb2df130cdc7cp+5), interior_chi = T(0x1.da6d446763a2ap+7);
        const T xi = T(0x1.02a725cde2cb9p+0);
        const T incidence = T(0x1.921fb54442d18p-1), scattering = T(0x1.2d97c7f3321d2p+1);
        const T incident_azimuth = T(0x1.921fb54442d18p+1), scattered_azimuth = T(0x1.921fb54442d18p+2);
        const T density_body = T(0x1.3d70a3d70a3d7p+0), density_water = T(1026.8);
        const std::vector<T> empty;
        const auto optimized = psms_public_field<T>(18, 185, chi, interior_chi, xi,
            incidence, scattering, incident_azimuth, scattered_azimuth, density_body, density_water,
            empty, empty, "Amn_fluid_gas", true, true, false);
        // Fixed support M=48,N=320 agrees with independent M=60 and N=400 controls.
        const std::complex<T> reference(T(0.125805993390489), T(-0.000718870932272577));
        const T reference_norm = psms_direct_refinement::modulus(reference);
        expect_true(psms_direct_refinement::finite(optimized));
        const T difference_db = T(20) * precabs(log10q(
            psms_direct_refinement::modulus(optimized) / reference_norm));
        expect_true(bool(difference_db <= T(1e-8)));
        const T complex_difference = psms_direct_refinement::modulus(optimized - reference);
        expect_true(bool(complex_difference <= psms_modal_abs_tol<T>() +
            psms_modal_rel_tol<T>() * reference_norm));
#endif
    }

    test_that("adaptive quad liquid expands an undersized seed to match enlarged support") {
#if ACOUSTICTS_HAVE_QUADMATH
        using T = acousticts_quad_t;
        const T chi = T(0.4), interior_chi = T(0.7), xi = T(1.5), density_ratio = T(1.1);
        const double pi = std::acos(-1.0), angle = 0.61;
        const std::vector<T> empty;
        auto evaluate = [&](int maximum_order, int maximum_degree, bool adaptive) {
            return psms_public_field<T>(maximum_order, maximum_degree, chi, interior_chi, xi,
                T(angle), T(pi - angle), T(0), T(pi), density_ratio, T(1),
                empty, empty, "Amn_fluid_coefficients", true, adaptive, false);
        };
        const auto undersized = evaluate(0, 0, false);
        const auto optimized = evaluate(0, 0, true);
        // Independently retain every order and degree in the enlarged coefficient solve.
        const auto reference = evaluate(12, 24, false);
        expect_true(psms_direct_refinement::finite(optimized));
        expect_true(psms_direct_refinement::finite(reference));
        expect_true(bool(optimized != undersized));
        const T reference_norm = psms_direct_refinement::modulus(reference);
        expect_true(bool(reference_norm > T(0)));
        if (reference_norm > T(0)) {
            const T difference_db = T(20) * precabs(log10q(
                psms_direct_refinement::modulus(optimized) / reference_norm));
            expect_true(bool(difference_db <= T(1e-8)));
            const T complex_difference = psms_direct_refinement::modulus(optimized - reference);
            expect_true(bool(complex_difference <= psms_modal_abs_tol<T>() +
                psms_modal_rel_tol<T>() * reference_norm));
        }
#endif
    }

    test_that("adaptive quad liquid recovers a degree allowance invalidated by cancellation") {
#if ACOUSTICTS_HAVE_QUADMATH
        using T = acousticts_quad_t;
        // Preserve the public double inputs for the 128 kHz, 60 degree liquid null.
        const T chi = T(0x1.2dbd3cb9bd207p+5), interior_chi = T(0x1.2d20b0d584431p+5);
        const T xi = T(0x1.02a725cde2cb9p+0);
        const T incidence = T(0x1.0c152382d7365p+0), scattering = T(0x1.0c152382d7366p+1);
        const T incident_azimuth = T(0x1.921fb54442d18p+1), scattered_azimuth = T(0x1.921fb54442d18p+2);
        const T density_body = T(0x1.013999999999ap+10), density_water = T(1026.8);
        const std::vector<T> empty;
        auto evaluate = [&](int maximum_order, int maximum_degree, bool adaptive) {
            return psms_public_field<T>(maximum_order, maximum_degree, chi, interior_chi, xi,
                incidence, scattering, incident_azimuth, scattered_azimuth, density_body, density_water,
                empty, empty, "Amn_fluid_coefficients", true, adaptive, false);
        };
        const auto optimized = evaluate(14, 35, true);
        // Separate order and degree enlargements check the fixed-support reference.
        const auto reference = evaluate(32, 80, false);
        const auto order_control = evaluate(40, 80, false);
        const auto degree_control = evaluate(32, 100, false);
        auto check_agreement = [&](const std::complex<T>& actual, const std::complex<T>& expected) {
            expect_true(psms_direct_refinement::finite(actual));
            expect_true(psms_direct_refinement::finite(expected));
            const T reference_norm = psms_direct_refinement::modulus(expected);
            expect_true(bool(reference_norm > T(0)));
            if (reference_norm > T(0)) {
                const T difference_db = T(20) * precabs(log10q(
                    psms_direct_refinement::modulus(actual) / reference_norm));
                expect_true(bool(difference_db <= T(1e-8)));
                const T complex_difference = psms_direct_refinement::modulus(actual - expected);
                expect_true(bool(complex_difference <= psms_modal_abs_tol<T>() +
                    psms_modal_rel_tol<T>() * reference_norm));
            }
        };
        check_agreement(reference, order_control);
        check_agreement(reference, degree_control);
        check_agreement(optimized, reference);
#endif
    }

    test_that("both embedded deletions preserve distinct reduced forcing and complete pivots") {
        auto check = [](auto scalar) {
            using T = decltype(scalar);
            using Complex = std::complex<T>;
            const int dimension = 4;
            const std::vector<Complex> matrix{
                {0, 0}, {2, 0}, {0, 1}, {1, 0},
                {1, 0}, {0, 0}, {3, 0}, {0, 2},
                {2, 0}, {0, 1}, {4, 0}, {1, 0},
                {1, 0}, {3, 0}, {0, 2}, {5, 0}
            };
            const std::array<std::vector<Complex>, 2> expected{{
                {{1, 1}, {-2, T(0.5)}, {T(0.75), -1}},
                {{0, -1}, {2, 1}}
            }};
            std::array<std::vector<Complex>, 2> rhs;
            for (int comparison = 0; comparison < 2; ++comparison) {
                rhs[comparison].resize(dimension);
                const int retained = dimension - comparison - 1;
                for (int row = 0; row < retained; ++row)
                    for (int column = 0; column < retained; ++column)
                        rhs[comparison][row] += matrix[row * dimension + column] * expected[comparison][column];
            }
            const auto factors = lup_decompose<T>(matrix, dimension);
            expect_false(factors.singular);
            expect_true(bool(factors.row_piv[0] != 0 && factors.col_piv[0] != 0));
            const auto result = psms_embedded_support::reduced_solutions<T>(matrix, factors, dimension, rhs);
            for (int comparison = 0; comparison < 2; ++comparison) {
                CATCH_INFO("embedded deletion=" << comparison + 1);
                expect_true(result.valid[comparison]);
                if (!result.valid[comparison]) continue;
                const int retained = dimension - comparison - 1;
                std::vector<Complex> lower_matrix(retained * retained), lower_rhs(retained);
                for (int row = 0; row < retained; ++row) {
                    lower_rhs[row] = rhs[comparison][row];
                    for (int column = 0; column < retained; ++column)
                        lower_matrix[row * retained + column] = matrix[row * dimension + column];
                }
                const auto direct = solve_linear_system_lup_refined<T>(lower_matrix, lower_rhs, retained, 2);
                for (int row = 0; row < retained; ++row) {
                    const T allowance = T(2048) * preceps<T>() *
                        psms_direct_refinement::abs1(expected[comparison][row]);
                    expect_true(bool(psms_direct_refinement::abs1(
                        result.solutions[comparison][row] - expected[comparison][row]) <= allowance));
                    expect_true(bool(psms_direct_refinement::abs1(
                        result.solutions[comparison][row] - direct[row]) <= allowance));
                }
                expect_true(bool(result.backward_error[comparison] <= T(64) * preceps<T>()));
            }
            // Deleted equations must not contribute to either supplied reduced forcing.
            rhs[0][dimension - 1] = Complex(1, 0);
            const auto invalid = psms_embedded_support::reduced_solutions<T>(matrix, factors, dimension, rhs);
            expect_false(invalid.valid[0]);
        };
        check(double(0));
#if ACOUSTICTS_HAVE_QUADMATH
        check(acousticts_quad_t(0));
#endif
    }

    test_that("singular reduced systems are unavailable even when the full system is invertible") {
        auto check = [](auto scalar) {
            using T = decltype(scalar);
            using Complex = std::complex<T>;
            const std::vector<Complex> permutation{
                {0, 0}, {1, 0}, {0, 0}, {0, 0}, {0, 0}, {1, 0}, {1, 0}, {0, 0}, {0, 0}
            };
            std::array<std::vector<Complex>, 2> rhs{{
                std::vector<Complex>(3, Complex(0, 0)), std::vector<Complex>(3, Complex(0, 0))
            }};
            rhs[0][0] = rhs[1][0] = Complex(1, 0);
            const auto factors = lup_decompose<T>(permutation, 3);
            expect_false(factors.singular);
            const auto reduced = psms_embedded_support::reduced_solutions<T>(permutation, factors, 3, rhs);
            expect_false(reduced.valid[0]);
            expect_false(reduced.valid[1]);
            const std::vector<Complex> singular(9, Complex(1, 0));
            const auto singular_factors = lup_decompose<T>(singular, 3);
            expect_true(singular_factors.singular);
            const auto invalid = psms_embedded_support::reduced_solutions<T>(singular, singular_factors, 3, rhs);
            expect_false(invalid.valid[0]);
            expect_false(invalid.valid[1]);
        };
        check(double(0));
#if ACOUSTICTS_HAVE_QUADMATH
        check(acousticts_quad_t(0));
#endif
    }

    test_that("selective refinement recovers cancellation-prone reduced coefficients") {
        auto check = [](auto scalar) {
            using T = decltype(scalar);
            using Complex = std::complex<T>;
            const T diagonal = precsqrt(preceps<T>()), small = precsqrt(diagonal) / T(3);
            const std::vector<Complex> matrix{
                {1, 0}, {0, 0}, {1, 0}, {0, 0}, {diagonal, 0}, {1, 0}, {1, 0}, {1, 0}, {1, 0}
            };
            std::array<std::vector<Complex>, 2> rhs{{
                {{1, 0}, {diagonal * small, 0}, {0, 0}}, {{1, 0}, {0, 0}, {0, 0}}
            }};
            const auto factors = lup_decompose<T>(matrix, 3);
            const auto result = psms_embedded_support::reduced_solutions<T>(matrix, factors, 3, rhs);
            expect_true(result.valid[0]);
            expect_true(result.valid[1]);
            if (result.valid[0]) {
                expect_true(bool(psms_direct_refinement::abs1(result.solutions[0][0] - Complex(1, 0))
                    <= T(512) * preceps<T>()));
                expect_true(bool(psms_direct_refinement::abs1(result.solutions[0][1] - Complex(small, 0))
                    <= T(512) * preceps<T>() * small));
                expect_true(bool(result.backward_error[0] <= T(64) * preceps<T>()));
                expect_true(result.correction_count[0] <= 2);
            }
        };
        check(double(0));
#if ACOUSTICTS_HAVE_QUADMATH
        check(acousticts_quad_t(0));
#endif
    }

    test_that("embedded liquid fields match independently evaluated lower degree supports") {
#if ACOUSTICTS_HAVE_QUADMATH
        using T = acousticts_quad_t;
        using namespace psms_automatic;
        const int m = 1, degree = 12;
        const T chi = T(2), interior_chi = T(1.98), xi = T(1.2), density_ratio = T(1.01);
        const auto incidence = psms_direction_from_angle(T(0.37));
        const std::vector<T> arguments{incidence.cosine};
        const std::vector<SpheroidalDirection<T>> directions{incidence};
        LegendreWorkspace<T> workspace;
        auto evaluate = [&](int maximum_degree, bool embedded) {
            auto exterior = cprofcn_mblock_coefficients<T>(chi, m, 1,
                maximum_degree - m + 1, arguments, 2, 0, xi - T(1), true);
            evaluate_coefficient_angular_values(exterior, m, directions, workspace);
            return solve_order<T>(m, maximum_degree, maximum_degree, interior_chi, xi,
                incidence.cosine, precpi<T>(), density_ratio, true, false, exterior,
                nullptr, nullptr, false, false, false, false, embedded);
        };
        const auto embedded = evaluate(degree, true);
        expect_false(embedded.singular);
        for (int comparison = 0; comparison < 2; ++comparison) {
            const auto reference = evaluate(degree - 2 * (comparison + 1), false);
            expect_false(reference.singular);
            expect_true(embedded.embedded_valid[comparison]);
            const auto difference = embedded.embedded_fields[comparison] - reference.field;
            const T relative_error = psms_direct_refinement::modulus(difference) /
                psms_direct_refinement::modulus(reference.field);
            expect_true(bool(relative_error < T(1e-24)));
        }
#endif
    }
}

#if ACOUSTICTS_HAVE_QUADMATH
context("Complete-coefficient Robin operator") {
    test_that("implicit projections preserve full dense actions and residuals") {
        using T = acousticts_quad_t;
        using Complex = std::complex<T>;
        int exact_comparisons = 0, rejected_inputs = 0;
        auto require = [](bool condition, const char* message) {
            if (!condition) Rcpp::stop(message);
        };
        auto equal_vectors = [&](const std::vector<Complex>& actual, const std::vector<Complex>& expected) {
            require(actual == expected, "Implicit operator dyadic identity failed.");
            exact_comparisons += int(actual.size());
        };
        auto batch = [](int capacity, std::initializer_list<int> lengths) {
            ProfcnCoefficientBatchResult<T> result;
            result.capacity = capacity;
            result.functions.lnum = int(lengths.size());
            result.functions.m_count = 1;
            result.lengths.assign(lengths);
            result.coefficients.assign(std::size_t(capacity) * lengths.size(), T(0));
            for (int degree = 0; degree < result.functions.lnum; ++degree) {
                const int length = result.lengths[degree];
                for (int position = degree % 2; position < length; position += 2) {
                    const int numerator = (3 * degree + 2 * position) % 7 + 1;
                    const T sign = ((degree + position / 2) % 2) ? T(-1) : T(1);
                    result.coefficients[std::size_t(capacity) * degree + position] = sign * T(numerator) / T(8);
                }
                // Nonzero padding detects any accidental use beyond the native column length.
                for (int position = length; position < capacity; ++position)
                    result.coefficients[std::size_t(capacity) * degree + position] = T(13) / T(2);
            }
            return result;
        };
        const auto exterior = batch(7, {5, 4, 7, 6, 3});
        const auto interior = batch(9, {9, 2, 5, 8, 3, 4, 7});
        for (int parity = 0; parity < 2; ++parity) {
            const int columns = (exterior.functions.lnum - parity + 1) / 2;
            const int rows = (interior.functions.lnum - parity + 1) / 2;
            std::vector<T> overlap(std::size_t(rows) * columns, T(0));
            std::vector<Complex> response(rows), diagonal(columns), input(columns), rhs(columns);
            for (int row = 0; row < rows; ++row)
                response[row] = Complex(T(row + 1) / T(8), T(1 - row) / T(16));
            for (int column = 0; column < columns; ++column) {
                diagonal[column] = Complex(T(2) + T(column) / T(4), -T(column + 1) / T(8));
                input[column] = Complex(T(column + 1) / T(8), T(2 - column) / T(16));
                rhs[column] = Complex(-T(column + 1) / T(4), T(column + 2) / T(8));
            }
    
            // Independently form C from only the shared retained Legendre positions.
            for (int row = 0; row < rows; ++row) {
                const int internal_degree = parity + 2 * row;
                for (int column = 0; column < columns; ++column) {
                    const int external_degree = parity + 2 * column;
                    const int length = std::min(interior.lengths[internal_degree], exterior.lengths[external_degree]);
                    for (int position = parity; position < length; position += 2)
                        overlap[std::size_t(row) * columns + column] +=
                            interior.coefficients[std::size_t(interior.capacity) * internal_degree + position] *
                            exterior.coefficients[std::size_t(exterior.capacity) * external_degree + position];
                }
            }
            std::vector<Complex> expected_projection(rows), weighted(rows), expected_map(columns);
            for (int row = 0; row < rows; ++row) {
                for (int column = 0; column < columns; ++column)
                    expected_projection[row] += overlap[std::size_t(row) * columns + column] * input[column];
                weighted[row] = response[row] * expected_projection[row];
            }
            for (int column = 0; column < columns; ++column)
                for (int row = 0; row < rows; ++row)
                    expected_map[column] += overlap[std::size_t(row) * columns + column] * weighted[row];
            equal_vectors(psms_implicit_overlap::overlap_projection(exterior, interior, parity, input), expected_projection);
            equal_vectors(psms_implicit_overlap::apply_map(exterior, interior, parity, response, input), expected_map);
    
            // Compare both affine signs and the actual dense infinity norm exactly.
            std::vector<Complex> expected_action(columns), expected_residual(columns), actual;
            T dense_norm = 0;
            for (int column = 0; column < columns; ++column) {
                expected_action[column] = diagonal[column] * input[column] - expected_map[column];
                expected_residual[column] = rhs[column] - expected_action[column];
                T row_norm = 0;
                for (int other = 0; other < columns; ++other) {
                    Complex entry = column == other ? diagonal[column] : Complex(0, 0);
                    for (int row = 0; row < rows; ++row)
                        entry -= overlap[std::size_t(row) * columns + column] * response[row] *
                            overlap[std::size_t(row) * columns + other];
                    row_norm += precabs(entry.real()) + precabs(entry.imag());
                }
                dense_norm = std::max(dense_norm, row_norm);
            }
            psms_implicit_overlap::Operator<T> operation(exterior, interior, parity, response, diagonal);
            require(operation.initialize(), "Implicit operator initialization failed.");
            require(operation.diagonal == diagonal, "Implicit preconditioner changed the exterior diagonal.");
            require(operation.norm_bound >= dense_norm, "Implicit norm bound is below the explicit matrix norm.");
            require(operation.apply(input, actual), "Implicit operator action failed.");
            equal_vectors(actual, expected_action);
            require(operation.action(input, actual, &rhs), "Implicit affine residual failed.");
            equal_vectors(actual, expected_residual);
    
            // Recover an independent manufactured solution with the implicit GMRES solve.
            const auto solution = psms_robin_solver::solve_implicit(operation, expected_action);
            require(solution.success, "Implicit GMRES rejected a manufactured system.");
            for (int column = 0; column < columns; ++column)
                require(psms_automatic::combine_norms(
                    precabs(solution.solution[column].real() - input[column].real()),
                    precabs(solution.solution[column].imag() - input[column].imag())) < T(1e-28),
                    "Implicit GMRES disagrees with its manufactured solution.");
    
            // Reject malformed dimensions and nonfinite inputs before entering native arithmetic.
            auto invalid = exterior;
            invalid.lengths[parity] = invalid.capacity + 1;
            psms_implicit_overlap::Operator<T> bad_length(invalid, interior, parity, response, diagonal);
            require(!bad_length.initialize(), "Implicit operator accepted an invalid coefficient length.");
            ++rejected_inputs;
            invalid = exterior;
            invalid.coefficients.pop_back();
            psms_implicit_overlap::Operator<T> bad_rectangle(invalid, interior, parity, response, diagonal);
            require(!bad_rectangle.initialize(), "Implicit operator accepted an invalid coefficient rectangle.");
            ++rejected_inputs;
            auto bad_input = input;
            bad_input.pop_back();
            require(!operation.apply(bad_input, actual), "Implicit operator accepted an invalid input length.");
            ++rejected_inputs;
            bad_input = input;
            bad_input[0] = Complex(precnan<T>(), T(0));
            require(!operation.apply(bad_input, actual), "Implicit operator accepted a nonfinite input.");
            ++rejected_inputs;
            auto bad_rhs = rhs;
            bad_rhs.pop_back();
            require(!operation.action(input, actual, &bad_rhs), "Implicit operator accepted an invalid RHS length.");
            ++rejected_inputs;
        }
        expect_true(exact_comparisons > 0);
        expect_true(rejected_inputs == 10);
    }
}
#endif
