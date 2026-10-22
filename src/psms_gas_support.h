#pragma once

#if ACOUSTICTS_HAVE_QUADMATH
namespace psms_gas_support {
using T = acousticts_quad_t;
using Complex = std::complex<T>;
using namespace psms_automatic;

// Embedded principal-system comparisons reuse the checked direct factorization.
inline bool direct_field(T chi_sw, T chi_body, T xi, int degree_seed,
        T density_ratio, const SpheroidalDirection<T>& incidence,
        T phi_difference, Complex& field) {
    // Native iopang=0 skips eta entirely, while coefficient angular values are rebuilt below.
    const std::vector<T> arguments{T(0)};
    const std::vector<SpheroidalDirection<T>> directions{incidence};
    LegendreWorkspace<T> workspace;
    auto norm = [](const Complex& value) { return combine_norms(precabs(value.real()), precabs(value.imag())); };
    // Require both the TS target and the native complex-amplitude tolerance.
    const T ts_amplitude_allowance = -expm1q(-logq(T(10)) * T(1e-8) / T(20));
    const T amplitude_allowance = std::min(ts_amplitude_allowance, psms_modal_rel_tol<T>());
    struct DegreeContribution {
        int degree = 0, attempts = 0;
        Complex field{0, 0};
        T estimate = 0;
        bool arithmetic_floor = false;
    };
    // Retain scalar summaries so cancellation can tighten previously accepted orders.
    std::vector<DegreeContribution> contributions;
    auto resolve_degree = [&](int m, int degree, int previous_attempts,
            const Complex& other_orders, bool cancellation_recovery,
            DegreeContribution& accepted) {
        SpheroidalCoefficientOverlapCache<T> overlap_cache;
        OrderSolution<T> order;
        T estimate = 0, allowance = 0;
        for (int refinement = previous_attempts; refinement < 12; ++refinement) {
            Rcpp::checkUserInterrupt();
            auto exterior = cprofcn_mblock_coefficients<T>(T(chi_sw), m, 1,
                degree - m + 1, arguments, 2, 0, T(xi) - T(1), true);
            evaluate_coefficient_angular_values(exterior, m, directions, workspace);
            order = solve_order<T>(m, degree, degree, T(chi_body), T(xi), incidence.cosine,
                phi_difference, density_ratio, true, false, exterior, nullptr,
                &overlap_cache, false, false, true, false, true);
            if (order.singular) return false;
            if (!order.embedded_valid[0] || !order.embedded_valid[1])
                return false;
            const T scale = norm(other_orders + order.field);
            allowance = amplitude_allowance * scale / (T(2) * T(m + 1) * T(m + 2));
            const T last = norm(order.field - order.embedded_fields[0]);
            const T preceding = norm(order.embedded_fields[0] - order.embedded_fields[1]);
            // Account for cancellation between modal terms when identifying quad roundoff.
            const T arithmetic_floor = T(64) * preceps<T>() * order.embedded_roundoff_scale;
            bool degree_pass = false;
            accepted.arithmetic_floor = false;
            // Geometric continuation is conditional on contraction, not a universal bound.
            if (order.embedded_valid[0] && order.embedded_valid[1]) {
                if (last <= arithmetic_floor && preceding <= arithmetic_floor) {
                    estimate = arithmetic_floor;
                    degree_pass = estimate <= allowance;
                    accepted.arithmetic_floor = true;
                } else if (preceding > T(0) && last < preceding) {
                    const T ratio = last / preceding;
                    estimate = last / (T(1) - ratio);
                    degree_pass = estimate <= allowance && preceding <= allowance;
                }
            }
            // A floor-limited recovery still has to pass the complete error sum below.
            if (degree_pass || (cancellation_recovery && accepted.arithmetic_floor)) {
                accepted.degree = degree;
                accepted.attempts = refinement + 1;
                accepted.field = order.field;
                accepted.estimate = estimate;
                return true;
            }
            degree = checked_profcn_integer((5LL * (degree + 1LL) + 3LL) / 4LL - 1LL,
                "optimized degree refinement");
        }
        Rcpp::stop("Optimized degree tail did not converge within twelve refinements.");
        return false;
    };
    field = Complex(0, 0);
    Complex field_correction(0, 0);
    T previous_order = -1, penultimate_order = -1, total_degree_estimate = 0;
    int small_orders = 0;
    // Independent order blocks permit refinement without discarding accepted work.
    for (int m = 0; !(chi_sw == chi_body && density_ratio == T(1)); ++m) {
        Rcpp::checkUserInterrupt();
        const int degree = std::max(degree_seed, m + 5);
        DegreeContribution order;
        if (!resolve_degree(m, degree, 0, field, false, order)) return false;
        contributions.push_back(order);
        total_degree_estimate += order.estimate;
        add_field(field, field_correction, order.field);
        T magnitude = norm(order.field);
        const T order_allowance = amplitude_allowance * norm(field) / T(2);
        small_orders = m >= 2 && magnitude <= order_allowance ? small_orders + 1 : 0;
        bool order_pass = incidence.exact_axis;
        if (small_orders >= 2 && previous_order > T(0) && penultimate_order > T(0)) {
            const T ratio = std::max(magnitude / previous_order, previous_order / penultimate_order);
            order_pass = ratio < T(1) && magnitude * ratio / (T(1) - ratio) <= order_allowance;
        }
        if (order_pass) {
            bool recovered = false;
            // Refine the largest remaining error against the complete cancelled field.
            while (total_degree_estimate > amplitude_allowance * norm(field) / T(2)) {
                Rcpp::checkUserInterrupt();
                T fixed_estimate = 0, largest_estimate = -1;
                int selected = -1;
                for (int index = 0; index <= m; ++index) {
                    const auto& retained = contributions[index];
                    if (retained.arithmetic_floor || retained.attempts >= 12) {
                        fixed_estimate += retained.estimate;
                    } else if (retained.estimate > largest_estimate) {
                        largest_estimate = retained.estimate;
                        selected = index;
                    }
                }
                if (selected < 0 || fixed_estimate > amplitude_allowance * norm(field) / T(2))
                    Rcpp::stop("Adaptive gas cancellation exhausts the degree arithmetic or refinement allowance.");
                auto& retained = contributions[selected];
                // Rebuild other contributions instead of subtracting nearly cancelling terms.
                Complex other_orders(0, 0), other_correction(0, 0);
                for (int index = 0; index <= m; ++index)
                    if (index != selected)
                        add_field(other_orders, other_correction, contributions[index].field);
                const int enlarged_degree = checked_profcn_integer(
                    (5LL * (retained.degree + 1LL) + 3LL) / 4LL - 1LL,
                    "adaptive gas cancellation refinement");
                if (!resolve_degree(selected, enlarged_degree, retained.attempts,
                        other_orders, true, retained)) return false;
                field = field_correction = Complex(0, 0);
                total_degree_estimate = 0;
                for (const auto& contribution : contributions) {
                    add_field(field, field_correction, contribution.field);
                    total_degree_estimate += contribution.estimate;
                }
                recovered = true;
            }
            if (recovered) {
                // Renew the azimuthal-tail check after replacing any degree contribution.
                const T updated_allowance = amplitude_allowance * norm(field) / T(2);
                small_orders = 0;
                for (int index = m; index >= 2 &&
                        norm(contributions[index].field) <= updated_allowance; --index)
                    ++small_orders;
                magnitude = norm(contributions[m].field);
                previous_order = m > 0 ? norm(contributions[m - 1].field) : T(-1);
                penultimate_order = m > 1 ? norm(contributions[m - 2].field) : T(-1);
                order_pass = incidence.exact_axis;
                if (small_orders >= 2 && previous_order > T(0) && penultimate_order > T(0)) {
                    const T ratio = std::max(magnitude / previous_order, previous_order / penultimate_order);
                    order_pass = ratio < T(1) && magnitude * ratio / (T(1) - ratio) <= updated_allowance;
                }
            }
            if (order_pass) break;
        }
        penultimate_order = previous_order;
        previous_order = magnitude;
    }
    return true;
}

struct DirectContribution {
    int degree = 0, attempts = 0;
    Complex field{0, 0};
    T estimate = 0;
    bool arithmetic_floor = false;
};

struct DirectState {
    SpheroidalDirection<T> incidence;
    int degree_seed;
    std::vector<DirectContribution> contributions;
    Complex field{0, 0}, correction{0, 0};
    T previous_order = -1, penultimate_order = -1, total_degree_estimate = 0;
    int small_orders = 0;
    bool finished = false, rejected = false;
    DirectState(const SpheroidalDirection<T>& direction, int seed)
        : incidence(direction), degree_seed(seed) {}
};

inline T direct_norm(const Complex& value) {
    return combine_norms(precabs(value.real()), precabs(value.imag()));
}

inline bool resolve_direct_degree(DirectState& state,
        psms_angle_prepared::Workspace<T>& prepared, LegendreWorkspace<T>& angular_workspace,
        int m, int degree, int previous_attempts, const Complex& other_orders,
        bool cancellation_recovery, T chi_sw, T chi_body, T xi, T density_ratio,
        T phi_difference, T amplitude_allowance, DirectContribution& accepted) {
    SpheroidalCoefficientOverlapCache<T> overlap_cache;
    OrderSolution<T> order;
    for (int refinement = previous_attempts; refinement < 12; ++refinement) {
        Rcpp::checkUserInterrupt();
        // Refresh this angle's angular values on the shared native degree block.
        auto& exterior = prepared.get(true, chi_sw, m, degree, xi - T(1)).values;
        evaluate_coefficient_angular_values(exterior, m,
            std::vector<SpheroidalDirection<T>>{state.incidence}, angular_workspace);
        const auto& interior = prepared.get(false, chi_body, m, degree, xi - T(1));
        order = solve_order<T>(m, degree, degree, chi_body, xi, state.incidence.cosine,
            phi_difference, density_ratio, true, false, exterior, nullptr,
            &overlap_cache, false, false, true, false, true, &interior);
        if (order.singular || !order.embedded_valid[0] || !order.embedded_valid[1])
            return false;
        const T scale = direct_norm(other_orders + order.field);
        const T allowance = amplitude_allowance * scale / (T(2) * T(m + 1) * T(m + 2));
        const T last = direct_norm(order.field - order.embedded_fields[0]);
        const T preceding = direct_norm(order.embedded_fields[0] - order.embedded_fields[1]);
        const T arithmetic_floor = T(64) * preceps<T>() * order.embedded_roundoff_scale;
        bool degree_pass = false;
        accepted.arithmetic_floor = false;
        if (last <= arithmetic_floor && preceding <= arithmetic_floor) {
            accepted.estimate = arithmetic_floor;
            degree_pass = accepted.estimate <= allowance;
            accepted.arithmetic_floor = true;
        } else if (preceding > T(0) && last < preceding) {
            const T ratio = last / preceding;
            accepted.estimate = last / (T(1) - ratio);
            degree_pass = accepted.estimate <= allowance && preceding <= allowance;
        }
        if (degree_pass || (cancellation_recovery && accepted.arithmetic_floor)) {
            accepted.degree = degree;
            accepted.attempts = refinement + 1;
            accepted.field = order.field;
            return true;
        }
        degree = checked_profcn_integer((5LL * (degree + 1LL) + 3LL) / 4LL - 1LL,
            "optimized degree refinement");
    }
    Rcpp::stop("Optimized degree tail did not converge within twelve refinements.");
    return false;
}

inline void advance_direct(DirectState& state,
        psms_angle_prepared::Workspace<T>& prepared, LegendreWorkspace<T>& angular_workspace,
        int m, T chi_sw, T chi_body, T xi, T density_ratio,
        T phi_difference, T amplitude_allowance) {
    Rcpp::checkUserInterrupt();
    const int degree = std::max(state.degree_seed, m + 5);
    DirectContribution order;
    if (!resolve_direct_degree(state, prepared, angular_workspace, m, degree, 0,
            state.field, false, chi_sw, chi_body, xi, density_ratio,
            phi_difference, amplitude_allowance, order)) {
        state.rejected = true;
        return;
    }
    state.contributions.push_back(order);
    state.total_degree_estimate += order.estimate;
    add_field(state.field, state.correction, order.field);
    T magnitude = direct_norm(order.field);
    const T order_allowance = amplitude_allowance * direct_norm(state.field) / T(2);
    state.small_orders = m >= 2 && magnitude <= order_allowance ? state.small_orders + 1 : 0;
    bool order_pass = state.incidence.exact_axis;
    if (state.small_orders >= 2 && state.previous_order > T(0) && state.penultimate_order > T(0)) {
        const T ratio = std::max(magnitude / state.previous_order,
            state.previous_order / state.penultimate_order);
        order_pass = ratio < T(1) && magnitude * ratio / (T(1) - ratio) <= order_allowance;
    }
    if (order_pass) {
        bool recovered = false;
        // Tighten only this angle's previously accepted degree contributions.
        while (state.total_degree_estimate > amplitude_allowance * direct_norm(state.field) / T(2)) {
            Rcpp::checkUserInterrupt();
            T fixed_estimate = 0, largest_estimate = -1;
            int selected = -1;
            for (int index = 0; index <= m; ++index) {
                const auto& retained = state.contributions[index];
                if (retained.arithmetic_floor || retained.attempts >= 12) {
                    fixed_estimate += retained.estimate;
                } else if (retained.estimate > largest_estimate) {
                    largest_estimate = retained.estimate;
                    selected = index;
                }
            }
            if (selected < 0 || fixed_estimate > amplitude_allowance * direct_norm(state.field) / T(2))
                Rcpp::stop("Adaptive gas cancellation exhausts the degree arithmetic or refinement allowance.");
            auto& retained = state.contributions[selected];
            Complex other_orders(0, 0), other_correction(0, 0);
            for (int index = 0; index <= m; ++index)
                if (index != selected)
                    add_field(other_orders, other_correction, state.contributions[index].field);
            const int enlarged_degree = checked_profcn_integer(
                (5LL * (retained.degree + 1LL) + 3LL) / 4LL - 1LL,
                "adaptive gas cancellation refinement");
            if (!resolve_direct_degree(state, prepared, angular_workspace, selected,
                    enlarged_degree, retained.attempts, other_orders, true, chi_sw,
                    chi_body, xi, density_ratio, phi_difference, amplitude_allowance,
                    retained)) {
                state.rejected = true;
                return;
            }
            state.field = state.correction = Complex(0, 0);
            state.total_degree_estimate = 0;
            for (const auto& contribution : state.contributions) {
                add_field(state.field, state.correction, contribution.field);
                state.total_degree_estimate += contribution.estimate;
            }
            recovered = true;
        }
        if (recovered) {
            const T updated_allowance = amplitude_allowance * direct_norm(state.field) / T(2);
            state.small_orders = 0;
            for (int index = m; index >= 2 &&
                    direct_norm(state.contributions[index].field) <= updated_allowance; --index)
                ++state.small_orders;
            magnitude = direct_norm(state.contributions[m].field);
            state.previous_order = m > 0 ? direct_norm(state.contributions[m - 1].field) : T(-1);
            state.penultimate_order = m > 1 ? direct_norm(state.contributions[m - 2].field) : T(-1);
            order_pass = state.incidence.exact_axis;
            if (state.small_orders >= 2 && state.previous_order > T(0) &&
                    state.penultimate_order > T(0)) {
                const T ratio = std::max(magnitude / state.previous_order,
                    state.previous_order / state.penultimate_order);
                order_pass = ratio < T(1) && magnitude * ratio / (T(1) - ratio) <= updated_allowance;
            }
        }
        if (order_pass) {
            state.finished = true;
            return;
        }
    }
    state.penultimate_order = state.previous_order;
    state.previous_order = magnitude;
}

struct State {
    SpheroidalDirection<T> incidence;
    int previous_degree, previous_interior_degree;
    Complex field{0, 0}, correction{0, 0};
    T previous_order = -1, penultimate_order = -1;
    int small_orders = 0;
    bool finished = false;
    Quadrature<T> boundary_workspace;
    LegendreWorkspace<T> angular_workspace;
    State(const SpheroidalDirection<T>& direction, int seed, int exterior_seed)
        : incidence(direction), previous_degree(exterior_seed),
          previous_interior_degree(std::max(1, seed)) {
        boundary_workspace.incident.reset_start_on_extension = true;
    }
};

// Robin traces test pressure and normal-velocity continuity in separate spaces.
inline void advance(State& state, psms_angle_prepared::Workspace<T>& prepared,
        int m,
        T chi_sw, T chi_body, T xi, T density_ratio, T phi_difference) {
    const auto& incidence = state.incidence;
    const T pressure_reference = T(1) / precsqrt(T(2));
    const T derivative_reference = T(chi_sw) * precsqrt((incidence.cosine * incidence.cosine +
        T(xi) * T(xi) * incidence.sine * incidence.sine / ((T(xi) - T(1)) * (T(xi) + T(1)))) / T(6));
    const T amplitude_allowance = std::min(-expm1q(-logq(T(10)) * T(1e-8) / T(20)), psms_modal_rel_tol<T>());
    auto field_norm = [](const Complex& value) { return combine_norms(precabs(value.real()), precabs(value.imag())); };
    auto projection_defect = [](const OrderSolution<T>& order) {
        TailNorm<T> defect;
        for (size_t index = 0; index < order.traces[0].size(); ++index) {
            const auto pressure = order.traces[0][index] - order.traces[2][index];
            const auto derivative = order.traces[1][index] - order.traces[3][index];
            defect.append(derivative + Complex(0, 1) * pressure);
        }
        return defect.value();
    };
    auto& field = state.field;
    auto& correction = state.correction;
    auto& previous_order = state.previous_order;
    auto& penultimate_order = state.penultimate_order;
    auto& small_orders = state.small_orders;
    auto& previous_degree = state.previous_degree;
    auto& previous_interior_degree = state.previous_interior_degree;
    auto& boundary_workspace = state.boundary_workspace;
    auto& angular_workspace = state.angular_workspace;

    Rcpp::checkUserInterrupt();
    int degree = std::max(previous_degree, m + 1);
    int interior_degree = std::max(previous_interior_degree, m + 1);
    SpheroidalCoefficientOverlapCache<T> overlap_cache;
    int exterior_degree = -1;
    OrderSolution<T> order;
    BoundaryCheck<T> previous_boundary;
    std::vector<Complex> previous_coefficients;
    bool have_previous = false, accepted = false;
    for (int attempt = 0; attempt < 24; ++attempt) {
        Rcpp::checkUserInterrupt();
        // Consecutive angles reference the same current-order native arrays.
        auto& exterior = prepared.get(true, chi_sw, m, degree, xi - T(1)).values;
        // Interior-only refinement retains every exterior special-function result.
        if (exterior_degree != degree) {
            evaluate_coefficient_angular_values(exterior, m,
                std::vector<SpheroidalDirection<T>>{incidence}, angular_workspace);
            exterior_degree = degree;
        }
        order = solve_order<T>(m, degree, interior_degree, T(chi_body), T(xi), incidence.cosine,
            phi_difference, density_ratio, true, true, exterior, nullptr, &overlap_cache,
            true, false, false, true, false,
            &prepared.get(false, chi_body, m, interior_degree, xi - T(1)));
        if (order.singular)
            Rcpp::stop("Gas Robin system failed its native iterative and dense arithmetic safeguards.");
        check_boundary(order, T(chi_sw), T(xi), incidence,
            pressure_reference, derivative_reference, boundary_workspace);
        const auto budgets = boundary_budgets(order.boundary, m, pressure_reference, derivative_reference);
        const T projection = projection_defect(order);
        accepted = order.boundary.quadrature_converged && boundary_passes(order.boundary, budgets);
        const bool incident_unresolved = order.boundary.incident_pressure_error > budgets.first ||
            order.boundary.incident_derivative_error > budgets.second;
        // Each omitted Robin component is bounded by the orthogonal trace defect.
        const bool refine_interior = !incident_unresolved && projection > std::min(budgets.first, budgets.second);
        if (accepted) break;
        if (have_previous && coefficient_difference(previous_coefficients, order.coefficients) == T(0) &&
                order.boundary.pressure_error == previous_boundary.pressure_error &&
                order.boundary.derivative_error == previous_boundary.derivative_error &&
                order.boundary.incident_pressure_error == previous_boundary.incident_pressure_error &&
                order.boundary.incident_derivative_error == previous_boundary.incident_derivative_error)
            Rcpp::stop("Gas Robin coefficients and boundary residuals stagnated before acceptance.");
        previous_coefficients = order.coefficients;
        previous_boundary = order.boundary;
        have_previous = true;
        if (refine_interior) interior_degree = grow_degree(interior_degree, m);
        else degree = grow_degree(degree, m);
    }
    if (!accepted) Rcpp::stop("Gas Robin support exceeded its 24-trial execution budget.");
    previous_degree = degree;
    previous_interior_degree = interior_degree;
    add_field(field, correction, order.field);
    const T magnitude = field_norm(order.field);
    const T order_allowance = amplitude_allowance * field_norm(field) / T(2);
    small_orders = m >= 2 && magnitude <= order_allowance ? small_orders + 1 : 0;
    bool order_pass = incidence.exact_axis;
    if (small_orders >= 2 && previous_order > T(0) && penultimate_order > T(0)) {
        const T ratio = std::max(magnitude / previous_order, previous_order / penultimate_order);
        order_pass = ratio < T(1) && magnitude * ratio / (T(1) - ratio) <= order_allowance;
    }
    if (order_pass) state.finished = true;
    penultimate_order = previous_order;
    previous_order = magnitude;
}

// Each order shares native coefficients across angles with identical physical inputs.
inline std::vector<Complex> fields(T chi_sw, T chi_body, T xi,
        const std::vector<int>& degree_seeds, T density_ratio,
        const std::vector<SpheroidalDirection<T>>& directions, T phi_difference) {
    if (directions.empty() || directions.size() != degree_seeds.size() ||
            !(chi_sw > 0) || !(chi_body > 0) || !(xi > 1) || !(density_ratio > 0) ||
            !is_finite_real(chi_sw) || !is_finite_real(chi_body) ||
            !is_finite_real(xi) || !is_finite_real(density_ratio) ||
            !is_finite_real(phi_difference))
        throw std::invalid_argument("Invalid gas support inputs.");
    // Furusawa's exterior estimate seeds refinement rather than limiting support.
    const T exterior_seed = ceilq(T(2) * chi_sw * precsqrt((xi - T(1)) * (xi + T(1)))) +
        ceilq(chi_sw / T(2));
    if (!(exterior_seed >= 0) || exterior_seed > T(INT_MAX - 2))
        throw std::overflow_error("Gas support exceeds native indexing.");
    std::vector<State> states;
    states.reserve(directions.size());
    for (size_t index = 0; index < directions.size(); ++index) {
        if (degree_seeds[index] < 0)
            throw std::invalid_argument("Invalid gas degree seed.");
        states.emplace_back(directions[index], degree_seeds[index], std::max(1, int(exterior_seed)));
    }
    if (directions.size() == 1) {
        auto& state = states.front();
        state.finished = direct_field(chi_sw, chi_body, xi, degree_seeds.front(),
            density_ratio, state.incidence, phi_difference, state.field);
        // A rejected direct sum cannot seed the Robin solve.
        if (!state.finished) state.field = Complex(0, 0);
    } else {
        std::vector<DirectState> direct_states;
        direct_states.reserve(directions.size());
        for (size_t index = 0; index < directions.size(); ++index)
            direct_states.emplace_back(directions[index], degree_seeds[index]);
        if (chi_sw == chi_body && density_ratio == T(1)) {
            for (auto& state : direct_states) state.finished = true;
        } else {
            const T ts_allowance = -expm1q(-logq(T(10)) * T(1e-8) / T(20));
            const T amplitude_allowance = std::min(ts_allowance, psms_modal_rel_tol<T>());
            psms_angle_prepared::Workspace<T> direct_prepared;
            LegendreWorkspace<T> angular_workspace;
            // Interleave independent direct solves to reuse exact native degree blocks.
            for (int m = 0; ; ++m) {
                bool active = false;
                for (auto& state : direct_states) {
                    if (state.finished || state.rejected) continue;
                    active = true;
                    advance_direct(state, direct_prepared, angular_workspace, m,
                        chi_sw, chi_body, xi, density_ratio, phi_difference,
                        amplitude_allowance);
                }
                if (!active) break;
            }
        }
        for (size_t index = 0; index < states.size(); ++index) {
            states[index].finished = direct_states[index].finished;
            // A rejected direct sum never contributes to the Robin fallback.
            if (states[index].finished) states[index].field = direct_states[index].field;
        }
    }
    psms_angle_prepared::Workspace<T> prepared;
    for (int m = 0; ; ++m) {
        bool active = false;
        for (auto& state : states) {
            if (state.finished) continue;
            active = true;
            advance(state, prepared, m,
                chi_sw, chi_body, xi, density_ratio, phi_difference);
        }
        if (!active) break;
    }
    std::vector<Complex> result;
    result.reserve(states.size());
    for (const auto& state : states) result.push_back(state.field);
    return result;
}
}
#endif
