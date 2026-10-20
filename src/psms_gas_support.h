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
    field = Complex(0, 0);
    Complex field_correction(0, 0);
    T previous_order = -1, penultimate_order = -1, total_degree_estimate = 0;
    int small_orders = 0;
    // Independent order blocks permit refinement without discarding accepted work.
    for (int m = 0; !(chi_sw == chi_body && density_ratio == T(1)); ++m) {
        Rcpp::checkUserInterrupt();
        int degree = std::max(degree_seed, m + 5);
        SpheroidalCoefficientOverlapCache<T> overlap_cache;
        OrderSolution<T> order;
        T estimate = 0, allowance = 0;
        bool degree_pass = false;
        for (int refinement = 0; refinement < 12; ++refinement) {
            auto exterior = cprofcn_mblock_coefficients<T>(T(chi_sw), m, 1,
                degree - m + 1, arguments, 2, 0, T(xi) - T(1), true);
            evaluate_coefficient_angular_values(exterior, m, directions, workspace);
            order = solve_order<T>(m, degree, degree, T(chi_body), T(xi), incidence.cosine,
                phi_difference, density_ratio, true, false, exterior, nullptr,
                &overlap_cache, false, false, true, false, true);
            if (order.singular) return false;
            if (!order.embedded_valid[0] || !order.embedded_valid[1])
                return false;
            const T scale = norm(field + order.field);
            allowance = amplitude_allowance * scale / (T(2) * T(m + 1) * T(m + 2));
            const T last = norm(order.field - order.embedded_fields[0]);
            const T preceding = norm(order.embedded_fields[0] - order.embedded_fields[1]);
            // Account for cancellation between modal terms when identifying quad roundoff.
            const T arithmetic_floor = T(64) * preceps<T>() * order.embedded_roundoff_scale;
            degree_pass = false;
            // Geometric continuation is conditional on contraction, not a universal bound.
            if (order.embedded_valid[0] && order.embedded_valid[1]) {
                if (last <= arithmetic_floor && preceding <= arithmetic_floor) {
                    estimate = arithmetic_floor;
                    degree_pass = estimate <= allowance;
                } else if (preceding > T(0) && last < preceding) {
                    const T ratio = last / preceding;
                    estimate = last / (T(1) - ratio);
                    degree_pass = estimate <= allowance && preceding <= allowance;
                }
            }
            if (degree_pass) break;
            degree = checked_profcn_integer((5LL * (degree + 1LL) + 3LL) / 4LL - 1LL,
                "optimized degree refinement");
        }
        if (!degree_pass) Rcpp::stop("Optimized degree tail did not converge within twelve refinements.");
        total_degree_estimate += estimate;
        add_field(field, field_correction, order.field);
        const T magnitude = norm(order.field);
        const T order_allowance = amplitude_allowance * norm(field) / T(2);
        small_orders = m >= 2 && magnitude <= order_allowance ? small_orders + 1 : 0;
        bool order_pass = incidence.exact_axis;
        if (small_orders >= 2 && previous_order > T(0) && penultimate_order > T(0)) {
            const T ratio = std::max(magnitude / previous_order, previous_order / penultimate_order);
            order_pass = ratio < T(1) && magnitude * ratio / (T(1) - ratio) <= order_allowance;
        }
        if (order_pass) {
            if (total_degree_estimate > amplitude_allowance * norm(field) / T(2))
                Rcpp::stop("Cancellation invalidated the accumulated degree-error allowance.");
            break;
        }
        penultimate_order = previous_order;
        previous_order = magnitude;
    }
    return true;
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
        int m, int shared_degree, int shared_interior_degree,
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
    int degree = std::max(std::max(previous_degree, shared_degree), m + 1);
    int interior_degree = std::max(std::max(previous_interior_degree, shared_interior_degree), m + 1);
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
            !(chi_sw > 0) || !(chi_body > 0) || !(xi > 1) || !(density_ratio > 0))
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
        auto& state = states.back();
        state.finished = direct_field(chi_sw, chi_body, xi, degree_seeds[index],
            density_ratio, state.incidence, phi_difference, state.field);
        // Discard the entire rejected direct sum before changing formulation.
        if (!state.finished) state.field = Complex(0, 0);
    }
    psms_angle_prepared::Workspace<T> prepared;
    int shared_degree = 0, shared_interior_degree = 0;
    for (int m = 0; ; ++m) {
        bool active = false;
        for (auto& state : states) {
            if (state.finished) continue;
            active = true;
            advance(state, prepared, m, shared_degree, shared_interior_degree,
                chi_sw, chi_body, xi, density_ratio, phi_difference);
            shared_degree = std::max(shared_degree, state.previous_degree);
            shared_interior_degree = std::max(shared_interior_degree, state.previous_interior_degree);
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
