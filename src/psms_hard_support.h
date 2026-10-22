#ifndef ACOUSTICTS_PSMS_HARD_SUPPORT_H
#define ACOUSTICTS_PSMS_HARD_SUPPORT_H

#include <algorithm>
#include <array>
#include <climits>
#include <complex>
#include <map>
#include <stdexcept>
#include <utility>
#include <vector>

#if ACOUSTICTS_HAVE_QUADMATH

// These helpers use the normalized expansions and packed radial values in psms_backscatter.h.
namespace psms_hard_support {
using Real = acousticts_quad_t;
using Complex = std::complex<Real>;

inline Real magnitude(const Complex& value) {
    return psms_automatic::combine_norms(precabs(value.real()), precabs(value.imag()));
}

struct DegreeTail {
    bool resolved = false;
    Real estimate = 0, log_estimate = precnan<Real>();
    std::array<Real, 2> parity_estimates = {0, 0};
    std::array<Real, 2> log_ratios = {precnan<Real>(), precnan<Real>()};
    const char* reason = "transition_gate";
    int unresolved_parity = -1, unresolved_degree = -1;
};

inline bool coefficient_envelope_log(const ProfcnCoefficientBatchResult<Real>& exterior,
        int order, int index, bool rigid, std::vector<Real>& weights,
        Real& logarithmic_envelope, const char*& reason) {
    const int length = exterior.lengths[index], parity = index % 2;
    const Real* coefficients = exterior.coefficients.data() +
        size_t(exterior.capacity) * index;
    // DLMF 14.30.1 and 14.30.9 bound each unit-norm Ferrers function uniformly.
    while (weights.size() < size_t(order) + size_t(length)) {
        const Real degree = Real(weights.size());
        weights.push_back(precsqrt((Real(2) * degree + Real(1)) / Real(2)));
    }
    Real scale = 0;
    for (int offset = parity; offset < length; offset += 2) {
        if (!is_finite_real(coefficients[offset])) {
            reason = "nonfinite_coefficient";
            return false;
        }
        scale = std::max(scale, precabs(coefficients[offset]));
    }
    if (!(scale > 0)) {
        reason = "zero_coefficient_column";
        return false;
    }
    // Sum positive scaled coefficients without angular cancellation or squaring.
    Real weighted_sum = 0, correction = 0;
    for (int offset = parity; offset < length; offset += 2)
        psms_automatic::add(weighted_sum, correction,
            (precabs(coefficients[offset]) / scale) * weights[size_t(order) + offset]);
    if (!(weighted_sum > 0) || !is_finite_real(weighted_sum)) {
        reason = "invalid_coefficient_bound";
        return false;
    }
    const Real logarithmic_angular_bound = logq(scale) + logq(weighted_sum);
    const auto& functions = exterior.functions;
    const auto first = psms_packed::component(
        rigid ? functions.r1dc : functions.r1c,
        rigid ? functions.ir1de : functions.ir1e, size_t(index));
    const auto second = psms_packed::component(
        rigid ? functions.r2dc : functions.r2c,
        rigid ? functions.ir2de : functions.ir2e, size_t(index));
    if (first.mantissa == Real(0)) {
        reason = "zero_first_kind_radial";
        return false;
    }
    // Packed radial exponents preserve a positive ratio below the native range.
    Real logarithmic_radial_ratio = 0;
    if (second.mantissa != Real(0)) {
        const Real difference = logq(precabs(second.mantissa) / precabs(first.mantissa)) +
            Real(second.exponent - first.exponent) * logq(Real(10));
        logarithmic_radial_ratio = difference > Real(0) ?
            -difference - log1pq(expq(-Real(2) * difference)) / Real(2) :
            -log1pq(expq(Real(2) * difference)) / Real(2);
    }
    logarithmic_envelope = logq(order == 0 ? Real(1) : Real(2)) +
        Real(2) * logarithmic_angular_bound + logarithmic_radial_ratio;
    if (!is_finite_real(logarithmic_envelope)) {
        reason = "nonfinite_envelope_log";
        return false;
    }
    return true;
}

inline Real positive_tail_estimate(Real logarithmic_estimate) {
    const Real estimate = psms_automatic::finite_norm(expq(logarithmic_estimate));
    return estimate == Real(0) ? nextafterq(Real(0), Real(1)) : estimate;
}

inline DegreeTail degree_tail(const ProfcnCoefficientBatchResult<Real>& exterior, int order,
        int transition_degree, bool broadside, bool rigid, std::vector<Real>& weights) {
    DegreeTail result;
    bool have_tail = false;
    // Only three coefficient-envelope columns per active parity are inspected.
    for (int parity = 0; parity < 2; ++parity) {
        if (broadside && parity == 1) continue;
        result.unresolved_parity = parity;
        int last = exterior.functions.lnum - 1;
        if (last % 2 != parity) --last;
        result.unresolved_degree = order + last;
        if (last < 4 || order + last - 4 <= transition_degree) return result;
        std::array<Real, 3> logarithmic_envelopes;
        for (int sample = 0; sample < 3; ++sample) {
            const int index = last - 4 + 2 * sample;
            result.unresolved_degree = order + index;
            if (!coefficient_envelope_log(exterior, order, index, rigid, weights,
                    logarithmic_envelopes[sample], result.reason)) return result;
        }
        const Real logarithmic_ratio = std::max(
            logarithmic_envelopes[1] - logarithmic_envelopes[0],
            logarithmic_envelopes[2] - logarithmic_envelopes[1]);
        result.log_ratios[parity] = logarithmic_ratio;
        if (!is_finite_real(logarithmic_ratio) || !(logarithmic_ratio < Real(0))) {
            result.reason = "noncontracting_envelope";
            return result;
        }
        // The geometric extrapolation remains conditional on continued contraction.
        const Real logarithmic_tail = logarithmic_envelopes[2] + logarithmic_ratio -
            logq(-expm1q(logarithmic_ratio));
        if (!is_finite_real(logarithmic_tail)) {
            result.reason = "nonfinite_tail_log";
            return result;
        }
        result.parity_estimates[parity] = positive_tail_estimate(logarithmic_tail);
        if (have_tail) {
            const Real larger = std::max(result.log_estimate, logarithmic_tail);
            const Real smaller = std::min(result.log_estimate, logarithmic_tail);
            result.log_estimate = larger + log1pq(expq(smaller - larger));
        } else {
            result.log_estimate = logarithmic_tail;
            have_tail = true;
        }
    }
    result.estimate = positive_tail_estimate(result.log_estimate);
    result.reason = "resolved";
    result.unresolved_parity = result.unresolved_degree = -1;
    result.resolved = true;
    return result;
}

inline Real term_mass(const std::vector<Complex>& terms) {
    Real total = 0, correction = 0;
    for (const auto& term : terms)
        psms_automatic::add(total, correction, magnitude(term));
    return psms_automatic::finite_norm(total);
}
}

inline std::complex<acousticts_quad_t> psms_adaptive_hard_field(
        acousticts_quad_t chi, acousticts_quad_t xi,
        const SpheroidalDirection<acousticts_quad_t>& incidence,
        acousticts_quad_t phi_difference, bool rigid) {
    using namespace psms_automatic;
    using T = acousticts_quad_t;
    using Complex = std::complex<T>;
    if (!(chi > T(0) && xi > T(1)) || !is_finite_real(chi) || !is_finite_real(xi) ||
            !is_finite_real(phi_difference) || !is_finite_real(incidence.cosine) ||
            precabs(incidence.cosine) > T(1) || !is_finite_real(incidence.sine) ||
            incidence.sine < T(0) || incidence.sine > T(1))
        throw std::invalid_argument("Adaptive hard-boundary inputs are outside their native physical domain.");
    const bool broadside = incidence.cosine == T(0);
    // DLMF 30.3.3-4 and 30.13.9 justify these transition gates, not final cutoffs.
    const T ka = chi * xi, kb = chi * precsqrt((xi - T(1)) * (xi + T(1)));
    if (!is_finite_real(ka) || !is_finite_real(kb) || ka > T(INT_MAX - 1))
        throw std::length_error("The radial transition seed exceeds the native degree range.");
    const int transition_degree = static_cast<int>(ceilq(ka));
    const T absolute_allowance = psms_modal_abs_tol<T>();
    const T order_allowance = absolute_allowance / T(2);
    Complex field(0, 0), field_correction(0, 0);
    T degree_error = 0, degree_error_correction = 0, order_error = 0;
    std::array<T, 3> order_envelopes = {0, 0, 0};
    bool previous_order_pass = false;
    LegendreWorkspace<T> workspace;
    std::vector<T> angular_bound_weights;
    const std::vector<T> arguments = {incidence.cosine};
    const std::vector<SpheroidalDirection<T>> directions = {incidence};
    // Summable absolute allocations remain valid if later orders cancel the field.
    for (int m = 0;; ++m) {
        Rcpp::checkUserInterrupt();
        checked_profcn_integer(static_cast<long long>(m) + 1, "joint order count");
        const T degree_allowance = absolute_allowance /
            (T(2) * (T(m) + T(1)) * (T(m) + T(2)));
        // Request three tail samples per parity beyond the transition before testing.
        int degree = checked_profcn_integer(
            static_cast<long long>(std::max(m, transition_degree)) + 6,
            "joint initial degree");
        bool have_candidate = false;
        std::vector<Complex> previous_terms;
        OrderSolution<T> accepted;
        psms_hard_support::DegreeTail accepted_tail;
        // Only the unresolved order regenerates its native degree block.
        for (;;) {
            Rcpp::checkUserInterrupt();
            const int count = checked_profcn_integer(static_cast<long long>(degree) - m + 1,
                "joint inclusive degree count");
            auto exterior = cprofcn_mblock_coefficients<T>(chi, m, 1, count,
                arguments, 2, 0, xi - T(1), true);
            evaluate_coefficient_angular_values(exterior, m, directions, workspace);
            auto trial = solve_hard_order(m, degree, incidence.cosine,
                phi_difference, true, rigid, exterior, false, true);
            const auto tail = psms_hard_support::degree_tail(
                exterior, m, transition_degree, broadside, rigid, angular_bound_weights);
            T prefix_change = 0, prefix_correction = 0, added_mass = 0, added_correction = 0;
            if (have_candidate) {
                for (size_t index = 0; index < previous_terms.size(); ++index)
                    add(prefix_change, prefix_correction, psms_hard_support::magnitude(
                        trial.outgoing_terms[index] - previous_terms[index]));
                for (size_t index = previous_terms.size(); index < trial.outgoing_terms.size(); ++index)
                    add(added_mass, added_correction, psms_hard_support::magnitude(trial.outgoing_terms[index]));
            }
            const bool tail_pass = tail.resolved && tail.log_estimate <= logq(degree_allowance);
            const bool confirmation_pass = have_candidate && tail_pass &&
                finite_norm(prefix_change + added_mass) <= degree_allowance;
            if (confirmation_pass) {
                accepted = std::move(trial);
                accepted_tail = tail;
                break;
            }
            have_candidate = tail_pass;
            if (have_candidate) previous_terms = std::move(trial.outgoing_terms);
            else previous_terms.clear();
            degree = grow_degree(degree, m);
        }
        add_field(field, field_correction, accepted.field);
        add(degree_error, degree_error_correction, accepted_tail.estimate);
        const T envelope = finite_norm(psms_hard_support::term_mass(accepted.outgoing_terms) +
            accepted_tail.estimate);
        order_envelopes[0] = order_envelopes[1];
        order_envelopes[1] = order_envelopes[2];
        order_envelopes[2] = envelope;
        bool order_pass = false, order_confirmed = false;
        order_error = T(0);
        // The outgoing azimuthal estimate is conditional on continued contraction.
        if (!incidence.exact_axis && m >= 2 && T(m - 2) > kb &&
                order_envelopes[0] > 0 && order_envelopes[1] > 0 && order_envelopes[2] > 0) {
            const T ratio = std::max(order_envelopes[1] / order_envelopes[0],
                order_envelopes[2] / order_envelopes[1]);
            if (ratio < 1) {
                order_error = finite_norm(envelope * ratio / (T(1) - ratio));
                order_pass = order_error <= order_allowance;
                order_confirmed = previous_order_pass && order_pass &&
                    finite_norm(envelope + order_error) <= order_allowance;
            }
        }
        previous_order_pass = order_pass;
        if (incidence.exact_axis || order_confirmed) break;
        if (m == INT_MAX - 1) throw std::length_error("Hard-boundary support exceeds the native order range.");
    }
    const T estimated_error = finite_norm(degree_error + order_error);
    if (estimated_error > absolute_allowance)
        throw std::runtime_error("The combined hard-boundary tail estimate exceeds its absolute allowance.");
    finite(field);
    return field;
}

namespace psms_hard_support {
struct BatchState {
    SpheroidalDirection<Real> incidence;
    bool rigid = false, finished = false, previous_order_pass = false;
    Complex field = {0, 0}, field_correction = {0, 0};
    Real degree_error = 0, degree_error_correction = 0, order_error = 0;
    std::array<Real, 3> order_envelopes = {0, 0, 0};
};

struct DegreeState {
    int degree = 0;
    bool accepted = false, have_candidate = false;
    std::vector<Complex> previous_terms;
    psms_automatic::OrderSolution<Real> solution;
    DegreeTail tail;
};
}

inline std::vector<std::complex<acousticts_quad_t>> psms_adaptive_rigid_release_fields(
        acousticts_quad_t chi, acousticts_quad_t xi,
        const std::vector<SpheroidalDirection<acousticts_quad_t>>& incidences,
        const std::vector<bool>& rigid, acousticts_quad_t phi_difference) {
    using namespace psms_automatic;
    using namespace psms_hard_support;
    using T = acousticts_quad_t;
    if (incidences.size() != rigid.size() || incidences.empty())
        throw std::invalid_argument("Fixed-rigid/pressure-release batch cases need matching dimensions.");
    if (!(chi > T(0) && xi > T(1)) || !is_finite_real(chi) || !is_finite_real(xi) ||
            !is_finite_real(phi_difference))
        throw std::invalid_argument("Fixed-rigid/pressure-release batch inputs are outside their native physical domain.");
    for (const auto& incidence : incidences)
        if (!is_finite_real(incidence.cosine) || precabs(incidence.cosine) > T(1) ||
                !is_finite_real(incidence.sine) || incidence.sine < T(0) || incidence.sine > T(1))
            throw std::invalid_argument("Fixed-rigid/pressure-release batch inputs are outside their native physical domain.");
    // Share only complete native blocks; every case retains its own modal test.
    const T ka = chi * xi, kb = chi * precsqrt((xi - T(1)) * (xi + T(1)));
    if (!is_finite_real(ka) || !is_finite_real(kb) || ka > T(INT_MAX - 1))
        throw std::length_error("The radial transition seed exceeds the native degree range.");
    const int transition_degree = static_cast<int>(ceilq(ka));
    const T absolute_allowance = psms_modal_abs_tol<T>();
    const T order_allowance = absolute_allowance / T(2);
    std::vector<BatchState> states(incidences.size());
    for (size_t case_index = 0; case_index < states.size(); ++case_index) {
        states[case_index].incidence = incidences[case_index];
        states[case_index].rigid = rigid[case_index];
    }
    LegendreWorkspace<T> workspace;
    std::vector<T> angular_bound_weights;
    for (int m = 0;; ++m) {
        Rcpp::checkUserInterrupt();
        checked_profcn_integer(static_cast<long long>(m) + 1, "joint order count");
        const T degree_allowance = absolute_allowance /
            (T(2) * (T(m) + T(1)) * (T(m) + T(2)));
        const int initial_degree = checked_profcn_integer(
            static_cast<long long>(std::max(m, transition_degree)) + 6,
            "joint initial degree");
        std::vector<DegreeState> degree_states(states.size());
        size_t pending = 0;
        for (size_t case_index = 0; case_index < states.size(); ++case_index) {
            if (states[case_index].finished) continue;
            degree_states[case_index].degree = initial_degree;
            ++pending;
        }
        // Cases requesting the same degree share one native expansion.
        while (pending > 0) {
            Rcpp::checkUserInterrupt();
            std::map<int, std::vector<size_t>> groups;
            for (size_t case_index = 0; case_index < states.size(); ++case_index)
                if (!states[case_index].finished && !degree_states[case_index].accepted)
                    groups[degree_states[case_index].degree].push_back(case_index);
            const auto& group = *groups.begin();
            const int degree = group.first;
            const int count = checked_profcn_integer(
                static_cast<long long>(degree) - m + 1, "joint inclusive degree count");
            std::vector<T> arguments;
            std::vector<SpheroidalDirection<T>> directions;
            arguments.reserve(group.second.size());
            directions.reserve(group.second.size());
            for (size_t case_index : group.second) {
                arguments.push_back(states[case_index].incidence.cosine);
                directions.push_back(states[case_index].incidence);
            }
            auto exterior = cprofcn_mblock_coefficients<T>(
                chi, m, 1, count, arguments, 2, 0, xi - T(1), true);
            evaluate_coefficient_angular_values(exterior, m, directions, workspace);
            for (size_t angle_index = 0; angle_index < group.second.size(); ++angle_index) {
                const size_t case_index = group.second[angle_index];
                auto& state = states[case_index];
                auto& step = degree_states[case_index];
                auto trial = solve_hard_order(
                    m, degree, state.incidence.cosine, phi_difference, true,
                    state.rigid, exterior, false, true, static_cast<int>(angle_index));
                const auto tail = degree_tail(exterior, m, transition_degree,
                    state.incidence.cosine == T(0), state.rigid, angular_bound_weights);
                T prefix_change = 0, prefix_correction = 0;
                T added_mass = 0, added_correction = 0;
                if (step.have_candidate) {
                    for (size_t index = 0; index < step.previous_terms.size(); ++index)
                        add(prefix_change, prefix_correction, magnitude(
                            trial.outgoing_terms[index] - step.previous_terms[index]));
                    for (size_t index = step.previous_terms.size();
                            index < trial.outgoing_terms.size(); ++index)
                        add(added_mass, added_correction,
                            magnitude(trial.outgoing_terms[index]));
                }
                const bool tail_pass = tail.resolved &&
                    tail.log_estimate <= logq(degree_allowance);
                const bool confirmation_pass = step.have_candidate && tail_pass &&
                    finite_norm(prefix_change + added_mass) <= degree_allowance;
                if (confirmation_pass) {
                    step.solution = std::move(trial);
                    step.tail = tail;
                    step.accepted = true;
                    --pending;
                    continue;
                }
                step.have_candidate = tail_pass;
                if (step.have_candidate) step.previous_terms = std::move(trial.outgoing_terms);
                else step.previous_terms.clear();
                step.degree = grow_degree(degree, m);
            }
        }
        // Advance each independent azimuthal tail after its degree is accepted.
        size_t unfinished = 0;
        for (size_t case_index = 0; case_index < states.size(); ++case_index) {
            auto& state = states[case_index];
            if (state.finished) continue;
            auto& step = degree_states[case_index];
            add_field(state.field, state.field_correction, step.solution.field);
            add(state.degree_error, state.degree_error_correction, step.tail.estimate);
            const T envelope = finite_norm(term_mass(step.solution.outgoing_terms) +
                step.tail.estimate);
            state.order_envelopes[0] = state.order_envelopes[1];
            state.order_envelopes[1] = state.order_envelopes[2];
            state.order_envelopes[2] = envelope;
            bool order_pass = false, order_confirmed = false;
            state.order_error = T(0);
            if (!state.incidence.exact_axis && m >= 2 && T(m - 2) > kb &&
                    state.order_envelopes[0] > 0 && state.order_envelopes[1] > 0 &&
                    state.order_envelopes[2] > 0) {
                const T ratio = std::max(
                    state.order_envelopes[1] / state.order_envelopes[0],
                    state.order_envelopes[2] / state.order_envelopes[1]);
                if (ratio < 1) {
                    state.order_error = finite_norm(envelope * ratio / (T(1) - ratio));
                    order_pass = state.order_error <= order_allowance;
                    order_confirmed = state.previous_order_pass && order_pass &&
                        finite_norm(envelope + state.order_error) <= order_allowance;
                }
            }
            state.previous_order_pass = order_pass;
            state.finished = state.incidence.exact_axis || order_confirmed;
            if (!state.finished) ++unfinished;
        }
        if (unfinished == 0) break;
        if (m == INT_MAX - 1)
            throw std::length_error("Fixed-rigid/pressure-release support exceeds the native order range.");
    }
    std::vector<Complex> fields;
    fields.reserve(states.size());
    for (const auto& state : states) {
        if (finite_norm(state.degree_error + state.order_error) > absolute_allowance)
            throw std::runtime_error(
                "The combined fixed-rigid/pressure-release tail estimate exceeds its absolute allowance.");
        finite(state.field);
        fields.push_back(state.field);
    }
    return fields;
}

#endif
#endif
