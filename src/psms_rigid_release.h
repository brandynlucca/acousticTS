#ifndef ACOUSTICTS_PSMS_RIGID_RELEASE_H
#define ACOUSTICTS_PSMS_RIGID_RELEASE_H

#include <algorithm>
#include <array>
#include <climits>
#include <complex>
#include <map>
#include <stdexcept>
#include <utility>
#include <vector>

#if ACOUSTICTS_HAVE_QUADMATH

extern "C" void profcn_retain_geometry_order_quad(int order);

// These helpers use the normalized expansions and packed radial values in psms_backscatter.h.
namespace psms_rigid_release {
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

inline int refinement_start_from_hint(int degree, int order, int hint) {
    // Advance only on the existing refinement grid without using the hint as a cutoff.
    while (degree < hint) {
        const long long increment = std::max(2LL,
            (static_cast<long long>(degree) - order + 1) / 4);
        const long long next = static_cast<long long>(degree) + increment;
        if (next > static_cast<long long>(hint)) break;
        degree = psms_automatic::grow_degree(degree, order);
    }
    return degree;
}

inline void validate_rejected_radials(const ProfcnCoefficientBatchResult<Real>& exterior,
        bool broadside, bool rigid) {
    psms_packed::ScalingDiagnostics scaling;
    // Preserve native radial failures before discarding a rejected trial.
    for (int index = 0; index < exterior.functions.lnum; ++index) {
        if (broadside && index % 2 == 1) continue;
        const auto radial = psms_packed::exterior_pair(exterior.functions, 0, index, scaling);
        const Complex divisor = rigid ? radial.derivative : radial.value;
        if (divisor == Complex(0, 0))
            throw std::runtime_error("The outgoing fixed-rigid/pressure-release radial divisor vanished in native precision.");
    }
}

inline Real term_mass(const std::vector<Complex>& terms) {
    Real total = 0, correction = 0;
    for (const auto& term : terms)
        psms_automatic::add(total, correction, magnitude(term));
    return psms_automatic::finite_norm(total);
}
}


extern "C" void profcn_set_radial_quadrature_quad(int);
extern "C" int profcn_get_radial_quadrature_quad();
extern "C" void profcn_get_radial_quadrature_errors_quad(int,acousticts_quad_t*,acousticts_quad_t*,int*);
extern "C" void profcn_set_radial_quadrature_even_quad(int);
extern "C" int profcn_get_radial_quadrature_even_quad();
extern "C" int profcn_radial_quadrature_available_quad();
struct PsmsRadialQuadratureScope {
    bool enabled=false;
    int previous=profcn_get_radial_quadrature_quad();
    int previous_even=profcn_get_radial_quadrature_even_quad();
    explicit PsmsRadialQuadratureScope(bool broadside) {
        enabled=profcn_radial_quadrature_available_quad()!=0;
        profcn_set_radial_quadrature_quad(enabled?1:0);
        profcn_set_radial_quadrature_even_quad(broadside?1:0);
    }
    ~PsmsRadialQuadratureScope() {
        profcn_set_radial_quadrature_quad(previous);
        profcn_set_radial_quadrature_even_quad(previous_even);
    }
};
struct PsmsRadialQuadratureError {
    bool valid=false;
    bool angular_evaluated=false;
    acousticts_quad_t field_error=0;
    std::vector<acousticts_quad_t> relative_errors;
};
enum class PsmsRadialUncertaintyRejection {
    none=0, invalid_input=1, shifted_underflow=2, packed_underflow=3,
    invalid_divisor=4, relative_underflow=5, first_ratio_underflow=6,
    multiplier_underflow=7, invalid_output=8
};
struct PsmsScaledRadialError {
    bool valid=false;
    PsmsRadialUncertaintyRejection rejection=PsmsRadialUncertaintyRejection::invalid_input;
    acousticts_quad_t relative_error=0, multiplier_error=0;
};
inline PsmsScaledRadialError psms_scaled_radial_uncertainty(
        const psms_packed::Decimal<acousticts_quad_t>& first,
        const psms_packed::Decimal<acousticts_quad_t>& second,
        acousticts_quad_t packed_error,long long error_exponent) {
    using T=acousticts_quad_t;
    using Rejection=PsmsRadialUncertaintyRejection;
    PsmsScaledRadialError result;
    if(!(packed_error>=0) || !is_finite_real(packed_error) ||
            !is_finite_real(first.mantissa) || !is_finite_real(second.mantissa)) return result;
    const auto exponent=std::max(first.exponent,second.exponent);
    psms_packed::ScalingDiagnostics scaling;
    const T a=psms_packed::shifted(first,exponent,scaling);
    const T b=psms_packed::shifted(second,exponent,scaling);
    if(scaling.scaled_underflows>0) {
        result.rejection=Rejection::shifted_underflow;
        return result;
    }
    const T uncertainty=psms_packed::scale_decimal(packed_error,error_exponent-exponent);
    if(packed_error>0 && uncertainty==0) {
        result.rejection=Rejection::packed_underflow;
        return result;
    }
    const T high_modulus=hypotq(a,b);
    const T lower_modulus=hypotq(a,std::max(T(0),precabs(b)-uncertainty));
    if(!(lower_modulus>0) || !is_finite_real(uncertainty) || !is_finite_real(high_modulus)) {
        result.rejection=Rejection::invalid_divisor;
        return result;
    }
    result.relative_error=uncertainty/lower_modulus;
    if(uncertainty>0 && result.relative_error==0) {
        result.rejection=Rejection::relative_underflow;
        return result;
    }
    const T first_ratio=precabs(a)/high_modulus;
    if(a!=0 && first_ratio==0) {
        result.rejection=Rejection::first_ratio_underflow;
        return result;
    }
    result.multiplier_error=first_ratio*result.relative_error;
    if(first_ratio>0 && result.relative_error>0 && result.multiplier_error==0) {
        result.rejection=Rejection::multiplier_underflow;
        return result;
    }
    if(!is_finite_real(result.relative_error) || !is_finite_real(result.multiplier_error)) {
        result.rejection=Rejection::invalid_output;
        return result;
    }
    result.valid=true;
    result.rejection=Rejection::none;
    return result;
}
inline PsmsRadialQuadratureError psms_radial_quadrature_error(
        ProfcnCoefficientBatchResult<acousticts_quad_t>& exterior,
        int order,bool broadside,bool rigid,
        const std::vector<SpheroidalDirection<acousticts_quad_t>>& directions,
        psms_automatic::LegendreWorkspace<acousticts_quad_t>& workspace) {
    using T=acousticts_quad_t;
    PsmsRadialQuadratureError result;
    const int count=exterior.functions.lnum;
    std::vector<int> valid(count,0);
    std::vector<T> radial_error(count,0),derivative_error(count,0);
    profcn_get_radial_quadrature_errors_quad(count,radial_error.data(),derivative_error.data(),valid.data());
    std::vector<T> multiplier_errors(count,0);
    result.relative_errors.resize(count,0);
    bool have_error=false;
    for(int index=0;index<count;++index) {
        if(broadside && index%2) continue;
        if(!valid[index]) return result;
        const auto& data=exterior.functions;
        const T packed_error=rigid?derivative_error[index]:radial_error[index];
        if(!(packed_error>=0) || !is_finite_real(packed_error)) return result;
        if(packed_error==0) continue;
        const auto first=psms_packed::component(rigid?data.r1dc:data.r1c,
            rigid?data.ir1de:data.ir1e,index);
        const auto second=psms_packed::component(rigid?data.r2dc:data.r2c,
            rigid?data.ir2de:data.ir2e,index);
        const auto scaled=psms_scaled_radial_uncertainty(first,second,packed_error,
            static_cast<long long>(rigid?data.ir2de[index]:data.ir2e[index]));
        if(!scaled.valid) return result;
        result.relative_errors[index]=scaled.relative_error;
        const T error=scaled.multiplier_error;
        if(!is_finite_real(error) || !is_finite_real(result.relative_errors[index])) return result;
        multiplier_errors[index]=error;
        have_error=have_error || error>0;
    }
    T correction=0;
    if(have_error) {
        psms_automatic::evaluate_coefficient_angular_values(exterior,order,directions,workspace);
        result.angular_evaluated=true;
    }
    for(int index=0;have_error && index<count;++index) {
        if(multiplier_errors[index]==0) continue;
        const T angular=extract_angular_value_from_mblock<T>(exterior.functions,0,index,0);
        const T contribution=(order==0?T(1):T(2))*angular*angular*multiplier_errors[index];
        if(angular!=0 && contribution==0) return result;
        psms_automatic::add(result.field_error,correction,contribution);
    }
    result.valid=is_finite_real(result.field_error);
    return result;
}
inline psms_rigid_release::DegreeTail psms_quadrature_degree_tail(
        const ProfcnCoefficientBatchResult<acousticts_quad_t>& exterior,int order,
        int transition_degree,bool broadside,bool rigid,std::vector<acousticts_quad_t>& weights,
        const std::vector<acousticts_quad_t>& relative_errors) {
    using T=acousticts_quad_t;
    using namespace psms_rigid_release;
    DegreeTail result;
    bool have_tail=false;
    for(int parity=0;parity<2;++parity) {
        if(broadside && parity==1) continue;
        int last=exterior.functions.lnum-1;
        if(last%2!=parity) --last;
        if(last<4 || order+last-4<=transition_degree) return result;
        std::array<T,3> upper,lower;
        for(int sample=0;sample<3;++sample) {
            const int index=last-4+2*sample;
            T logarithm=0;
            if(!coefficient_envelope_log(exterior,order,index,rigid,weights,logarithm,result.reason)) return result;
            const T uncertainty=relative_errors[index];
            if(!(uncertainty>=0 && uncertainty<1)) return result;
            upper[sample]=logarithm+log1pq(uncertainty);
            lower[sample]=logarithm+log1pq(-uncertainty);
        }
        const T ratio=std::max(upper[1]-lower[0],upper[2]-lower[1]);
        if(!is_finite_real(ratio) || !(ratio<0)) return result;
        const T tail=upper[2]+ratio-logq(-expm1q(ratio));
        if(!is_finite_real(tail)) return result;
        result.parity_estimates[parity]=positive_tail_estimate(tail);
        result.log_ratios[parity]=ratio;
        if(have_tail) {
            const T larger=std::max(result.log_estimate,tail);
            const T smaller=std::min(result.log_estimate,tail);
            result.log_estimate=larger+log1pq(expq(smaller-larger));
        } else { result.log_estimate=tail; have_tail=true; }
    }
    result.estimate=positive_tail_estimate(result.log_estimate);
    result.resolved=true;
    result.reason="resolved_with_quadrature_estimate";
    return result;
}
inline std::complex<acousticts_quad_t> psms_adaptive_rigid_release_field(
        acousticts_quad_t chi, acousticts_quad_t xi,
        const SpheroidalDirection<acousticts_quad_t>& incidence,
        acousticts_quad_t phi_difference, bool rigid, bool retain_geometry_orders = true) {
    using namespace psms_automatic;
    using T = acousticts_quad_t;
    using Complex = std::complex<T>;
    PsmsRadialQuadratureScope quadrature_scope(incidence.cosine==T(0));
    if (!(chi > T(0) && xi > T(1)) || !is_finite_real(chi) || !is_finite_real(xi) ||
            !is_finite_real(phi_difference) || !is_finite_real(incidence.cosine) ||
            precabs(incidence.cosine) > T(1) || !is_finite_real(incidence.sine) ||
            incidence.sine < T(0) || incidence.sine > T(1))
        throw std::invalid_argument("Adaptive fixed-rigid/pressure-release inputs are outside their native physical domain.");
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
    int previous_order_seed = 0;
    bool original_quadrature_for_point=!quadrature_scope.enabled;
    // Summable absolute allocations remain valid if later orders cancel the field.
    for (int m = 0;; ++m) {
        Rcpp::checkUserInterrupt();
        checked_profcn_integer(static_cast<long long>(m) + 1, "joint order count");
        // Independent point solves need geometry only for the current order.
        if (!retain_geometry_orders) profcn_retain_geometry_order_quad(m);
        const T degree_allowance = absolute_allowance /
            (T(2) * (T(m) + T(1)) * (T(m) + T(2)));
        // Request three tail samples per parity beyond the transition before testing.
        int degree = checked_profcn_integer(
            static_cast<long long>(std::max(m, transition_degree)) + 6,
            "joint initial degree");
        degree = psms_rigid_release::refinement_start_from_hint(degree, m, previous_order_seed);
        int previous_degree = degree;
        bool have_candidate = false;
        bool original_quadrature = original_quadrature_for_point;
        const int original_initial_degree = degree;
        T previous_quadrature_error=0, accepted_quadrature_error=0;
        std::vector<Complex> previous_terms;
        OrderSolution<T> accepted;
        psms_rigid_release::DegreeTail accepted_tail;
        // Only the unresolved order regenerates its native degree block.
        for (;;) {
            Rcpp::checkUserInterrupt();
            const int count = checked_profcn_integer(static_cast<long long>(degree) - m + 1,
                "joint inclusive degree count");
            auto fallback = [&]() {
                original_quadrature=true;
                original_quadrature_for_point=true;
                degree=original_initial_degree;
                previous_degree=degree;
                have_candidate=false;
                previous_terms.clear();
                previous_quadrature_error=0;
                profcn_retain_geometry_order_quad(-1);
                profcn_set_radial_quadrature_quad(0);

            };
            profcn_set_radial_quadrature_quad(original_quadrature?0:1);
            if(!original_quadrature) profcn_set_radial_quadrature_even_quad(broadside?1:0);
            ProfcnCoefficientBatchResult<T> exterior;
            try {
                exterior=cprofcn_mblock_coefficients<T>(chi,m,1,count,
                    arguments,2,0,xi-T(1),true);
            } catch(const std::runtime_error&) {
                if(original_quadrature) throw;
                fallback(); continue;
            }

            psms_rigid_release::DegreeTail ordinary_tail;
            try {
                ordinary_tail=psms_rigid_release::degree_tail(
                    exterior,m,transition_degree,broadside,rigid,angular_bound_weights);
            } catch(const std::runtime_error&) {
                if(original_quadrature) throw;
                fallback(); continue;
            }
            auto tail=ordinary_tail;
            T quadrature_error=0;
            bool angular_evaluated=false;
            if(!original_quadrature) {
                PsmsRadialQuadratureError estimate;
                try {
                    estimate=psms_radial_quadrature_error(exterior,m,broadside,rigid,directions,workspace);
                    angular_evaluated=estimate.angular_evaluated;
                } catch(const std::runtime_error&) {
                    estimate.valid=false;
                }
                if(!estimate.valid) {
                    fallback(); continue;
                }
                quadrature_error=estimate.field_error;
                if(!(quadrature_error<degree_allowance)) {
                    fallback(); continue;
                }
                const bool radial_uncertainty=std::any_of(estimate.relative_errors.begin(),
                    estimate.relative_errors.end(),[](T value) { return value>T(0); });
                if(radial_uncertainty)
                    tail=psms_quadrature_degree_tail(exterior,m,transition_degree,broadside,
                        rigid,angular_bound_weights,estimate.relative_errors);
                if(ordinary_tail.resolved && ordinary_tail.estimate<=degree_allowance &&
                        (!tail.resolved || tail.estimate+quadrature_error>degree_allowance)) {
                    fallback(); continue;
                }
            }
            const bool tail_pass=tail.resolved && (quadrature_error==0
                ? tail.log_estimate<=logq(degree_allowance)
                : tail.estimate+quadrature_error<=degree_allowance);
            // Unresolved tails cannot pass confirmation and need no angular field sum.
            if (!tail_pass) {
                try { psms_rigid_release::validate_rejected_radials(exterior,broadside,rigid); }
                catch(const std::runtime_error&) {
                    if(original_quadrature) throw;
                    fallback(); continue;
                }
                have_candidate = false;
                previous_terms.clear();
                previous_degree = degree;
                degree = grow_degree(degree, m);
                continue;
            }
            if(!angular_evaluated) evaluate_coefficient_angular_values(exterior, m, directions, workspace);
            OrderSolution<T> trial;
            try {
                trial=solve_rigid_release_order(m,degree,incidence.cosine,
                    phi_difference,true,rigid,exterior,false,true);
            } catch(const std::runtime_error&) {
                if(original_quadrature) throw;
                fallback(); continue;
            }
            T prefix_change = 0, prefix_correction = 0, added_mass = 0, added_correction = 0;
            if (have_candidate) {
                for (size_t index = 0; index < previous_terms.size(); ++index)
                    add(prefix_change, prefix_correction, psms_rigid_release::magnitude(
                        trial.outgoing_terms[index] - previous_terms[index]));
                for (size_t index = previous_terms.size(); index < trial.outgoing_terms.size(); ++index)
                    add(added_mass, added_correction, psms_rigid_release::magnitude(trial.outgoing_terms[index]));
            }
            const bool confirmation_pass = have_candidate && tail_pass &&
                finite_norm(prefix_change + added_mass + previous_quadrature_error + quadrature_error) <= degree_allowance;
            if (confirmation_pass) {
                accepted = std::move(trial);
                accepted_tail = tail;
                accepted_quadrature_error=quadrature_error;

                previous_order_seed = previous_degree;
                break;
            }
            if(have_candidate && prefix_change+added_mass<=degree_allowance &&
                    prefix_change+added_mass+previous_quadrature_error+quadrature_error>degree_allowance &&
                    !original_quadrature) {
                fallback(); continue;
            }
            previous_quadrature_error=quadrature_error;
            have_candidate = tail_pass;
            if (have_candidate) previous_terms = std::move(trial.outgoing_terms);
            else previous_terms.clear();
            previous_degree = degree;
            degree = grow_degree(degree, m);
        }
        add_field(field, field_correction, accepted.field);
        add(degree_error, degree_error_correction, accepted_tail.estimate + accepted_quadrature_error);
        const T envelope = finite_norm(psms_rigid_release::term_mass(accepted.outgoing_terms) +
            accepted_tail.estimate + accepted_quadrature_error);

        order_envelopes[0] = order_envelopes[1];
        order_envelopes[1] = order_envelopes[2];
        order_envelopes[2] = envelope;
        bool order_pass = false, order_confirmed = false;
        order_error = T(0);
        // The outgoing azimuthal estimate is conditional on continued contraction.
        if (!incidence.exact_axis && m >= 2 && T(m - 2) > kb &&
                order_envelopes[0] > 0 && order_envelopes[1] > 0 && order_envelopes[2] > 0) {
            const T ratio=std::max(order_envelopes[1]/order_envelopes[0],
                order_envelopes[2]/order_envelopes[1]);
            if (ratio < 1) {
                order_error = finite_norm(envelope * ratio / (T(1) - ratio));
                order_pass = order_error <= order_allowance;
                order_confirmed = previous_order_pass && order_pass &&
                    finite_norm(envelope + order_error) <= order_allowance;
            }
        }
        previous_order_pass = order_pass;
        if (incidence.exact_axis || order_confirmed) break;
        if (m == INT_MAX - 1) throw std::length_error("Fixed-rigid/pressure-release support exceeds the native order range.");
    }
    const T estimated_error = finite_norm(degree_error + order_error);
    if (estimated_error > absolute_allowance)
        throw std::runtime_error("The combined fixed-rigid/pressure-release tail estimate exceeds its absolute allowance.");
    finite(field);
    return field;
}
namespace psms_rigid_release {
struct BatchState {
    SpheroidalDirection<Real> incidence;
    bool rigid = false, finished = false, previous_order_pass = false;
    int previous_order_seed = 0;
    Complex field = {0, 0}, field_correction = {0, 0};
    Real degree_error = 0, degree_error_correction = 0, order_error = 0;
    std::array<Real, 3> order_envelopes = {0, 0, 0};
};

struct DegreeState {
    int degree = 0, previous_degree = 0;
    bool accepted = false, have_candidate = false;
    std::vector<Complex> previous_terms;
    psms_automatic::OrderSolution<Real> solution;
    DegreeTail tail;
};
}

inline std::vector<std::complex<acousticts_quad_t>> psms_adaptive_rigid_release_fields(
        acousticts_quad_t chi, acousticts_quad_t xi,
        const std::vector<SpheroidalDirection<acousticts_quad_t>>& incidences,
        const std::vector<bool>& rigid, acousticts_quad_t phi_difference,
        bool retain_geometry_orders = true) {
    using namespace psms_automatic;
    using namespace psms_rigid_release;
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
        // Keep shared degree refinement without retaining completed orders.
        if (!retain_geometry_orders) profcn_retain_geometry_order_quad(m);
        const T degree_allowance = absolute_allowance /
            (T(2) * (T(m) + T(1)) * (T(m) + T(2)));
        const int initial_degree = checked_profcn_integer(
            static_cast<long long>(std::max(m, transition_degree)) + 6,
            "joint initial degree");
        std::vector<DegreeState> degree_states(states.size());
        size_t pending = 0;
        for (size_t case_index = 0; case_index < states.size(); ++case_index) {
            if (states[case_index].finished) continue;
            auto& step = degree_states[case_index];
            step.degree = refinement_start_from_hint(initial_degree, m,
                states[case_index].previous_order_seed);
            step.previous_degree = step.degree;
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
            std::vector<DegreeTail> tails;
            std::vector<bool> tail_passes;
            tails.reserve(group.second.size());
            tail_passes.reserve(group.second.size());
            bool any_tail_pass = false;
            // Test each case before evaluating angular values for the shared block.
            for (size_t case_index : group.second) {
                const auto& state = states[case_index];
                tails.push_back(degree_tail(exterior, m, transition_degree,
                    state.incidence.cosine == T(0), state.rigid, angular_bound_weights));
                const auto& tail = tails.back();
                const bool tail_pass = tail.resolved &&
                    tail.log_estimate <= logq(degree_allowance);
                tail_passes.push_back(tail_pass);
                any_tail_pass = any_tail_pass || tail_pass;
            }
            // Preserve every angular index when any case needs the shared values.
            if (any_tail_pass)
                evaluate_coefficient_angular_values(exterior, m, directions, workspace);
            for (size_t angle_index = 0; angle_index < group.second.size(); ++angle_index) {
                const size_t case_index = group.second[angle_index];
                auto& state = states[case_index];
                auto& step = degree_states[case_index];
                const auto& tail = tails[angle_index];
                const bool tail_pass = tail_passes[angle_index];
                if (!tail_pass) {
                    validate_rejected_radials(exterior, state.incidence.cosine == T(0), state.rigid);
                    step.have_candidate = false;
                    step.previous_terms.clear();
                    step.previous_degree = degree;
                    step.degree = grow_degree(degree, m);
                    continue;
                }
                auto trial = solve_rigid_release_order(
                    m, degree, state.incidence.cosine, phi_difference, true,
                    state.rigid, exterior, false, true, static_cast<int>(angle_index));
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
                const bool confirmation_pass = step.have_candidate && tail_pass &&
                    finite_norm(prefix_change + added_mass) <= degree_allowance;
                if (confirmation_pass) {
                    step.solution = std::move(trial);
                    step.tail = tail;
                    state.previous_order_seed = step.previous_degree;
                    step.accepted = true;
                    --pending;
                    continue;
                }
                step.have_candidate = tail_pass;
                if (step.have_candidate) step.previous_terms = std::move(trial.outgoing_terms);
                else step.previous_terms.clear();
                step.previous_degree = degree;
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
