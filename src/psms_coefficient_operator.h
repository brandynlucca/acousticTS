#pragma once

#include <atomic>
#include <cfenv>
#include <cstdint>
#include <memory>
#include <string>

// Complete Legendre projections implement the Robin map without truncating overlaps.
// DLMF 30.8.1 and 14.17.6 supply the expansion and orthogonality identities.
namespace psms_coefficient_operator_cache {
constexpr std::size_t limit_bytes = std::size_t(16) * 1024 * 1024;
struct EntryBase {
    std::string key;
    std::size_t bytes = 0;
    std::unique_ptr<EntryBase> next;
    explicit EntryBase(std::string value) : key(std::move(value)) {}
    virtual ~EntryBase() = default;
};
template<typename Value> void append_exact(std::string& key, const Value& value) {
    key.append(reinterpret_cast<const char*>(&value), sizeof(value));
}

template<typename T> struct Entry : EntryBase {
    const T norm_bound;
    const std::vector<std::complex<T>> diagonal;
    Entry(std::string key, T bound, const std::vector<std::complex<T>>& values)
        : EntryBase(std::move(key)), norm_bound(bound), diagonal(values) {}
};
struct State {
    bool enabled = true;
    std::unique_ptr<EntryBase> first;
    EntryBase* last = nullptr;
    std::size_t retained_bytes = 0, peak_retained_bytes = 0, entries = 0;
    std::size_t requests = 0, hits = 0, misses = 0, bypasses = 0;
    std::size_t evictions = 0, oversized = 0, allocation_bypasses = 0;

    void remove_first(bool eviction) {
        auto removed = std::move(first);
        first = std::move(removed->next);
        retained_bytes -= removed->bytes;
        --entries;
        if (!first) last = nullptr;
        if (eviction) ++evictions;
    }
    void clear() {
        while (first) remove_first(false);
    }
    void reset(bool clear_entries) {
        if (clear_entries) clear();
        requests = hits = misses = bypasses = 0;
        evictions = oversized = allocation_bypasses = 0;
        peak_retained_bytes = retained_bytes;
    }
    ~State() { clear(); }
};
inline thread_local State state;

template<typename T> std::string input_key(const ProfcnCoefficientBatchResult<T>& exterior,
        const ProfcnCoefficientBatchResult<T>& interior, int exterior_parity, int interior_parity,
        const std::vector<std::complex<T>>& response, const std::vector<std::complex<T>>& ratio) {
    ++state.requests;
    if (!state.enabled || !exterior.cache_identity || !interior.cache_identity) {
        ++state.bypasses;
        return {};
    }
    try {
        std::string key;
        const unsigned precision = std::is_same_v<T, double> ? 64 : 128;
        append_exact(key, precision);
        for (const auto* batch : {&exterior, &interior}) {
            append_exact(key, batch->cache_identity);
            append_exact(key, batch->functions.lnum);
            append_exact(key, batch->functions.m_count);
            append_exact(key, batch->capacity);
            const std::size_t length_count = batch->lengths.size();
            append_exact(key, length_count);
            for (int length : batch->lengths) append_exact(key, length);
        }
        append_exact(key, exterior_parity);
        append_exact(key, interior_parity);
        for (const auto* values : {&response, &ratio}) {
            const std::size_t count = values->size();
            append_exact(key, count);
            for (const auto& value : *values) {
                append_exact(key, value.real());
                append_exact(key, value.imag());
            }
        }
        const int rounding = std::fegetround();
        append_exact(key, rounding);
        return key;
    } catch (const std::bad_alloc&) {
        ++state.allocation_bypasses;
        return {};
    }
}

template<typename T> bool restore(const std::string& key, T& norm_bound,
        std::vector<std::complex<T>>& diagonal) {
    if (key.empty()) return false;
    for (const EntryBase* entry = state.first.get(); entry; entry = entry->next.get()) {
        if (entry->key != key) continue;
        try {
            const auto& stored = *static_cast<const Entry<T>*>(entry);
            diagonal = stored.diagonal;
            norm_bound = stored.norm_bound;
            ++state.hits;
            return true;
        } catch (const std::bad_alloc&) {
            state.clear();
            ++state.allocation_bypasses;
            break;
        }
    }
    ++state.misses;
    return false;
}

template<typename T> std::size_t storage_bytes(const std::string& key,
        const std::vector<std::complex<T>>& diagonal) {
    std::size_t bytes = sizeof(Entry<T>);
    if (key.capacity() >= limit_bytes - bytes) return limit_bytes + 1;
    bytes += key.capacity() + 1;
    if (diagonal.capacity() > (limit_bytes - bytes) / sizeof(std::complex<T>)) return limit_bytes + 1;
    return bytes + diagonal.capacity() * sizeof(std::complex<T>);
}

template<typename T> void retain(std::string key, T norm_bound,
        const std::vector<std::complex<T>>& diagonal) {
    if (key.empty()) return;
    const std::size_t estimate = storage_bytes(key, diagonal);
    if (estimate > limit_bytes) { ++state.oversized; return; }
    while (state.first && estimate > limit_bytes - state.retained_bytes) state.remove_first(true);
    try {
        auto stored = std::make_unique<Entry<T>>(std::move(key), norm_bound, diagonal);
        stored->bytes = storage_bytes(stored->key, stored->diagonal);
        if (stored->bytes > limit_bytes) { ++state.oversized; return; }
        while (state.first && stored->bytes > limit_bytes - state.retained_bytes) state.remove_first(true);
        state.retained_bytes += stored->bytes;
        state.peak_retained_bytes = std::max(state.peak_retained_bytes, state.retained_bytes);
        ++state.entries;
        EntryBase* next_last = stored.get();
        if (state.last) state.last->next = std::move(stored);
        else state.first = std::move(stored);
        state.last = next_last;
    } catch (const std::bad_alloc&) {
        ++state.allocation_bypasses;
    }
}
} // namespace psms_coefficient_operator_cache

namespace psms_implicit_overlap {
template<typename T> using Complex = std::complex<T>;

template<typename T> T upward(T value) {
    if (!is_finite_real(value) || value < T(0)) return precnan<T>();
    const T infinity = T(std::numeric_limits<double>::infinity());
#if ACOUSTICTS_HAVE_QUADMATH
    if constexpr (std::is_same_v<T, acousticts_quad_t>) return nextafterq(value, infinity);
    else
#endif
    return std::nextafter(value, infinity);
}

template<typename T> T positive_sum(T first, T second) {
    if (!is_finite_real(first) || !is_finite_real(second) || first < T(0) || second < T(0))
        return precnan<T>();
    if (first == T(0)) return second;
    if (second == T(0)) return first;
    return upward(first + second);
}

template<typename T> T positive_product(T first, T second) {
    if (!is_finite_real(first) || !is_finite_real(second) || first < T(0) || second < T(0))
        return precnan<T>();
    if (first == T(0) || second == T(0)) return T(0);
    return upward(first * second);
}

// Use the existing exact native accumulator for every coefficient projection.
template<typename T, typename Coefficient, typename Value>
Complex<T> projection(std::size_t count, Coefficient coefficient, Value value) {
#if ACOUSTICTS_HAVE_EXACT_QUAD_DOT
    if constexpr (std::is_same_v<T, acousticts_quad_t>) {
        if (std::fegetround() == FE_TONEAREST) {
            // Share each coefficient load across independent exact component sums.
            static thread_local quad_dot_product::Accumulator real_accumulator;
            static thread_local quad_dot_product::Accumulator imaginary_accumulator;
            real_accumulator.reset();
            imaginary_accumulator.reset();
            bool valid = true;
            for (std::size_t index = 0; index < count; ++index) {
                const T factor = coefficient(index);
                const Complex<T> sample = value(index);
                valid = real_accumulator.add(factor, sample.real()) && valid;
                valid = imaginary_accumulator.add(factor, sample.imag()) && valid;
            }
            return valid ? Complex<T>(real_accumulator.rounded(), imaginary_accumulator.rounded()) :
                Complex<T>(precnan<T>(), precnan<T>());
        }
    }
#endif
    Complex<T> sum(0, 0), correction(0, 0);
    for (std::size_t index = 0; index < count; ++index)
        psms_robin_solver::add(coefficient(index) * value(index), sum, correction);
    return sum;
}

// Outward rounding bounds the stored finite coefficients without orthogonality assumptions.
template<typename T, typename Coefficient, typename Value>
T positive_projection(std::size_t count, Coefficient coefficient, Value value) {
#if ACOUSTICTS_HAVE_EXACT_QUAD_DOT
    if constexpr (std::is_same_v<T, acousticts_quad_t>) {
        if (std::fegetround() == FE_TONEAREST) {
            static thread_local quad_dot_product::Accumulator accumulator;
            accumulator.reset();
            bool valid = true, positive = false;
            for (std::size_t index = 0; index < count; ++index) {
                const T first = precabs(coefficient(index)), second = value(index);
                if (!is_finite_real(second) || second < T(0)) return precnan<T>();
                positive = positive || (first > T(0) && second > T(0));
                valid = accumulator.add(first, second) && valid;
            }
            if (!valid) return precnan<T>();
            return positive ? upward(accumulator.rounded()) : T(0);
        }
    }
#endif
    T sum = 0;
    for (std::size_t index = 0; index < count; ++index) {
        const T first = precabs(coefficient(index)), second = value(index);
        if (!is_finite_real(first) || !is_finite_real(second) || second < T(0))
            return precnan<T>();
        sum = positive_sum(sum, positive_product(first, second));
    }
    return sum;
}

// A view retains every native coefficient and pads only beyond each recorded length.
template<typename T> struct CoefficientBlock {
    const ProfcnCoefficientBatchResult<T>& batch;
    int parity, columns = 0;
    std::size_t support = 0;

    CoefficientBlock(const ProfcnCoefficientBatchResult<T>& values, int selected_parity)
        : batch(values), parity(selected_parity) {
        if (batch.functions.lnum >= 0 && (parity == 0 || parity == 1))
            columns = batch.functions.lnum / 2 + (parity == 0 ? batch.functions.lnum % 2 : 0);
    }

    bool initialize() {
        support = 0;
        if ((parity != 0 && parity != 1) || batch.functions.m_count != 1 ||
                batch.functions.lnum < 0 || batch.capacity < 0 ||
                batch.lengths.size() != std::size_t(batch.functions.lnum) ||
                !psms_robin_solver::valid_dimensions(batch.coefficients.size(),
                    batch.capacity, batch.functions.lnum)) return false;
        for (int column = 0; column < columns; ++column) {
            const std::size_t degree = std::size_t(parity) + 2 * std::size_t(column);
            const int length = batch.lengths[degree];
            if (length < 1 || length > batch.capacity) return false;
            const std::size_t count = length > parity ?
                (std::size_t(length - parity) + 1) / 2 : 0;
            support = std::max(support, count);
            for (std::size_t position = 0; position < count; ++position)
                if (!is_finite_real(coefficient(column, position))) return false;
        }
        return true;
    }

    T coefficient(std::size_t column, std::size_t position) const {
        const std::size_t degree = std::size_t(parity) + 2 * column;
        const std::size_t legendre = std::size_t(parity) + 2 * position;
        if (legendre >= std::size_t(batch.lengths[degree])) return T(0);
        return batch.coefficients[std::size_t(batch.capacity) * degree + legendre];
    }

    bool forward(const std::vector<Complex<T>>& input, std::vector<Complex<T>>& output,
            std::size_t common_support) const {
        if (input.size() != std::size_t(columns) || common_support < support) return false;
        for (const auto& value : input) if (!psms_robin_solver::finite(value)) return false;
        output.assign(common_support, Complex<T>(0, 0));
        for (std::size_t position = 0; position < support; ++position) {
            output[position] = projection<T>(input.size(),
                [&](std::size_t column) { return coefficient(column, position); },
                [&](std::size_t column) { return input[column]; });
            if (!psms_robin_solver::finite(output[position])) return false;
        }
        return true;
    }

    bool transpose(const std::vector<Complex<T>>& input, std::vector<Complex<T>>& output) const {
        if (input.size() < support) return false;
        for (const auto& value : input) if (!psms_robin_solver::finite(value)) return false;
        output.resize(columns);
        for (int column = 0; column < columns; ++column) {
            const int length = batch.lengths[parity + 2 * column];
            const std::size_t count = length > parity ? (std::size_t(length - parity) + 1) / 2 : 0;
            output[column] = projection<T>(count,
                [&](std::size_t position) { return coefficient(column, position); },
                [&](std::size_t position) { return input[position]; });
            if (!psms_robin_solver::finite(output[column])) return false;
        }
        return true;
    }

    bool absolute_forward(const std::vector<T>& input, std::vector<T>& output,
            std::size_t common_support) const {
        if (input.size() != std::size_t(columns) || common_support < support) return false;
        output.assign(common_support, T(0));
        for (std::size_t position = 0; position < support; ++position) {
            output[position] = positive_projection<T>(input.size(),
                [&](std::size_t column) { return coefficient(column, position); },
                [&](std::size_t column) { return input[column]; });
            if (!is_finite_real(output[position])) return false;
        }
        return true;
    }

    bool absolute_transpose(const std::vector<T>& input, std::vector<T>& output) const {
        if (input.size() < support) return false;
        output.resize(columns);
        for (int column = 0; column < columns; ++column) {
            const int length = batch.lengths[parity + 2 * column];
            const std::size_t count = length > parity ? (std::size_t(length - parity) + 1) / 2 : 0;
            output[column] = positive_projection<T>(count,
                [&](std::size_t position) { return coefficient(column, position); },
                [&](std::size_t position) { return input[position]; });
            if (!is_finite_real(output[column])) return false;
        }
        return true;
    }
};

template<typename T> struct Operator {
    CoefficientBlock<T> exterior, interior;
    int rows, size;
    const std::vector<Complex<T>>& interior_response;
    const std::vector<Complex<T>>& exterior_ratio;
    std::vector<Complex<T>> diagonal;
    T norm_bound = 0;
    std::size_t common_support = 0;
    bool ready = false;

    Operator(const ProfcnCoefficientBatchResult<T>& external_coefficients,
            const ProfcnCoefficientBatchResult<T>& internal_coefficients, int parity,
            const std::vector<Complex<T>>& response, const std::vector<Complex<T>>& ratio)
        : exterior(external_coefficients, parity), interior(internal_coefficients, parity),
          rows(interior.columns), size(exterior.columns), interior_response(response), exterior_ratio(ratio) {}

    bool initialize() {
        ready = false;
        norm_bound = 0;
        diagonal.clear();
        if (!exterior.initialize() || !interior.initialize() ||
                interior_response.size() != std::size_t(rows) ||
                exterior_ratio.size() != std::size_t(size)) return false;
        for (const auto& value : interior_response) if (!psms_robin_solver::finite(value)) return false;
        for (const auto& value : exterior_ratio) if (!psms_robin_solver::finite(value)) return false;
        common_support = std::max(exterior.support, interior.support);

        // Reuse only immutable bounds after the original dimension and finite checks.
        auto cache_key = psms_coefficient_operator_cache::input_key(
            exterior.batch, interior.batch, exterior.parity, interior.parity, interior_response, exterior_ratio);
        if (psms_coefficient_operator_cache::restore(cache_key, norm_bound, diagonal)) {
            ready = true;
            return true;
        }

        // Bound |De|^T |Di| |v| |Di|^T |De| 1 using only positive coefficient maps.
        std::vector<T> ones(size, T(1)), legendre, weighted, bound;
        if (!exterior.absolute_forward(ones, legendre, common_support) ||
                !interior.absolute_transpose(legendre, weighted)) return false;
        for (int row = 0; row < rows; ++row) {
            const T response_norm = positive_sum(precabs(interior_response[row].real()),
                precabs(interior_response[row].imag()));
            weighted[row] = positive_product(response_norm, weighted[row]);
            if (!is_finite_real(weighted[row])) return false;
        }
        if (!interior.absolute_forward(weighted, legendre, common_support) ||
                !exterior.absolute_transpose(legendre, bound)) return false;
        diagonal = exterior_ratio;
        for (int column = 0; column < size; ++column) {
            const T diagonal_norm = positive_sum(precabs(diagonal[column].real()),
                precabs(diagonal[column].imag()));
            const T total = positive_sum(diagonal_norm, bound[column]);
            if (!is_finite_real(total)) return false;
            norm_bound = std::max(norm_bound, total);
        }
        psms_coefficient_operator_cache::retain(std::move(cache_key), norm_bound, diagonal);
        ready = true;
        return true;
    }

    bool weighted_legendre(const std::vector<Complex<T>>& input,
            std::vector<Complex<T>>& output) const {
        std::vector<Complex<T>> exterior_legendre, interior_modes;
        if (!exterior.forward(input, exterior_legendre, common_support) ||
                !interior.transpose(exterior_legendre, interior_modes)) return false;
        for (int row = 0; row < rows; ++row) {
            interior_modes[row] *= interior_response[row];
            if (!psms_robin_solver::finite(interior_modes[row])) return false;
        }
        return interior.forward(interior_modes, output, common_support);
    }

    // Fuse the last coefficient projection with the affine terms in the exact sum.
    bool action(const std::vector<Complex<T>>& input, std::vector<Complex<T>>& output,
            const std::vector<Complex<T>>* rhs = nullptr) const {
        if (!ready || input.size() != std::size_t(size) ||
                (rhs && rhs->size() != std::size_t(size))) return false;
        if (rhs) for (const auto& value : *rhs) if (!psms_robin_solver::finite(value)) return false;
        std::vector<Complex<T>> mapped;
        if (!weighted_legendre(input, mapped)) return false;
        output.assign(size, Complex<T>(0, 0));
        const T sign = rhs ? T(-1) : T(1);
        for (int column = 0; column < size; ++column) {
            const int length = exterior.batch.lengths[exterior.parity + 2 * column];
            const std::size_t count = length > exterior.parity ?
                (std::size_t(length - exterior.parity) + 1) / 2 : 0;
#if ACOUSTICTS_HAVE_EXACT_QUAD_DOT
            if constexpr (std::is_same_v<T, acousticts_quad_t>) {
                if (std::fegetround() == FE_TONEAREST) {
                    // Keep each affine component in its original exact summation order.
                    static thread_local quad_dot_product::Accumulator real_accumulator;
                    static thread_local quad_dot_product::Accumulator imaginary_accumulator;
                    real_accumulator.reset();
                    imaginary_accumulator.reset();
                    bool valid = !rhs || real_accumulator.add((*rhs)[column].real(), T(1));
                    valid = real_accumulator.add(sign * exterior_ratio[column].real(), input[column].real()) && valid;
                    valid = real_accumulator.add(-sign * exterior_ratio[column].imag(), input[column].imag()) && valid;
                    valid = (!rhs || imaginary_accumulator.add((*rhs)[column].imag(), T(1))) && valid;
                    valid = imaginary_accumulator.add(sign * exterior_ratio[column].real(), input[column].imag()) && valid;
                    valid = imaginary_accumulator.add(sign * exterior_ratio[column].imag(), input[column].real()) && valid;
                    for (std::size_t position = 0; position < count; ++position) {
                        const T factor = -sign * exterior.coefficient(column, position);
                        const Complex<T> sample = mapped[position];
                        valid = real_accumulator.add(factor, sample.real()) && valid;
                        valid = imaginary_accumulator.add(factor, sample.imag()) && valid;
                    }
                    if (!valid) return false;
                    output[column] = {real_accumulator.rounded(), imaginary_accumulator.rounded()};
                    if (!psms_robin_solver::finite(output[column])) return false;
                    continue;
                }
            }
#endif
            Complex<T> correction(0, 0);
            output[column] = rhs ? (*rhs)[column] : Complex<T>(0, 0);
            psms_robin_solver::add(sign * exterior_ratio[column] * input[column], output[column], correction);
            for (std::size_t position = 0; position < count; ++position)
                psms_robin_solver::add(-sign * exterior.coefficient(column, position) * mapped[position],
                    output[column], correction);
            if (!psms_robin_solver::finite(output[column])) return false;
        }
        return true;
    }

    bool apply(const std::vector<Complex<T>>& input, std::vector<Complex<T>>& output) const {
        return action(input, output);
    }
};

// Form Di^T De x for interior boundary traces without storing the overlap matrix.
template<typename T> std::vector<Complex<T>> overlap_projection(
        const ProfcnCoefficientBatchResult<T>& exterior, const ProfcnCoefficientBatchResult<T>& interior,
        int parity, const std::vector<Complex<T>>& input) {
    CoefficientBlock<T> external(exterior, parity), internal(interior, parity);
    std::vector<Complex<T>> legendre, output;
    if (!external.initialize() || !internal.initialize() ||
            !external.forward(input, legendre, std::max(external.support, internal.support)) ||
            !internal.transpose(legendre, output))
        throw std::runtime_error("Implicit Robin overlap projection dimension or range.");
    return output;
}

// Form De^T Di diag(v) Di^T De x for the complete coupled right-hand side.
template<typename T> std::vector<Complex<T>> apply_map(
        const ProfcnCoefficientBatchResult<T>& exterior, const ProfcnCoefficientBatchResult<T>& interior,
        int parity, const std::vector<Complex<T>>& response, const std::vector<Complex<T>>& input) {
    CoefficientBlock<T> external(exterior, parity), internal(interior, parity);
    std::vector<Complex<T>> legendre, weighted, output;
    if (!external.initialize() || !internal.initialize() || response.size() != std::size_t(internal.columns))
        throw std::runtime_error("Implicit Robin map dimension or coefficient range.");
    const std::size_t support = std::max(external.support, internal.support);
    if (!external.forward(input, legendre, support) || !internal.transpose(legendre, weighted))
        throw std::runtime_error("Implicit Robin map projection range.");
    for (std::size_t row = 0; row < weighted.size(); ++row) {
        if (!psms_robin_solver::finite(response[row]))
            throw std::runtime_error("Nonfinite implicit Robin response.");
        weighted[row] *= response[row];
        if (!psms_robin_solver::finite(weighted[row]))
            throw std::runtime_error("Implicit Robin response product range.");
    }
    if (!internal.forward(weighted, legendre, support) || !external.transpose(legendre, output))
        throw std::runtime_error("Implicit Robin map projection range.");
    return output;
}
} // namespace psms_implicit_overlap

namespace psms_robin_solver {
template<typename T> Result<T> solve_implicit(
        psms_implicit_overlap::Operator<T>& operation, const std::vector<Complex<T>>& rhs,
        const Options<T>& options = Options<T>()) {
    Result<T> result;
    if (rhs.size() > std::size_t(std::numeric_limits<int>::max())) {
        result.status = "dimension_range";
        return result;
    }
    const int size = int(rhs.size());
    if (operation.size != size) {
        result.status = "dimension_mismatch";
        return result;
    }
    if (options.restart <= 0 || options.maximum_iterations < 0 ||
            !is_finite_real(options.relative_tolerance) || !(options.relative_tolerance > T(0)) ||
            !is_finite_real(options.absolute_tolerance) || options.absolute_tolerance < T(0)) {
        result.status = "invalid_controls";
        return result;
    }
    for (std::size_t index = 0; index < rhs.size(); ++index) {
        if (!finite(rhs[index])) {
            result.status = "nonfinite_input";
            return result;
        }
    }
    result.solution.assign(size, Complex<T>(0, 0));


    if (!operation.initialize()) { result.status = "operator_bound_range"; return result; }
    result.operator_norm = operation.norm_bound;
    const T diagonal_guard = T(64) * preceps<T>() * std::max(T(1), result.operator_norm);
    for (const auto& value : operation.diagonal) {
        if (!(magnitude(value) > diagonal_guard)) {
            result.status = "diagonal_preconditioner_unavailable";
            return result;
        }
    }
    std::vector<Complex<T>> residual(rhs), candidate;
    const T rhs_norm = norm_infinity(rhs);
    if (!is_finite_real(rhs_norm)) { result.status = "rhs_norm_range"; return result; }
    auto inspect = [&](const std::vector<Complex<T>>& value) {
        ++result.true_residual_checks;
        if (!operation.action(value, residual, &rhs)) { result.status = "true_residual_range"; return -1; }
        result.residual_norm = norm_infinity(residual);
        const T reference = result.operator_norm * norm_infinity(value) + rhs_norm;
        const T allowance = options.absolute_tolerance + options.relative_tolerance * reference;
        if (!is_finite_real(result.residual_norm) || !is_finite_real(reference) || !is_finite_real(allowance)) {
            result.status = "residual_budget_range";
            return -1;
        }
        result.backward_error = reference == T(0) ? result.residual_norm : result.residual_norm / reference;
        result.rhs_relative_residual = rhs_norm == T(0) ? result.residual_norm : result.residual_norm / rhs_norm;
        const T rhs_allowance = T(64) * preceps<T>() * rhs_norm;
        return result.residual_norm <= allowance && result.residual_norm <= rhs_allowance ? 1 : 0;
    };
    auto accept = [&]() {
        result.success = true;
        result.needs_lu_fallback = false;
        result.status = "converged_true_residual";
    };
    int checked = inspect(result.solution);
    if (checked < 0) return result;
    if (checked > 0) { accept(); return result; }
    const int restart = std::min(size, options.restart);
    // Retain each actual native diagonal-preconditioned Arnoldi vector.
    while (result.iterations < options.maximum_iterations) {
        const int width = std::min(restart, options.maximum_iterations - result.iterations);
        const T beta = norm_two(residual);
        if (!is_finite_real(beta) || !(beta > T(0))) {
            result.status = "residual_norm_range";
            return result;
        }
        std::vector<std::vector<Complex<T>>> basis(std::size_t(width) + 1, std::vector<Complex<T>>(size));
        std::vector<std::vector<Complex<T>>> preconditioned(width, std::vector<Complex<T>>(size));
        std::vector<Complex<T>> hessenberg((std::size_t(width) + 1) * width), transformed(std::size_t(width) + 1);
        std::vector<Complex<T>> sines(width), coefficients;
        std::vector<T> cosines(width);
        const auto base = result.solution;
        for (int row = 0; row < size; ++row) basis[0][row] = residual[row] / beta;
        transformed[0] = beta;
        bool restart_ready = false;
        for (int column = 0; column < width; ++column) {

            for (int row = 0; row < size; ++row) {
                if (!divide(basis[column][row], operation.diagonal[row],
                        preconditioned[column][row])) {
                    result.status = "diagonal_preconditioner_range";
                    return result;
                }
            }
            std::vector<Complex<T>> next;
            if (!operation.apply(preconditioned[column], next)) { result.status = "arnoldi_product_range"; return result; }
            // Two modified Gram-Schmidt passes preserve native complex orthogonality.
            for (int pass = 0; pass < 2; ++pass) {
                for (int row = 0; row <= column; ++row) {
                    const auto projection = inner_product(basis[row], next);
                    if (!finite(projection)) { result.status = "arnoldi_inner_product_range"; return result; }
                    hessenberg[std::size_t(row) * width + column] += projection;
                    for (int index = 0; index < size; ++index) next[index] -= projection * basis[row][index];
                }
            }
            const T next_norm = norm_two(next);
            if (!is_finite_real(next_norm)) { result.status = "arnoldi_norm_range"; return result; }
            hessenberg[std::size_t(column + 1) * width + column] = next_norm;
            if (next_norm > T(0))
                for (int row = 0; row < size; ++row) basis[column + 1][row] = next[row] / next_norm;
            for (int previous = 0; previous < column; ++previous) {
                auto& upper = hessenberg[std::size_t(previous) * width + column];
                auto& lower = hessenberg[std::size_t(previous + 1) * width + column];
                const auto value = cosines[previous] * upper + sines[previous] * lower;
                lower = -std::conj(sines[previous]) * upper + cosines[previous] * lower;
                upper = value;
            }
            auto& diagonal = hessenberg[std::size_t(column) * width + column];
            if (!rotation(diagonal, hessenberg[std::size_t(column + 1) * width + column],
                    cosines[column], sines[column], diagonal)) {
                result.status = "arnoldi_rotation_range";
                return result;
            }
            hessenberg[std::size_t(column + 1) * width + column] = 0;
            transformed[column + 1] = -std::conj(sines[column]) * transformed[column];
            transformed[column] *= cosines[column];
            result.estimated_residual = magnitude(transformed[column + 1]);
            ++result.iterations;
            const T trigger = options.absolute_tolerance + options.relative_tolerance * rhs_norm;
            const bool at_end = column + 1 == width;
            if (!(result.estimated_residual <= trigger || next_norm == T(0) || at_end)) continue;
            coefficients.assign(column + 1, Complex<T>(0, 0));
            T triangular_norm = 0;
            for (int row = 0; row <= column; ++row)
                for (int col = row; col <= column; ++col)
                    triangular_norm = std::max(triangular_norm, magnitude(hessenberg[std::size_t(row) * width + col]));
            const T triangular_guard = T(64) * preceps<T>() * std::max(T(1), triangular_norm);
            for (int row = column; row >= 0; --row) {
                Complex<T> value = transformed[row], correction(0, 0);
                for (int col = row + 1; col <= column; ++col)
                    add(-hessenberg[std::size_t(row) * width + col] * coefficients[col], value, correction);
                const auto pivot = hessenberg[std::size_t(row) * width + row];
                if (!(magnitude(pivot) > triangular_guard) || !divide(value, pivot, coefficients[row])) {
                    result.status = "projected_pivot_guard";
                    return result;
                }
            }
            candidate = base;
            for (int row = 0; row < size; ++row) {
                Complex<T> correction(0, 0);
                for (int col = 0; col <= column; ++col)
                    add(preconditioned[col][row] * coefficients[col], candidate[row], correction);
            }
            checked = inspect(candidate);
            if (checked < 0) return result;
            // Accept only after evaluating the original factorized residual.
            if (checked > 0) {
                result.solution = std::move(candidate);
                accept();
                return result;
            }
            if (next_norm == T(0)) { result.status = "arnoldi_breakdown_without_convergence"; return result; }
            if (at_end) {
                if (candidate == base) { result.status = "native_iteration_stagnation"; return result; }
                result.solution = std::move(candidate);
                restart_ready = true;
            }
        }
        if (!restart_ready) { result.status = "restart_state_invalid"; return result; }
        ++result.restarts;
    }
    result.status = "iteration_budget_exhausted";
    return result;
}


}
