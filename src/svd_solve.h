#ifndef ACOUSTICTS_SVD_SOLVE_H
#define ACOUSTICTS_SVD_SOLVE_H

#include <RcppArmadillo.h>
#include <algorithm>
#include <limits>

// Compute the minimum-norm least-squares solution from an economical SVD.
// This is the same numerical fallback formerly obtained by explicitly forming
// a pseudoinverse, but avoids that construction and Armadillo's sort_index path.
template <typename ResultType, typename RightHandSideType>
inline bool acousticts_svd_solve(ResultType& result,
                                 const arma::cx_mat& matrix,
                                 const RightHandSideType& right_hand_side) {
    if (matrix.n_rows == 0 || matrix.n_cols == 0 ||
        matrix.n_rows != right_hand_side.n_rows) {
        return false;
    }

    arma::cx_mat left_vectors;
    arma::vec singular_values;
    arma::cx_mat right_vectors;
    bool ok = arma::svd_econ(
        left_vectors,
        singular_values,
        right_vectors,
        matrix,
        "both",
        "dc"
    );
    if (!ok || singular_values.n_elem == 0 || !singular_values.is_finite()) {
        return false;
    }

    const double tolerance =
        static_cast<double>((std::max)(matrix.n_rows, matrix.n_cols)) *
        singular_values[0] * std::numeric_limits<double>::epsilon();

    arma::uword retained = 0;
    for (arma::uword i = 0; i < singular_values.n_elem; ++i) {
        retained += singular_values[i] >= tolerance ? 1u : 0u;
    }

    if (retained == 0) {
        result = matrix.t() * right_hand_side;
        result.zeros();
        return true;
    }

    arma::vec inverse_singular_values(retained);
    for (arma::uword i = 0; i < retained; ++i) {
        inverse_singular_values[i] = 1.0 / singular_values[i];
    }

    const arma::cx_mat left_retained = left_vectors.cols(0, retained - 1);
    const arma::cx_mat right_retained = right_vectors.cols(0, retained - 1);
    const arma::cx_vec inverse_complex =
        arma::conv_to<arma::cx_vec>::from(inverse_singular_values);

    result = right_retained *
        (arma::diagmat(inverse_complex) *
         (left_retained.t() * right_hand_side));

    return result.is_finite();
}

#endif
