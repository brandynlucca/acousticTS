################################################################################
################################################################################
# Mathematical functions
################################################################################
################################################################################
#' Along-matrix summing function
#' @param rpos Position vector
#' @param iterations Number of iterations
#' @rdname along_sum
#' @keywords internal
#' @noRd
along_sum <- function(rpos, iterations) {
  if (iterations < 2) {
    stop("along_sum requires at least 2 columns to sum", call. = FALSE)
  }
  lhs <- rpos[, 1:(iterations - 1), drop = FALSE]
  rhs <- rpos[, 2:iterations, drop = FALSE]
  lhs + rhs
}
################################################################################
#' Numerical integration via adaptive quadrature of complex values
#' @param integral Input integration function that is indexed so that it can be
#' used within `apply`.
#' @param x Indexing argument for multi-row objects
#' @param y Indexing argument for multi-column objects
#' @rdname contour_integrate
#' @keywords internal
#' @noRd
contour_integrate <- function(integral, x, y) {
  complex(
    real = stats::integrate(function(s) Re(integral(s, x, y)),
      lower = 0, upper = 1
    )$value,
    imaginary = stats::integrate(function(s) Im(integral(s, x, y)),
      lower = 0, upper = 1
    )$value
  )
}
################################################################################
#' Wrapper function incorporating phase deviation into contour integration
#' @param integral Input integration function that is indexed so that it can be
#' used within `apply`
#' @param x Indexing argument for multi-row objects
#' @param y Indexing argument for multi-column objects
#' @param n_iterations Number of phase deviations to average and summarize
#' @param integral Integral function used for numerical integration via adaptive
#' quadrature
#' @param phase_sd Phase standard deviation
#' @rdname phase_integrate
#' @keywords internal
#' @noRd
phase_integrate <- function(x, y, n_iterations, integral, phase_sd) {
  rng <- stats::rnorm(n_iterations, 0, 1)
  phase <- exp(1i * rng * phase_sd)
  contour_integrate(integral, x, y) * phase
}
################################################################################
#' Convert angular measurements from radians to degrees
#' @param x A real value in radians
#' @examples
#' orientation <- pi / 2 # radians
#' degrees(orientation) # this should return a value equal to 90 degrees
#' @return
#' Angle in degrees
#' @rdname degrees
#' @export
degrees <- function(x) x * 180.0 / pi
################################################################################
#' Convert angular measurements from degrees to radians.
#' @param x A real value in degrees
#' @examples
#' orientation <- 90 # degrees
#' radians(orientation) # this should return a value equal to pi / 2 radians
#' @return
#' Angle in radians.
#' @rdname radians
#' @export
radians <- function(x) x * pi / 180.0
################################################################################
#' Calculates the Euclidean norm across each row of a given matrix.
#' @param x A matrix with numeric, real values.
#' @usage
#' vecnorm( x )
#' @examples
#' values <- matrix(c(1, 2, 3), ncol = 3)
#' vecnorm(values) # should yield 3.741657
#' @return
#' Calculates the Euclidean norm of a vector.
#' @rdname vecnorm
#' @export
vecnorm <- function(x) sqrt(rowSums(x * x))
################################################################################
#' Compute the Neumann factor \eqn{\nu_{n}}
#'
#' The Neumann factor, denoted \eqn{\nu_{n}}, is a simple multiplicative
#' constant commonly used in spherical or spheroidal function theory. It
#' accounts for the symmetry of cosine terms or the duplication of even-order
#' contributions in integrals.
#'
#' Formally, the Neumann factor is defined as:
#' \deqn{
#'   \eta_n = \begin{cases}
#'     1, & n = 0, \\
#'     2, & n > 0.
#'   \end{cases}
#' }
#'
#' @param x An integer iterator.
#'
#' @details
#' This factor frequently appears in the normalization of spherical Bessel
#' functions, Legendre expansions, and spheroidal wave functions. It is not
#' related to the "von Neumann ordinals" used in set theory.
#'
#' @return A numeric vector of the same length as `x`, containing values of
#' \eqn{\eta_x}, each equal to 1 or 2.
#'
#' @examples
#' neumann(0) # should return 1
#' neumann(1) # should return 2
#' neumann(2) # should return 2
#'
#' # Vectorized use:
#' neumann(0:4)
#'
#' @export
neumann <- function(x) {
  # Validation =================================================================
  if (any(x < 0 | !x %% 1 == 0)) {
    stop(
      paste0(
        "Value 'x', ", x, ", must be an integer greater than or equal to 0."
      )
    )
  }
  # Compute ====================================================================
  ifelse(x == 0, 1, 2)
}
################################################################################
#' Gauss–Legendre nodes and weights
#'
#' Compute Gauss–Legendre quadrature nodes and weights on an interval
#' \eqn{[a,~b]}.
#'
#' @param n Number of quadrature nodes (n >= 1).
#' @param a Left endpoint of the integration interval.
#' @param b Right endpoint of the integration interval (b > a).
#'
#' @return A list with components:
#' \describe{
#'   \item{nodes}{Quadrature abscissae \eqn{x_i} in \eqn{[a,~b]}.}
#'   \item{weights}{Quadrature weights \eqn{w_i} such that
#'     \eqn{\int_a^b f(x)\,dx \approx \sum_{i=1}^n w_i\,f(x_i).}}
#' }
#'
#' @details
#' Gauss–Legendre quadrature provides exact integration for polynomials of
#' degree up to \eqn{2n-1} using n nodes and weights chosen as the roots of the
#' Legendre polynomial \eqn{P_n(x)} on the canonical interval \eqn{[-1,1]}. For
#' a general interval \eqn{[a,b]} the mapping
#' \deqn{x = \tfrac{a+b}{2} + \tfrac{b-a}{2}\,t,\quad t\in[-1,1],}
#' transforms canonical nodes \eqn{t_i} to \eqn{x_i} and scales weights by
#' \deqn{w_i = \tfrac{b-a}{2}\,w_i^{(0)},}
#' where \eqn{w_i^{(0)}} are the standard weights on \eqn{[-1,1]}.
#'
#' This wrapper performs basic argument validation and calls the C++ routine
#' to obtain nodes and weights with high accuracy for moderate \code{n}.
#'
#' @references
#' Davis, P. J., & Rabinowitz, P. (2007). Methods of Numerical Integration
#' (2nd ed.).
#'
#' @export
gauss_legendre <- function(n, a = -1, b = 1) {
  # Validation =================================================================
  if (!is.numeric(n) || length(n) != 1 || n < 1 && n%%1 != 0) {
    stop("n must be a positive integer")
  }
  if (!is.numeric(a) || !is.numeric(b) || length(a) != 1 || length(b) != 1) {
    stop("a and b must be numeric scalars")
  }
  if (b <= a) stop("b must be greater than a")
  # Get nodes and their weights ================================================
  gauss_legendre_cpp(n = n, a = a, b = b)
}
