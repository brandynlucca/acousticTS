################################################################################
# Harmonics functions
################################################################################
################################################################################
#' Legendre Polynomial of the First Kind, \eqn{P_\nu(x)}
#'
#' @description
#' Computes the Legendre polynomial of the first kind, \eqn{P_\nu(x)}, for
#' real order \eqn{\nu} (integer or fractional) and real argument \eqn{x}.
#'
#' @details
#' The Legendre polynomial of the first kind satisfies the differential
#' equation:
#' \deqn{(1 - x^2) \frac{d^2 P_\nu}{dx^2} - 2x \frac{dP_\nu}{dx} +
#' \nu(\nu + 1) P_\nu = 0}
#'
#' For integer order \eqn{n}, the function uses the recurrence relation:
#' \deqn{P_0(x) = 1}
#' \deqn{P_1(x) = x}
#' \deqn{P_n(x) = \frac{(2n-1) x P_{n-1}(x) - (n-1) P_{n-2}(x)}{n}}
#'
#' For fractional order \eqn{\nu}, the function uses:
#' \itemize{
#'   \item For \eqn{|x| \leq 1}: The Ferrers function via the hypergeometric
#'   series
#'         \eqn{P_\nu(x) = {}_2F_1(-\nu, \nu+1; 1; \frac{1-x}{2})}
#'   \item For \eqn{|x| > 1}: A numerical contour integral representation
#' }
#'
#' @param n Numeric vector. Degree (order) of the Legendre polynomial. Can be
#'   integer or fractional (e.g., 0, 1, 2.5, 3.7). Can also be negative.
#' @param x Numeric vector. Real argument(s) at which to evaluate the
#' polynomial. Valid for all real values including \eqn{|x| > 1}.
#'
#' @return
#' A numeric matrix of dimension \code{length(n)} by \code{length(x)}, where
#' element \code{[i, j]} contains \eqn{P_{n_i}(x_j)}.
#'
#' @examples
#' # Single values
#' Pn(1, 1)
#'
#' # Multiple orders, single argument
#' Pn(c(1, 2, 3), 1)
#'
#' # Single order, multiple arguments
#' Pn(1, c(1, 2, 3))
#'
#' # Multiple orders and arguments (returns a matrix)
#' Pn(c(1, 2, 3), c(1, 2, 3))
#'
#' # Fractional orders
#' Pn(c(0.5, 2.2, 3), c(1, 2, 3))
#'
#' # Negative arguments
#' Pn(c(0.5, 2, -3), c(1, -2.5, 3))
#'
#' @note
#' This function calls underlying \eqn{C++} code via \code{Rcpp} for
#' computational efficiency and  to support different cases for both order and
#' argument that are not readily available in \code{R}.
#'
#' @references
#' Abramowitz, M. and Stegun, I. A. (1972). \emph{Handbook of Mathematical
#' Functions with Formulas, Graphs, and Mathematical Tables}. Dover
#' Publications.
#' Chapter 8: Legendre Functions.
#'
#' NIST Digital Library of Mathematical Functions.
#' \url{https://dlmf.nist.gov/14}
#'
#' @seealso
#' \code{\link{Pndk}} for the k<sup>th</sup> derivative of the Legendre
#' polynomial of the first kind,
#' \code{\link{Qn}} for Legendre functions of the second kind.
#'
#' @useDynLib acousticTS, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @keywords legendre
#' @rdname Pn
#' @export
Pn <- function(n, x) {
  # Validation =================================================================
  if (is.complex(n) | is.complex(x)) {
    stop("'n' and 'x' must be real numbers.")
  }
  if (any(!is.finite(n))) {
    stop("'n' must be finite.")
  }
  if (any(!is.finite(x))) {
    stop("'x' must be a real, finite numeric.")
  }
  # Compute Legendre polynomial ================================================
  Pn_cpp(n, x)
}

#' Derivative of the Legendre Polynomial of the First Kind
#'
#' @description
#' Computes the \eqn{k}-th derivative of the Legendre polynomial of the first
#' kind, \eqn{\frac{d^k}{dx^k} P_\nu(x)}, with respect to the argument \eqn{x}.
#'
#' @details
#' **For integer order \eqn{n}:**
#'
#' The derivative is computed using the relationship with associated Legendre
#' polynomials:
#' \deqn{\frac{d^m P_n(x)}{dx^m} = \frac{1}{(1-x^2)^{m/2}} P_n^m(x)}
#'
#' where \eqn{P_n^m(x)} is the associated Legendre polynomial (the \code{Boost}
#' \code{C++} uses Condon-Shortley phase convention).
#'
#' At the endpoints \eqn{x = \pm 1}, the known closed-form expressions are used:
#' \deqn{P_n^{(k)}(1) = \frac{1}{2^k k!} \prod_{j=0}^{k-1} (n-j)(n+1+j)}
#' \deqn{P_n^{(k)}(-1) = (-1)^{n+k} P_n^{(k)}(1)}
#'
#' If \eqn{k > n} for integer \eqn{n}, the result is 0 (derivative of a
#' polynomial of degree \eqn{n} taken more than \eqn{n} times).
#'
#' **For fractional order \eqn{\nu}:**
#'
#' Derivatives are computed using central finite differences:
#' \deqn{\frac{dP_\nu}{dx} \approx \frac{P_\nu(x+h) - P_\nu(x-h)}{2h}}
#'
#' Higher-order derivatives use the generalized finite difference stencil.
#' Note that accuracy may be limited for fractional orders due to the
#' numerical integration underlying \code{\link{Pn}}.
#'
#' @param n Numeric vector. Degree (order) of the Legendre polynomial. Can be
#'   integer or fractional.
#' @param x Numeric vector. Real argument(s) at which to evaluate the
#'   derivative.
#' @param k Integer. Order of the derivative (\eqn{k \geq 0}). Default is 1
#'   for the first derivative.
#'
#' @return
#' A numeric matrix of dimension \code{length(n)} by \code{length(x)}, where
#' element \code{[i, j]} contains \eqn{\frac{d^k}{dx^k} P_{n_i}(x_j)}.
#'
#' @examples
#' # First derivative of P_2(x) at x = 0.5

#' # P_2(x) = (3x^2 - 1)/2, so P'_2(x) = 3x, P'_2(0.5) = 1.5
#' Pndk(2, 0.5, 1)
#'
#' # Second derivative of P_2(x)
#' # P''_2(x) = 3
#' Pndk(2, 0.5, 2)
#'
#' # Multiple orders
#' Pndk(c(1, 2, 3), 0.5, 1)
#'
#' # First derivative at multiple points
#' Pndk(2, c(-0.5, 0, 0.5), 1)
#'
#' # Fractional order (uses finite differences)
#' Pndk(0.5, 0.5, 1)
#'
#' # Derivative at endpoint
#' # P'_n(1) = n(n+1)/2
#' Pndk(3, 1, 1)  # Should be 3*4/2 = 6
#'
#' @references
#' Abramowitz, M. and Stegun, I. A. (1972). \emph{Handbook of Mathematical
#' Functions}. Dover Publications. Section 8.5: Associated Legendre Functions.
#'
#' NIST Digital Library of Mathematical Functions.
#' \url{https://dlmf.nist.gov/14.10}
#'
#' @note
#' For fractional orders, the finite difference approximation may have reduced
#' accuracy (typically 4-6 significant digits) compared to integer orders.
#' A warning is issued for higher-order derivatives (\eqn{k > 1}) of
#' fractional orders.
#'
#' This function calls underlying \eqn{C++} code via \code{Rcpp} for
#' computational efficiency and  to support different cases for both order and
#' argument that are not readily available in \code{R}.
#'
#' @seealso
#' \code{\link{Pn}} for the Legendre polynomial of the first kind.
#'
#' @useDynLib acousticTS, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @keywords legendre derivative
#' @rdname Pndk
#' @export
Pndk <- function(n, x, k = 1L) {
  # Validation =================================================================
  if (is.complex(n) | is.complex(x)) {
    stop("'n' and 'x' must be real numbers.")
  }
  if (any(!is.finite(n))) {
    stop("'n' must be finite.")
  }
  if (any(!is.finite(x))) {
    stop("'x' must be a real, finite numeric.")
  }
  if (!is.numeric(k) || length(k) != 1 || k < 0 || k != floor(k)) {
    stop("'k' must be a non-negative integer.")
  }
  # Compute derivative =========================================================
  Pn_deriv_cpp(as.numeric(n), as.numeric(x), as.integer(k))
}

################################################################################
#' Legendre Function of the Second Kind, \eqn{Q_\nu(x)}
#'
#' @description
#' Computes the Legendre function of the second kind, \eqn{Q_\nu(x)}, for
#' real order \eqn{\nu} (integer or fractional) and real argument \eqn{x}.
#' Returns complex values when \eqn{|x| > 1}.
#'
#' @details
#' The Legendre function of the second kind satisfies the same differential
#' equation as \eqn{P_\nu(x)}:
#' \deqn{(1 - x^2) \frac{d^2 Q_\nu}{dx^2} - 2x \frac{dQ_\nu}{dx} + \nu(\nu + 1)
#' Q_\nu = 0}
#'
#' but represents the linearly independent second solution.
#'
#' **For \eqn{|x| < 1} (real result):**
#' \itemize{
#'   \item For integer order \eqn{n}: Uses the \code{C++} \code{Boost} library
#'   implementation
#'   \item For fractional order \eqn{\nu}: Uses the Ferrers identity
#'         \deqn{Q_\nu(x) = \frac{\pi}{2 \sin(\pi \nu)} \left[ \cos(\pi \nu)
#'         P_\nu(x) - P_\nu(-x) \right]}
#' }
#'
#' **For \eqn{x = \pm 1}:**
#' Returns infinity (singularity).
#'
#' **For \eqn{|x| > 1} (complex result):**
#' \itemize{
#'   \item Real part: Computed via the integral representation
#'         \deqn{\text{Re}[Q_\nu(x)] = \int_0^\infty \frac{dt}{(x +
#'         \sqrt{x^2-1} \cosh t)^{\nu+1}}}
#'   \item Imaginary part: \deqn{\text{Im}[Q_\nu(x)] = -\frac{\pi}{2} P_\nu(x)}
#' }
#'
#' @param n Numeric vector. Degree (order) of the Legendre function. Can be
#'   integer or fractional (e.g., 0, 1, 2.5, 3.7).
#' @param x Numeric vector. Real argument(s) at which to evaluate the function.
#'   Valid for all real values. Note that \eqn{x = \pm 1} are singularities.
#'
#' @return
#' A complex matrix of dimension \code{length(n)} by \code{length(x)}, where
#' element \code{[i, j]} contains \eqn{Q_{n_i}(x_j)}. For \eqn{|x| < 1}, the
#' imaginary part is zero.
#'
#' @examples
#' # Single values
#' Qn(1, 0.5)
#'
#' # Multiple orders, single argument
#' Qn(c(1, 2, 3), 0.5)
#'
#' # Single order, multiple arguments
#' Qn(1, c(-0.5, 0, 0.5))
#'
#' # Multiple orders and arguments (returns a matrix)
#' Qn(c(1, 2, 3), c(0.25, 0.5, 0.75))
#'
#' # Fractional orders
#' Qn(c(0.5, 1.5, 2.7), c(0.5, 0.9))
#'
#' # Arguments |x| > 1 return complex values
#' Qn(c(1, 2, 3.5), c(2, 3.5))
#'
#' # Singularity at x = 1
#' Qn(1, 1)
#'
#' @note
#' This function calls underlying \eqn{C++} code via \code{Rcpp} for
#' computational efficiency and  to support different cases for both order and
#' argument that are not readily available in \code{R}.
#'
#' @references
#' Abramowitz, M. and Stegun, I. A. (1972). \emph{Handbook of Mathematical
#' Functions with Formulas, Graphs, and Mathematical Tables}. Dover
#' Publications.
#' Chapter 8: Legendre Functions.
#'
#' NIST Digital Library of Mathematical Functions.
#' \url{https://dlmf.nist.gov/14}
#'
#' @seealso
#' \code{\link{Qndk}} for the k<sup>th</sup> derivative of the Legendre
#' polynomial of the second kind,
#' \code{\link{Pn}} for Legendre functions of the first kind.
#'
#' @useDynLib acousticTS, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @keywords legendre
#' @rdname Qn
#' @export
Qn <- function(n, x) {
  # Validation =================================================================
  if (is.complex(n) | is.complex(x)) {
    stop("'n' and 'x' must be real numbers.")
  }
  if (any(!is.finite(n))) {
    stop("'n' must be finite.")
  }
  if (any(!is.finite(x))) {
    stop("'x' must be a real, finite numeric.")
  }
  # Compute Legendre function ==================================================
  Qn_cpp(n, x)
}

#' Derivative of the Legendre Function of the Second Kind
#'
#' @description
#' Computes the \eqn{k}-th derivative of the Legendre function of the second
#' kind, \eqn{\frac{d^k}{dx^k} Q_\nu(x)}, with respect to the argument \eqn{x}.
#'
#' @details
#' Derivatives are computed using central finite differences:
#' \deqn{\frac{dQ_\nu}{dx} \approx \frac{Q_\nu(x+h) - Q_\nu(x-h)}{2h}}
#'
#' Higher-order derivatives use the generalized finite difference stencil with
#' binomial coefficients.
#'
#' For \eqn{|x| > 1}, the result is complex since \eqn{Q_\nu(x)} is complex
#' in that domain.
#'
#' **Singularities:** The function \eqn{Q_\nu(x)} has logarithmic singularities
#' at \eqn{x = \pm 1}. Derivatives near these points will have reduced accuracy
#' or may be undefined.
#'
#' @param n Numeric vector. Degree (order) of the Legendre function.
#' @param x Numeric vector. Real argument(s) at which to evaluate the
#'   derivative. Avoid \eqn{x = \pm 1}.
#' @param k Integer. Order of the derivative (\eqn{k \geq 0}). Default is 1.
#'
#' @return
#' A complex matrix of dimension \code{length(n)} by \code{length(x)}, where
#' element \code{[i, j]} contains \eqn{\frac{d^k}{dx^k} Q_{n_i}(x_j)}.
#'
#' @examples
#' # First derivative of Q_1(x) at x = 0.5
#' Qndk(1, 0.5, 1)
#'
#' # Compare with numerical derivative
#' h <- 1e-6
#' (Qn(1, 0.5 + h) - Qn(1, 0.5 - h)) / (2 * h)
#'
#' # Second derivative
#' Qndk(2, 0.5, 2)
#'
#' # Multiple orders
#' Qndk(c(1, 2, 3), 0.5, 1)
#'
#' # Complex result for |x| > 1
#' Qndk(1, 2.0, 1)
#'
#' @references
#' Abramowitz, M. and Stegun, I. A. (1972). \emph{Handbook of Mathematical
#' Functions}. Dover Publications. Chapter 8: Legendre Functions.
#'
#' @seealso
#' \code{\link{Qn}} for Legendre functions of the second kind.
#'
#' @note
#' Derivatives are computed via finite differences with step size
#' \eqn{h = 10^{-6}}. Accuracy is typically 4-6 significant digits for
#' first derivatives, less for higher orders.
#'
#' This function calls underlying \eqn{C++} code via \code{Rcpp} for
#' computational efficiency and  to support different cases for both order and
#' argument that are not readily available in \code{R}.
#'
#' @useDynLib acousticTS, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @keywords legendre derivative
#' @rdname Qndk
#' @export
Qndk <- function(n, x, k = 1L) {
  # Validation =================================================================
  if (is.complex(n) | is.complex(x)) {
    stop("'n' and 'x' must be real numbers.")
  }
  if (any(!is.finite(n))) {
    stop("'n' must be finite.")
  }
  if (any(!is.finite(x))) {
    stop("'x' must be a real, finite numeric.")
  }
  if (!is.numeric(k) || length(k) != 1 || k < 0 || k != floor(k)) {
    stop("'k' must be a non-negative integer.")
  }
  if (any(abs(abs(x) - 1) < 1e-10)) {
    warning("Derivatives near x = +/-1 may be inaccurate due to singularity.")
  }
  # Compute derivative =========================================================
  Qn_deriv_cpp(as.numeric(n), as.numeric(x), as.integer(k))
}
