#' Prolate Spheroidal Angular Function of the First Kind, \eqn{S^{1}_{mn}(c,
#' \eta)}
#'
#' @description
#' Computes the prolate spheroidal angular function of the first kind,
#' \eqn{S_{mn}^{(1)}(c, \eta)}, and its first derivative with respect to
#' \eqn{\eta} for given order \eqn{m}, degree \eqn{n}, size parameter \eqn{c},
#' and angular coordinate \eqn{\eta}.
#'
#' This function is an R wrapper for compiled C++ code, which in turn calls the
#' underlying Fortran library (\code{prolate_swf}) for high-performance
#' numerical computation. All heavy computation is performed in compiled code
#' for speed and
#' accuracy.
#'
#' @details
#' The prolate spheroidal angular functions are solutions to the angular part
#' of the scalar Helmholtz equation in prolate spheroidal coordinates. They
#' satisfy the differential equation:
#' \deqn{(1 - \eta^2) \frac{d^2 S_{mn}}{d\eta^2} - 2\eta \frac{dS_{mn}}{d\eta}
#' + \left(\lambda_{mn} - c^2 \eta^2 + \frac{m^2}{1 - \eta^2}\right) S_{mn} = 0}
#'
#' where \eqn{\lambda_{mn}} is the separation constant (eigenvalue).
#'
#' **Domain restrictions:**
#' \itemize{
#'   \item The angular coordinate must satisfy \eqn{|\eta| \leq 1}.
#'   \item The order must satisfy \eqn{m \geq 0}.
#'   \item The degree must satisfy \eqn{n \geq m}.
#' }
#'
#' **Normalization:**
#' When \code{normalize = FALSE} (default), the functions use the
#' Meixner-Schäfke normalization, which matches the normalization of the
#' corresponding associated Legendre functions. When \code{normalize = TRUE},
#' the functions are scaled to have unity norm.
#'
#' **Implementation:**
#' This function is an R wrapper for a compiled C++ interface (\code{Smn_cpp}),
#' which itself wraps the Fortran subroutine \code{profcn} from the
#' \code{prolate_swf} library developed by Arnie Lee Van Buren and Jeffrey
#' Boisvert. The underlying algorithm uses a combination of forward and
#' backward recursion with the Bouwkamp eigenvalue method for high accuracy
#' across wide parameter ranges. The C++ layer manages memory, precision
#' selection, and data conversion between R and Fortran for robust and
#' efficient computation.
#'
#' @param m Non-negative integer. The order of the spheroidal function (\eqn{m
#' \geq 0}).
#' @param n Non-negative integer. The degree of the spheroidal function (\eqn{n
#' \geq m}).
#' @param c Numeric. The scalar size parameter.
#' @param eta Numeric vector. The angular coordinate(s) at which to evaluate the
#'   function. Must satisfy \eqn{|\eta| \leq 1}.
#' @param normalize Logical. If \code{TRUE}, the angular functions are
#' normalized to have unity norm. If \code{FALSE} (default), the
#' Meixner-Schäfke normalization is used.
#' @param precision Character. Either \code{"double"} (default) or
#' \code{"quad"}. Controls the  floating-point precision used in the underlying
#' Fortran computation.
#'   \describe{
#'     \item{\code{"double"}}{Uses standard double-precision (64-bit)
#'     arithmetic. Fastest and  sufficient for most applications.}
#'     \item{\code{"quad"}}{Uses quadruple-precision (128-bit) arithmetic for
#'     higher numerical  accuracy in challenging parameter regimes (e.g., large
#'     \eqn{m}, \eqn{n}, or near  singularities). Computation is significantly
#'     slower.}
#'   }

#' @return A list containing:
#' \describe{
#'   \item{\code{value}}{Numeric vector of function values \eqn{S_{mn}^{(1)}(c,
#'   \eta)} at each input \code{eta}.}
#'   \item{\code{derivative}}{Numeric vector of first derivatives
#'     \eqn{\frac{d}{d\eta}S_{mn}^{(1)}(c, \eta)} at each input \code{eta}.}
#' }
#'
#' @examples
#' # Single evaluation
#' Smn(m = 2, n = 3, c = 1, eta = 0.5)
#'
#' # Multiple eta values
#' Smn(m = 0, n = 2, c = 5, eta = c(-0.5, 0, 0.5))
#'
#' # With unity normalization
#' Smn(m = 1, n = 1, c = 2, eta = 0.3, normalize = TRUE)
#'
#' # Double precision (default)
#' Smn(m = 2, n = 3, c = 1, eta = 0.5, precision = "double")
#'
#' # Quad precision
#' Smn(m = 2, n = 3, c = 1, eta = 0.5, precision = "quad")
#'
#' @references
#' Van Buren, A. L. and Boisvert, J. E. "Prolate Spheroidal Wave Functions."
#' GitHub repository:
#' \url{https://github.com/MathieuandSpheroidalWaveFunctions/prolate_swf}
#'
#' Meixner, J. and Schäfke, F. W. (1954). \emph{Mathieusche Funktionen und
#' Sphäroidfunktionen}. Springer-Verlag, Berlin.
#'
#' Flammer, C. (1957). \emph{Spheroidal Wave Functions}. Stanford University
#' Press.
#'
#' NIST Digital Library of Mathematical Functions. Chapter 30: Spheroidal Wave
#' Functions. \url{https://dlmf.nist.gov/30}
#'
#' @seealso \code{\link{Rmn}} for prolate spheroidal radial functions.
#'
#' @useDynLib acousticTS, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @export
Smn <- function(m, n, c, eta, normalize = FALSE, precision = "double") {
  # Validation =================================================================
  if (!is.numeric(n) || !all(is.finite(n)) || !all(n %% 1 == 0)) {
    stop("'n' must be a real integer, or a vector of real integers.")
  }
  if (!is.numeric(m) || !all(is.finite(m)) || !all(m %% 1 == 0)) {
    stop("'m' must be a real integer, or a vector of real integers.")
  }
  if (!is.numeric(eta) || !all(is.finite(eta))) {
    stop("'eta' must be a real number, or a vector of real numbers.")
  }
  if (!is.numeric(c) || length(c) > 1 || !is.finite(c)) {
    stop("'c' must be a single, real number.")
  }
  if (!precision %in% c("double", "quad")) {
    stop("'precision' must either be 'double' (default) or 'quad'.")
  }
  # Run compiled function ======================================================
  Smn_cpp(m, n, c, eta, normalize, precision)
}

#' Prolate Spheroidal Radial Functions
#'
#' @description
#' Computes the prolate spheroidal radial functions of the first
#' (\eqn{R_{mn}^{(1)}}), second (\eqn{R_{mn}^{(2)}}), third
#' (\eqn{R_{mn}^{(3)}}), and fourth (\eqn{R_{mn}^{(4)}}) kinds and their first
#' derivatives with respect to \eqn{\xi} for given order \eqn{m}, degree
#' \eqn{n}, size parameter \eqn{c}, and radial coordinate \eqn{\xi}.
#'
#' This function is an R wrapper for compiled C++ code, which in turn calls the
#' underlying Fortran library (\code{prolate_swf}) for high-performance
#' numerical computation. All heavy computation is performed in compiled code
#' for speed and accuracy.
#'
#' @details
#' The prolate spheroidal radial functions are solutions to the radial part
#' of the scalar Helmholtz equation in prolate spheroidal coordinates. They
#' satisfy the differential equation:
#' \deqn{(\xi^2 - 1) \frac{d^2 R_{mn}}{d\xi^2} + 2\xi \frac{dR_{mn}}{d\xi} -
#' \left(\lambda_{mn} - c^2 \xi^2 + \frac{m^2}{\xi^2 - 1}\right) R_{mn} = 0}
#'
#' where \eqn{\lambda_{mn}} is the separation constant (eigenvalue).
#'
#' **Function kinds:**
#' \itemize{
#'   \item \code{kind = 1}: Radial function of the first kind
#'   \eqn{R_{mn}^{(1)}(c, \xi)}. Regular at the origin; analogous to spherical
#'   Bessel functions \eqn{j_n}.
#'   \item \code{kind = 2}: Radial function of the second kind
#'   \eqn{R_{mn}^{(2)}(c, \xi)}. Singular at the focal points (\eqn{\xi = 1});
#'   analogous to spherical Neumann functions \eqn{y_n}.
#'   \item \code{kind = 3}: Radial function of the third kind (outgoing
#'   Hankel-type) \eqn{R_{mn}^{(3)}(c, \xi) = R_{mn}^{(1)}(c, \xi) + i
#'   R_{mn}^{(2)}(c, \xi)}. Used for outgoing wave solutions.
#'   \item \code{kind = 4}: Radial function of the fourth kind (incoming
#'   Hankel-type) \eqn{R_{mn}^{(4)}(c, \xi) = R_{mn}^{(1)}(c, \xi) - i
#'   R_{mn}^{(2)}(c, \xi)}. Used for incoming wave solutions.
#' }
#'
#' **Domain restrictions:**
#' \itemize{
#'   \item The radial coordinate must satisfy \eqn{\xi \geq 1} (prolate domain).
#'   \item The order must satisfy \eqn{m \geq 0}.
#'   \item The degree must satisfy \eqn{n \geq m}.
#' }
#'
#' **Normalization:**
#' The functions use the Morse-Feshbach normalization by default, where the
#' radial functions reduce to spherical Bessel/Neumann functions as
#' \eqn{c \rightarrow 0}.
#'
#' **Implementation:**
#' This function is an R wrapper for a compiled C++ interface (\code{Rmn_cpp}),
#' which itself wraps the Fortran subroutine \code{profcn} from the
#' \code{prolate_swf} library developed by Arnie Lee Van Buren and Jeffrey
#' Boisvert. The underlying algorithm uses a combination of forward and
#' backward recursion with the Bouwkamp eigenvalue method for high accuracy
#' across wide parameter ranges. The C++ layer manages memory, precision
#' selection, and data conversion between R and Fortran for robust and
#' efficient computation.
#'
#' @param m Non-negative integer. The order of the spheroidal function (\eqn{m
#' \geq 0}).
#' @param n Non-negative integer. The degree of the spheroidal function (\eqn{n
#' \geq m}).
#' @param c Numeric. The scalar size parameter (also denoted \eqn{\gamma} in
#' some references).
#' @param xi Numeric. The radial coordinate at which to evaluate the function.
#'   Must satisfy \eqn{\xi \geq 1}.
#' @param kind Integer. Specifies which kind of radial function to compute:
#'   \code{1} (first kind), \code{2} (second kind), \code{3} (third
#'   kind/outgoing), or \code{4} (fourth kind/incoming). Default is \code{1}.
#' @param precision Character. Either \code{"double"} (default) or
#' \code{"quad"}. Controls the  floating-point precision used in the underlying
#' Fortran computation.
#'   \describe{
#'     \item{\code{"double"}}{Uses standard double-precision (64-bit)
#'     arithmetic. Fastest and  sufficient for most applications.}
#'     \item{\code{"quad"}}{Uses quadruple-precision (128-bit) arithmetic for
#'     higher numerical  accuracy in challenging parameter regimes (e.g., large
#'     \eqn{m}, \eqn{n}, \eqn{c}, or near  singularities). Computation is
#'     significantly slower.}
#'   }

#' @return A list containing:
#' \describe{
#'   \item{\code{value}}{The function value \eqn{R_{mn}^{(k)}(c, \xi)}. Real for
#'     \code{kind = 1} or \code{2}; complex for \code{kind = 3} or \code{4}.}
#'   \item{\code{derivative}}{The first derivative
#'     \eqn{\frac{d}{d\xi}R_{mn}^{(k)}(c, \xi)}. Real for \code{kind = 1} or
#'     \code{2}; complex for \code{kind = 3} or \code{4}.}
#' }
#'
#' @examples
#' # First kind radial function
#' Rmn(m = 2, n = 3, c = 1, xi = 1.5, kind = 1)
#'
#' # Second kind radial function
#' Rmn(m = 0, n = 2, c = 5, xi = 2.0, kind = 2)
#'
#' # Third kind (outgoing) radial function
#' Rmn(m = 1, n = 1, c = 2, xi = 1.2, kind = 3)
#'
#' # Fourth kind (incoming) radial function
#' Rmn(m = 1, n = 1, c = 2, xi = 1.2, kind = 4)
#'
#' # Double precision (default)
#' Rmn(m = 2, n = 3, c = 1, xi = 1.5)
#'
#' # Quad precision
#' Rmn(m = 2, n = 3, c = 1, xi = 1.5, precision = "quad")
#'
#' @references
#' Van Buren, A. L. and Boisvert, J. E. "Prolate Spheroidal Wave Functions."
#' GitHub repository:
#' \url{https://github.com/MathieuandSpheroidalWaveFunctions/prolate_swf}
#'
#' Morse, P. M. and Feshbach, H. (1953). \emph{Methods of Theoretical Physics}.
#' McGraw-Hill, New York. Chapter 21.
#'
#' Flammer, C. (1957). \emph{Spheroidal Wave Functions}. Stanford University
#' Press.
#'
#' NIST Digital Library of Mathematical Functions. Chapter 30: Spheroidal Wave
#' Functions. \url{https://dlmf.nist.gov/30}
#'
#' @seealso \code{\link{Smn}} for prolate spheroidal angular functions.
#'
#' @useDynLib acousticTS, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @export
Rmn <- function(m, n, c, xi, kind = 1, precision = "double") {
  # Validation =================================================================
  if (!is.numeric(n) || !all(is.finite(n)) || !all(n %% 1 == 0)) {
    stop("'n' must be a real integer, or a vector of real integers.")
  }
  if (!is.numeric(m) || !all(is.finite(m)) || !all(m %% 1 == 0)) {
    stop("'m' must be a real integer, or a vector of real integers.")
  }
  if (!is.numeric(c) || length(c) > 1 || !is.finite(c)) {
    stop("'c' must be a single, real number.")
  }
  if (!is.numeric(xi) || length(xi) > 1 || !is.finite(xi)) {
    stop("'xi' must be a single, real number.")
  }
  if (!precision %in% c("double", "quad")) {
    stop("'precision' must either be 'double' (default) or 'quad'.")
  }
  # Run compiled function ======================================================
  Rmn_cpp(m, n, c, xi, kind)
}
