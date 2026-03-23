################################################################################
# Spherical Bessel functions
################################################################################

################################################################################
#' Spherical Bessel function of the first kind, \eqn{j_\nu(z)}, and its
#' respective derivatives
#'
#' @description
#' Computes the spherical Bessel function of the first kind (\eqn{j_\nu(z)}) 
#' and its k-th derivative.
#'
#' @details
#' The spherical Bessel function of the first kind is related to the cylindrical
#' Bessel function by:
#' \deqn{j_\nu(z) = \sqrt{\frac{\pi}{2z}} J_{\nu+1/2}(z)}
#'
#' where \eqn{J_\nu(z)} is the cylindrical Bessel function of the first kind.
#'
#' The spherical Bessel functions satisfy the differential equation:
#' \deqn{
#'  z^2 \frac{d^2 j_\nu}{dz^2} + 
#'  2z \frac{dj_\nu}{dz} + [z^2 - l(l+1)] j_\nu = 0
#' }
#'
#' **Special cases:**
#' \itemize{
#'   \item \eqn{j_\nu(0) = 0} for all \eqn{\nu}.
#'   \item \eqn{j_0(z) = \frac{\sin(z)}{z}}
#'   \item \eqn{j_1(z) = \frac{\sin(z)}{z^2} - \frac{\cos(z)}{z}}
#' }
#'
#' **Derivatives:**
#' \itemize{
#'   \item First derivative: \eqn{
#'    j'_\nu(z) = j_{\nu-1}(z) - \frac{\nu+1}{z} j_\nu(z)
#'   }
#'   \item Second derivative: \eqn{j''_\nu(z) = \frac{(\nu+1)(\nu+2) - z^2}{z^2}
#'   j_\nu(z) - \frac{2}{z} j_{\nu-1}(z)}
#' }
#'
#' @param l Numeric. The order of the spherical Bessel function. Can be integer
#'   or fractional.
#' @param n Numeric or matrix. The argument (\eqn{z}) at which to evaluate the
#'   function. If a matrix is provided, the function is applied column-wise.
#'
#' @return
#' A numeric vector or matrix (matching the input structure) containing:
#' \itemize{
#'   \item \code{js}: \eqn{j_\nu(z)}
#'   \item \code{jsdk}: \eqn{j^{(k)'}_l(z)} (k-th derivative)
#' }
#'
#' @examples
#' # Spherical Bessel function
#' js(0, 1)
#' js(1, 2.5)
#'
#' # Fractional order
#' js(0.5, 3)
#'
#' # Vector input
#' js(0, c(1, 2, 3))
#'
#' # Matrix input (applied column-wise)
#' js(1, matrix(1:6, nrow = 2))
#'
#' # First derivative
#' jsdk(1, 2, 1)
#'
#' # Second derivative
#' jsdk(1, 2, 2)
#'
#' @references
#' Abramowitz, M. and Stegun, I.A. (Eds.). (1964). \emph{Handbook of
#' Mathematical Functions with Formulas, Graphs, and Mathematical Tables}.
#' National Bureau of Standards, Applied Mathematics Series 55. Chapter 10.
#'
#' NIST Digital Library of Mathematical Functions. \url{https://dlmf.nist.gov/}
#' \itemize{
#'   \item Spherical Bessel functions: \url{https://dlmf.nist.gov/10.47}
#'   \item Relation to cylindrical Bessel functions: Eq. 10.47.3 at
#'         \url{https://dlmf.nist.gov/10.47}
#' }
#'
#' @seealso
#' \code{\link{jc}} for cylindrical Bessel functions of the first kind,
#' \code{\link{ys}} for spherical Bessel functions of the second kind,
#' \code{\link{hs}} for spherical Hankel functions.
#'
#' @keywords bessel
#' @rdname js
#' @export
js <- function(l, n) {
  # Function check =============================================================
  if (!is.numeric(l) || !(is.numeric(n) || is.complex(n))) {
    stop("Inputs must be numeric or complex vectors.")
  }
  if (is.complex(n)) {
    js_complex_cpp(as.integer(l), as.complex(n))
  } else {
    js_cpp(l, n)
  }
}
#' @rdname js
#' @keywords internal
#' @noRd
jsd <- function(l, n) {
  if (is.complex(n)) {
    js_complex_deriv_cpp(as.integer(l), as.complex(n), 1)
  } else {
    js_deriv_cpp(l, n, 1)
  }
}
#' @rdname js
#' @keywords internal
#' @noRd
jsdd <- function(l, n) {
  if (is.complex(n)) {
    js_complex_deriv_cpp(as.integer(l), as.complex(n), 2)
  } else {
    js_deriv_cpp(l, n, 2)
  }
}
#' @rdname js
#' @export
jsdk <- function(l, n, k) {
  if (is.complex(n)) {
    js_complex_deriv_cpp(as.integer(l), as.complex(n), k)
  } else {
    js_deriv_cpp(l, n, k)
  }
}
################################################################################
#' Spherical Bessel function of the second kind, \eqn{y_\nu(z)}, and its
#' respective derivatives
#'
#' @description
#' Computes the spherical Bessel function of the second kind (\eqn{y_\nu(z)}),
#' also known as the spherical Neumann function, and its k-th derivatives.
#'
#' @details
#' The spherical Bessel function of the second kind is related to the
#' cylindrical Bessel function by:
#' 
#' \deqn{y_\nu(z) = \sqrt{\frac{\pi}{2z}} Y_{\nu+1/2}(z)}
#'
#' where \eqn{Y_\nu(z)} is the cylindrical Bessel function of the second kind.
#'
#' The spherical Bessel functions satisfy the same differential equation as
#' \eqn{j_\nu(z)}:
#' \deqn{
#'  z^2 \frac{d^2 y_\nu}{dz^2} + 2z \frac{dy_\nu}{dz} + 
#'  [z^2 - \nu(\nu+1)] y_\nu = 0
#' }
#'
#' **Special cases:**
#' \itemize{
#'   \item \eqn{y_\nu(0) = -\infty} (singularity at the origin).
#'   \item \eqn{y_0(z) = -\frac{\cos(z)}{z}}
#'   \item \eqn{y_1(z) = -\frac{\cos(z)}{z^2} - \frac{\sin(z)}{z}}
#' }
#'
#' **Derivatives:**
#' \itemize{
#'   \item First derivative: \eqn{
#'    y'_\nu(z) = \frac{\nu}{z} y_\nu(z) - y_{\nu+1}(z)
#'  }
#'   \item Second derivative: \eqn{y''_\nu(z) = -\frac{\nu}{z^2} y_\nu(z) +
#'   \frac{\nu}{z} y'_\nu(z) - y'_{\nu+1}(z)}
#' }
#'
#' @param l Numeric. The order of the spherical Bessel function. Can be integer
#'   or fractional.
#' @param n Numeric or matrix. The argument (\eqn{z}) at which to evaluate the
#'   function. If a matrix is provided, the function is applied column-wise.
#'
#' @return
#' A numeric vector or matrix (matching the input structure) containing:
#' \itemize{
#'   \item \code{ys}: \eqn{y_\nu(z)}
#'   \item \code{ysd}: \eqn{y^{(k)'}_l(z)} (k-th derivative)
#' }
#'
#' @examples
#' # Spherical Bessel function of the second kind
#' ys(0, 1)
#' ys(1, 2.5)
#'
#' # Fractional order
#' ys(0.5, 3)
#'
#' # Vector input
#' ys(0, c(1, 2, 3))
#'
#' # Singularity at origin
#' ys(0, 0)  # Returns -Inf
#'
#' # First derivative
#' ysdk(1, 2, 1)
#'
#' # Second derivative
#' ysdk(1, 2, 2)
#'
#' # 3rd derivative
#' ysdk(1, 1, 3)
#' 
#' # 4th derivative
#' ysdk(1, 1, 4)
#'
#' @references
#' Abramowitz, M. and Stegun, I.A. (Eds.). (1964). \emph{Handbook of
#' Mathematical Functions with Formulas, Graphs, and Mathematical Tables}.
#' National Bureau of Standards, Applied Mathematics Series 55. Chapter 10.
#'
#' NIST Digital Library of Mathematical Functions. \url{https://dlmf.nist.gov/}
#' \itemize{
#'   \item Spherical Bessel functions: \url{https://dlmf.nist.gov/10.47}
#'   \item Relation to cylindrical Bessel functions: Eq. 10.47.4 at
#'         \url{https://dlmf.nist.gov/10.47}
#' }
#'
#' @seealso
#' \code{\link{yc}} for cylindrical Bessel functions of the second kind,
#' \code{\link{js}} for spherical Bessel functions of the first kind,
#' \code{\link{hs}} for spherical Hankel functions.
#'
#' @keywords bessel
#' @rdname ys
#' @export
ys <- function(l, n) {
  # Function check =============================================================
  if (!is.numeric(l) || !(is.numeric(n) || is.complex(n))) {
    stop("Inputs must be numeric or complex vectors.")
  }
  if (is.complex(n)) {
    ys_complex_cpp(as.integer(l), as.complex(n))
  } else {
    ys_cpp(l, n)
  }
}
#' @rdname ys
#' @keywords internal
#' @noRd
ysd <- function(l, n) {
  if (is.complex(n)) {
    ys_complex_deriv_cpp(as.integer(l), as.complex(n), 1)
  } else {
    ys_deriv_cpp(l, n, 1)
  }
}
#' @rdname ys
#' @keywords internal
#' @noRd
ysdd <- function(l, n) {
  if (is.complex(n)) {
    ys_complex_deriv_cpp(as.integer(l), as.complex(n), 2)
  } else {
    ys_deriv_cpp(l, n, 2)
  }
}
#' @rdname ys
#' @export
ysdk <- function(l, n, k) {
  if (is.complex(n)) {
    ys_complex_deriv_cpp(as.integer(l), as.complex(n), k)
  } else {
    ys_deriv_cpp(l, n, k)
  }
}
################################################################################
#' Spherical Bessel function of the third kind (Hankel), \eqn{h_\nu(x)}, and its
#' respective derivatives
#'
#' @description
#' Computes the spherical Hankel function of the first kind 
#' (\eqn{h^{(1)}_\nu(z)}) and its first-th derivative.
#'
#' @details
#' The spherical Hankel function of the first kind is defined as:
#' \deqn{h^{(1)}_\nu(z) = j_\nu(z) + i y_\nu(z)}
#'
#' where \eqn{j_\nu(z)} is the spherical Bessel function of the first kind and
#' \eqn{y_\nu(z)} is the spherical Bessel function of the second kind.
#'
#' It is related to the cylindrical Hankel function by:
#' \deqn{h^{(1)}_\nu(z) = \sqrt{\frac{\pi}{2z}} H^{(1)}_{\nu+1/2}(z)}
#'
#' The spherical Hankel functions are used extensively in scattering theory
#' to represent outgoing spherical waves.
#'
#' **Derivative:**
#' \deqn{
#'  \frac{d}{dz}h^{(1)}_\nu(z) = \frac{\nu}{z} h^{(1)}_\nu(z) 
#'  - h^{(1)}_{\nu+1}(z)
#' }
#'
#' @param l Numeric. The order of the spherical Hankel function. Can be integer
#'   or fractional.
#' @param n Numeric. The argument (\eqn{z}) at which to evaluate the function.
#'
#' @return
#' A complex vector containing:
#' \itemize{
#'   \item \code{hs}: \eqn{h^{(1)}_\nu(z)}
#'   \item \code{hsdK}: \eqn{\frac{d}{dz^k}h^{(1)}_l(z)} (k-th derivative)
#' }
#'
#' @examples
#' # Spherical Hankel function
#' hs(0, 1)
#' hs(1, 2.5)
#'
#' # Fractional order
#' hs(0.5, 3)
#'
#' # Vector input
#' hs(0, c(1, 2, 3))
#'
#' # First derivative
#' hsdk(1, 2, 1)
#'
#' @references
#' Abramowitz, M. and Stegun, I.A. (Eds.). (1964). \emph{Handbook of
#' Mathematical Functions with Formulas, Graphs, and Mathematical Tables}.
#' National Bureau of Standards, Applied Mathematics Series 55. Chapter 10.
#'
#' NIST Digital Library of Mathematical Functions. \url{https://dlmf.nist.gov/}
#' \itemize{
#'   \item Spherical Bessel functions: \url{https://dlmf.nist.gov/10.47}
#'   \item Relation to cylindrical Hankel functions: Eq. 10.47.5 at
#'         \url{https://dlmf.nist.gov/10.47}
#' }
#'
#' @seealso
#' \code{\link{hc}} for cylindrical Hankel functions,
#' \code{\link{js}} for spherical Bessel functions of the first kind,
#' \code{\link{ys}} for spherical Bessel functions of the second kind.
#'
#' @keywords bessel
#' @rdname hs
#' @export
hs <- function(l, n) {
  if (is.complex(n)) {
    hs_complex_cpp(as.integer(l), as.complex(n))
  } else {
    hs_cpp(l, n)
  }
}
#' @rdname hs
#' @keywords internal
#' @noRd
hsd <- function(l, n) {
  if (is.complex(n)) {
    hs_complex_deriv_cpp(as.integer(l), as.complex(n), 1)
  } else {
    hs_deriv_cpp(l, n, 1)
  }
}
#' @rdname hs
#' @keywords internal
#' @noRd
hsdd <- function(l, n) {
  if (is.complex(n)) {
    hs_complex_deriv_cpp(as.integer(l), as.complex(n), 2)
  } else {
    hs_deriv_cpp(l, n, 2)
  }
}
#' @rdname hs
#' @export
hsdk <- function(l, n, k) {
  if (is.complex(n)) {
    hs_complex_deriv_cpp(as.integer(l), as.complex(n), k)
  } else {
    hs_deriv_cpp(l, n, k)
  }
}
