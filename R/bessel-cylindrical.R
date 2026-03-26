################################################################################
# Cyndrical Bessel functions
################################################################################
################################################################################
#' Cylindrical Bessel function of the first kind, \eqn{J_\nu(z)}, and its
#' respective derivatives
#'
#' @description
#' Computes the cylindrical Bessel function of the first kind (\eqn{J_\nu(z)})
#' and its k-th derivatives.
#'
#' @details
#' The cylindrical Bessel function of the first kind satisfies Bessel's
#' differential equation:
#'
#' \deqn{z^2 \frac{d^2 J_\nu}{dz^2} + z \frac{dJ_\nu}{dz} + (z^2 - \nu^2) J
#' _\nu = 0}
#'
#' **Supported argument types:**
#' \itemize{
#'   \item **Purely real arguments** (\eqn{z = x}, where \eqn{x \in
#'   \mathbb{R}}):
#'         Fully supported for both positive and negative values.
#'   \item **Purely imaginary arguments** (\eqn{z = iy}, where \eqn{y \in
#'   \mathbb{R}}):
#'         Computed using the identity \eqn{J_\nu(iy) = e^{i\pi\nu/2} I_\nu(y)}
#'         where \eqn{I_\nu} is the modified Bessel function of the first kind.
#'   \item **General complex arguments** (\eqn{z = x + iy}, where \eqn{x \neq 0}
#'         and \eqn{y \neq 0}): **Not supported**.
#' }
#'
#' **Special cases:**
#' \itemize{
#'   \item \eqn{J_\nu(0) = 1} if \eqn{\nu = 0}, otherwise \eqn{J_\nu(0) = 0}.
#'   \item For negative real arguments with integer order \eqn{n}:
#'         \eqn{J_n(-x) = (-1)^n J_n(x)}.
#'   \item For negative real arguments with non-integer order \eqn{\nu}:
#'         \eqn{J_\nu(-x) = e^{i\pi\nu} J_\nu(x)} (complex result).
#' }
#'
#' **Derivatives:**
#' \itemize{
#'   \item First derivative: \eqn{J'_\nu(z) = J_{\nu-1}(z) - \frac{\nu}{z}
#'   J_\nu(z)}
#'   \item Second derivative: \eqn{J''_\nu(z) = \frac{1}{4}\left[J_{\nu-2}(z) -
#'   2J_\nu(z) + J_{\nu+2}(z)\right]}
#'   \item k-th derivative (DLMF 10.6.1):
#'         \deqn{
#'           \frac{d^k}{dz^k} J_\nu(z) = \frac{1}{2^k} \sum_{j=0}^{k} (-1)^j
#'            \binom{k}{j} J_{\nu - k + 2j}(z)
#'         }
#' }
#'
#' @param l Numeric. The order (\eqn{\nu}) of the Bessel function. Must be
#' purely real; complex orders are not supported.
#' @param n Numeric or complex. The argument (\eqn{z}) at which to evaluate the
#'   function. Supports purely real or purely imaginary values. General complex
#'   arguments (\eqn{x + iy} with \eqn{x \neq 0} and \eqn{y \neq 0}) are not
#'   supported.
#' @param k Non-negative integer. The order of the derivative for \code{jcdk}.
#'
#' @return
#' A complex vector containing:
#' \itemize{
#'   \item \code{jc}: \eqn{J_\nu(z)}
#'   \item \code{jcdk(..., k = 1)}: \eqn{J'_\nu(z)} (first derivative)
#'   \item \code{jcdk(..., k = 2)}: \eqn{J''_\nu(z)} (second derivative)
#'   \item \code{jcdk}: \eqn{J_\nu^{(k)}(z)} (k-th derivative)
#' }
#'
#' @examples
#' # Real argument, integer order
#' jc(0, 1)
#' jc(1, 2.5)
#'
#' # Real argument, fractional order
#' jc(0.5, 3)
#'
#' # Negative real argument (integer order gives real result)
#' jc(2, -1.5)
#'
#' # Negative real argument (fractional order gives complex result)
#' jc(0.5, -2)
#'
#' # Purely imaginary argument
#' jc(1, 1i)
#' jc(2, 3i)
#'
#' # First derivative
#' jcdk(1, 2, 1)
#'
#' # Second derivative
#' jcdk(1, 2, 2)
#'
#' @references
#' Abramowitz, M. and Stegun, I.A. (Eds.). (1964). \emph{Handbook of
#' Mathematical  Functions with Formulas, Graphs, and Mathematical Tables}.
#' National Bureau of Standards, Applied Mathematics Series 55. Chapter 9.
#'
#' NIST Digital Library of Mathematical Functions. \url{https://dlmf.nist.gov/}
#' \itemize{
#'   \item Bessel's equation: \url{https://dlmf.nist.gov/10.2}
#'   \item Negative argument identity (\eqn{J_\nu(-z)}): Eq. 10.4.1 at
#'         \url{https://dlmf.nist.gov/10.4}
#'   \item Imaginary argument identity (\eqn{J_\nu(iz)}): Eq. 10.27.6 at
#'         \url{https://dlmf.nist.gov/10.27}
#' }
#'
#' @seealso
#' \code{\link{yc}} for Bessel functions of the second kind,
#' \code{\link{hc}} for Hankel functions (third kind),
#' \code{\link{js}} for spherical Bessel functions of the first kind.
#'
#' @keywords bessel
#' @rdname jc
#' @export
jc <- function(l, n) {
  jc_cpp(n, l)
}
#' @rdname jc
#' @keywords internal
#' @noRd
jcd <- function(l, n) {
  # Internal helper function ===================================================
  jc_deriv_cpp(n, l, 1)
}
#' @rdname jc
#' @keywords internal
#' @noRd
jcdd <- function(l, n) {
  # Internal helper function ===================================================
  jc_deriv_cpp(n, l, 2)
}
#' @rdname jc
#' @export
jcdk <- function(l, n, k) {
  jc_deriv_cpp(n, l, k)
}

################################################################################
#' Cylindrical Bessel function of the second kind, \eqn{Y_\nu(x)}, and its
#' respective derivatives
#'
#' @description
#' Computes the cylindrical Bessel function of the second kind
#' (\eqn{Y_\nu(z)}), also known as the Neumann function or Weber function, and
#' its derivatives through \code{ycdk()}.
#'
#' @details
#' The cylindrical Bessel function of the second kind satisfies the same
#' differential equation as \eqn{J_\nu(z)}:
#' \deqn{z^2 \frac{d^2 Y_\nu}{dz^2} + z \frac{dY_\nu}{dz} + (z^2 - \nu^2)
#' Y_\nu = 0}
#'
#' but represents the linearly independent second solution.
#'
#' **Supported argument types:**
#' \itemize{
#'   \item **Purely real arguments** (\eqn{z = x}, where
#'   \eqn{x \in \mathbb{R}}):
#'         Fully supported for both positive and negative values.
#'   \item **Purely imaginary arguments** (\eqn{z = iy}, where
#'   \eqn{y \in \mathbb{R}}):
#'         Computed using the identity
#'         \eqn{Y_\nu(iy) = i e^{-i\pi\nu/2} I_\nu(y) - \frac{2}{\pi}
#'         e^{i\pi\nu/2} K_\nu(y)}
#'         where \eqn{I_\nu} and \eqn{K_\nu} are modified Bessel functions.
#'   \item **General complex arguments** (\eqn{z = x + iy}, where \eqn{x \neq 0}
#'         and \eqn{y \neq 0}): **Not supported**.
#' }
#'
#' **Special cases:**
#' \itemize{
#'   \item \eqn{Y_\nu(0) = -\infty} (singularity at the origin).
#'   \item For negative real arguments:
#'         \eqn{Y_\nu(-x) = \cos(\pi\nu) Y_\nu(x) + \sin(\pi\nu) J_\nu(x)}.
#'   \item For integer order \eqn{n}: \eqn{Y_n(-x) = (-1)^n Y_n(x)}.
#'   \item k-th derivative (DLMF 10.6.1):
#'         \deqn{
#'           \frac{d^k}{dz^k} Y_\nu(z) = \frac{1}{2^k} \sum_{j=0}^{k} (-1)^j
#'           \binom{k}{j}
#'            Y_{\nu - k + 2j}(z)
#'         }
#' }
#'
#' **Derivative:**
#' \deqn{Y'_\nu(z) = Y_{\nu-1}(z) - \frac{\nu}{z} Y_\nu(z)}
#'
#' @param l Numeric. The order (\eqn{\nu}) of the Bessel function. Must be
#' purely real; complex orders are not supported.
#' @param n Numeric or complex. The argument (\eqn{z}) at which to evaluate the
#'   function. Supports purely real or purely imaginary values. General complex
#'   arguments (\eqn{x + iy} with \eqn{x \neq 0} and \eqn{y \neq 0}) are not
#'   supported.
#' @param k Non-negative integer. The order of the derivative for \code{ycdk}.
#'
#' @return
#' A complex vector containing:
#' \itemize{
#'   \item \code{yc}: \eqn{Y_\nu(z)}
#'   \item \code{ycdk(..., k = 1)}: \eqn{Y'_\nu(z)} (first derivative)
#'   \item \code{ycdk(..., k = 2)}: \eqn{Y''_\nu(z)} (second derivative)
#'   \item \code{ycdk}: \eqn{Y_\nu^{(k)}(z)} (k-th derivative)
#' }
#'
#' @examples
#' # Real argument, integer order
#' yc(0, 1)
#' yc(1, 2.5)
#'
#' # Real argument, fractional order
#' yc(0.5, 3)
#'
#' # Negative real argument
#' yc(2, -1.5)
#'
#' # Purely imaginary argument
#' yc(1, 1i)
#' yc(2, 3i)
#'
#' # Singularity at origin
#' yc(0, 0)  # Returns -Inf
#'
#' # First derivative
#' ycdk(1, 2, 1)
#'
#' @references
#' Abramowitz, M. and Stegun, I.A. (Eds.). (1964).
#' \emph{Handbook of Mathematical Functions with Formulas, Graphs, and
#' Mathematical Tables}. National Bureau of Standards, Applied Mathematics
#' Series 55. Chapter 9.
#'
#' NIST Digital Library of Mathematical Functions. \url{https://dlmf.nist.gov/}
#' \itemize{
#'   \item Bessel's equation: \url{https://dlmf.nist.gov/10.2}
#'   \item Negative argument identity (\eqn{Y_\nu(-z)}): Eq. 10.4.1 at
#'         \url{https://dlmf.nist.gov/10.4}
#'   \item Imaginary argument identity (\eqn{Y_\nu(iz)}): Eq. 10.27.8 at
#'         \url{https://dlmf.nist.gov/10.27}
#' }
#'
#' @seealso
#' \code{\link{jc}} for Bessel functions of the first kind,
#' \code{\link{hc}} for Hankel functions (third kind),
#' \code{\link{ys}} for spherical Bessel functions of the second kind.
#'
#' @keywords bessel
#' @rdname yc
#' @export
yc <- function(l, n) {
  yc_cpp(n, l)
}
#' @rdname yc
#' @keywords internal
#' @noRd
ycd <- function(l, n) {
  # Internal helper function ===================================================
  yc_deriv_cpp(n, l, 1)
}
#' @rdname yc
#' @keywords internal
#' @noRd
ycdd <- function(l, n) {
  # Internal helper function ===================================================
  yc_deriv_cpp(n, l, 2)
}
#' @rdname yc
#' @export
ycdk <- function(l, n, k) {
  yc_deriv_cpp(n, l, k)
}

################################################################################
#' Cylindrical Bessel function of the third kind (Hankel), \eqn{H_\nu(x)}, and
#' its respective derivatives
#'
#' @description
#' Computes the cylindrical Hankel function of the first kind
#' (\eqn{H^{(1)}_\nu(z)}) and its derivatives through \code{hcdk()}.
#'
#' @details
#' The Hankel function of the first kind is defined as:
#' \deqn{H^{(1)}_\nu(z) = J_\nu(z) + i Y_\nu(z)}
#'
#' where \eqn{J_\nu(z)} is the Bessel function of the first kind and
#' \eqn{Y_\nu(z)} is the Bessel function of the second kind.
#'
#' **Supported argument types:**
#' Since \eqn{H^{(1)}_\nu(z)} is computed from \eqn{J_\nu(z)} and
#' \eqn{Y_\nu(z)}, the same restrictions apply:
#' \itemize{
#'   \item **Purely real arguments** (\eqn{z = x}): Fully supported.
#'   \item **Purely imaginary arguments** (\eqn{z = iy}): Supported.
#'   \item **General complex arguments** (\eqn{z = x + iy}): **Not supported**.
#' }
#'
#' **Derivatives:**
#' \itemize{
#'   \item First derivative:
#'         \eqn{\frac{d}{dz}H^{(1)}_\nu(z) = \frac{\nu}{z} H^{(1)}_\nu(z) -
#'         H^{(1)}_{\nu+1}(z)}
#'   \item Second derivative:
#'         \eqn{\frac{d^2}{dz^2}H^{(1)}_\nu(z) = H^{(1)}_{\nu-2}(z) -
#'         \frac{2\nu-1}{z} H^{(1)}_{\nu-1}(z) + \frac{\nu^2+\nu}{z^2}
#'         H^{(1)}_\nu(z)}
#'   \item k-th derivative (DLMF 10.6.1):
#'         \eqn{\frac{d^k}{dz^k}H^{(1)}_\nu(z) = \frac{1}{2^k} \sum_{j=0}^{k}
#'         (-1)^j \binom{k}{j} H^{(1)}_{\nu-k+2j}(z)}
#' }
#'
#' @param l Numeric. The order (\eqn{\nu}) of the Hankel function. Must be
#' purely real; complex orders are not supported.
#' @param n Numeric or complex. The argument (\eqn{z}) at which to evaluate the
#'   function. Supports purely real or purely imaginary values. General complex
#'   arguments are not supported.
#' @param k Non-negative integer. The order of the derivative for \code{hcdk}.
#'
#' @return
#' A complex vector containing:
#' \itemize{
#'   \item \code{hc}: \eqn{H^{(1)}_\nu(z)}
#'   \item \code{hcdk(..., k = 1)}: \eqn{\frac{d}{dz}H^{(1)}_\nu(z)}
#'   (first derivative)
#'   \item \code{hcdk(..., k = 2)}: \eqn{\frac{d^2}{dz^2}H^{(1)}_\nu(z)}
#'   (second derivative)
#'   \item \code{hcdk}: \eqn{\frac{d^k}{dz^k}H^{(1)}_\nu(z)} (k-th derivative)
#' }
#'
#' @examples
#' # Hankel function
#' hc(0, 1)
#' hc(1, 2.5)
#'
#' # Fractional order
#' hc(0.5, 3)
#'
#' # Purely imaginary argument
#' hc(1, 1i)
#'
#' # First derivative
#' hcdk(1, 2, 1)
#'
#' # Second derivative
#' hcdk(1, 2, 2)
#'
#' # k-th derivative
#' hcdk(1, 2, 3)  # Third derivative
#'
#' @references
#' Abramowitz, M. and Stegun, I.A. (Eds.). (1964). \emph{Handbook of
#' Mathematical Functions with Formulas, Graphs, and Mathematical Tables}.
#' National Bureau of Standards, Applied Mathematics Series 55. Chapter 9.
#'
#' NIST Digital Library of Mathematical Functions.
#' \url{https://dlmf.nist.gov/}
#' \itemize{
#'   \item Hankel function definition: \url{https://dlmf.nist.gov/10.2}
#'   \item k-th derivative formula: Eq. 10.6.1 at
#'   \url{https://dlmf.nist.gov/10.6}
#' }
#'
#' @seealso
#' \code{\link{jc}} for Bessel functions of the first kind,
#' \code{\link{yc}} for Bessel functions of the second kind,
#' \code{\link{hs}} for spherical Hankel functions.
#'
#' @keywords bessel
#' @rdname hc
#' @export
hc <- function(l, n) {
  hc_cpp(n, l)
}
#' @rdname hc
#' @keywords internal
#' @noRd
hcd <- function(l, n) {
  # Internal helper function ===================================================
  hc_deriv_cpp(n, l, 1)
}
#' @rdname hc
#' @keywords internal
#' @noRd
hcdd <- function(l, n) {
  # Internal helper function ===================================================
  hc_deriv_cpp(n, l, 2)
}
#' @rdname hc
#' @export
hcdk <- function(l, n, k) {
  hc_deriv_cpp(n, l, k)
}
