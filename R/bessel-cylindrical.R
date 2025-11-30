################################################################################
# Cyndrical Bessel functions
################################################################################
################################################################################
#' Cylindrical Bessel function of the first kind, \eqn{J_\nu(x)} and its
#' respective derivatives
#'
#' @description
#' Computes the cylindrical Bessel function of the first kind (\eqn{J_\nu(z)})
#' and its first (\code{jcd}) and second (\code{jcdd}) derivatives.
#'
#' @details
#' The cylindrical Bessel function of the first kind satisfies Bessel's
#' differential equation:
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
#' }
#'
#' @param l Numeric. The order (\eqn{\nu}) of the Bessel function. Must be
#' purely real; complex orders are not supported.
#' @param n Numeric or complex. The argument (\eqn{z}) at which to evaluate the
#'   function. Supports purely real or purely imaginary values. General complex
#'   arguments (\eqn{x + iy} with \eqn{x \neq 0} and \eqn{y \neq 0}) are not
#'   supported.
#'
#' @return
#' A complex vector containing:
#' \itemize{
#'   \item \code{jc}: \eqn{J_\nu(z)}
#'   \item \code{jcd}: \eqn{J'_\nu(z)} (first derivative)
#'   \item \code{jcdd}: \eqn{J''_\nu(z)} (second derivative)
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
#' jcd(1, 2)
#'
#' # Second derivative
#' jcdd(1, 2)
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
  # besselJ(l, n)
  jc_cpp(n, l)
  # jc_cpp(n, l)
}

jc_old <- function(l, n) {
  besselJ(n, l)

}


#' @rdname jcdk
#' @export
jcdk <- function(l, n, k) {
  jc_deriv_cpp(l, n, k)
}

#' @rdname jc
#' @export
jcd <- function(l, n) {
  jc_internal <- function(l, n) {
    if (n == 0) {
      if (l == 1) {
        0.5
      } else {
        0.0
      }
    } else {
      jc(l - 1, n) - (l / n) * jc(l, n)
    }
  }
  jc_vec <- Vectorize(jc_internal)
  # Return based on input class ================================================
  switch(class(n)[1],
         numeric = jc_vec(l, n),
         matrix = apply(n, 2, FUN = function(x) {
           jc_vec(l, x)
         })
  )
}

#' @rdname jc
#' @export
jcdd <- function(l, n) {
  0.25 * (jc(l - 2, n) - 2 * jc(l, n) + jc(l + 2, n))
}

################################################################################
#' Cylindrical Bessel (Neumann) Function of the Second Kind and Its Derivatives
#'
#' @description
#' Computes the cylindrical Bessel function of the second kind
#' (\eqn{Y_\nu(z)}), also known as the Neumann function or Weber function, and
#' its first derivative (\code{ycd}).
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
#'
#' @return
#' A complex vector containing:
#' \itemize{
#'   \item \code{yc}: \eqn{Y_\nu(z)}
#'   \item \code{ycd}: \eqn{Y'_\nu(z)} (first derivative)
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
#' ycd(1, 2)
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
  # besselY(n, l)
}

yc_old <- function(l, n) besselY(n, l)

#' @rdname yc
#' @export
ycdk <- function(l, n, k) {
  yc_deriv_cpp(l, n, k)
}

#' @rdname yc
#' @export
ycd <- function(l, n) {
  yc_internal <- function(l, n) {
    if (n == 0.0) {
      if (l == 1) {
        0.5
      } else {
        0.0
      }
    } else {
      yc(l - 1, n) - (l / n) * yc(l, n)
    }
  }
  yc_vec <- Vectorize(yc_internal)
  # Return based on input class ================================================
  switch(class(n)[1],
         numeric = yc_vec(l, n),
         matrix = apply(n, 2, FUN = function(x) {
           yc_vec(l, x)
         })
  )
}

################################################################################
#' Cylindrical Hankel Function of the First Kind and Its Derivatives
#'
#' @description
#' Computes the cylindrical Hankel function of the first kind
#' (\eqn{H^{(1)}_\nu(z)}) and its first (\code{hcd}), second (\code{hcdd}),
#' and k-th (\code{hcdk}) derivatives.
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
#'         \frac{2\nu-1}{z} H^{(1)}_{\nu-1}(z) + \frac{\nu^2+\nu}{z^2} H^{(1)}_\nu(z)}
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
#'   \item \code{hcd}: \eqn{\frac{d}{dz}H^{(1)}_\nu(z)} (first derivative)
#'   \item \code{hcdd}: \eqn{\frac{d^2}{dz^2}H^{(1)}_\nu(z)} (second derivative)
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
#' hcd(1, 2)
#'
#' # Second derivative
#' hcdd(1, 2)
#'
#' # k-th derivative
#' hcdk(1, 2, 3)  # Third derivative
#'
#' @references
#' Abramowitz, M. and Stegun, I.A. (Eds.). (1964). \emph{Handbook of Mathematical
#' Functions with Formulas, Graphs, and Mathematical Tables}. National Bureau
#' of Standards, Applied Mathematics Series 55. Chapter 9.
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
  # jc(l, n) + 1i * yc(l, n)
}

hc_old <- function(l, n) {
  jc_old(l, n) + 1i * yc_old(l, n)
}

#' @rdname hc
#' @export
hcd <- function(l, n) {
  (l * hc(l, n) / n - hc(l + 1, n))
}

#' @rdname hc
#' @export
hcdd <- function(l, n) {
  hcdd_internal <- function(l, n) {
    if (n == 0) {
      stop(
        paste0(
          "Second derivative for the cylindrical Hankel function of the ",
          "first kind is not defined at n = 0."
        )
      )
    }
    hc(l - 2, n) -
      ((2 * l - 1) / n) * hc(l - 1, n) +
      ((l^2 + l) / (n^2)) * hc(l, n)
  }
  hcdd_vec <- Vectorize(hcdd_internal)
  switch(class(n)[1],
         numeric = hcdd_vec(l, n),
         matrix = apply(n, 2, function(x) hcdd_vec(l, x))
  )
}

#' @rdname hc
#' @export
hcdk <- function(l, n, k) {
  # hc_deriv_cp(n, l, k)
  hcdk_internal <- function(l, n, k) {
    sum <- 0
    for (j in 0:k) {
      sum <- sum + (-1)^j * choose(k, j) * hc(l - k + 2 * j, n)
    }
    sum / (2^k)
  }
  hcdk_vec <- Vectorize(function(l, n) hcdk_internal(l, n, k))
  switch(class(n)[1],
         numeric = hcdk_vec(l, n),
         matrix = apply(n, 2, function(x) hcdk_vec(l, x))
  )
}
