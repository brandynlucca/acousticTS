################################################################################
# Spherical Bessel functions
################################################################################

################################################################################
#' Spherical Bessel Function of the First Kind and Its Derivatives
#'
#' @description
#' Computes the spherical Bessel function of the first kind (\eqn{j_l(z)}) and
#' its first (\code{jsd}) and second (\code{jsdd}) derivatives.
#'
#' @details
#' The spherical Bessel function of the first kind is related to the cylindrical
#' Bessel function by:
#' \deqn{j_l(z) = \sqrt{\frac{\pi}{2z}} J_{l+1/2}(z)}
#'
#' where \eqn{J_\nu(z)} is the cylindrical Bessel function of the first kind.
#'
#' The spherical Bessel functions satisfy the differential equation:
#' \deqn{z^2 \frac{d^2 j_l}{dz^2} + 2z \frac{dj_l}{dz} + [z^2 - l(l+1)] j_l = 0}
#'
#' **Special cases:**
#' \itemize{
#'   \item \eqn{j_l(0) = 0} for all \eqn{l}.
#'   \item \eqn{j_0(z) = \frac{\sin(z)}{z}}
#'   \item \eqn{j_1(z) = \frac{\sin(z)}{z^2} - \frac{\cos(z)}{z}}
#' }
#'
#' **Derivatives:**
#' \itemize{
#'   \item First derivative: \eqn{j'_l(z) = j_{l-1}(z) - \frac{l+1}{z} j_l(z)}
#'   \item Second derivative: \eqn{j''_l(z) = \frac{(l+1)(l+2) - z^2}{z^2}
#'   j_l(z) - \frac{2}{z} j_{l-1}(z)}
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
#'   \item \code{js}: \eqn{j_l(z)}
#'   \item \code{jsd}: \eqn{j'_l(z)} (first derivative)
#'   \item \code{jsdd}: \eqn{j''_l(z)} (second derivative)
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
#' jsd(1, 2)
#'
#' # Second derivative
#' jsdd(1, 2)
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
  if (!is.numeric(l) || !is.numeric(n)) {
    stop("Inputs must be numeric vectors.")
  }
  js_cpp(l, n)
  # Internal function ==========================================================
  # js_internal <- function(l, n) {
  #   if (n == 0) {
  #     0
  #   } else {
  #     jc(l + 0.5, n) * sqrt(pi / (2 * n))
  #   }
  # }
  # js_vec <- Vectorize(js_internal)
  # # Return based on input class ================================================
  # switch(class(n)[1],
  #        numeric = js_vec(l, n),
  #        matrix = apply(n, 2, FUN = function(x) {
  #          js_vec(l, x)
  #        })
  # )
}

js_old <- function(l, n) {
  # Internal function ==========================================================
  js_internal <- function(l, n) {
    if (n == 0) {
      0
    } else {
      jc_old(l + 0.5, n) * sqrt(pi / (2 * n))
    }
  }
  js_vec <- Vectorize(js_internal)
  # Return based on input class ================================================
  switch(class(n)[1],
         numeric = js_vec(l, n),
         matrix = apply(n, 2, FUN = function(x) {
           js_vec(l, x)
         })
  )
}

#' @rdname js
#' @export
jsd <- function(l, n) {
  js_deriv_cpp(l, n, 1)
}

jsd_old <- function(l, n) {
  ifelse(n == 0,
         0,
         js_old(l - 1, n) - (l + 1) / n * js_old(l, n)
  )
}


#' @rdname js
#' @export
jsdd <- function(l, n) {
  js_deriv_cpp(l, n, 2)
  # (1 / (n^2) *
  #    ((l + 1) * (l + 2) - n^2) * js(l, n) -
  #    2 / n * js(l - 1, n))
}

jsdd_old <- function(l, n) {
  (1 / (n^2) *
     ((l + 1) * (l + 2) - n^2) * js_old(l, n) -
     2 / n * js_old(l - 1, n))
}

jsdk <- function(l, n, k) {
  js_deriv_cpp(l, n, k)
}
################################################################################
#' Spherical Bessel Function of the Second Kind and Its Derivatives
#'
#' @description
#' Computes the spherical Bessel function of the second kind (\eqn{y_l(z)}),
#' also known as the spherical Neumann function, and its first (\code{ysd}) and
#' second (\code{ysdd}) derivatives.
#'
#' @details
#' The spherical Bessel function of the second kind is related to the
#' cylindrical Bessel function by:
#' \deqn{y_l(z) = \sqrt{\frac{\pi}{2z}} Y_{l+1/2}(z)}
#'
#' where \eqn{Y_\nu(z)} is the cylindrical Bessel function of the second kind.
#'
#' The spherical Bessel functions satisfy the same differential equation as
#' \eqn{j_l(z)}:
#' \deqn{z^2 \frac{d^2 y_l}{dz^2} + 2z \frac{dy_l}{dz} + [z^2 - l(l+1)] y_l = 0}
#'
#' **Special cases:**
#' \itemize{
#'   \item \eqn{y_l(0) = -\infty} (singularity at the origin).
#'   \item \eqn{y_0(z) = -\frac{\cos(z)}{z}}
#'   \item \eqn{y_1(z) = -\frac{\cos(z)}{z^2} - \frac{\sin(z)}{z}}
#' }
#'
#' **Derivatives:**
#' \itemize{
#'   \item First derivative: \eqn{y'_l(z) = \frac{l}{z} y_l(z) - y_{l+1}(z)}
#'   \item Second derivative: \eqn{y''_l(z) = -\frac{l}{z^2} y_l(z) +
#'   \frac{l}{z} y'_l(z) - y'_{l+1}(z)}
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
#'   \item \code{ys}: \eqn{y_l(z)}
#'   \item \code{ysd}: \eqn{y'_l(z)} (first derivative)
#'   \item \code{ysdd}: \eqn{y''_l(z)} (second derivative)
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
#' ysd(1, 2)
#'
#' # Second derivative
#' ysdd(1, 2)
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
  if (!is.numeric(l) || !is.numeric(n)) {
    stop("Inputs must be numeric vectors.")
  }
  ys_cpp(l, n)
  # Internal function ==========================================================
  # ys_internal <- function(l, n) {
  #   if (n > 0) {
  #     yc(l + 0.5, n) * sqrt(pi / (2 * n))
  #   } else if (n < 0) {
  #     (yc(l + 0.5, -n) * sqrt(pi / (2 * -n)))
  #   } else if (n == 0) {
  #     -Inf
  #   }
  # }
  # ys_vec <- Vectorize(ys_internal)
  # # Return based on input class ================================================
  # switch(class(n)[1],
  #        numeric = ys_vec(l, n),
  #        matrix = apply(n, 2, FUN = function(x) {
  #          ys_vec(l, x)
  #        })
  # )
}

ys_old <- function(l, n) {
  # Internal function ==========================================================
  ys_internal <- function(l, n) {
    if (n > 0) {
      yc_old(l + 0.5, n) * sqrt(pi / (2 * n))
    } else if (n < 0) {
      (yc_old(l + 0.5, -n) * sqrt(pi / (2 * -n)))
    } else if (n == 0) {
      -Inf
    }
  }
  ys_vec <- Vectorize(ys_internal)
  # Return based on input class ================================================
  switch(class(n)[1],
         numeric = ys_vec(l, n),
         matrix = apply(n, 2, FUN = function(x) {
           ys_vec(l, x)
         })
  )
}

#' @rdname ys
#' @export
ysd <- function(l, n) {
  ys_deriv_cpp(l, n, 1)
  # ifelse(n == 0,
  #        -Inf,
  #        l / n * ys(l, n) - ys(l + 1, n)
  # )
}

ysd_old <- function(l, n) {
  ifelse(n == 0,
         -Inf,
         l / n * ys_old(l, n) - ys_old(l + 1, n)
  )
}

#' @rdname ys
#' @export
ysdd <- function(l, n) {
  ys_deriv_cpp(l, n, 2)
  # (-l / (n^2) * ys(l, n) +
  #    l / n * (l / n * ys(l, n) - ys(l + 1, n)) -
  #    ((l + 1) / n * ys(l + 1, n) - ys(l + 2, n)))
}

ysdd_old <- function(l, n) {
  (-l / (n^2) * ys_old(l, n) +
     l / n * (l / n * ys_old(l, n) - ys_old(l + 1, n)) -
     ((l + 1) / n * ys_old(l + 1, n) - ys_old(l + 2, n)))
}
################################################################################
#' Spherical Hankel Function of the First Kind and Its Derivative
#'
#' @description
#' Computes the spherical Hankel function of the first kind (\eqn{h^{(1)}_l(z)})
#' and its first derivative (\code{hsd}).
#'
#' @details
#' The spherical Hankel function of the first kind is defined as:
#' \deqn{h^{(1)}_l(z) = j_l(z) + i y_l(z)}
#'
#' where \eqn{j_l(z)} is the spherical Bessel function of the first kind and
#' \eqn{y_l(z)} is the spherical Bessel function of the second kind.
#'
#' It is related to the cylindrical Hankel function by:
#' \deqn{h^{(1)}_l(z) = \sqrt{\frac{\pi}{2z}} H^{(1)}_{l+1/2}(z)}
#'
#' The spherical Hankel functions are used extensively in scattering theory
#' to represent outgoing spherical waves.
#'
#' **Derivative:**
#' \deqn{\frac{d}{dz}h^{(1)}_l(z) = \frac{l}{z} h^{(1)}_l(z) - h^{(1)}_{l+1}(z)}
#'
#' @param l Numeric. The order of the spherical Hankel function. Can be integer
#'   or fractional.
#' @param n Numeric. The argument (\eqn{z}) at which to evaluate the function.
#'
#' @return
#' A complex vector containing:
#' \itemize{
#'   \item \code{hs}: \eqn{h^{(1)}_l(z)}
#'   \item \code{hsd}: \eqn{\frac{d}{dz}h^{(1)}_l(z)} (first derivative)
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
#' hsd(1, 2)
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
  # sqrt(pi / (2 * n)) * hc(l + 0.5, n)
  hs_cpp(l, n)
}

hs_old <- function(l, n) {
  sqrt(pi / (2 * n)) * hc_old(l + 0.5, n)
}

#' @rdname hs
#' @export
hsd <- function(l, n) {
  # -hs(l + 1, n) + (l / n) * hs(l, n)
  hs_deriv_cpp(l, n, 1)
}

hsd_old <- function(l, n) {
  -hs_old(l + 1, n) + (l / n) * hs_old(l, n)
}
