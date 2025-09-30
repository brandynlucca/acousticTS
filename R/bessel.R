################################################################################
# DIFFERENTIAL EQUATION SOLUTION FUNCTIONS
################################################################################
################################################################################
# Cyndrical Bessel functions
################################################################################
################################################################################
#' Cylindrical Bessel function of the first kind and its respective derivatives
#' @description
#' Calculate the cylindrical Bessel function of the first kind (jc) and both its
#' first (jcd) and second (jcdd) derivatives.
#' @param l A real integer or fractional order
#' @param n A real argument
#' @return
#' Calculates the cylindrical Bessel function of the first kind (\eqn{J_v}) and
#' its respective derivatives.
#' @references
#' Amos D.E., "AMOS, A Portable Package for Bessel Functions of a Complex
#' Argument and Nonnegative Order", http://netlib.org/amos
#'
#' @keywords bessel1
#' @rdname jc
#' @export
jc <- function(l, n) {
  besselJ(n, l)
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
#' Cylindrical Bessel (Neumann) function of the second kind and its derivatives
#' @description
#' Calculate the cylindrical Bessel function of the first kind (yc) and its
#' first (ycd) derivatives.
#' @param l A real integer or fractional order
#' @param n A real argument
#' @return
#' Calculates the cylindrical Bessel function of the first kind (\eqn{Y_v}) and
#' its respective derivatives.
#' @references
#' Amos D.E., "AMOS, A Portable Package for Bessel Functions of a Complex
#' Argument and Nonnegative Order", http://netlib.org/amos
#'
#' @keywords bessel2
#' @rdname yc
#' @export
yc <- function(l, n) besselY(n, l)
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
#' Cylindrical Bessel (Hankel) function of the third kind and its derivatives
#' @param l A real integer or fractional order
#' @param n A real argument
#' @return
#' Calculates the cylindrical Bessel function of the first kind (\eqn{H_v}) and
#' its respective derivatives.
#' @references
#' Amos D.E., "AMOS, A Portable Package for Bessel Functions of a Complex
#' Argument and Nonnegative Order", http://netlib.org/amos
#'
#' @keywords bessel3
#' @rdname hc
#' @export
hc <- function(l, n) {
  jc(l, n) + 1i * yc(l, n)
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
  # Return based on input class ================================================
  switch(class(n)[1],
    numeric = hcdd_vec(l, n),
    matrix = apply(n, 2, function(x) hcdd_vec(l, x))
  )
}
#' k-th derivative of the cylindrical Hankel function of the first kind
#'
#' @param l A real integer or fractional order
#' @param n A real argument
#' @param k Order of the derivative (non-negative integer)
#' @return
#' Calculates the k-th derivative of the cylindrical Hankel function of the
#' first kind (\eqn{H_v^{(1)}}) with respect to its order.
#' @references
#' DLMF 10.6.1, https://dlmf.nist.gov/10.6.E1
#' @keywords bessel3
#' @rdname hcdk
#' @export
hcdk <- function(l, n, k) {
  hcdk_internal <- function(l, n, k) {
    sum <- 0
    for (j in 0:k) {
      sum <- sum + (-1)^j * choose(k, j) * hc(l - k + 2 * j, n)
    }
    sum / (2^k)
  }
  hcdk_vec <- Vectorize(function(l, n) hcdk_internal(l, n, k))
  # Return based on input class
  switch(class(n)[1],
    numeric = hcdk_vec(l, n),
    matrix = apply(n, 2, function(x) hcdk_vec(l, x))
  )
}
################################################################################
# Spherical Bessel functions
################################################################################
################################################################################
#' Spherical Bessel function of the first kind and its respective derivatives
#' @param l A real integer or fractional order
#' @param n A real argument
#' @return
#' Calculates the spherical Bessel function of the first kind (js) and both its
#' first (jsd) and second (jsdd) derivatives.
#' @references
#' Amos D.E., "AMOS, A Portable Package for Bessel Functions of a Complex
#' Argument and Nonnegative Order", http://netlib.org/amos
#'
#' @keywords bessel1
#' @rdname js
#' @export
js <- function(l, n) {
  # Function check =============================================================
  if (!is.numeric(l) || !is.numeric(n)) {
    stop("Inputs must be numeric vectors.")
  }
  # Internal function ==========================================================
  js_internal <- function(l, n) {
    if (n == 0) {
      0
    } else {
      jc(l + 0.5, n) * sqrt(pi / (2 * n))
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
  js(l - 1, n) - (l + 1) / n * js(l, n)
}
#' @rdname js
#' @export
jsdd <- function(l, n) {
  (1 / (n^2) *
    ((l + 1) * (l + 2) - n^2) * js(l, n) -
    2 / n * js(l - 1, n))
}
################################################################################
#' Spherical Bessel function of the second kind and its respective derivative
#' @param l A real integer or fractional order
#' @param n A real argument
#' @return
#' Calculates the spherical Bessel function of the second kind (ys) and its
#' first (ysd) and second (ysdd) derivatives.
#' @references
#' Amos D.E., "AMOS, A Portable Package for Bessel Functions of a Complex
#' Argument and Nonnegative Order", http://netlib.org/amos
#'
#' @keywords bessel2
#' @rdname ys
#' @export
ys <- function(l, n) {
  # Function check =============================================================
  if (!is.numeric(l) || !is.numeric(n)) {
    stop("Inputs must be numeric vectors.")
  }
  # Internal function ==========================================================
  ys_internal <- function(l, n) {
    if (n > 0) {
      yc(l + 0.5, n) * sqrt(pi / (2 * n))
    } else if (n < 0) {
      (yc(l + 0.5, -n) * sqrt(pi / (2 * -n)))
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
  l / n * ys(l, n) - ys(l + 1, n)
}
#' @rdname ys
#' @export
ysdd <- function(l, n) {
  (-l / (n^2) * ys(l, n) +
    l / n * (l / n * ys(l, n) - ys(l + 1, n)) -
    ((l + 1) / n * ys(l + 1, n) - ys(l + 2, n)))
}
################################################################################
#' Spherical Bessel function of the third kind and its respective derivative
#' @param l A real integer or fractional order
#' @param n A real argument
#' @return
#' Calculates the spherical Bessel function of the third kind (hs) and its
#' first (hsd) derivative.
#' @references
#' Amos D.E., "AMOS, A Portable Package for Bessel Functions of a Complex
#' Argument and Nonnegative Order", http://netlib.org/amos
#'
#' @keywords bessel3
#' @rdname hs
#' @export
hs <- function(l, n) {
  sqrt(pi / (2 * n)) * hc(l + 0.5, n)
}
#' @rdname hs
#' @export
hsd <- function(l, n) {
  -hs(l + 1, n) + (l / n) * hs(l, n)
}
################################################################################
# Legendre Polynomials
################################################################################
################################################################################
#' Legendre Polynomial function (Pn) of the first kind.
#' @param n Degree of Legendre polynomial
#' @param x Real value
#' @return
#' Returns a matrix array calculated by computing the outer product that
#' generates a series of values when evaluated over a vector.
#'
#' @keywords legendre
#' @rdname Pn
#' @export
Pn <- function(n, x) {
  # Define internal recursive function =========================================
  Pn_internal <- function(n, x) {
    # Boolean statement ========================================================
    ## Pre-allocate result vector ==============================================
    result <- rep(
      x = 0,
      times = length(x)
    )
    ## Iterate through modal series constant vector ============================
    if (n == 0) {
      result <- rep(
        x = 1,
        times = length(x)
      )
    } else if (n == 1) {
      result <- x
    } else {
      for (i in seq_along(x)) {
        x_val <- x[i]
        P_n_1 <- x_val
        P_n_2 <- 1
        P_n <- 0
        for (j in 2:n) {
          P_n <- ((2 * j - 1) * x_val * P_n_1 - (j - 1) * P_n_2) / j
          P_n_2 <- P_n_1
          P_n_1 <- P_n
        }
        result[i] <- P_n
      }
    }
    result
  }
  # Vectorize the internal function ============================================
  Pn_vec <- Vectorize(Pn_internal)
  # Calculate outer product to generate vector of values =======================
  outer(n, x, Pn_vec)
}
#' Calculate Bessel function cache for the Goodman and Stern (1962) model
#' @param ka_matrix_m Modal ka matrix
#' @param m Modal vector
#' @return Cached Bessel function values
#' @keywords internal
#' @noRd
.calculate_bessel_cache <- function(ka_matrix_m, m) {
  # Define all Bessel functions to apply =======================================
  bessel_functions <- list(
    js = js, jsd = jsd, jsdd = jsdd,
    ys = ys, ysd = ysd, ysdd = ysdd,
    hs = hs, hsd = hsd
  )
  # Map out the function assignment ============================================
  bessel_map <- list(
    k1a_shell = c("js", "jsd", "hs", "hsd"),
    kLa_shell = c("js", "jsd", "jsdd", "ys", "ysd", "ysdd"),
    kTa_shell = c("js", "jsd", "jsdd", "ys", "ysd", "ysdd"),
    kTa_fluid = c("js", "jsd", "jsdd", "ys", "ysd", "ysdd"),
    kLa_fluid = c("js", "jsd", "jsdd", "ys", "ysd", "ysdd"),
    k3a_fluid = c("js", "jsd")
  )
  # Pre-calculate the Bessel functions and their respective derivatives ========
  bessel_cache <- lapply(names(ka_matrix_m), function(ka_m) {
    ka_series <- ka_matrix_m[[ka_m]]
    bessel_match <- bessel_map[[ka_m]]

    # Calculate for only matching functions
    lapply(bessel_functions[bessel_match], function(func) {
      func(m, ka_series)
    })
  })
  # Add the names ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  names(bessel_cache) <- names(ka_matrix_m)
  bessel_cache
}
