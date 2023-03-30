################################################################################
# DIFFERENTIAL EQUATION SOLUTION FUNCTIONS
################################################################################
################################################################################
# Cyndrical Bessel functions
################################################################################
################################################################################
#' Wrapper function for the Cylindrical Bessel function of the first kind (Jc)
#' @param l An integer or fractional order
#' @param n A complex or real argument
#' @usage
#' jc( l , n )
#' @return
#' Calculates the Cylindrical Bessel function of the first kind (J_a).
#' @references
#' Amost D.E., "AMOS, A Portable Package for Bessel Functions of a Complex
#' Argument and Nonnegative Order", http://netlib.org/amos
#' @rdname jc
#' @export
jc <- function( l , n ) return( base::besselJ( n , l ) )
################################################################################
#' Wrapper function for the Cylindrical Bessel (Neumann) function of the second kind (Yc)
#' @param l An integer or fractional order
#' @param n A complex or real argument
#' @usage
#' yc( l , n )
#' @return
#' Calculates the Bessel function of the second kind (Y_v). Functions are
#' referenced from Amost D.E., "AMOS, A Portable Package
#' for Bessel Functions of a Complex Argument and Nonnegative Order",
#' http://netlib.org/amos
#' @rdname yc
#' @export
yc <- function( l , n ) return( base::besselY( n , l ) )
################################################################################
#' Cylindrical Bessel (Hankel) function of the third kind (Hc)
#' @param l An integer or fractional order
#' @param n A complex or real argument
#' @usage
#' hc( l , n )
#' @return
#' Calculates the Bessel function of the second kind (H_v). Functions are
#' referenced from Amost D.E., "AMOS, A Portable Package
#' for Bessel Functions of a Complex Argument and Nonnegative Order",
#' http://netlib.org/amos
#' @rdname hc
#' @export
hc <- function( l , n ) return( jc( l , n ) + 1i * yc( l , n ) )
################################################################################
# Cyndrical Bessel function derivatives (first and second)
################################################################################
################################################################################
#' First derivative of the Cylindrical Bessel function of the first kind (Jcd)
#' @param l An integer or fractional order
#' @param n A complex or real argument
#' @usage
#' jcd( l , n )
#' @rdname jc
#' @export
jcd <- function( l , n ) {
  ## Iterate through modal series constant vector ==============================
  if ( n == 0.0 ) {
    if ( l == 1 ) {
      return( 0.5 )
    } else {
      return( 0.0 )
    }
  } else {
    return( jc( l - 1 , n ) - ( l / n ) * jc( l , n ) )
  }
}
################################################################################
#' Second derivative of the Cylindrical Bessel function of the first kind (Jcdd)
#' @param l An integer or fractional order
#' @param n A complex or real argument
#' @usage
#' jcdd( l , n )
#' @rdname jc
#' @export
jcdd <- function( l , n ) {
  return( 0.25 * ( jc( l - 2 , n ) - 2 * jc( l , n ) + jc( l + 2 , n ) ) )
}
################################################################################
#' First derivative of Cylindrical Bessel function of the second kind (Ycd)
#' @param l An integer or fractional order
#' @param n A complex or real argument
#' @usage
#' ycd( l , n )
#' @rdname yc
#' @export
ycd <- function( l , n ) {
  ## Iterate through modal series constant vector ==============================
  if ( n == 0.0 ) {
    if ( l == 1 ) {
      return( 0.5 )
    } else {
      return( 0.0 )
    }
  } else {
    return( yc( l - 1 , n ) - ( l / n ) * yc( l , n ) )
  }
}
################################################################################
#' First derivative of the Bessel function of the third kind (Hcd)
#' @param l An integer or fractional order
#' @param n A complex or real argument
#' @usage
#' hcd( l , n )
#' @rdname hc
#' @export
hcd <- function( l, n ) return( ( l * hc( l , n ) / n - hc( l + 1 , n ) ) )
################################################################################
# Spherical Bessel functions
################################################################################
################################################################################
#' Spherical Bessel function of the first kind (ja)
#' @param l An integer or fractional order
#' @param n A complex or real argument
#' @usage
#' js( l , n )
#' @return
#' js( l , n ) calculates the spherical Bessel function of the first kind
#' (zbesj AMOS routine, Iv).
#' @rdname js
#' @export
#Spherical Bessel function of first kind
js <- function( l , n ) {
  # Function check =============================================================
  if ( !base::is.numeric( l ) && !base::numeric( n ) )
    stop( "Inputs must be numeric vectors." )
  # Internal function ==========================================================
  js_internal <- function( l , n ) {
    if( n == 0 ) {
      return( 0 )
    } else {
      return( jc( l + 0.5 , n ) * base::sqrt( pi / ( 2 * n ) ) )
    }
  }
  js_vec <- Vectorize( js_internal )
  # Return based on input class ================================================
  base::switch( class( n )[1] ,
                numeric = js_vec( l , n ) ,
                matrix = apply( n , 2 , FUN = function( x ) {
                  js_vec( l , x )
                } ) ) -> result
  return( result )
}
################################################################################
#' Calculate spherical Bessel function of the second kind (ya)
#' @description
#' This function calculates the spherical Bessel functions of the second kind.
#' Referenced from Amost D.E.,
#' "AMOS, A Portable Package for Bessel Functions of a Complex Argument and
#' Nonnegative Order", http://netlib.org/amos
#'
#' @param l An integer or fractional order
#' @param n A complex or real argument
#' @usage
#' ys( l , n )
#' @return
#' ys( l , n ) calculates the spherical Bessel function of the second kind (zbesy
#' AMOS routine, y_v)
#' @rdname ys
#' @export
#Spherical Bessel function of second kind
ys <- function( l , n ) {
  # Function check =============================================================
  if ( !base::is.numeric( l ) && !base::numeric( n ) )
    stop( "Inputs must be numeric vectors." )
  # Internal function ==========================================================
  ys_internal <- function( l , n ) {
    if( n > 0 ) {
      return( yc( l + 0.5 , n ) * base::sqrt( pi / ( 2 * n ) ) )
    }else if( n < 0 ){
      return( - ( yc( l + 0.5 , n ) * base::sqrt (pi / ( 2 * n ) ) ) )
    }else if( n == 0 ){
      return( -Inf )
    }
  }
  ys_vec <- Vectorize( ys_internal )
  # Return based on input class ================================================
  base::switch( class( n )[1] ,
                numeric = ys_vec( l , n ) ,
                matrix = apply( n , 2 , FUN = function( x ) {
                  ys_vec( l , x )
                } ) ) -> result
  return( result )
}
################################################################################
#' Spherical Bessel function of the third kind (ha)
#' @param l An integer or fractional order
#' @param n A complex or real argument
#' @usage
#' hs( l , n )
#' @rdname hs
#' @export
hs <- function( l , n) return( base::sqrt( pi / ( 2 * n ) ) * hc( l + 0.5 , n ) )
################################################################################
# Spherical Bessel function derivatives (first and second)
################################################################################
################################################################################
#' Calculate the first derivative of the spherical Bessel function of the
#' first kind (jad)
#'
#' @param l An integer or fractional order
#' @param n A complex or real argument
#' @usage
#' jsd( l , n )
#' @rdname js
#' @export
jsd <- function( l , n ) {
  return( js( l - 1 , n ) - ( l + 1 ) / n * js( l , n ) )
}
################################################################################
#' Calculate the second derivative of the spherical Bessel function of the
#' first kind (jadd)
#' @param l An integer or fractional order
#' @param n A complex or real argument
#' @usage
#' jsdd( l , n )
#' @rdname js
#' @export
jsdd <- function( l , n ) {
  return( 1 / ( n ^ 2 ) *
            ( ( l + 1 ) * ( l + 2 ) - n ^ 2 ) * js( l , n ) -
            2 / n * js( l - 1 , n ) )
}
################################################################################
#' First derivatives of spherical Bessel functions (yad)
#' @param l An integer or fractional order
#' @param n A complex or real argument
#' @usage
#' ysd( l , n )
#' @rdname ys
#' @export
ysd <- function( l , n ) {
  return( l / n * ys( l , n ) - ys( l + 1 , n ) )
}
################################################################################
#' First derivative of the spherical Bessel function of the third kind (had)
#' @param l An integer or fractional order
#' @param n A complex or real argument
#' @usage
#' hsd( l , n )
#' @rdname hs
#' @export
hsd <- function( l , n ) return( - hs( l + 1 , n ) + ( l / n ) * hs( l , n ) )
################################################################################
# Legendre Polynomials
################################################################################
################################################################################
#' Legendre Polynomial function (Pn) of the first kind.
#' @param n Polynomial degree.
#' @param x Interval bounded by -1 and 1
#' @rdname Pn
#' @export
Pn <- function( n , x ) {
  # Define internal recursive function =========================================
  Pn_internal <- function( n , x ) {
    # Boolean statement ========================================================
    ## Pre-allocate result vector ==============================================
    result <- base::rep( x = 0 ,
                         times = base::length( x ) )
    ## Iterate through modal series constant vector ============================
    if ( n == 0 ) {
      result <- base::rep( x = 1 ,
                           times = base::length( x ) )
    } else if ( n == 1 ) {
      result <- x
    } else {
      for ( i in 1 : base::length( x ) ) {
        x_val <- x[i]
        P_n_1 <- x_val
        P_n_2 <- 1
        P_n <- 0
        for ( j in 2 : n ) {
          P_n <- ( ( 2 * j - 1 ) * x_val * P_n_1 - ( j - 1 ) * P_n_2 ) / j
          P_n_2 <- P_n_1
          P_n_1 <- P_n
        }
        result[i] <- P_n
      }
    }
    return( result )
  }
  # Vectorize the internal function ============================================
  Pn_vec <- Vectorize( Pn_internal )
  # Calculate outer product to generate vector of values =======================
  return( base::outer( n , x , Pn_vec ) )
}
