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
#' @param l An integer or fractional order
#' @param n A complex or real argument
#' @return
#' Calculates the cylindrical Bessel function of the first kind (J_v) and its
#' respective derivatives.
#' @references
#' Amos D.E., "AMOS, A Portable Package for Bessel Functions of a Complex
#' Argument and Nonnegative Order", http://netlib.org/amos
#' @rdname jc
#' @export
jc <- function( l , n ) besselJ( n , l )
#' @rdname jc
#' @export
jcd <- function( l , n ) {
  jc_internal <- function( l , n ) {
    if ( n == 0 ) {
      if ( l == 1 ) {
        return( 0.5 )
      } else {
        return( 0.0 )
      }
    } else {
      return( jc( l - 1 , n ) - ( l / n ) * jc( l , n ) )
    }
  }
  jc_vec <- Vectorize( jc_internal )
  # Return based on input class ================================================
  switch( class( n )[ 1 ] ,
          numeric = jc_vec( l , n ) ,
          matrix = apply( n , 2 , FUN = function( x ) {
            jc_vec( l , x )
          } ) ) -> result
  return( result )
}
#' @rdname jc
#' @export
jcdd <- function( l , n ) {
  return( 0.25 * ( jc( l - 2 , n ) - 2 * jc( l , n ) + jc( l + 2 , n ) ) )
}
################################################################################
#' Cylindrical Bessel (Neumann) function of the second kind and its derivative
#' @description
#' Calculate the cylindrical Bessel function of the first kind (yc) and its
#' first (ycd) derivatives. 
#' @param l An integer or fractional order
#' @param n A complex or real argument
#' @return
#' Calculates the cylindrical Bessel function of the first kind (Y_v) and its
#' respective derivatives.
#' @references
#' Amos D.E., "AMOS, A Portable Package for Bessel Functions of a Complex
#' Argument and Nonnegative Order", http://netlib.org/amos
#' @rdname yc
#' @export
yc <- function( l , n ) besselY( n , l )
#' @rdname yc
#' @export 
ycd <- function( l , n ) {
  yc_internal <- function( l , n ) {
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
  yc_vec <- Vectorize( yc_internal )
  # Return based on input class ================================================
  switch( class( n )[ 1 ] ,
          numeric = yc_vec( l , n ) ,
          matrix = apply( n , 2 , FUN = function( x ) {
            yc_vec( l , x )
          } ) ) -> result
  return( result )
}
################################################################################
#' Cylindrical Bessel (Hankel) function of the third kind and its derivative
#' @param l An integer or fractional order
#' @param n A complex or real argument
#' @return
#' Calculates the cylindrical Bessel function of the first kind (H_v) and its
#' respective derivatives.
#' @references
#' Amos D.E., "AMOS, A Portable Package for Bessel Functions of a Complex
#' Argument and Nonnegative Order", http://netlib.org/amos
#' @rdname hc
#' @export
hc <- function( l , n ) return( jc( l , n ) + 1i * yc( l , n ) )
#' @rdname hc
#' @export
hcd <- function( l, n ) return( ( l * hc( l , n ) / n - hc( l + 1 , n ) ) )
################################################################################
# Spherical Bessel functions
################################################################################
################################################################################
#' Spherical Bessel function of the first kind and its respective derivatives
#' @param l An integer or fractional order
#' @param n A complex or real argument
#' @return
#' Calculates the spherical Bessel function of the first kind (js) and both its
#' first (jsd) and second (jsdd) derivatives.
#' @references
#' Amos D.E., "AMOS, A Portable Package for Bessel Functions of a Complex
#' Argument and Nonnegative Order", http://netlib.org/amos
#' @rdname js
#' @export
js <- function( l , n ) {
  # Function check =============================================================
  if ( ! is.numeric( l ) && ! numeric( n ) )
    stop( "Inputs must be numeric vectors." )
  # Internal function ==========================================================
  js_internal <- function( l , n ) {
    if( n == 0 ) {
      return( 0 )
    } else {
      return( jc( l + 0.5 , n ) * sqrt( pi / ( 2 * n ) ) )
    }
  }
  js_vec <- Vectorize( js_internal )
  # Return based on input class ================================================
  switch( class( n )[1] ,
          numeric = js_vec( l , n ) ,
          matrix = apply( n , 2 , FUN = function( x ) {
            js_vec( l , x )
          } ) ) -> result
  return( result )
}
#' @rdname js
#' @export
jsd <- function( l , n ) {
  return( js( l - 1 , n ) - ( l + 1 ) / n * js( l , n ) )
} 
#' @rdname js
#' @export
jsdd <- function( l , n ) {
  return( 1 / ( n ^ 2 ) *
            ( ( l + 1 ) * ( l + 2 ) - n ^ 2 ) * js( l , n ) -
            2 / n * js( l - 1 , n ) )
}
################################################################################
#' Spherical Bessel function of the second kind and its respective derivative
#' @param l An integer or fractional order
#' @param n A complex or real argument
#' @return
#' Calculates the spherical Bessel function of the first kind (js) and its
#' first (ysd) derivative.
#' @references
#' Amos D.E., "AMOS, A Portable Package for Bessel Functions of a Complex
#' Argument and Nonnegative Order", http://netlib.org/amos
#' @rdname ys
#' @export
ys <- function( l , n ) {
  # Function check =============================================================
  if ( ! is.numeric( l ) && ! numeric( n ) )
    stop( "Inputs must be numeric vectors." )
  # Internal function ==========================================================
  ys_internal <- function( l , n ) {
    if( n > 0 ) {
      return( yc( l + 0.5 , n ) * sqrt( pi / ( 2 * n ) ) )
    }else if( n < 0 ){
      return( - ( yc( l + 0.5 , n ) * sqrt (pi / ( 2 * n ) ) ) )
    }else if( n == 0 ){
      return( -Inf )
    }
  }
  ys_vec <- Vectorize( ys_internal )
  # Return based on input class ================================================
  switch( class( n )[1] ,
          numeric = ys_vec( l , n ) ,
          matrix = apply( n , 2 , FUN = function( x ) {
            ys_vec( l , x )
          } ) ) -> result
  return( result )
}
#' @rdname ys
#' @export
ysd <- function( l , n ) {
  return( l / n * ys( l , n ) - ys( l + 1 , n ) )
}
################################################################################
#' Spherical Bessel function of the third kind and its respective derivative
#' @param l An integer or fractional order
#' @param n A complex or real argument
#' @return
#' Calculates the spherical Bessel function of the third kind (hs) and its
#' first (hsd) derivative.
#' @references
#' Amos D.E., "AMOS, A Portable Package for Bessel Functions of a Complex
#' Argument and Nonnegative Order", http://netlib.org/amos
#' @rdname hs
#' @export
hs <- function( l , n) return( sqrt( pi / ( 2 * n ) ) * hc( l + 0.5 , n ) )
#' @rdname hs
#' @export
hsd <- function( l , n ) return( - hs( l + 1 , n ) + ( l / n ) * hs( l , n ) )
################################################################################
# Legendre Polynomials
################################################################################
################################################################################
#' Legendre Polynomial function (Pn) of the first kind.
#' @param n Degree of Legendre polynomial
#' @param x Real value
#' @return
#' Returns a matrix array calculated by computing the outer product that generates
#' a series of values when evaluated over a vector.
#' @rdname Pn
#' @export
Pn <- function( n , x ) {
  # Define internal recursive function =========================================
  Pn_internal <- function( n , x ) {
    # Boolean statement ========================================================
    ## Pre-allocate result vector ==============================================
    result <- rep( x = 0 ,
                   times = length( x ) )
    ## Iterate through modal series constant vector ============================
    if ( n == 0 ) {
      result <- rep( x = 1 ,
                     times = length( x ) )
    } else if ( n == 1 ) {
      result <- x
    } else {
      for ( i in 1 : length( x ) ) {
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
  return( outer( n , x , Pn_vec ) )
}