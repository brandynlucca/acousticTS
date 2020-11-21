#' Calculate Bessel function of the first kind.
#'
#' @param l An integer or fractional order
#' @param n A complex or real argument
#' @usage
#' ja(l,n)
#' @examples
#' l <- 1.5
#' n <- 1
#' ja(l,n)
#' [1] 0.2402978
#' @return
#' Calculates the Bessel function of the first kind (J_a).
#' @references
#' Amost D.E., "AMOS, A Portable Package for Bessel Functions of a Complex Argument and Nonnegative Order", http://netlib.org/amos
#' @export

#Bessel function of first or second kind calculation
ja <- function(l,n){
  if(n == 0){
    s <- 0
  }else{
    if(n > 0){
      nn <- n
    }else{
      nn <- abs(n)
    }

    if(nn > 10){
      s <- sqrt(2/(pi*nn)) * cos(nn - (l/2 + 0.25)*pi)
    }else if(nn < 0.01){
      s <- (0.5*nn)^l / gamma(l+1)
    }else{
      if(l >= 0){
        ll <- l
      }else if(l < 0 & l%%1==0){
        ll <- abs(l)
      }else{
        ll <- l
      }

      f <- 0; t <- 0; s <- 0; i <- 0
      while(f == 0){
        s <- s + (-1)^i * (nn/2)^(ll+2*i) / (factorial(i) * gamma(i+ll+1))

        if(abs(s - t) < 1e-30){f <- 1}
        t <- s; i <- i + 1

      }

      if(n < 0){
        if(l > 0){
          s <- -1^l * s
        }else if(l%%1!=0){
          if(s < 0){
            s <- 1i*s
          }else{
            s <- -1i*s
          }
        }
      }else{
        if(l < 0 & l%%1==0){
          s <- -1^l * s
        }
      }
    }
  }
  return(s)
}

#' Calculate Bessel function of the second kind.
#'
#' @param l An integer or fractional order
#' @param n A complex or real argument
#' @usage
#' ya(l,n)
#' @examples
#' l <- 1.5
#' n <- 1
#' ya(l,n)
#' [1] -1.102496
#' @return
#' Calculates the Bessel function of the second kind (Y_a). Functions are referenced from Amost D.E., "AMOS, A Portable Package
#' for Bessel Functions of a Complex Argument and Nonnegative Order", http://netlib.org/amos
#' @export

#Bessel function of first or second kind calculation
ya <- function(l,n){

  nn <- abs(n)
  ll <- abs(l)

  if(ll%%1==0){
    p1 <- -(nn/2)^(-ll)/pi
    p2 <- 0
    for (i in seq(0, (ll - 1), 1)) {
      p2 <- p2 + (factorial(ll - i - 1)/factorial(i)) *
        (nn^2/4)^i
    }
    p3 <- 2/pi * log(nn/2) * ja(ll, nn)
    f <- 0; t <- 0; i <- 0
    p4 <- 0

    while (f == 0) {
      p4 <- p4 + (digamma(i + 1) + digamma(ll + i + 1)) *
        ((-nn^2/4)^i)/(factorial(i) * factorial(ll + i))
      if (abs(p4 - t) < 1e-30) {
        f <- 1
      }
      t <- p4
      i <- i + 1
    }
    s <- p1 * p2 + p3 - (nn/2)^ll/pi * p4
  }else{
    s <- (ja(ll, nn) * cos(ll * pi) - ja(-ll, nn))/sin(ll * pi)
  }

  if(n < 0){
    if(n%%1 == 0){
      if(l%%1==0){
        s <- -1^ll * s + 1i*(-1)^ll*2*ja(ll,nn)
      }else{
        s <- 1i*-1^(ll+1) * s
      }
    }
  }

  if(l < 0){
    if(l%%1==0){
      s <- -1^l * s
    }else{
      s <- -1^l * ja(ll+0.5, n)
    }
  }

  return(s)
}

#' Calculate spherical Bessel function of the first kind.
#'
#' @description
#' This function calculates the spherical Bessel functions of the first kind. Referenced from Amost D.E.,
#' "AMOS, A Portable Package for Bessel Functions of a Complex Argument and Nonnegative Order", http://netlib.org/amos
#'
#' @param l An integer or fractional order
#' @param n A complex or real argument
#' @usage
#' jl(l,n)
#' @examples
#' l <- 1
#' n < 1
#' jl(l,n)
#' [1] 0.3011687
#' @return
#' jl(l,n) calculates the spherical Bessel function of the first kind (zbesj AMOS routine, Iv).
#' @export

#Spherical Bessel function of first kind
jl <- function(l,n){
  if(n == 0){
    return(0)
  }else if(n > 0){
    return(ja(l+0.5, n) * sqrt(pi/2/n))
  }else{
    return(-1^l * (ja(abs(l)+0.5, abs(n)) * sqrt(pi/2/abs(n))))
  }
}

#' Calculate spherical Bessel function of the second kind.
#'
#' @description
#' This function calculates the spherical Bessel functions of the second kind. Referenced from Amost D.E.,
#' "AMOS, A Portable Package for Bessel Functions of a Complex Argument and Nonnegative Order", http://netlib.org/amos
#'
#' @param l An integer or fractional order
#' @param n A complex or real argument
#' @usage
#' yl(l,n)
#' @examples
#' l <- 1
#' n < 1
#' yl(l,n)
#' [1] -1.381773
#' @return
#' yl(l,n) calculates the spherical Bessel function of the second kind (zbesy AMOS routine, Yv)
#' @export

#Spherical Bessel function of second kind
yl <- function(l,n){
  if(n > 0){
    return(ya(l+0.5,n) * sqrt(pi/2/n))
  }else if(n < 0){
    return(-(ya(l+0.5,n) * sqrt(pi/2/n)))
  }else if(n == 0){
    return(Inf)
  }
}

#' Calculate the first derivative of the spherical Bessel function of the first kind.
#'
#' @description
#' These functions alculate the first and second derivatives of the spherical Bessel functions required for the elastic sphere
#' TS model. These comprise the first derivatives of the spherical Bessel functions of the first and second kind,
#' as well as the second derivative of the spherical Bessel function of the first kind.
#' Functions are referenced from Amost D.E., "AMOS, A Portable Package for Bessel Functions of a Complex Argument
#' and Nonnegative Order", http://netlib.org/amos
#'
#' @param l An integer or fractional order
#' @param n A complex or real argument
#' @usage
#' jsd(l,n)
#' @examples
#' l <- 1
#' n < 1
#' jsd(l,n)
#' [1] 0.2391336
#' @return
#' jsd(l,n) calculates the first derivative of the spherical Bessel function of the first kind (zbesj AMOS routine, Iv).
#' @export
#'
#First derivatives of spherical Bessel functions
jsd <- function(l,n){
  return(-jl(l+1, n) + (l/n)*jl(l,n))
} #first kind

#'
#' @export
jd <- function(l,n){
  if(n == 0){
    if(l == 0){
      s <- 0
    }else if(l == 1){
      s <- 0.5
    }else{
      s <- 0
    }
  }else{
    s <- ja(l-1,n)-(l/n)*ja(l,n)
  }
  return(s)
}

#' @export
yd <- function(l,n){
  sign <- 1
  if(l < 0){
    sign <- -1
    l <- abs(l)
  }
  s <- ya(l-1,n)-(l/n)*ya(l,n)
  return(ifelse(sign<1,(-1^l)*s,s))
}

#' Bessel function of the third kind
#'@export
ha <- function(l,n){return(ja(l,n) + 1i*ya(l,n))}

#' Spherical Bessel function of the third kind
#' @export
hl <- function(l,n){return(sqrt(pi/(2*n))*ha(l+0.5,n))}

#' First derivative of the spherical Bessel function of the third kind
#' @export
hsd <- function(l,n){return(-hl(l+1,n)+(l/n)*hl(l,n))}

#' First derivative of the Bessel function of the third kind
#' @export
had <- function(l, n){
return((l*ha(l,n)/n - ha(l+1,n)))
}

#' Calculate the first derivative of the spherical Bessel function of the second kind.
#'
#' @description
#' These functions alculate the first and second derivatives of the spherical Bessel functions required for the elastic sphere
#' TS model. These comprise the first derivatives of the spherical Bessel functions of the first and second kind,
#' as well as the second derivative of the spherical Bessel function of the first kind.
#' Functions are referenced from Amost D.E., "AMOS, A Portable Package for Bessel Functions of a Complex Argument
#' and Nonnegative Order", http://netlib.org/amos
#'
#' @param l An integer or fractional order
#' @param n A complex or real argument
#' @usage
#' ysd(l,n)
#' @examples
#' l <- 1
#' n < 1
#' ysd(l,n)
#' [1] 2.223244
#' @return
#' yd(l,n) calculates the first derivative of the spherical Bessel function of the second kind (zbesy AMOS routine, Yv).
#' @export

#First derivatives of spherical Bessel functions
ysd <- function(l,n){l / n * yl(l,n) - yl(l+1)} #second kind

#' Calculate the second derivative of the spherical Bessel function of the first kind.
#'
#' @description
#' These functions alculate the first and second derivatives of the spherical Bessel functions required for the elastic sphere
#' TS model. These comprise the first derivatives of the spherical Bessel functions of the first and second kind,
#' as well as the second derivative of the spherical Bessel function of the first kind.
#' Functions are referenced from Amost D.E., "AMOS, A Portable Package for Bessel Functions of a Complex Argument
#' and Nonnegative Order", http://netlib.org/amos
#'
#' @param l An integer or fractional order
#' @param n A complex or real argument
#' @usage
#' jsdd(l,n)
#' @examples
#' l <- 1
#' n < 1
#' jsdd(l,n)
#' [1] -0.1770986
#' @return
#' jsdd(l,n) calculates the second derivative of the spherical Bessel function of the first kind (zbesj, AMOS routine, Iv).
#' @export

#Second derivative of spherical Bessel function of first kind
jsdd <- function(l,n){1 / (n^2) * ((l+1)*(l+2) - n^2) * jl(l,n) - 2 / n * jl(l,n)}

#' Second derivative of the Bessel function of first kind
#' @export
jdd <- function(l,n){
  s <- (0.25)*(ja(l-2,n) - 2*ja(l,n)+ja(l+2,n))
  return(s)
}

