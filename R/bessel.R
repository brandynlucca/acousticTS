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
  f <- 0; t <- 0; s <- 0; i <- 0; sign <- 1 #pre-allocate terms for while loop
  if(l < 0 & l%%1==0){
    sign <- -1
    l <- abs(l)
  }
  while(f == 0){
    s <- s + (-1)^i * (n/2)^(l+2*i) / (factorial(i) * gamma(i+l+1))
    if(abs(s - t) < 1e-30){f <- 1} #thresholding term
    t <- s; i <- i + 1 #iterate forward
  }
  return(ifelse(sign<1,(-1^l)*s,s))
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
  if(l%%1==0 & l != 0){
    sign <- 1
    if(l < 0){
      sign <- -1
      l <- abs(l)
    }
    p1 <- -(n/2)^(-l)/pi
    p2 <- 0
    for(i in seq(0,(l-1),1)){
      p2 <- p2 +
        (factorial(l-i-1)/factorial(i))*(n^2/4)^i
    }
    p3 <- 2/pi*log(n/2)*ja(l,n)
    f <- 0; t <- 0; p4 <- 0; i <- 0
    while(f == 0){
      p4 <- p4 +
        (digamma(i+1)+digamma(l+i+1))*((-n^2/4)^i)/(factorial(i)*factorial(l+i))
      if(abs(p4-t) < 1e-30){f <- 1}
      t <- p4; i <- i + 1
    }
    quant <- p1*p2+p3-(n/2)^l/pi*p4
    return(ifelse(sign<1,(-1^l)*quant,quant))
  }else if(l == 0){
    f <- 0; i <- 0; q <- 0; t <- 0
    p1 <- 2/pi * besselJ(n,l)*log(n/2)
    while(f == 0){
      q <- q + (-1)^i/factorial(i)^2*(n/2)^(2*i)*digamma(i+1)
      if(abs(q-t) < 1e-30){f <- 1}
      t <- q; i <- i + 1
    }
    qualt <- p1 - 2/pi*q
    return(qualt)
  }else{return((ja(l,n)*cos(l*pi)-ja(-l,n))/sin(l*pi))}
}

#' Calculate spherical Bessel function of the first kind.
#'
#' @description
#' This function calculates the spherical Bessel functions of the first kind. Referenced from Amost D.E.,
#' "AMOS, A Portable Package for Bessel Functions of a Complex Argument and Nonnegative Order", http://netlib.org/amos
#'
#' @param l An integer or fractional order
#' @param n A complex or real argument
#' @param sign Flag used to to determine whether to add or subtract 0.5 from the function's order (l); only used
#' for spherical Bessel function of first kind.
#' @usage
#' jl(l,n,sign)
#' @examples
#' l <- 1
#' n < 1
#' jl(l,n,sign=1)
#' [1] 0.3011687
#' @return
#' jl(l,n,sign) calculates the spherical Bessel function of the first kind (zbesj AMOS routine, Iv).
#' @export

#Spherical Bessel function of first kind
jl <- function(l,n,sign){ifelse(sign==1,ja(l+0.5,n) * sqrt(pi/2/n),ja(l-0.5,n) * sqrt(pi/2/n))}

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
yl <- function(l,n){ya(l+0.5,n) * sqrt(pi/2/n)}

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
jsd <- function(l,n){jl(l,n,-1) - (l+1) / n * jl(l,n,1)} #first kind

#' @export
jd <- function(l,n){
  s <- ja(l-1,n)-(l/n)*ja(l,n)
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
ysd <- function(l,n){l / n * yl(l,n) - yl(l+1,n)} #second kind

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
jsdd <- function(l,n){1 / (n^2) * ((l+1)*(l+2) - n^2) * jl(l,n,1) - 2 / n * jl(l,n,-1)}
