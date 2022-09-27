#' Calculate the acoustic wavenumber based on the sound speed of water.
#'
#' @param sound_speed Sound speed (c, m/s)
#' @param frequency Frequency (f, Hz)
#' @usage
#' kcalc(frequency, sound_speed)
#' @examples
#' c <- 1500 #m/s
#' f <- 120e3 #Hz
#' kcalc(f,c)
#' # 502.6547
#' @return
#' Calculates the acoustic wavenumber based on the sound speed of water
#' @export
kcalc <- function(frequency, sound_speed){
  return(2 * pi * frequency / sound_speed)
}

#' Primary access function.
#' @param object Scatterer-class object.
#' @param feature Feature of interest (e.g. body).
#' @export
setGeneric("extract", function(object, feature)
  standardGeneric("extract"))

#' Method for what is printed for objects.
#' @param object Scatterer-class object.
#' @param feature Feature of interest (e.g. body).
#' @export
setMethod("extract",
          signature(object = "scatterer"),
          function(object,
                   feature) {
            if(feature == "model" & class(object) == "CAL") {
              data_TS <- slot(object, feature)$calibration
              data_x <- slot(object, "model_parameters")$calibration$parameters$acoustics
              output <- data.frame(frequency = data_x$frequency,
                                   k_sw = data_x$k_sw,
                                   k_l = data_x$k_l,
                                   k_t = data_x$k_t,
                                   f_bs = data_TS$f_bs,
                                   sigma_bs = data_TS$sigma_bs,
                                   TS = data_TS$TS)
              return(output)
            } else {
              return(slot(object, feature))
            }
          })

#' Calculate reflection coefficient
#' @param material_properties_1 Interface 1
#' @param material_properties_2 Interface 2
#' @export
reflection_coefficient <- function(material_properties_1,
                                   material_properties_2) {
  mat1_condensed <- material_properties_1$density*material_properties_1$sound_speed
  mat2_condensed <- material_properties_2$density*material_properties_2$sound_speed
  R <- (mat2_condensed - mat1_condensed) / (mat2_condensed + mat1_condensed)
  return(R)
}

#' Calculate reflection coefficient
#' @param body Body material properties (i.e. g, h)
#' @param medium Seawater material properties (i.e., density, sound speed)
#' @export
R12 <- function(body, medium) {
  # Animal sound speed ==========================
  c1 <- body$h * medium$sound_speed
  # Animal density ==============================
  rho1 <- body$g * medium$density
  # Calculate Reflection coefficient
  R12 <- ((c1 * rho1) / (medium$sound_speed * medium$density) - 1) /
    ((c1 * rho1) / (medium$sound_speed * medium$density) + 1)
  return(R12)
}
#' Sums along the position vector
#' @param rpos Position vector
#' @param iterations Number of iterations
#' @export
along_sum <- function(rpos, iterations) {
  output <- rpos[, 1:(iterations-1)] + rpos[, 2:iterations]
  return(output)
  }

#' Calculates the Euclidean norm of a matrix input.
#'
#' @param x A matrix with numeric, real values.
#' @usage
#' vecnorm(x)
#' @examples
#' values <- matrix(c(1,2,3), ncol=3)
#' vecnorm(values)
#' # 3.741657
#' @return
#' Calculates the Euclidean norm of a vector.
#' @export
vecnorm <- function(x){
  return(sqrt(rowSums(x^2)))
}

#' Convert backscatter values from log- to linear-domain
#'
#' @param value Logarithmic (e.g. TS or \eqn{S_V}) or linear (\eqn{\sigma_bs}) value
#' @param coefficient Optional. Numeric coefficient preceding the logarithm. Default is 10.
#' @examples
#' \dontrun{
#' TS <- -50 #dB re. 1 m^2
#' sigma_bs <- linear(TS) #convert to sigma_bs, 10^(TS/10)
#' print(sigma_bs)
#' TS_new <- dB(sigma_bs) #convert back to TS
#' print(TS_new)
#' }
#' @return
#' Calculates the acoustic wavenumber based on the sound speed of water
#' @export
linear <- function(value, coefficient = 10){
  return(coefficient^(value / coefficient))
}
#' @rdname linear
dB <- function(value, coefficient=10){
  return(coefficient * log10(value))
}

#' Support rotatin function for KRM (swimbladder)
#' @inheritParams body_rotation
#' @export
bladder_rotation <- function(sum_rpos, rpos, theta, k_length){
  v <- (sum_rpos[1, ]*cos(theta) + sum_rpos[3, ]*sin(theta)) / 2
  v <- matrix(data = rep(v, each = k_length),
              ncol = length(v),
              nrow = k_length)
  delta_u <- diff(rpos[1, ]) * sin(theta)
  return(list(v = v, delta_u = delta_u))
}

#' Support rotating function for KRM (body)
#' @param sum_rpos Summed position matrix
#' @param rpos Position matrix
#' @param theta Orientation angle
#' @param k_length Length of wavenumber vector
#'
#' @export
body_rotation <- function(sum_rpos, rpos, theta, k_length){
  dorsal_sum <- matrix(data = rep(sum_rpos[3, ], each = k_length),
                       ncol = length(sum_rpos[3, ]),
                       nrow = k_length)
  ventral_sum <- matrix(data = rep(sum_rpos[4, ], each = k_length),
                        ncol=length(sum_rpos[4, ]),
                        nrow = k_length)

  vbU <- (dorsal_sum * cos(theta) + dorsal_sum * sin(theta)) / 2
  vbL <- (ventral_sum * cos(theta) + ventral_sum * sin(theta)) / 2
  delta_u <- diff(rpos[1, ]) * sin(theta)
  return(list(vbU = vbU, vbL = vbL, delta_u = delta_u))
}


#' Bessel functions of the first kind
#'
#' @param l An integer or fractional order
#' @param n A complex or real argument
#' @usage
#' ja(l,n)
#' @examples
#' l <- 1.5
#' n <- 1
#' ja(l,n)
#' # 0.2402978
#' @return
#' Calculates the Bessel function of the first kind (J_a).
#' @references
#' Amost D.E., "AMOS, A Portable Package for Bessel Functions of a Complex
#' Argument and Nonnegative Order", http://netlib.org/amos
#' @rdname ja
#' @export

#Bessel function of first or second kind calculation
ja <- function(l, n){
  return(ifelse(l >= 0 & n >= 0,
                besselJ(n,l),
                ifelse(l < 0 & n > 0,
                       besselJ(abs(n), l),
                       -1^abs(l)*besselJ(abs(n), abs(l)))))

}

#'First derivative of Bessel function of the first kind
#' @param l An integer or fractional order
#' @param n A complex or real argument
#' @rdname ja
#' @export
jd <- function(l, n){
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

#' Second derivative of the Bessel function of first kind
#' @param l An integer or fractional order
#' @param n A complex or real argument
#' @rdname ja
#' @export
jdd <- function(l, n){
  s <- (0.25)*(ja(l-2,n) - 2*ja(l,n)+ja(l+2,n))
  return(s)
}

#' Spherical Bessel function of the first kind.
#'
#' @param l An integer or fractional order
#' @param n A complex or real argument
#' @usage
#' js(l,n)
#' @examples
#' l <- 1
#' n <- 1
#' js(l,n)
#' # 0.3011687
#' @return
#' js(l,n) calculates the spherical Bessel function of the first kind
#' (zbesj AMOS routine, Iv).
#' @rdname js
#' @export

#Spherical Bessel function of first kind
js <- function(l, n) {
  s <- function(l, n) {
    if(n == 0){
      return(0)
    }else if(n > 0){
      return(ja(l+0.5, n) * sqrt(pi/(2*n)))
    }else{
      return(-1^l * (ja(abs(l)+0.5, abs(n)) * sqrt(pi/(2*abs(n)))))
    }
  }

  s <- Vectorize(s)

  switch(class(n)[1],
         numeric = s(l, n),
         matrix = apply(n, 2, FUN=function(x) s(l, x))) -> out
  return(out)
}

#' Calculate the first derivative of the spherical Bessel function of the
#' first kind.
#'
#' @param l An integer or fractional order
#' @param n A complex or real argument
#'
#' @rdname js
#' @export
jsd <- function(l, n){
  return((js(l-1,n)-(l+1)/n*js(l,n)))
}

#' Calculate the second derivative of the spherical Bessel function of the
#' first kind.
#'
#' @param l An integer or fractional order
#' @param n A complex or real argument
#'
#' @rdname js
#' @export
jsdd <- function(l,n){return(1/(n^2)*((l+1)*(l+2)-n^2)*js(l,n)-2/n*js(l-1,n))}

#' Calculate Bessel function of the second kind.
#'()
#' @param l An integer or fractional order
#' @param n A complex or real argument
#' @usage
#' ya(l,n)
#' @examples
#' l <- 1.5
#' n <- 1
#' ya(l,n)
#' # -1.102496
#' @return
#' Calculates the Bessel function of the second kind (Y_a). Functions are
#' referenced from Amost D.E., "AMOS, A Portable Package
#' for Bessel Functions of a Complex Argument and Nonnegative Order",
#' http://netlib.org/amos
#' @rdname ya
#' @export

#Bessel function of first or second kind calculation
ya <- function(l,n){
  s <- ifelse(n >= 0,
              besselY(n, l),
              ifelse(l >= 0,
                     -1^l*besselY(abs(n), abs(l)),
                     -besselY(abs(n), abs(l))))
  return(s)
}

#'First derivative of Bessel function of the second kind
#' @param l An integer or fractional order
#' @param n A complex or real argument
#'
#' @rdname ya
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

#' Calculate spherical Bessel function of the second kind.
#'
#' @description
#' This function calculates the spherical Bessel functions of the second kind.
#' Referenced from Amost D.E.,
#' "AMOS, A Portable Package for Bessel Functions of a Complex Argument and
#' Nonnegative Order", http://netlib.org/amos
#'
#' @param l An integer or fractional order
#' @param n A complex or real argument
#' @usage
#' ys(l,n)
#' @examples
#' l <- 1
#' n <- 1
#' ys(l,n)
#' # -1.381773
#' @return
#' ys(l,n) calculates the spherical Bessel function of the second kind (zbesy
#' AMOS routine, Yv)
#' @rdname ys
#' @export
#Spherical Bessel function of second kind
ys <- function(l,n){
  s <- function(l, n){
    if(n > 0){
      return(ya(l+0.5,n) * sqrt(pi/(2*n)))
    }else if(n < 0){
      return(-(ya(l+0.5,n) * sqrt(pi/(2*n))))
    }else if(n == 0){
      return(Inf)
    }
  }

  s <- Vectorize(s)

  switch(class(n)[1],
         numeric = s(l, n),
         matrix = apply(n, 2, FUN=function(x) s(l, x))) -> out
  return(out)
}

#' First derivatives of spherical Bessel functions
#'
#' @param l An integer or fractional order
#' @param n A complex or real argument
#' @rdname ys
#' @export
ysd <- function(l,n){l / n * ys(l,n) - ys(l+1,n)} #second kind

#' Bessel function of the third kind
#' @param l An integer or fractional order
#' @param n A complex or real argument
#' @rdname ha
#' @export
ha <- function(l,n){return(ja(l,n) + 1i*ya(l,n))}

#' First derivative of the Bessel function of the third kind
#' @param l An integer or fractional order
#' @param n A complex or real argument
#' @rdname ha
#' @export
had <- function(l, n){
  return((l*ha(l,n)/n - ha(l+1,n)))
}


#' Spherical Bessel function of the third kind
#' @param l An integer or fractional order
#' @param n A complex or real argument
#' @rdname hs
#' @export
hs <- function(l,n){return(sqrt(pi/(2*n))*ha(l+0.5,n))}

#' First derivative of the spherical Bessel function of the third kind
#' @param l An integer or fractional order
#' @param n A complex or real argument
#' @rdname hs
#' @export
hsd <- function(l,n){return(-hs(l+1,n)+(l/n)*hs(l,n))}
