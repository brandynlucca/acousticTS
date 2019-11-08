# All formulas for calculation of TS derived from:
#

#' Calculates the Euclidean norm of a vector.
#'
#' @param x A vector with numeric, real values.
#' @usage
#' vecnorm(x)
#' @examples
#' values <- c(1,2,3)
#' vecnorm(values)
#' [1] 3.741657
#' @return
#' Calculates the Euclidean norm of a vector.
#' @export

#Euclidean vector norm
vecnorm <- function(x){sqrt(sum(x**2))} #Calculates Euclidean norm of a vecto

#' Toggle between radians and degrees.
#'
#' @param x A real value in degrees or radians
#' @param d The value input-type. Two input types: "deg" for degrees and "rad" for radians.
#' @usage
#' degrad(x,d)
#' @examples
#' x <- 180 #degrees
#' degrad(x, "deg")
#' [1] 3.141593
#' x <- pi #radians
#' degrad(x, "rad")
#' [1] 180
#' @return
#' Converts degrees to radians or radians to degrees
#' @export

#Toggle between radians and dgrees
degrad <- function(x,d){
  if(d == "deg"){
    value <- x*pi/180.0
  }else if(d == "rad"){
    value <- x*180.0/pi
  }
  return(value)
}


#' Calculate the acoustic wavenumber based on the sound speed of water.
#'
#' @param c Sound speed (m/s)
#' @param f Frequency (Hz)
#' @usage
#' kcalc(f,c)
#' @examples
#' c <- 1500 #m/s
#' f <- 120e3 #Hz
#' kcalc(f,c)
#' [1] 502.6547
#' @return
#' Calculates the acoustic wavenumber based on the sound speed of water
#' @export

#Calculate acoustic wavenumber based on the sound speed of water
kcalc <- function(f,c){2*pi*f/c}

#' Fluid-filled scatterer (FFS) object/class.
#'
#' @description
#' A S4 class that provides slots to contain relevant animal metadata for parameterizing models for fluid-filled scatterers (FFS) partitioned
#' into discretized cylinders. This, specifically, includes a position matrix, radius, material properties (g, h), orientation,
#' animal shape, and body curvature. This class is used within the DWBA and DFCM model functions. In the future, this will also allow for
#' converting one class of scatterer into another for seemless usage for model comparisons.
#' @export

#Create S4 class object to contain all animal metadata
FFS <- setClass("FFS",slots=c(rpos="matrix", a="numeric", g="numeric", h="numeric", theta="numeric", shape="character", pc="numeric"))

#' Calculates the theoretical TS of a fluid-filled scatterer at a given frequency using the distorted Born wave approximation (DWBA) model.
#'
#'
#' @param shape Desired object/animal shape. Must be class "FFS".
#' @param c Sound speed of surrounding medium (m/s). Default value is 1500 m/s.
#' @param f Frequency (Hz).
#' @param phi Phase deviation (\eqn{\phi}), or phase variability. Accounts for complexities in animal shape and stochasticity of noise in scattering field.
#' Default value is 0.0.
#' @param tilt Orientation of the target relative to the transmit source (\eqn{\theta}). Broadside incidence is considered 90 degrees, or pi/2.
#' Default value is pi/2; input should be in radians.
#' @usage
#' SDWBA(shape, c, f, phi, tilt)
#' @details
#' Calculates the theoretical TS of a fluid-filled scatterer at a given frequency using the distorted Born wave approximation (DWBA) model.
#' @return
#' Target strength (TS, dB re: 1 m^2)
#' @references
#' Stanton, T.K., Chu, D., and Wiebe, P.H. 1998. Sound scattering by several zooplankton groups. II. Scattering models. Journal of the Acoustical Society of America, 103(1), 236-253.
#'
#' Hankin, R.K.S. 2006. Introducing elliptic, an R package for elliptic and modular functions. Journal of Statistical Software, 15(7).
#' @export

SDWBA <- function(shape, c=1500, f, phi=0.0, tilt=pi/2){
  require(elliptic)
  k_1 <- cbind(cos(tilt), rep(0,length(tilt)), sin(tilt))
  k1 <- kcalc(f,c) * k_1
  k2 <- vecnorm(k1) / animal@h
  n <- length(animal@a)
  f.bs <- 0 + 0i

  for(j in 1:(n-1)){
    r1 <- c(animal@rpos[1,j], animal@rpos[2,j], animal@rpos[3,j])
    r2 <- c(animal@rpos[1,j+1], animal@rpos[2,j+1], animal@rpos[3,j+1])
    a1 <- animal@a[j]
    a2 <- animal@a[j+1]
    alphatilt <- acos((k1%*%(r2-r1)) / (vecnorm(k1)*vecnorm(r2-r1)))
    betatilt <- abs(alphatilt - pi/2)

    integrand <- function(s){
      rx <- s * (r2[1] - r1[1]) + r1[1]
      ry <- s * (r2[2] - r1[2]) + r1[2]
      rz <- s * (r2[3] - r1[3]) + r1[3]
      r <- c(rx, ry, rz)
      a <- s * (a2 - a1) + a1
      gamgam <- 1/(animal@g*animal@h^2)+1/animal@g-2

      if(abs(abs(betatilt) - (pi/2)) < 1e-10){
        bessy <- k2 * a
      }else{
        bessy <- ja(1,2*k2*a*cos(betatilt))/cos(betatilt)
      }

      if(animal@shape == "straight"){
        return(vecnorm(k1)/4*gamgam*a*exp(2i*k1%*%r/animal@h)*bessy*vecnorm(r2-r1))
      }else{
        pc <- animal@pc
        return(vecnorm(k1)*pc/4*gamgam*a*exp(1i*2.0*k2*pc)*exp(-1i*2.0*k2*pc*cos(betatilt))*bessy*(vecnorm(r2-r1)/pc))
      }
    }

    integrand <- Vectorize(integrand)

    f <- myintegrate(integrand, 0, 1)
    f.bs <- f.bs + f * exp(1i * rnorm(1,0,phi))
  }
  TS <- 20*log10(abs(f.bs))
  return(TS)
}

#' Calculates the theoretical TS of a fluid-filled scatterer at a given frequency using the deformed finite cylinder model (DFCM).
#'
#'
#' @param shape Desired object/animal shape. Must be class "FFS".
#' @param L Maximum length (m) of the scatterer. Can either be provided by shape input, or a manual value.
#' @param a Maximum radius (m) of the scatterer. Can either be provided by shape input, or a manual value.
#' @param g Density contrast (g) of the scatterer. Can either be provided by shape input, or a manual value.
#' @param h Sound speed contrast (h) of the scatterer. Can either be provided by shape input, or a manual value.
#' @param pc Radius of curvature (m) of the scatterer. Can either be provided by shape input, or a manual value.
#' @param numdiv Number of divisions to partition the animal shape into.
#' Can either be provided by the length of the shape position vector, or a manual value.
#' @param pm Density of surrounding medium (kg/m^3). Default value 1.025 kg/m^3.
#' @param cm Sound speed of surrounding medium (m/s). Default value is 1500 m/s.
#' @param f Frequency (Hz).
#' @param alpha Numerically derived coefficient (\eqn{\alpha_\beta}). Default is 0.8, and is the suggested value.
#' @param method Currently only the two-ray method is available, but the six-ray and other formulations will be available in the future.
#' @usage
#' SDWBA(shape, f)
#' SDWBA(L, a, g, h, pc, numdiv, tilt, f)
#' SDWBA(shape, L, a, g, h, pc, numdiv, tilt, cm, pm, f, alpha, method)
#' @details
#' Calculates the theoretical TS of a fluid-filled scatterer at a given frequency using the deformed finite cylinder model (DFCM).
#' @return
#' Target strength (TS, dB re: 1 m^2)
#' @references
#' Stanton, T.K., Chu, D., and Wiebe, P.H. 1998. Sound scattering by several zooplankton groups. II. Scattering models. Journal of the Acoustical Society of America, 103(1), 236-253.
#' @export

DFCM <- function(shape=NULL, L=max(shape@rpos[,1]), a=max(shape@a), g=shape@g, h=shape@h, pc=shape@pc, numdiv=length(shape@rpos[1,]),
                 tilt=shape@theta, cm=1500, pm=1.025, f, alpha=0.8, method="two-ray"){

  if(method == "two-ray"){
    ca <- h*cm #animal soundspeed
    pa <- g*pm #animal density
    a_new <- L/numdiv
    R12 <- ((pa*ca)/(cm*pm)-1) / (pa*ca/(pm*cm)+1)
    km <- kcalc(f,cm)
    ka <- kcalc(f,ca)
    T12 <- 2*(pa*ca/(pm*cm))/(1+(pa*ca/(pm*cm)))
    T21 <- 2*(pm*cm/(pa*ca))/(1+(pm*cm/(pa*ca)))
    mu <- -pi/2*km*a_new / (km*a_new+0.4)
    Io <- 1 - T12 * T21 * exp(1i*4*ka*a_new) * exp(1i*mu)
    fbs <- 0.5*sqrt(pc*a_new) * R12 * exp(-1i*2*km*a_new) * Io * exp(-alpha*(2*(tilt-pi/2)*pc/L)^2)
    TS <- 20*log10(abs(fbs))
    return(TS)
  }
}

