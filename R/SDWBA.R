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
#' @param niterations Number of times/iterations the model will be ran.
#' @param summary Will output the
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

SDWBA <- function(animal=NULL, x=animal@rpos[1,], y=animal@rpos[2,], z=animal@rpos[3,],
                  c=1500, frequency, phase=0.0, tilt=animal@theta, a=animal@a, h=animal@h, g=animal@g, pc=animal@pc){
  require(elliptic)
  rpos <- as.matrix(rbind(x,y,z))
  kt <- cbind(cos(tilt),rep(0,length(tilt)),sin(tilt))
  k1 <- kcalc(frequency,c)*kt; k2 <- vecnorm(k1) / h
  fbs <- 0 + 0i

  for(j in 1:(animal@ncyl-1)){
    r1 <- rpos[,j]; r2 <- rpos[,j+1]
    a1 <- a[j]; a2 <- a[j+1]
    beta <- abs(acos((k1%*%(r2-r1))/(vecnorm(k1)*vecnorm(r2-r1))) - pi/2)

    SDWBAint <- function(s){
      rint <- s * (r2-r1)+r1
      aint <- s * (a2-a1)+a1
      gamma <- 1/(g*h^2)+1/g-2

      if(abs(abs(beta)-pi/2)<1e-10){
        bessel <- k2*aint
      }else{
        bessel <- ja(1,2*k2*aint*cos(beta))/cos(beta)
      }

      if(animal@shape == "straight"){
        return(vecnorm(k1)/4*gamma*aint*exp(2i*k1%*%rint/h)*bessel*vecnorm(r2-r1))
      }else{
        return(vecnorm(k1)*pc/4*gamma*aint*exp(2i*k2*pc)*exp(-2i*k2*pc*cos(beta))*bessel*(vecnorm(r2-r1)/pc))
      }
    }

    SDWBAint <- Vectorize(SDWBAint)
    integral <- myintegrate(SDWBAint,0,1)
    fbs <- fbs + integral * exp(1i * rnorm(1,0,phase))
  }
  return(20*log10(abs(fbs)))
}
