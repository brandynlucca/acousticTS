#' Calculates the theoretical TS of a fluid-like scatterer at a given frequency 
#' using the distorted Born wave approximation (DWBA) model.
#'
#' @param shape Desired object/animal shape. Must be class "FLS".
#' @param c Sound speed of surrounding medium (m/s). Default value is 1500 m/s.
#' @param frequency Frequency (Hz).
#' @param x,y,z The x-, y-, and z-axis coordinates that make up the position 
#' matrix, \eqn{r_0}.
#' @param a Radius vector of an animal (m).
#' @param h Sound speed contrast.
#' @param g Density contrast.
#' @param pc Radius of curvature. Default is 3.3.
#' @param curve A boolean value that dictates whether an animal is curved or 
#' not.
#' @param phase Phase deviation (\eqn{\phi}), or phase variability. Accounts for
#'  complexities in animal shape and stochasticity of noise in scattering field.
#' Default value is 0.0.
#' @param theta Orientation of the target relative to the transmit source 
#' (\eqn{\theta}). Broadside incidence is considered 90 degrees, or pi/2.
#' Default value is pi/2; input should be in radians.
#' @param ncyl Number of segments comprising the scatterer shape.
#' @usage
#' SDWBA(shape, c, frequency, phase)
#'
#' SDWBA(shape=NULL, c, frequency, x, y, z, a, h, g, pc, curve, phase, theta, 
#' ncyl)
#' @details
#' Calculates the theoretical TS of a fluid-filled scatterer at a given 
#' frequency using the distorted Born wave approximation (DWBA) model.
#' @return
#' Target strength (TS, dB re: 1 m^2)
#' @references
#' Stanton, T.K., Chu, D., and Wiebe, P.H. 1998. Sound scattering by several 
#' zooplankton groups. II. Scattering models. Journal of the Acoustical Society 
#' of America, 103(1), 236-253.
#'
#' Hankin, R.K.S. 2006. Introducing elliptic, an R package for elliptic and 
#' modular functions. Journal of Statistical Software, 15(7).
#' @importFrom elliptic myintegrate
#' @export

SDWBA <- function(shape=NULL, 
                  x=shape@rpos[1,], y=shape@rpos[2,], z=shape@rpos[3,],
                  c=1500, frequency, phase=0.0, a=shape@a, h=shape@h, g=shape@g,
                  curve=ifelse(is.null(shape),F,shape@curve),
                  pc=ifelse(is.null(shape),
                            ifelse(curve==T, 3.0, NA),
                            ifelse(curve==T,
                                   ifelse(shape@pc == 0, 3.0, shape@pc),
                                   NA)),
                  theta=ifelse(is.null(shape),pi/2,shape@theta),
                  ncyl=ifelse(is.null(shape),length(x),shape@ncyl),
                  L=shape@L,
                  progress=TRUE){
  
  rpos <- as.matrix(rbind(x,y,z))
  kt <- cbind(cos(theta),rep(0,length(theta)),sin(theta))
  k1 <- kcalc(frequency,c)*kt; k2 <- vecnorm(k1) / h
  fbs <- 0 + 0i
  
  for(j in 1:(ncyl-1)){
    r1 <- rpos[,j]; r2 <- rpos[,j+1]
    a1 <- a[j]; a2 <- a[j+1]
    alpha <- acos((k1 %*% (r2-r1)) / (vecnorm(k1)*vecnorm(r2-r1)))
    beta <- abs(alpha - pi/2)
    
    SDWBAint <- function(s){
      rint <- s * (r2-r1)+r1
      aint <- s * (a2-a1)+a1
      gamma <- 1/(g*h^2)+1/g-2
      
      if(abs(abs(beta) - pi/2) < 1e-10){
        bessel <- k2*aint
      }else{
        bessel <- ja(1,2*k2*aint*cos(beta))/cos(beta)
      }
      
      if(curve == F){
        ts <- vecnorm(k1)/4*gamma*aint*exp(2i*k1%*%rint/h)*bessel*vecnorm(r2-r1)
        return(ts)
      }else{
        pct <- pc*L
        ts <- vecnorm(k1)*pct/4*gamma*aint*exp(2i*k2*pct)*
          exp(-2i*k2*pct*cos(beta))*bessel*(vecnorm(r2-r1)/pct)
        return(ts)
      }
    }
    
    SDWBAint <- Vectorize(SDWBAint)
    integral <- elliptic::myintegrate(SDWBAint,0,1)
    fbs <- fbs + integral * exp(1i * rnorm(1,0,phase))
  }
  return(20*log10(abs(fbs)))
}

#' Calculates the theoretical TS of a target using the deformed cysinder 
#' model (DCM). 
#'
#' @param shape Scattering-class object.
#' @param frequency The acoustic frequency (Hz)
#' @param c The ambient seawater sound speed (m/s).
#' @param rho The ambient seawater density (kg/m^3).
#' @param method Specific DCM formulation. See 'Details' for options. 
#' 
#' @param alpha Conditional. Required for the two-ray path equation. Numerically 
#' derived coefficient, (\eqn{\alpha_\beta}). Default is 0.8, and is the 
#' suggested value.
#' 
#' @param L Optional. Scatterer length (m). 
#' @param a Optional. Scatterer radius (m). 
#' @param LAratio Optional. Length:radius ratio.
#' @param pc Optional. Ratio of the radius of curvature and length.
#' @param g Optional. Density contrast. 
#' @param h Optional. Sound speed contrast. 
#' @param theta Optional. Orientation (\eqn{\theta}, radians). 
#' @param N Number of discrete cysinders. 
#' 
#' @usage 
#' DCM(shape, frequency, c, rho, method, alpha, LAratio, L, a, g, h, pc, theta, N)
#' 
#' @details
#' Only the two-ray path model is afvailable at the moment.. 
#' @return
#' Target strength (TS, dB re: 1 m^2)
#' @references
#' Stanton, T.K., Chu, D., and Wiebe, P.H. 1998. Sound scattering by several 
#' zooplankton groups. II. Scattering models. Journal of the Acoustical Society 
#' of America, 103(1), 236-253.
#' @export
DCM <- function(shape=NULL, frequency, c=1500, rho=1026, method="two-ray",
                alpha=0.8, LAratio=NULL,
                L=max(shape@L), a=NULL, g=shape@g, h=shape@h, 
                pc=shape@pc,theta=shape@theta, N=shape@ncyl-1){
  if(method == "two-ray"){
    c1 <- h*c; rho1 <- g*rho
    a <- ifelse(is.null(a), L/LAratio, a)
    pcL <- pc*L
    R12 <- ((rho1*c1)/(c*rho)-1)/(rho1*c1/(rho*c)+1)
    k <- kcalc(frequency, c); k1 <- k/h
    T12 <- 2*(rho1*c1/(rho*c))/(1+(rho1*c1/(rho*c)))
    T21 <- 2*(rho*c/(rho1*c1))/(1+(rho*c/(rho1*c1)))
    mu <- -pi/2*k*a/(k*a+0.4)
    I0 <- 1-T12*T21*exp(4i*k1*a)*exp(1i*mu*a)
    fbs <- 0.5*sqrt(pcL*a)*R12*exp(-2i*k*a)*I0*
      exp(-alpha*(2*(theta-pi/2)*pcL/L)^2)
    return(20*log10(Mod(fbs)))
  }
}

#' Calculates theoretical TS of a solid sphere of a certain material at a given 
#' frequency.
#'
#' @description
#' This function is a wrapper around TS_calculate(...) that parametrizes the 
#' remainder of the model, while also doing simple calculations that do not 
#' need to be looped. This function provides a TS estimate at a given frequency.
#' @param shape CAL-class object.
#' @param frequency The acoustic frequency (Hz)
#' @param c The ambient seawater sound speed (m/s).
#' @param rho The ambient seawater density (kg/m^3).
#' 
#' @param a Optional. Radius (m).
#' @param c1 Optional. Longitudinal sound speed (m/s).
#' @param c2 Optional. Transversal sound speed (m/s).
#' @param rho1 Optional. Density of sphere (kg/m^3).
#' @param g Optional. Density contrast. 
#' @param h1 Optional. Longitudinal sound speed contrast. 
#' @param h2 Optional. Transversal sound speed contrast. 
#' 
#' @usage
#' ssphere(shape, frequency, c, rho, ...)
#' 
#' @return The theoretical acoustic target strength (TS, dB re: 1 m^2) of a 
#' solid sphere at a given frequency.
#' @references 
#' MacLennan D. N. (1981). The theory of solid spheres as sonar calibration 
#' targets. Scottish Fisheries Research No. 22, Department of Agriculture and 
#' Fisheries for Scotland.
#' @export
ssphere <- function(shape=NULL, frequency, c=1500, rho=1026,
                    a=shape@a, c1=shape@c1, c2=shape@c2, rho1=shape@rho1,
                    g=shape@rho1/rho, h1=shape@c1/c, h2=shape@c2/c){
  
  ka <- kcalc(frequency,c)*a; ka1 <- ka/h1; ka2 <- ka/h2  #Equations 6a
  alpha <- 2*g*h2^2; beta <- g*h1^2 - alpha #Equations 6d & 6e
  m <- 0; fsub <- 0; t <- 0; f <- TRUE #initialize
  #Iterate through while loop
  while(f == TRUE){
    A2 <- (m^2+m-2)*js(m,ka2)+ka2^2*jsdd(m,ka2)  #Equation 6b
    A1 <- 2*m*(m+1)*(ka1*jsd(m,ka1)-js(m,ka1)) #Equation 6c
    B2 <- A2*ka1^2*(beta*js(m,ka1)-alpha*jsdd(m,ka1))-
      A1*alpha*(js(m,ka2)-ka2*jsd(m,ka2)) #Equation 6f
    B1 <- ka*(A2*ka1*jsd(m,ka1)-A1*js(m,ka2)) #Equation 6g
    eta <- atan(-(B2*jsd(m,ka)-B1*js(m,ka))/
                  (B2*ysd(m,ka)-B1*ys(m,ka))) #Equation 6h
    t <- (-1)^m*(2*m+1)*sin(eta)*exp(1i*eta) #Equation 7a
    fsub <- fsub + t
    
    if(abs(t/fsub) < 1e-10){
      f <- FALSE
    }else{
      m <- m+1
    }
  }
  fsphere <- -2/ka*fsub #Equation 7b
  sigmabs <- pi*a^2*abs(Mod(fsphere))^2 #Equation 8
  return(10*log10(sigmabs/(4*pi)))
}

#' Calculates the theoretical TS of a fluid sphere using an exact modal series 
#' solution proposed by Andersen (1950).
#'
#' @param shape Desired animal object. Must be class FS.
#' @param c Sound speed of surrounding medium (m/s). Default value is 1500 m/s.
#' @param frequency Frequency (Hz).
#' @param a Radius vector of sphere (m).
#' @param h Sound speed contrast.
#' @param g Density contrast.
#' @param method Type of sphere being modeled. There are several options: 1 
#' (fluid-filled), 2 (fixed-rigid), 3 (pressure-release). Currently, fluid 
#' shells with fluid interiors are not available, but will be in future
#' versions.
#' @usage
#' fsphere(shape, frequency, c)
#'
#' fsphere(shape=NULL, frequency, c, a, g, h)
#' @details
#' Calculates the theoretical TS of a fluid sphere at a given frequency using 
#' an exact modal series solution.
#' @return
#' Target strength (TS, dB re: 1 m^2)
#' @references
#' Anderson, V.C. (1950). Sound scattering from a fluid sphere. Journal of the 
#' Acoustical Society of America, 22, 426-431.
#' @export
fsphere <- function(shape=NULL, frequency, c=1500, a=shape@a, 
                    g=shape@g, h=shape@h, method=1){
  k1 <- kcalc(frequency, c); k2 <- kcalc(frequency, c*h)
  m <- 0; fsub <- 0; t <- 0; f <- TRUE #initialize
  #Iterate through while loop
  while(f == TRUE){
    if(method == 1){
      cmn <- (jsd(m,k2*a)*ys(m,k1*a)) / (js(m,k2*a)*jsd(m,k1*a)) - 
        (g*h*ysd(m,k1*a)/jsd(m,k1*a)) #numerator term
      cmd <- (jsd(m,k2*a)*js(m,k1*a)) / (js(m,k2*a)*jsd(m,k1*a)) - (g*h)
      cm <- cmn/cmd
      bm <- -1 / (1 + 1i*cm)
    }else if(method == 2){
      bm <- -(jsd(m,k1*a)/hsd(m,k1*a))
    }else if(method == 3){
      bm <- -(js(m,k1*a)/hs(m,k1*a))
    }
    t <- ((2*m+1)*(-1)^m*bm)
    fsub <- fsub + t
    
    if(abs(t)/abs(fsub) < 1e-30){
      f <- FALSE
    }else{
      m <- m + 1
    }
  }
  fbubble <- -1i/k1 * fsub
  return(20*log10(abs(fbubble)))
}

#' Calculates the theoretical TS using Kirchoff-ray Mode approximation.
#'
#' @param shape Desired object/animal shape. Must be class "SBF".
#' @param c Sound speed of surrounding medium (m/s). Default value is 1500 m/s.
#' @param rho Density of medium (kg m^-1). Default value is 1030 kg m^-1.
#' @param frequency Frequency (Hz).
#' @param theta Orientation of the target relative to the transmit source 
#' (\eqn{\theta}). Broadside incidence is considered 90 degrees, or pi/2.
#' Default value is pi/2; input should be in radians.
#' @usage
#' KRM(shape, c, rho, frequency, theta)
#' @details
#' Calculates the theoretical TS using the Kirchoff-ray Mode model.
#' @return
#' Target strength (TS, dB re: 1 m^2)
#' @references
#' Clay C.S. and Horne J.K. (1994). Acoustic models of fish: The Atlantic cod 
#' (Gadus morhua). Journal of the Acoustical Society of AMerica, 96, 1661-1668.
#' @export
#' @import foreach
#' @import doSNOW
#' @importFrom parallel detectCores
#' @import snow

KRM <- function(shape=NULL, c=1500, rho=1030, frequency, theta=pi/2){
  shape@theta <- theta
  k <- kcalc(frequency, c); kb <- kcalc(frequency, shape@cb)
  RBC <- (shape@psb*shape@csb - shape@pb*shape@cb) / (shape@psb*shape@csb + 
                                                        shape@pb*shape@cb)
  RWB <- (shape@pb*shape@cb - rho*c) / (shape@pb*shape@cb + rho*c)
  TT <- 1 - RWB^2
  frequency.soft <- 0i; frequency.fluid <- 0i
  
  for(i in 1:(length(shape@bladder[1,])-1)){
    p1 <- shape@bladder[,i]; p2 <- shape@bladder[,i+1]
    as <- (p1[2]+p2[2])/4
    Asb <- k*as/(k*as+0.083)
    Psip <- k*as/(40+k*as)-1.05
    vs <- ((p1[1]+p2[1])*cos(shape@theta) + (p1[3]+p2[3])*sin(shape@theta))/2
    delus <- (p2[1]-p1[1])*sin(shape@theta)
    frequency.soft <- frequency.soft +
      -1i*(RBC*TT)/(2*sqrt(pi))*Asb*sqrt((k*as+1)*sin(shape@theta))*
      exp(-1i*(2*k*vs+Psip))*delus
  }
  
  for(i in 1:(length(shape@body[1,])-1)){
    b1 <- shape@body[,i]; b2 <- shape@body[,i+1]
    ab <- (b1[2]+b2[2])/4
    Psib <- -pi*kb*((b1[3]+b2[3])/2)/(2*(kb*((b1[3]+b2[3])/2)+0.4))
    delub <- (b2[1]-b1[1])*sin(shape@theta)
    vbU <- ((b1[3]+b2[3])*cos(shape@theta) + (b1[3]+b2[3])*sin(shape@theta))/2
    vbL <- ((b1[4]+b2[4])*cos(shape@theta) + (b1[4]+b2[4])*sin(shape@theta))/2
    
    frequency.fluid <- frequency.fluid +
      -1i*(RWB/(2*sqrt(pi)))*sqrt(k*ab)*delub*
      (exp(-2i*k*vbU)-TT*exp(-2i*k*vbU+2i*kb*(vbU-vbL)+1i*Psib))
    
  }
  sigma <- abs(frequency.fluid + frequency.soft)
  TS <- 20*log10(sigma)
  return(TS)
}

#' Calculates the theoretical TS of a shelled organism using the non-modal High Pass (HP) model
#' 
#' @param shape Desired animal object. Must be class FS.
#' @param c Sound speed of surrounding medium (m/s). Default value is 1500 m/s.
#' @param frequency Frequency (Hz).
#' 
#' @param a optional. Radius vector of sphere (m).
#' @param h Optional. Sound speed contrast.
#' @param g Optional. Density contrast.
#' 
#' @return
#' Target strength (TS, dB re: 1 m^2)
#' @references
#' Lavery, A.C., Wiebe, P.H., Stanton, T.K., Lawson, G.L., Benfield, M.C., 
#' Copley, N. 2007. Determining dominant scatterers of sound in mixed 
#' zooplankton popuilations. The Journal of the Acoustical Society of America,
#' 122(6): 3304-3326.
#' 
#' @export
HP <- function(shape=NULL, c=1500, frequency, 
               g=shape@g, h=shape@h, a=shape@a){
  k <- kcalc(frequency, c)
  api <- (1-g*h^2)/(3*g*h^2) + (1-g)/(1+2*g)
  R <- (g*h-1)/(g*h+1)
  fbs <- a^2*(k*a)^4*api^2/((1+4*(k*a)^4*api^2)/R^2)
  TS <- 10*log10(fbs)
  return(TS)
}
