# Theoretical TS for elastic spheres
# All formulas for calculation of TS derived from: MacLennan D.N. (1981) The Theory of Solid Spheres as Sonar
# Calibration Targets. Scottish Fisheries Research No. 22, Department of Agriculture and Fisheries for Scotland.

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
#' Calculates the Bessel function of the first kind (J_a). Functions are referenced from Amost D.E., "AMOS, A Portable Package 
#' for Bessel Functions of a Complex Argument and Nonnegative Order", http://netlib.org/amos
#' @export

#Bessel function of first or second kind calculation
ja <- function(l,n){
    f <- 0; t <- 0; s <- 0; i <- 0 #pre-allocate terms for while loop

    while(f == 0){
      s <- s + (-1)^i * (n/2)^(l+2*i) / (factorial(i) * gamma(i+l+1))
      
      if(abs(s - t) < 1e-30){f <- 1} #thresholding term
      t <- s; i <- i + 1 #iterate forward
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
ya <- function(l,n){return((ja(l,n) * cos(l*pi) - ja(-l,n)) / sin(l*pi))}

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
#' @param sign Flag used to to determine whether to add or subtract 0.5 from the function's order (l); only used 
#' for spherical Bessel function of first kind.
#' @usage 
#' jl(l,n,sign) 
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
#' jd(l,n)
#' @examples 
#' l <- 1
#' n < 1 
#' jd(l,n)
#' [1] 0.2391336
#' @return 
#' jd(l,n) calculates the first derivative of the spherical Bessel function of the first kind (zbesj AMOS routine, Iv).
#' @export

#First derivatives of spherical Bessel functions
jd <- function(l,n){jl(l,n,-1) - (l+1) / n * jl(l,n,1)} #first kind

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
#' yd(l,n)
#' @examples 
#' l <- 1
#' n < 1 
#' yd(l,n)
#' [1] 2.223244
#' @return 
#' yd(l,n) calculates the first derivative of the spherical Bessel function of the second kind (zbesy AMOS routine, Yv).
#' @export

#First derivatives of spherical Bessel functions
yd <- function(l,n){l / n * yl(l,n) - yl(l+1,n)} #second kind

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
#' jdd(l,n)
#' @examples 
#' l <- 1
#' n < 1 
#' jdd(l,n)
#' [1] -0.1770986
#' @return 
#' jdd(l,n) calculates the second derivative of the spherical Bessel function of the first kind (zbesj, AMOS routine, Iv).
#' @export

#Second derivative of spherical Bessel function of first kind
jdd <- function(l,n){1 / (n^2) * ((l+1)*(l+2) - n^2) * jl(l,n,1) - 2 / n * jl(l,n,-1)}

#' Parse the acoustic material properties of an elastic sphere with respect to the type of mateiral.
#'
#' @description 
#' This function parses the acoustic material properties of four typical types of calibration spheres: tungsten carbide, copper,
#' stainless steel, and aluminum. These material properties include the transerval and longitudinal sound speeds (m/s), as well as 
#' the density (kg/m^3). 
#' 
#' @param material Accepts one of four character arguments: Tungsten carbide (default), Copper, Stainless steel, and Aluminum.
#' @usage 
#' sphere_param(material="Tungsten carbide")
#' @return Returns a 1x4 dataframe which includes the material (Material, factor), longitudinal sound speed (c1, numeric, m/s), transversal 
#' sound speed (c2, numeric, m/s), and density (rho1, numeric, kg/m^3).
#' @export

#Sphere parameters function
sphere_param <- function(material="Tungsten carbide"){
  sphere_props <- data.frame(Material=c("Tungsten carbide","Copper","Stainless steel","Aluminum"),
                             c1=c(6853,4760,5610,6260),
                             c2=c(4171,2288.5,3120,3080),
                             rho1=c(14900,8947,7800,2700))
  return(sphere_props[which(tolower(material) == tolower(sphere_props$Material)),])
}

#' Form function for calculating the acoustic target strength (TS, dB re: 1 m^2) of an elastic sphere.
#'
#' @description 
#' This function calculates the acoustic target strength (TS, dB re: 1 m^2) of an elastic sphere. All equations are taken from 
#' MacLennan (1981). 
#' @param q The acoustic wavenumber (k) multiplied by the sphere radius (a) (Equation 6a).
#' @param q1 The relationship between q (or ka) and the ratio of seawater-to-transversal soundspeed (Equation 6a).
#' @param q2 The relationship between q (or ka) and the ratio of seawater-to-longitudinal soundspeed (Equation 6a).
#' @param alpha The relationship between the seawater-to-sphere soundspeed and density ratios (Equation 6d). 
#' @param beta The relationship with alpha subtracted from between the seawater-to-sphere soundspeed 
#' and density ratios (Equation 6e).
#' @param a The radius of the target sphere (mm). 
#' @return The theoretical acoustic target strength (TS, dB re: 1 m^2) of an elastic sphere. 
#' @export

#Calculate theoretical TS at a single frequency for a calibration sphere
TS_calculate <- function(q,q1,q2,alpha,beta,a){
  q <- as.numeric(q)
  k <- TRUE #While loop flag
  l <- 0 #Iteration value
  foo_ll <- 0
  foo_l_prev <- 0 #Previous iteration of summed term in Equation 7
  
  while(k == TRUE){
    A2 <- (l^2 + l - 2)* jl(l,q2,1) + q2^2 * jdd(l,q2) #Equation 6b
    A1 <- 2 * l * (l + 1) * (q1 * jd(l,q1) - jl(l,q1,1)) #Equation 6c
    B2 <- A2 * q1^2 * (beta*jl(l,q1,1) - alpha*jdd(l,q1)) - A1 * alpha * (jl(l,q2,1) - q2 * jd(l,q2)) #Equation 6f
    B1 <- q * (A2 * q1 * jd(l,q1) - A1 * jl(l,q2,1)) #Equation 6g
    eta_l <- atan(-(B2 * jd(l,q) - B1 * jl(l,q,1)) / (B2 * yd(l,q) - B1*yl(l,q))) #Equation 6h
    foo_l <- (-1)^l * (2*l+1) * sin(eta_l) * exp(1i*eta_l) #Equation 7; what is calculated and summed for "l"th iteration
    foo_ll <- foo_ll + foo_l #Sum up to current iteration
    
    if(abs(foo_l/foo_ll) < 1e-10){
      k <- FALSE
    }else{
      l <- l + 1
    }
  }
  
  foo_q <- -2 / q * foo_ll #Equation 7; full calculation
  sigma <- pi * a^2 * abs(foo_q)^2 #Equation 8
  TS <- 10*log10(sigma/(4*pi)) #Equation 9
  return(TS)
}

#' Calculate the acoustic target strength (TS, dB re: 1 m^2) of an elastic sphere of a certain material at a given frequency.
#' 
#' @description 
#' This function is a wrapper around TS_calculate(...) that parametrizes the remainder of the model, while also doing 
#' simple calculations that do not need to be looped. This function provides a TS estimate at a given frequency.
#' @param frequency The acoustic frequency (Hz) 
#' @param c The ambient seawater sound speed (m/s). 
#' @param rho The ambient seawater density (kg/m^3).
#' @param material The material of the calibration sphere (tungsten carbide, aluminum, stainless steel, copper)
#' @param diameter The diameter of the sphere (mm).
#' @examples 
#' frequency <- 120e3
#' c <- 1500
#' rho <- 1030
#' material <- "Tungsten carbide"
#' diameter <- 38.1
#' sphere.ts(frequency,c,rho,material,diameter)
#' [1] -39.52
#' @return The theoretical acoustic target strength (TS, dB re: 1 m^2) of an elastic sphere at a given frequency. 
#' @export

#Calculate theoretical TS frequency spectrum for a calibration sphere
sphere.ts <- function(frequency,c,rho,material="Tungsten carbide",diameter=38.1){
  sphere_mat <- sphere_param(material) #call sphere material properties
  a <- diameter*1e-3/2 #calculate radius; convert from diameter to radius; m
  alpha <- 2*(as.numeric(sphere_mat[["rho1"]])/rho)*(as.numeric(sphere_mat[["c2"]])/c)^2 #Equation 6d
  beta <- (as.numeric(sphere_mat[["rho1"]])/rho)*(as.numeric(sphere_mat[["c1"]])/c)^2 - alpha #Equation 6e
  k <- 2*pi*frequency/c #acoustic wavenumber
  q <- k*a #ka variable; Equation 6a
  q1 <- q * c / as.numeric(sphere_mat[["c1"]]) #Equation 6a
  q2 <- q * c / as.numeric(sphere_mat[["c2"]]) #Equation 6a
  return(TS_calculate(q,q1,q2,alpha,beta,a))
}

#' Calculate the acoustic target strength (TS, dB re: 1 m^2) of an elastic sphere of a certain material at a given frequency.
#' 
#' @description 
#' This function is a wrapper around TS_calculate(...) that parametrizes the remainder of the model, while also doing 
#' simple calculations that do not need to be looped. This function provides a mean TS which includes a FM sweep around
#' a center frequency. This is calculated by adding and subtracting half an estimated bandwidth (1/pulse length) and taking
#' the mean (calculated in the linear domain).
#' 
#' @param frequency The acoustic frequency (Hz) 
#' @param c The ambient seawater sound speed (m/s). 
#' @param rho The ambient seawater density (kg/m^3).
#' @param material The material of the calibration sphere (tungsten carbide, aluminum, stainless steel, copper)
#' @param diameter The diameter of the sphere (mm).
#' @param pulse_length The pulse length (microseconds) of a signal. 
#' @examples 
#' frequency <- 120e3
#' c <- 1500
#' rho <- 1030
#' material <- "Tungsten carbide"
#' diameter <- 38.1
#' pulse_length <- 256
#' sphere.cw(frequency,c,rho,material,diameter,pulse_length)
#' [1] -39.549
#' @return The theoretical acoustic target strength (TS, dB re: 1 m^2) of 
#'    an elastic sphere at a given frequency and pulse length. 
#' @export

#Calculate theoretical TS frequency spectrum for a calibration sphere
sphere.cw <- function(frequency,c,rho,material="Tungsten carbide",diameter=38.1,pulse_length){
  bandwidth <- 1/(pulse_length/1e6) #BW = 1/pulse length, pulse length here is converted to seconds
  fs <- frequency - bandwidth/2 #frequency start
  fe <- frequency + bandwidth/2 #frequency end
  cw_spec <- seq(fs,fe,length.out=100) #calculate frequency spectrum with 100 steps
  sphere_mat <- sphere_param(material) #call sphere material properties
  a <- diameter*1e-3/2 #calculate radius; convert from diameter to radius; m
  alpha <- 2*(as.numeric(sphere_mat[["rho1"]])/rho)*(as.numeric(sphere_mat[["c2"]])/c)^2 #Equation 6d
  beta <- (as.numeric(sphere_mat[["rho1"]])/rho)*(as.numeric(sphere_mat[["c1"]])/c)^2 - alpha #Equation 6e
  
  ts <- rep(NA,length(cw_spec)) #pre-allocate ts vector to calculate average
  
  pb <- winProgressBar(title="Processing frequencies...", label="0% completed", min=0, max=100, initial=0)
  for(i in 1:length(cw_spec)){
    k <- 2*pi*cw_spec[i]/c #acoustic wavenumber
    q <- k*a #ka variable; Equation 6a
    q1 <- q * c / as.numeric(sphere_mat[["c1"]]) #Equation 6a
    q2 <- q * c / as.numeric(sphere_mat[["c2"]]) #Equation 6a
    ts[i] <- TS_calculate(q,q1,q2,alpha,beta,a)
    load.info <- sprintf("%d%% completed", round((i/length(ts))*100))
    setWinProgressBar(pb, i/(length(ts))*100, label=load.info)
  }
  close(pb)
  
  return(10*log10(mean(10^(ts/10)))) #calculate average in linear domain
}

#' Calculate the frequency spectrum of the acoustic target strength (TS, dB re: 1 m^2) of an elastic sphere .
#' 
#' @description 
#' This function is a wrapper around TS_calculate(...) that parametrizes the remainder of the model, while also doing 
#' simple calculations that do not need to be looped. This function provides TS estimates along a frequency spectrum.
#' 
#' @param c The ambient seawater sound speed (m/s). 
#' @param rho The ambient seawater density (kg/m^3).
#' @param material The material of the calibration sphere (tungsten carbide, aluminum, stainless steel, copper)
#' @param diameter The diameter of the sphere (mm).
#' @param fs Frequency start (Hz).
#' @param fe Frequency end (Hz). 
#' @param fi Frequency interval (Hz).
#' @examples 
#' frequency <- 120e3
#' c <- 1500
#' rho <- 1030
#' material <- "Tungsten carbide"
#' diameter <- 38.1
#' 
#' sphere.spec(frequency,c,rho,material,diameter,fs=18e3,fe=230e3,fi=1e3)
#' @return Returns the TS-frequency spectrum of an elastic sphere.
#' @export

sphere.spec <- function(c,rho,material="Tungsten carbide",diameter,fs,fe,fi){
  freq_spec <- seq(from=fs, to=fe, by=fi) #frequency spectrum
  sphere_mat <- sphere_param(material) #call sphere material properties
  a <- diameter*1e-3/2 #calculate radius; convert from diameter to radius; m
  alpha <- 2*(as.numeric(sphere_mat[["rho1"]])/rho)*(as.numeric(sphere_mat[["c2"]])/c)^2 #Equation 6d
  beta <- (as.numeric(sphere_mat[["rho1"]])/rho)*(as.numeric(sphere_mat[["c1"]])/c)^2 - alpha #Equation 6e
  
  tmp.df <- data.frame(Frequency=freq_spec, TS=NA) #Generate dataframe
  pb <- winProgressBar(title="Processing frequencies...", label="0% completed", min=0, max=100, initial=0)
  
  for(i in 1:nrow(tmp.df)){
    k <- 2*pi*tmp.df$Frequency[i]/c #acoustic wavenumber 
    q <- k*a #ka variable; Equation 6a
    q1 <- q * c / as.numeric(sphere_mat[["c1"]]) #Equation 6a
    q2 <- q * c / as.numeric(sphere_mat[["c2"]]) #Equation 6a
    tmp.df$TS[i] <- round(TS_calculate(q,q1,q2,alpha,beta,a),2)
    load.info <- sprintf("%d%% completed", round((i/nrow(tmp.df))*100))
    setWinProgressBar(pb, i/(nrow(tmp.df))*100, label=load.info)
  }
  
  close(pb)
  
  return(tmp.df)
}

#' Plots the TS-frequency spectrum.
#' 
#' @description 
#' This function generates a plot of the TS-frequency spectrum using ggplot2. 
#' 
#' @param df A dataframe generated from the sphere.spec(...) function with two columns: Frequency and TS.
#' @examples 
#' ts.df <- sphere.spec(frequency=120e3,c=1500,rho=1030,material="Tungsten carbide",diameter=38.1,fs=18e3,fe=230e3,fi=1e3)
#' sphere.spec_plot(ts.df)
#'  
#' @return 
#' Generates TS-frequency spectrum for an elastic sphere. 
#' @export
#' @import ggplot2

sphere.spec_plot <- function(df){
  require(ggplot2)
  
  p1 <- ggplot(data=df, aes(x=Frequency/1000,y=TS)) + geom_path(size=1) + geom_point(size=3) +
    theme_bw() +
    theme(text=element_text(size=18), axis.text=element_text(size=18, colour="black"),
          panel.grid=element_blank()) +
    labs(x="Frequency (kHz)",
         y=expression(paste("Target strength (dB re: ",m^2," at 1m)")))
  
  return(p1)
}
