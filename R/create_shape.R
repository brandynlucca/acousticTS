#' A wrapper function that automatically creates generalized and/or canonical
#' shapes for TS modeling. 
#'
#' @param L Length or diameter of an object (m). 
#' @param N Number of cylinders. 
#' @param shape Shape. Details for shape specification are given under 'Details', 
#' including mandatory additional arguments.
#' 
#' @param g Optional density contrast input.
#' @param h Optional sound speed contrast input
#' @param theta Optional orientation input. 
#' @param curve Optional shape curvature input [T/F].
#' @param pc Optional radius of curvature ratio input. 
#' @param object Optional scatterer class object.  
#' 
#' @details 
#' The \strong{shape} argument specifies what shape for the function to generate into 
#' the desired shape for TS modeling. Options currently include:
#' \tabular{rlllll}{
#'  \tab \strong{Object Shape} \tab \strong{shape=...} \tab  \tab 
#'  \strong{Parameters} \tab \strong{Root function}\cr
#'  \tab \emph{Discrete Cylinder} \tab "cyl" \tab \tab LAratio \tab
#'  \code{\link[=cylinder]{cylinder(...)}}\cr
#'  \tab \emph{Tapered Cylinder} \tab "tcyl" \tab \tab LAratio, taper \tab 
#'  \code{\link[=tapered_cylinder]{tapered_cylinder(...)}}\cr
#'  \tab \emph{Polynomial Cylinder} \tab "pcyl" \tab \tab LAratio, poly \tab 
#'  \code{\link[=polynomial_cylinder]{polynomial_cylinder(...)}}\cr
#'  \tab \emph{Prolate Spheroid} \tab "psph" \tab \tab LAratio \tab 
#'  \code{\link[=prolate_spheroid]{prolate_spheroid(...)}}\cr 
#'  \tab \emph{Sphere} \tab "sph" \tab \tab None \tab 
#'  \code{\link[=sphere]{sphere(...)}}\cr
#' }
#' 
#' \subsection{Model Parameter Definitions}{
#' \itemize{
#'  \item \strong{LAratio}: the length-to-radius ratio (L/A), which specifically
#'  refers to the radius at the mid-point of the cylinder and should be 
#'  the maximum value. A typical L/A ratio in the literature is 16 for krill.
#'  \item \strong{taper}: the taper order (n), which parameterizes the tapering
#'  function reported by Chu \emph{et al.} (1993) to create a tapered cylinder.
#'  The tapering order will converge on a prolate and oblate spheroid when
#'  L > 2a and L < 2a, respectively, and n = 2. A typical taper order in the 
#'  literature is 10.
#'  \item \strong{poly}: the vector of arbitrary polynomial coefficients to 
#'  generate a deformed cylinder as reported by Smith \emph{et al.} (2013). 
#'  Although listed as a mandatory argument for the polynomial cylinder 
#'  function, it has a default setting that uses the sixth-degree polynomial 
#'  coefficients reported by Smith \emph{et al.} (2013). 
#' }
#' }
#'  
#' @return
#' Creates a generalized and/or canonical shape for TS modeling. 
#' 
#' @references
#' Chu, D., Foote, K.G., and Stanton, T.K. 1993. Further analysis of target 
#' strength measurements of Antarctic krill at 38 and 120 kHz: Comparison and 
#' deformed cylinder model and inference of orientation distribution. The 
#' Journal of the Acoustical Society of America, 93(5): 2985-2988. 
#' https://doi.org/10.1121/1.405818
#' 
#' Smith, J.N., Ressler, P.H., and Warren, J.D. 2013. A distorted wave Born 
#' approximation target strength model for Bering Sea euphausiids. ICES Journal
#' of Marine Science, 70(1): 204-214. https://doi.org/10.1093/icesjms/fss140
#' @export
create_shape <- function(L, N, shape, ...){
  params <- names(as.list(substitute(list(...))))
  if(shape == "cyl"){
    if(!all(c("LAratio") %in% params))
      stop("A discrete cylinder requires the 'LAratio' argument.")
    cylinder(L,N, ...)
  }
  else if(shape == "tcyl"){
    if(!all(c("LAratio","taper") %in% params))
      stop("A tapered cylinder requires both the 'taper' and 'LAratio' arguments.")
    taper_cylinder(L,N, ...)
  }else if(shape == "pcyl"){
    if(!all(c("LAratio") %in% params))
      stop("A polynomial cylinder requires the 'LAratio' argument.")
    polynomial_cylinder(L,N, ...)
  }else if(shape == "psph"){
    if(!all(c("LAratio") %in% params))
      stop("A polynomial cylinder requires the 'LAratio' argument.")
    prolate_spheroid(L,N, ...)
  }else if(shape == "sph"){
    sphere(L,N, ...)
  }else{
    stop("Defined shape does not exist. Please select from the available options ",
          "in the ?create_sphere help documentation.")
  }
}


#' Support function to generate tapered cylinder. 
#'
#' @param L Length of object (m).
#' @param a Maximum radius of object (m).
#' @param LAratio Optional. Length:radius ratio.
#' @param g Density contrast.
#' @param h Sound speed contrast. 
#' @param theta Orientation. 
#' @param curve Curvature [T/F].
#' @param pc Radius of curvature ratio. 
#' @param N Number of cylinders. 
#' @param object Object class. 
#' @return
#' Creates a tapered cylinder. 
#' @export
cylinder <- function(L, LAratio, N,
                     a=NULL,
                           g=1.0357, h=1.0279, 
                           theta=pi/2, pc=3.0, curve=F,
                           object="FLS"){
  x <- seq(0,L,length.out=N+1)
  
  if(is.null(a)){
    a <- rep(L/LAratio,N+1)
  }else{
    a <- rep(a,N+1)
  }

  if(object == "FLS"){
    shape <- FLSgenerate(x=x, y=rep(0,N+1), z=rep(L/LAratio*2,N+1),
                         a=a, g=g, h=h, theta=theta, pc=pc)
  }
  return(shape)
}

#' Support function to generate tapered cylinder. 
#'
#' @param L Length of object (m).
#' @param a Maximum radius of object (m).
#' @param LAratio Length:radius ratio [optional].
#' @param g Density contrast.
#' @param h Sound speed contrast. 
#' @param theta Orientation. 
#' @param curve Curvature [T/F].
#' @param pc Radius of curvature ratio. 
#' @param N Number of cylinders. 
#' @param taper Degree of taper. 
#' @param object Object class. 
#' @return
#' Creates a tapered cylinder. 
#' @export
taper_cylinder <- function(L, LAratio, N, taper,
                           g=1.0357, h=1.0279, 
                           theta=pi/2, pc=3.0, curve=F,
                           object="FLS"){
  x <- seq(-1, 1, length.out=N+1) #along-body axis
  t <- sqrt(1-x^taper) #tapered radii
  
  if(object == "FLS"){
    shape <- FLSgenerate(x=x*L/2+L/2, y=rep(0,N+1), z=rep(L/LAratio*2,N+1),
                         a=L/LAratio*taper, g=g, h=h, theta=theta, pc=pc)
  }
  return(shape)
}

#' Support function to generate polynomial deformed cylinder. 
#'
#' @param L Length of object (m).
#' @param LAratio Length:radius ratio.
#' @param g Density contrast.
#' @param h Sound speed contrast. 
#' @param theta Orientation. 
#' @param curve Curvature [T/F].
#' @param pc Radius of curvature ratio. 
#' @param N Number of cylinders. 
#' @param poly Vector of polynomial coefficients.
#' @param object Object class. 
#' @return
#' Creates a polynomial deformed cylinder. 
#' @export
polynomial_cylinder <- function(L, LAratio, N,
                      g=1.0357, h=1.0279, 
                      theta=pi/2, pc=3.0, curve=F,
                      poly=c(0.83,0.36,-2.1,-1.2,0.63,0.82,0.64),
                      object="FLS"){
  x <- seq(-1, 1, length.out=N+1) #along-body axis
  n_order <- rev(seq_len(length(poly)))-1
  poly_fun <- paste(poly, paste("*x^",n_order, sep=""),collapse="+",sep="")
  a <- abs(eval(parse(text=poly_fun)))

  if(object == "FLS"){
    shape <- FLSgenerate(x=x*L/2+L/2, y=rep(0,N+1), 
                         z=rep(L/LAratio*2,N+1),
                         a=a*L/LAratio, g=g, h=h, theta=theta, pc=pc)
  }
  return(shape)
}

#' Support function to generate prolate spheroid. 
#'
#' @param L Length of object (m).
#' @param LAratio Length:radius ratio.
#' @param g Density contrast.
#' @param h Sound speed contrast. 
#' @param theta Orientation. 
#' @param curve Curvature [T/F].
#' @param pc Radius of curvature ratio. 
#' @param N Number of cylinders. 
#' @param object Object class. 
#' @return
#' Creates a prolate spheroid. 
#' @export
prolate_spheroid <- function(L, LAratio, N, 
                             g=1.0357, h=1.0279, 
                             theta=pi/2, pc=3.0, curve=F,
                             object="FLS"){
  x <- seq(0,L,length.out=N+1)
  a <- (L/LAratio)*sqrt(1-((x-L/2)/(L/2))^2)
  
  if(object == "FLS"){
    shape <- FLSgenerate(x=x, y=rep(0,N+1), 
                         z=rep(L/LAratio*2,N+1),
                         a=a, g=g, h=h, theta=theta, pc=pc)
  }
  return(shape)
}

#' Support function to generate sphere. 
#'
#' @param L Length of object (m).
#' @param LAratio Length:radius ratio.
#' @param g Density contrast.
#' @param h Sound speed contrast. 
#' @param theta Orientation. 
#' @param curve Curvature [T/F].
#' @param pc Radius of curvature ratio. 
#' @param N Number of cylinders. 
#' @param object Object class. 
#' @return
#' Creates a sphere. 
#' @export
sphere <- function(L, N, 
                   g=1.0357, h=1.0279, 
                   theta=pi/2, pc=0.0, curve=F,
                   object="FLS"){
  x <- seq(0,L,length.out=N+1)
  a <- sqrt((L/2)^2-(x-L/2)^2)
  
  if(object == "FLS"){
    shape <- FLSgenerate(x=x, y=rep(0,N+1), 
                         z=rep(max(a)*2,N+1),
                         a=a, g=g, h=h, theta=theta, pc=pc,)
  }
  return(shape)
}

#' Generate SBF shape
#' @export
SBFgenerate <- function(xb,wb,zbU,zbL,xsb,wsb,zsbU,zsbL,pb,cb,psb,csb,theta=pi/2){
  return(new("SBF", body=as.matrix(rbind(xb,wb,zbU,zbL)), bladder=as.matrix(rbind(xsb,wsb,zsbU,zsbL)),
             pb=pb, cb=cb, psb=psb, csb=csb, theta=theta, L=max(xb), ncylb=length(xb), ncylsb=length(xsb)))}

#' Generate CAL shape
#' @param material Material-type for the soldi sphere. See 'Details' built-in
#' material options. 
#' @param a Spherical radius (m).
#' @param c1 Longitudinal sound speed (m/s).
#' @param c2 Transversal sound speed (m/s).
#' @param rho1 Density (kg/m^3)
#' @export
CALgenerate <- function(material="WC", a=38.1e-3, 
                        c1=NULL, c2=NULL, rho1=NULL){
  if(material == "Cu"){
    c1=4760; c2=2288.5; rho1=8947
  }else if(material == "steel"){
    c1=5980; c2=3297; rho1=7970
  }else if(material == "Al"){
    c1=6260; c2=3080; rho1=2700
  }else if(material == "WC"){
    c1=6853; c2=4171; rho1=14900
  }else if(material == "brass"){
    c1=4372; c2=2100; rho1=8360
  }
  return(new("CAL", material=material, 
             a=a, 
             c1=c1, c2=c2, rho1=rho1))
}

#' Manually generate a FLS object.
#'
#' @param x Vector containing x-axis body (m) shape data.
#' @param y Vector containing y-axis body (m) shape data.
#' @param z Vector containing z-axis body (m) shape data.
#' @param a Vector containing radii (m).
#' @param g Density contrast.
#' @param h Soundspeed contrast
#' @param theta Orientation of the target relative to the transmit source 
#' (\eqn{\theta}). Broadside incidence is considered 90 degrees, or pi/2.
#' Default value is pi/2; input should be in radians.
#' @usage
#' FLSgenerate(x,y,z,a,g,h,theta,curve,pc)
#' @examples
#' #Manually parameterize shape
#' x <- seq(1,10,1)*1e-3; y <- rep(0,10); z <- c(seq(1,5,1),rev(seq(1,5,1)))*1e-4
#' a <- z/2
#' g <- 1.036
#' h <- 1.0279
#' new_target <- FLSgenerate(x=x,y=y,z=z,a=a,g=g,h=h)
#' #Let's model where sound speed (c) is 1500 m/s, frequency is 120 kHz, with no 
#' #phase deviation
#' c <- 1500
#' freq <- 120e3
#' SDWBA(shape=new_target, c=c, frequency=freq)
#' #-114.0107
#' @return
#' Calls in an FLS-class object from a *.csv file
#' @export
#' @import methods
FLSgenerate <- function(x,y,z,a,g,h,theta=pi/2,curve=F,pc=0.0){
  return(new("FLS", rpos=as.matrix(rbind(x,y,z)),a=a,g=g,h=h,theta=theta,
             curve=curve,pc=pc,L=max(x),ncyl=length(x)))}

#' Generate ESS shape
#' @param a Radius (m).
#' @param g Density contrast. 
#' @param h Sound speed contrast.
#' @param theta Orientation (\eqn{\theta}, radians)
#' 
#' @param L Optional. Length of the longitudinal axis (m) for non-spherical 
#' shapes.
#' @param spherical Optional. Boolean for specifying whether scatterer is 
#' spherical or oblong. 
#' @export
ESSgenerate <- function(a,g,h,theta=pi/2,L=numeric(),spherical=TRUE){
  return(new("ESS", a=a, g=g, h=h, theta=theta, 
             L=ifelse(spherical==TRUE,numeric(),L),
             spherical=spherical))
}
