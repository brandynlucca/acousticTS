#' Calculates the theoretical TS of a fluid-filled scatterer at a given frequency using the deformed finite cylinder model (DFCM).
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
