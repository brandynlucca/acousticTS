#' Calculates the theoretical TS of a shelled organism using the non-modal High Pass (HP) model
#' @param shape Desired object/animal shape.
#' @param frequency Frequency (Hz).
#' @param c Sound speed of surrounding medium (m/s). Default value is 1500 m/s.
#' @param a Radius of object (m).
#' @return
#' Target strength (TS, dB re: 1 m^2)
#' @references
#' @export

HP <- function(shape=NULL, c=1500, frequency, g=shape@g, h=shape@h, a=max(shape@a)){
  k <- kcalc(frequency, c)
  api <- (1-g*h^2)/(3*g*h^2) + (1-g)/(1+2*g)
  R <- (g*h-1)/(g*h+1)
  fbs <- a^2*(k*a)^4*api^2/((1+4*(k*a)^4*api^2)/R^2)
  TS <- 10*log10(fbs)
  return(TS)
}

