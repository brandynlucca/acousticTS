#' Calculates the theoretical TS of a shelled organism using the non-modal High Pass (HP) model
#'
#' @usage
#' @details
#' @return
#' Target strength (TS, dB re: 1 m^2)
#' @references
#' @export

HP <- function(animal, c=1500, frequency, a=max(animal@a)){
  k <- kcalc(frequency, c)
  alphapi <- (1-animal@g*animal@h^2)/3*animal@g*animal@h^2+(1-animal@g)/(1+2*animal@g)
  R <- (animal@g*animal@h-1)/(animal@g*animal@h+1)
  fbs <- a^2*(k*a)^4*alphapi^2/(1+4*(k*a)^4*alphapi^2/R^2)
  TS <- 10*log10(fbs)
  return(TS)
}
