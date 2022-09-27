#' Toggle between radians and degrees.
#'
#' @param x A real value in degrees or radians
#' @param d The value input-type. Two input types: "deg" for degrees and "rad"
#' for radians.
#' @usage
#' degrad(x,d)
#' @examples
#' x <- 180 #degrees
#' degrad(x, "deg")
#' # 3.141593
#' x <- pi #radians
#' degrad(x, "rad")
#' # 180
#' @return
#' Converts degrees to radians or radians to degrees
#' @export
degrad <- function(x, d){
  if(d == "deg"){
    value <- x*pi/180.0
  }else if(d == "rad"){
    value <- x*180.0/pi
  }
  return(value)
}

#' standardize <- function(object,
#'                         frequency,
#'                         ncyl_init=14,
#'                         ncyl_min=10,
#'                         phase_sd_init=sqrt(2)/2,
#'                         length_init=38.35e-3,
#'                         frequency_init=120e3){
#'   length_shape <- pull(object, "L")
#'   ncyl_shape <- pull(object, "ncyl")
#'   Nfl <- round((ncyl_shape*frequency*length_shape) / (frequency_init*length_init))
#'   Nfl[Nfl < ncyl_min] <- ncyl_min
#'   sigma_ph <- phase_sd_init * ((ncyl_init*length_shape)/(Nfl*length_init))
#'   return(data.frame(ncyl_b=Nfl, sigma_sd=sigma_ph))
#' }

