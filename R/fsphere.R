#' Calculates the theoretical TS of a fluid sphere using an exact modal series solution proposed by Andersen (1950).
#'
#' @param shape Desired animal object. Must be class FS.
#' @param c Sound speed of surrounding medium (m/s). Default value is 1500 m/s.
#' @param frequency Frequency (Hz).
#' @param a Radius vector of sphere (m).
#' @param h Sound speed contrast.
#' @param g Density contrast.
#' @param type Type of sphere being modeled. There are several options: 1 (fluid-filled), 2 (fixed-rigid),
#' 3 (pressure-release). Currently, fluid shells with fluid interiors are not available, but will be in future
#' versions.
#' @usage
#' fsphere(shape, frequency, c)
#'
#' fsphere(shape=NULL, frequency, c, a, g, h)
#' @details
#' Calculates the theoretical TS of a fluid sphere at a given frequency using an exact modal series solution.
#' @return
#' Target strength (TS, dB re: 1 m^2)
#' @references
#' \href{https://doi.org/10.1121/1.1906621}(Anderson, V.C. (1950). Sound scattering from a fluid sphere. Journal of the Acoustical Society of America, 22, 426-431. )
#' @export
#'

fsphere <- function(shape=NULL, frequency, c=1500, a=shape@a, g=shape@g, h=shape@h, method=1){
  k1 <- kcalc(frequency, c); k2 <- kcalc(frequency, c*h)
  m <- 0; fsub <- 0; t <- 0; f <- TRUE #initialize
  #Iterate through while loop
  while(f == TRUE){
    if(method == 1){
      cmn <- (jsd(m,k2*a)*yl(m,k1*a)) / (jl(m,k2*a)*jsd(m,k1*a)) - (g*h*ysd(m,k1*a)/jsd(m,k1*a)) #numerator term
      cmd <- (jsd(m,k2*a)*jl(m,k1*a)) / (jl(m,k2*a)*jsd(m,k1*a)) - (g*h)
      cm <- cmn/cmd
      bm <- -1 / (1 + 1i*cm)
    }else if(method == 2){
      bm <- -(jsd(m,k1*a)/hsd(m,k1*a))
    }else if(method == 3){
      bm <- -(jl(m,k1*a)/hl(m,k1*a))
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


#' Wrapper function that can simulate over distributions of values
#'
#' @param shape Desired animal object. Must be class FS.
#' @param c Sound speed of surrounding medium (m/s). Default value is 1500 m/s.
#' @param frequency Frequency (Hz).
#' @param a Radius vector of sphere (m).
#' @param h Sound speed contrast.
#' @param g Density contrast.
#' @param type Type of sphere being modeled. There are several options: 1 (fluid-filled), 2 (fixed-rigid),
#' 3 (pressure-release). Currently, fluid shells with fluid interiors are not available, but will be in future
#' versions.
#' @param aggregate Options to aggregate dataframe output into a series of summary statistics. Options include "mean", "median", "minimum", and "maximum".
#' @param permute Calculates for every permutation of variables.
#' @param parallel Boolean value that sets whether multicore CPU parallelization will be used to speed up calculations.
#' @param n.cores Number of CPU cores that will be dedicated to parallelizing model calculations.
#' @usage
#' fsphere.sim(shape, c=1500, frequency, permute=F,
#' aggregate=NULL, parallel=F, n.cores=NULL)
#'
#' fsphere.sim(shape=NULL, c=1500, frequency, h, g, a, permute=F, aggregate=NULL, parallel=F, n.cores=NULL)
#'
#' @details
#' Calculates the theoretical TS of a fluid sphere at a given frequency using an exact modal series solution.
#' @return
#' Target strength (TS, dB re: 1 m^2)
#' @references
#' \href{https://doi.org/10.1121/1.1906621}(Anderson, V.C. (1950). Sound scattering from a fluid sphere. Journal of the Acoustical Society of America, 22, 426-431. )
#' @export
#' @import foreach
#' @import doParallel
#' @import parallel
