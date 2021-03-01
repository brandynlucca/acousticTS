#' Wrapper function that can simulate over distributions of values
#'
#' @param shape Desired animal object. Must be class FS.
#' @param c Sound speed of surrounding medium (m/s). Default value is 1500 m/s.
#' @param frequency Frequency (Hz).
#' @param a Radius vector of sphere (m).
#' @param h Sound speed contrast.
#' @param g Density contrast.
#' @param type Type of sphere being modeled. There are several options: 1 
#' (fluid-filled), 2 (fixed-rigid),
#' 3 (pressure-release). Currently, fluid shells with fluid interiors are not 
#' available, but will be in future
#' versions.
#' @param aggregate Options to aggregate dataframe output into a series of 
#' summary statistics. Options include "mean", "median", "minimum", and 
#' "maximum".
#' @param permute Calculates for every permutation of variables.
#' @param parallel Boolean value that sets whether multicore CPU parallelization 
#' will be used to speed up calculations.
#' @param n.cores Number of CPU cores that will be dedicated to parallelizing 
#' model calculations.
#' @usage
#' fsphere.sim(shape, c=1500, frequency, permute=F,
#' aggregate=NULL, parallel=F, n.cores=NULL)
#'
#' fsphere.sim(shape=NULL, c=1500, frequency, h, g, a, permute=F, 
#' aggregate=NULL, parallel=F, n.cores=NULL)
#'
#' @details
#' Calculates the theoretical TS of a fluid sphere at a given frequency using an 
#' exact modal series solution.
#' @return
#' Target strength (TS, dB re: 1 m^2)
#' @references
#' \href{https://doi.org/10.1121/1.1906621}(Anderson, V.C. (1950). 
#' Sound scattering from a fluid sphere. Journal of the Acoustical Society of 
#' America, 22, 426-431. )
#' @export

