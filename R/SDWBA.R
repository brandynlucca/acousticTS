#' Wrapper function that can simulate over distributions of values
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
#' @param phase Phase deviation (\eqn{\phi}), or phase variability. Accounts 
#' for complexities in animal shape and stochasticity of noise in scattering 
#' field.
#' Default value is 0.0.
#' @param theta Orientation of the target relative to the transmit source 
#' (\eqn{\theta}). Broadside incidence is considered 0 degrees, or 0 radians. 
#' Animals are
#' either tilted anterior-end down, which corresponds to (\eqn{-\pi/2}) 
#' and (\eqn{pi/2}) respectively. The model will automatically convert these 
#' orientations
#' to the orientation definition used by McGehee et al. (1998) whereby 
#' orientation is based on the angle of the incident soundwave such that
#' broadside incidence is \eqn{\pi/2} (dorsal-up) and \eqn{3\pi/2} 
#' (dorsal-down), or when the transmitted soundwave is completely perpendicular 
#' to the target.
#' Default value is pi/2; input should be in radians.
#' @param length Option to change the length of the scatterer shape.
#' @param nrep Number of repeated iterations to run the model.
#' @param scale Log or linear domain. Log is the default, and outputs TS 
#' (dB re: 1 m^2). Linear represent \eqn{\sigma_bs} (m^2), which is the
#' backscattering cross-section. Values must be in the linear domain as 
#' \eqn{\sigma_bs} if any averaging is going to be conducted on the output of
#' this model. Inputs allowed are "linear", "log", or "both".
#' @param aggregate Options to aggregate dataframe output into a series of 
#' summary statistics. Options include "mean", "median", "minimum", and "maximum".
#' @param permute Calculates for every permutation of variables.
#' @param parallel Boolean value that sets whether multicore CPU parallelization 
#' will be used to speed up calculations.
#' @param n.cores Number of CPU cores that will be dedicated to parallelizing 
#' model calculations.
#' @usage
#' SDWBA.sim(shape, c=1500, frequency, phase=0.0, curve=F, nrep=NULL, permute=F,
#' aggregate=NULL, parallel=F, n.cores=NULL)
#'
#' SDWBA.sim(shape, c=1500, frequency, phase=0.0, h, g, curve=F, pc, theta, 
#' length,
#' nrep=NULL, permute=F, aggregate=NULL, parallel=F, n.cores=NULL)
#'
#' SDWBA.sim(shape, c=1500, frequency, phase=0.0, x, y, z, a, h, g, curve=F, 
#' pc, theta, length,
#' nrep=NULL, permute=F, aggregate=NULL, parallel=F, n.cores=NULL)
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
#' @export
#' @import foreach
#' @import doSNOW
#' @import snow
#' @importFrom elliptic myintegrate
#' @importFrom parallel detectCores
#' @import stats
#' @import methods

SDWBA.sim <- function(shape=shape, 
                      x=shape@rpos[1,], y=shape@rpos[2,], z=shape@rpos[3,],
                      c=1500, frequency, phase=0.0, a=shape@a, 
                      h=shape@h, g=shape@g,
                      curve=shape@curve,
                      pc=shape@pc,
                      theta=shape@theta,
                      length=shape@L,
                      scale="log",
                      nrep=NULL, permute=F, aggregate=NULL, parallel=F, 
                      n.cores=NULL,
                      progress=T){

  if(!is.null(nrep)){
    repseq <- seq(1,nrep,1)
  }else{
    repseq <- 1
  }

  if(permute == T){
    simdf <- expand.grid(iteration=repseq, c=c, frequency=frequency,
                         g=g, h=h, theta=theta, pc=pc, curve=curve,
                         phase=phase, length=length, TS=NA)
  }else{
    baselen <- data.frame(param=c("repseq", "c", "frequency"),
                          len=c(length(repseq), length(c), length(frequency)))
    baseprod <- prod(baselen$len)
    paramlen <- data.frame(param=c("h","g","curve","pc",
                                   "theta","length","phase"),
                           len=c(length(h),length(g),length(curve),length(pc),
                                 length(theta),length(length),length(phase)))
    simdf <- data.frame(iterator=seq_len(baseprod), c=c, frequency=frequency)
    simdf <- cbind(simdf,
                   data.frame(g=rep(g, baseprod), h=rep(h, baseprod),
                              theta=rep(theta, baseprod), 
                              pc=rep(pc, baseprod),
                              curve=rep(curve, baseprod), 
                              length=rep(length, baseprod),
                              phase=rep(phase, baseprod),
                              TS=NA))
    simdf <- simdf[,-1]
  }

  if(progress == T){
    bar <- txtProgressBar(min=0, max=nrow(simdf), style=3, 
                          title="Calculating TS values...")
    progbar <- function(n) setTxtProgressBar(bar, n)
  }


  if(parallel==F){
    for(i in 1:nrow(simdf)){
      target_sim <- Shapely(shape,curve=simdf$curve[i],pc=simdf$pc[i],
                            theta=simdf$theta[i],length=simdf$length[i])
      simdf$TS[i] <- SDWBA(target_sim,c=simdf$c[i],frequency=simdf$frequency[i],
                           phase=simdf$phase[i],g=simdf$g[i],h=simdf$h[i],
                           curve=simdf$curve[i],pc=simdf$pc[i])

      if(progress == T){
        progbar(i)
      }
    }
  }else if(parallel==T){
    requireNamespace("doSNOW", quietly=T)
    requireNamespace("snow", quietly=T)

    if(!is.null(n.cores)){
      n.cores <- n.cores
    }else{
      n.cores <- parallel::detectCores() - 1
    }
    cl <- snow::makeCluster(n.cores)
    registerDoSNOW(cl)

    if(progress == T){
      bar <- txtProgressBar(min=0, max=nrow(simdf), style=3)
      opts <- list(progress=progbar)
    }else{
      opts <- list()
    }

    simdf$TS<- foreach(
      i=seq_len(nrow(simdf)), .combine=c,
      .packages=c("acousticTS","elliptic"),
      .options.snow=opts) %dopar% {
        target_sim <- Shapely(shape,curve=simdf$curve[i],pc=simdf$pc[i],
                              theta=simdf$theta[i],length=simdf$length[i])

        return(SDWBA(target_sim,c=simdf$c[i],frequency=simdf$frequency[i],
                     phase=simdf$phase[i],g=simdf$g[i],h=simdf$h[i],
                     curve=simdf$curve[i],theta=simdf$theta[i],pc=simdf$pc[i]))
      }

    stopCluster(cl)

    if(progress == T){
      close(bar)
    }
  }

  if(!is.null(aggregate)){
    dum <- data.frame(Stat=NA,TS=NA)
    if("mean" %in% aggregate){
      dum <- rbind(dum,data.frame(Stat="Mean", 
                                  TS=10*log10(mean(10^(simdf$TS/10)))))
    }
    if("median" %in% aggregate){
      dum <- rbind(dum,data.frame(Stat="Median", TS=median(simdf$TS)))
    }
    if("minimum" %in% aggregate){
      dum <- rbind(dum,data.frame(Stat="Minimum", TS=min(simdf$TS)))
    }
    if("maximum" %in% aggregate){
      dum <- rbind(dum,data.frame(Stat="Maximum", TS=max(simdf$TS)))
    }
    dum <- dum[-(is.na(dum$Stat)),]
    return(dum)
  }else{
    if(scale %in% c("linear","both")){
      simdf$sigmabs <- 10^(simdf$TS/10)
      if(scale == "linear"){
        simdf <- simdf[-(which(colnames(simdf) == "TS"))]
      }
    }
    return(simdf)
  }
}
