#' Calculates the theoretical TS using Kirchoff-ray Mode approximation.
#'
#' @param shape Desired object/animal shape. Must be class "SBF".
#' @param c Sound speed of surrounding medium (m/s). Default value is 1500 m/s.
#' @param rho Density of medium (kg m^-1). Default value is 1030 kg m^-1.
#' @param frequency Frequency (Hz).
#' @param theta Orientation of the target relative to the transmit source (\eqn{\theta}). Broadside incidence is considered 90 degrees, or pi/2.
#' Default value is pi/2; input should be in radians.
#' @usage
#' KRM(shape, c, rho, frequency, theta)
#' @details
#' Calculates the theoretical TS using the Kirchoff-ray Mode model.
#' @return
#' Target strength (TS, dB re: 1 m^2)
#' @references
#' @export
#' @import elliptic
#' @import foreach
#' @import doParallel
#' @import parallel
KRM <- function(shape=NULL, c=1500, rho=1030, frequency, theta=pi/2){
  shape@theta <- theta
  k <- kcalc(frequency, c); kb <- kcalc(frequency, shape@cb)
  RBC <- (shape@psb*shape@csb - shape@pb*shape@cb) / (shape@psb*shape@csb + shape@pb*shape@cb)
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
      -1i*(RBC*TT)/(2*sqrt(pi))*Asb*sqrt((k*as+1)*sin(shape@theta))*exp(-1i*(2*k*vs+Psip))*delus
  }

  for(i in 1:(length(shape@body[1,])-1)){
    b1 <- shape@body[,i]; b2 <- shape@body[,i+1]
    ab <- (b1[2]+b2[2])/4
    Psib <- -pi*kb*((b1[3]+b2[3])/2)/(2*(kb*((b1[3]+b2[3])/2)+0.4))
    delub <- (b2[1]-b1[1])*sin(shape@theta)
    vbU <- ((b1[3]+b2[3])*cos(shape@theta) + (b1[3]+b2[3])*sin(shape@theta))/2
    vbL <- ((b1[4]+b2[4])*cos(shape@theta) + (b1[4]+b2[4])*sin(shape@theta))/2

    frequency.fluid <- frequency.fluid +
      -1i*(RWB/(2*sqrt(pi)))*sqrt(k*ab)*delub*(exp(-2i*k*vbU)-TT*exp(-2i*k*vbU+2i*kb*(vbU-vbL)+1i*Psib))

  }
  sigma <- abs(frequency.fluid + frequency.soft)
  TS <- 20*log10(sigma)
  return(TS)
}

#' Wrapper function to simulate KRM runs.
#'
#' @export
KRM.sim <- function(shape=shape, c=1500, rho=1030, frequency,
                    theta=shape@theta,
                    pb=shape@pb, cb=shape@cb, length=shape@L,
                    permute=T, progress=T, parallel=F, n.cores=NULL){
  if(permute == T){
    simdf <- expand.grid(c=c, frequency=frequency, pb=pb, cb=cb,
                         theta=theta, length=length, TS=NA)
  }

  if(progress == T){
    bar <- txtProgressBar(min=0, max=nrow(simdf), style=3, title="Calculating TS values...")
    progbar <- function(n) setTxtProgressBar(bar, n)
  }

  if(parallel==F){
    for(i in 1:nrow(simdf)){
      lscale <- simdf$length[i]/shape@L
      mscale <- cbind(c(1,0,0,0), c(0,1,0,0), c(0,0,1,0), c(0,0,0,1)) * lscale
      target_sim <- shape
      target_sim@body <- t(t(target_sim@body) %*% mscale)
      target_sim@bladder <- t(t(target_sim@bladder) %*% mscale)
      target_sim@L <- max(target_sim@body[1,])
      target_sim@theta <- simdf$theta[i]
      target_sim@pb <- simdf$pb[i]
      target_sim@cb <- simdf$cb[i]
      simdf$TS[i] <- KRM(target_sim, c=simdf$c[i], frequency=simdf$frequency[i],
                         rho=rho, theta=simdf$theta[i])
       if(progress == T){
        progbar(i)
      }
    }
  }else if(parallel==T){
    requireNamespace("parallel", quietly=T)
    requireNamespace("doSNOW", quietly=T)
    requireNamespace("snow", quietly=T)

    if(!is.null(n.cores)){
      n.cores <- n.cores
    }else{
      n.cores <- detectCores() - 1
    }
    cl <- snow::makeCluster(n.cores)
    registerDoSNOW(cl)

    if(progress == T){
      bar <- txtProgressBar(min=0, max=nrow(simdf), style=3)
      opts <- list(progress=progbar)
    }else{
      opts <- list()
    }

    simdf$TS <- foreach(
      i=seq_len(nrow(simdf)), .combine=c,
      .packages=c("acousticTS"),
      .options.snow=opts) %dopar% {
        lscale <- simdf$length[i]/shape@L
        mscale <- cbind(c(1,0,0,0), c(0,1,0,0), c(0,0,1,0), c(0,0,0,1)) * lscale
        target_sim <- shape
        target_sim@body <- t(t(target_sim@body) %*% mscale)
        target_sim@bladder <- t(t(target_sim@bladder) %*% mscale)
        target_sim@L <- max(target_sim@body[1,])
        target_sim@theta <- simdf$theta[i]
        target_sim@pb <- simdf$pb[i]
        target_sim@cb <- simdf$cb[i]
        return(KRM(target_sim, c=simdf$c[i], frequency=simdf$frequency[i],
                           rho=rho, theta=simdf$theta[i]))
      }

    stopCluster(cl)

    if(progress == T){
      close(bar)
    }
  }
  return(simdf)

}
