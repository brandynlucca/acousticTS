#' Calculates the theoretical TS of a fluid-filled scatterer at a given frequency using the distorted Born wave approximation (DWBA) model.
#'
#'
#' @param shape Desired object/animal shape. Must be class "FFS".
#' @param c Sound speed of surrounding medium (m/s). Default value is 1500 m/s.
#' @param f Frequency (Hz).
#' @param phi Phase deviation (\eqn{\phi}), or phase variability. Accounts for complexities in animal shape and stochasticity of noise in scattering field.
#' Default value is 0.0.
#' @param theta Orientation of the target relative to the transmit source (\eqn{\theta}). Broadside incidence is considered 90 degrees, or pi/2.
#' Default value is pi/2; input should be in radians.
#' @param niterations Number of times/iterations the model will be ran.
#' @param summary Will output the
#' @usage
#' SDWBA(shape, c, f, phi, theta)
#' @details
#' Calculates the theoretical TS of a fluid-filled scatterer at a given frequency using the distorted Born wave approximation (DWBA) model.
#' @return
#' Target strength (TS, dB re: 1 m^2)
#' @references
#' Stanton, T.K., Chu, D., and Wiebe, P.H. 1998. Sound scattering by several zooplankton groups. II. Scattering models. Journal of the Acoustical Society of America, 103(1), 236-253.
#'
#' Hankin, R.K.S. 2006. Introducing elliptic, an R package for elliptic and modular functions. Journal of Statistical Software, 15(7).
#' @export

SDWBA <- function(shape=NULL, x=shape@rpos[1,], y=shape@rpos[2,], z=shape@rpos[3,],
                  c=1500, frequency, phase=0.0, a=shape@a, h=shape@h, g=shape@g,
                  pc=ifelse(curve == T, ifelse(is.null(shape),3.0,shape@pc),0.0),
                  theta=ifelse(is.null(shape),pi/2,shape@theta),
                  curve=ifelse(is.null(shape),F,shape@curve),
                  ncyl=ifelse(is.null(shape),length(x),shape@ncyl)){
  require(elliptic)
  rpos <- as.matrix(rbind(x,y,z))
  kt <- cbind(cos(theta),rep(0,length(theta)),sin(theta))
  k1 <- kcalc(frequency,c)*kt; k2 <- vecnorm(k1) / h
  fbs <- 0 + 0i

  for(j in 1:(ncyl-1)){
    r1 <- rpos[,j]; r2 <- rpos[,j+1]
    a1 <- a[j]; a2 <- a[j+1]
    beta <- abs(acos((k1%*%(r2-r1))/(vecnorm(k1)*vecnorm(r2-r1))) - pi/2)

    SDWBAint <- function(s){
      rint <- s * (r2-r1)+r1
      aint <- s * (a2-a1)+a1
      gamma <- 1/(g*h^2)+1/g-2

      if(abs(abs(beta)-pi/2)<1e-10){
        bessel <- k2*aint
      }else{
        bessel <- ja(1,2*k2*aint*cos(beta))/cos(beta)
      }

      if(curve == F){
        return(vecnorm(k1)/4*gamma*aint*exp(2i*k1%*%rint/h)*bessel*vecnorm(r2-r1))
      }else{
        pc <- pc*max(x)
        return(vecnorm(k1)*pc/4*gamma*aint*exp(2i*k2*pc)*exp(-2i*k2*pc*cos(beta))*bessel*(vecnorm(r2-r1)/pc))
      }
    }

    SDWBAint <- Vectorize(SDWBAint)
    integral <- myintegrate(SDWBAint,0,1)
    fbs <- fbs + integral * exp(1i * rnorm(1,0,phase))
  }
  return(20*log10(abs(fbs)))
}

#Wrapper function that can simulate over distributions of values
#' @export

SDWBA.sim <- function(shape=shape, x=shape@rpos[1,], y=shape@rpos[2,], z=shape@rpos[3,],
                      c=1500, frequency, phase=0.0, a=shape@a, h=shape@h, g=shape@g,
                      pc=ifelse(curve == T, ifelse(is.null(shape),3.0,shape@pc),0.0),
                      theta=ifelse(is.null(shape),pi/2,shape@theta),
                      curve=ifelse(is.null(shape),F,shape@curve),
                      length=ifelse(is.null(shape),max(x),shape@L),
                      nrep=NULL, aggregate=NULL, parallel=F, n.cores=NULL){
  if(!is.null(nrep)){
    repseq <- seq(1,nrep,1)
  }else{
    repseq <- 1
  }
  simdf <- expand.grid(iteation=repseq, c=c, frequency=frequency, g=g, h=h, pc=pc, theta=theta, curve=curve, phase=phase, length=length, TS=NA)

  if(parallel==F){
    for(i in 1:nrow(simdf)){
      target_sim <- Shapely(shape,curve=simdf$curve[i],pc=simdf$pc[i],theta=simdf$theta[i],length=simdf$length[i])
      simdf$TS[i] <- SDWBA(target_sim,c=c,frequency=frequency,phase=simdf$phase[i],g=simdf$g[i],h=simdf$h[i])
    }
  }else if(parallel==T){
    require(foreach)
    require(parallel)
    require(doParallel)
    if(!is.null(n.cores)){
      n.cores <- ncores
    }else{
      n.cores <- detectCores()
    }
    cl <- makeCluster(n.cores)
    registerDoParallel(cl)

    simdf$TS <- foreach(i=1:nrow(simdf), .combine=c) %dopar% {
      target_sim <- Shapely(shape,curve=simdf$curve,pc=simdf$pc,theta=simdf$theta,length=simdf$length)
      SDWBA(target_sim,c=simdf$c[i],frequency=simdf$frequency[i],phase=simdf$phase[i],g=simdf$g[i],h=simdf$h[i])
    }

    stopCluster(cl)
  }

  if(!is.null(aggregate)){
    dum <- data.frame(Stat=NA,TS=NA)
    if("mean" %in% aggregate){
      dum <- rbind(dum,data.frame(Stat="Mean", TS=10*log10(mean(10^(simdf$TS/10)))))
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
    return(simdf)
  }
}
