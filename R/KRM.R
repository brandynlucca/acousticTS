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
    bar <- txtProgressBar(min=0, max=nrow(simdf), style=3, 
                          title="Calculating TS values...")
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
