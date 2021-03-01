#' Calls in a *.csv file as a FLS object
#'
#' @param file A *.csv file formatted with the following columns: x, y, z, a, g,
#'  h, theta [optional], pc [optional]
#' @usage
#' FLSread(file)
#' @return
#' Calls in an FLS-class object from a *.csv file
#' @export
#' @import methods
#' @import utils

FLSread <- function(file){
  animal <- read.csv(file, header=T) #Call in *.csv file; assumes headers are present
  return(new("FLS", rpos=as.matrix(rbind(animal$x,animal$y,animal$z)),
             a=animal$a,
             g=animal$g[1],h=animal$h[1],
             theta=ifelse(length(animal$theta) > 0, animal$theta, pi/2),
             curve=F,pc=0.0,L=max(animal$x),ncyl=length(animal$x)))}

#' Write FLS object as a .csv file.
#' @export
FLSwrite <- function(shape, 
                     filename=paste(getwd(),"/target_shape_",
                                    Sys.Date(),".csv",sep="")){
  object <- cbind(pos_matrix(shape), 
                  h=rep(shape@h,shape@ncyl),
                  g=rep(shape@g,shape@ncyl))
  write.csv(object, file=filename)
}

#' Calls in a *.csv file as a SBF
#'
#' @param file A *.csv file formatted with the following columns: xb, wb, zbU, 
#' zbL, xsb, wsb, zsbU, zsbL, pb, cb, psb, csb, theta [optional].
#' @usage
#' SBFread(file)
#' @return
#' Calls in an SBF-class object from a *.csv file
#' @export
SBFread <- function(file){
  animal <- read.csv(file, header=T)
  return(new("SBF",
             body=as.matrix(rbind(animal$xb, animal$wb, animal$zbU, 
                                  animal$zbL)),
             bladder=as.matrix(rbind(animal$xsb[!is.na(animal$xsb)], 
                                     animal$wsb[!is.na(animal$xsb)],
                                     animal$zsbU[!is.na(animal$xsb)], 
                                     animal$zsbL[!is.na(animal$xsb)])),
             pb=animal$pb[1], cb=animal$cb[1], psb=animal$psb[1], 
             csb=animal$csb[1],
             theta=ifelse(length(animal$theta) > 0, animal$theta, pi/2), 
             L=max(animal$xb),
             ncylb=length(animal$xb[!is.na(animal$xb)]),
             ncylsb=length(animal$xsb[!is.na(animal$xsb)])))}