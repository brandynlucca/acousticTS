#Write KRM function
KRM <- function(animal, c=1490, rho=1030, f, tilt=pi/2){
  k <- kcalc(f, c); kb <- kcalc(f, animal@cb)
  RBC <- (animal@psb*animal@csb - animal@pb*animal@cb) / (animal@psb*animal@csb + animal@pb*animal@cb)
  RWB <- (animal@pb*animal@cb - rho*c) / (animal@pb*animal@cb + rho*c)
  TT <- 1 - RWB^2
  f.soft <- 0i; f.fluid <- 0i
  
  for(i in 1:(length(animal@bladder[1,])-1)){
    p1 <- animal@bladder[,i]; p2 <- animal@bladder[,i+1]
    as <- (p1[2]+p2[2])/4
    Asb <- k*as/(k*as+0.083)
    Psip <- k*as/(40+k*as)-1.05
    vs <- ((p1[1]+p2[1])*cos(animal@theta) + (p1[3]+p2[3])*sin(animal@theta))/2
    delus <- (p2[1]-p1[1])*sin(animal@theta)
    f.soft <- f.soft + 
      -1i*(RBC*TT)/(2*sqrt(pi))*Asb*sqrt((k*as+1)*sin(animal@theta))*exp(-1i*(2*k*vs+Psip))*delus
  }
  
  for(i in 1:(length(animal@body[1,])-1)){
    b1 <- animal@body[,i]; b2 <- animal@body[,i+1]
    ab <- (b1[2]+b2[2])/4
    Psib <- -pi*kb*((b1[3]+b2[3])/2)/(2*(kb*((b1[3]+b2[3])/2)+0.4))
    delub <- (b2[1]-b1[1])*sin(animal@theta)
    vbU <- ((b1[3]+b2[3])*cos(animal@theta) + (b1[3]+b2[3])*sin(animal@theta))/2
    vbL <- ((b1[4]+b2[4])*cos(animal@theta) + (b1[4]+b2[4])*sin(animal@theta))/2
    
    f.fluid <- f.fluid +
      -1i*(RWB/(2*sqrt(pi)))*sqrt(k*ab)*delub*(exp(-2i*k*vbU)-TT*exp(-2i*k*vbU+2i*kb*(vbU-vbL)+1i*Psib))
    
  }
  sigma <- abs(f.fluid + f.soft)
  TS <- 20*log10(sigma)
  return(TS)
}