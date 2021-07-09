# THERMAL or kappa SENSITIVITY OF LAMBDA ####

# Sensitivity analysis to study how change in Lambda can be
# attributed to change in different vital rates originating from 
# a change in temperature. I will then test for a few kappa 
# values at the same temperature to understand how 
# kappa affects lambda

# Sensitivity based on Vindenes et al. 2014 ####

# GROWTH SESITIVITY 
#Sensitivity of the growth density function to length (approximate) 
# dG.dx <- function(x, z=10.34){
#   (growthfun(x=x+0.001, y=y, z=z) - growthfun(x=x, y=y, z=z))/0.001
# }
dG.dm <- function(m, Pars){
  (DEBgrowthfun(m+0.001, Pars) - DEBgrowthfun(m, Pars))/0.001
}

#Sensitivity of the growth density function to temperature (approximate) 
# dG.dz <- function(y, x, z=10.34){
#   (growthfun(x=x, y=y, z=z+.001) - growthfun(x=x, y=y, z=z))/.001
# }

dG.dT <- function(m, Pars){
  (DEBgrowthfun(m, c(T=Pars[["T"]]+0.001, 
                     kappa =Pars[["kappa"]], 
                     Y = Pars[["Y"]])) 
   - DEBgrowthfun(m, Pars))/0.001
}

##SURVIVAL SENSITIVITY
#Sensitivity of the survival probability function to temperature (approximate) 
# dS.dz <- function(x, z=10.34){
#   (survfun(x=x, z=z+.001) - survfun(x=x, z=z))/0.001
# }

dS.dT <- function(m, Pars){
  (DEBsurvfun(m, c(T=Pars[["T"]]+0.001, 
                   kappa =Pars[["kappa"]], 
                   Y = Pars[["Y"]])) 
   - DEBsurvfun(m, Pars))/0.001
}

#Sensitivity of the survival probability function to length (approximate) 
# dS.dx <- function(x, z=10.34){
#   (survfun(x=x+0.001, z=z) - survfun(x=x, z=z))/0.001
# }
# 
dS.dm <- function(m, Pars){
  (DEBsurvfun(m+0.001, Pars) - DEBsurvfun(m, Pars))/0.001
}

# OFFSPRING SIZE SENSITIVITY
#Sensitivity of the offspring length density function to temperature (approximate)
# dO.dz <- function(y, z=10.34){
#   (off.length(y=y, z=z+.001) - off.length(y=y, z=z))/.001
# }

dO.dT <- function(Pars){
  (DEBoff.size(c(T=Pars[["T"]]+0.001, 
                 kappa =Pars[["kappa"]], 
                 Y = Pars[["Y"]]))
   - (DEBoff.size(Pars)))/0.001
}

# Sensitivity of the fecundity function to length (approximate) 
# dF.dx <- function(x, z=10.34){
#   (fecfun(x=x+0.001, z=z) - fecfun(x=x, z=z))/0.001
# }
# 
dR.dm <- function(m, Pars){
  (DEBrepfun(m+0.001, Pars) - DEBrepfun(m, Pars))/0.001
}

# Sensitivity of the fecundity function to temperature (approximate) 
# dF.dz <- function(x,z=10.34){
#   (fecfun(x=x, z=z+.001) - fecfun(x=x,z=z))/.001
# }

dR.dT <- function(m, Pars){
  (DEBrepfun(m, c(T=Pars[["T"]]+0.001, 
                  kappa =Pars[["kappa"]], 
                  Y = Pars[["Y"]])) 
   - DEBrepfun(m, Pars))/0.001
}

# Vindenes etal 2014 code
# thermal.sensitivity <- function(z=10.34){
#   res <- wvlambda(Kmat = K.matrix(z=z))
#   Fmat <- Smat <- Omat <- Gmat <- matrix(NA, n, n)
#   sx <- survfun(x, z=z)
#   dsx <- dS.dz(x, z=z)
#   fx <- fecfun(x, z=z)
#   dfx <- dF.dz(x, z=z)
#   offs <- off.length(y=x,z=z)
#   doff <- dO.dz(y=x,z=z)
#   gx <- dgx <- matrix(NA,n,n)
#   for(j in 1:n){
#     gx[,j] <- growthfun(x[j], y=x, z=z)
#     dgx[,j] <- dG.dz(y=x, x=x[j], z=z)
#   }
#   gx[,n] <- gx[,n-1]
#   dgx[,n] <- dgx[,n-1]
#   for(i in 1:n){			
#     Fmat[i,] <- res$v[i]*res$w*(dfx*offs[i])*dx
#     Omat[i,] <- res$v[i]*res$w*(doff[i]*fx)*dx
#     Smat[i,] <- res$v[i]*res$w*(dsx*gx[i,])*dx
#     Gmat[i,] <- res$v[i]*res$w*(dgx[i,]*sx)*dx
#   }
#   list("Fecundity contribution"= apply(Fmat,2,sum),"Survival contribution"=apply(Smat,2,sum),"Growth contribution"= apply(Gmat,2,sum),"Offspring length contribution"= apply(Omat,2,sum),"Total contribution"=apply(Fmat,2,sum)+apply(Smat,2,sum)+apply(Gmat,2,sum)+apply(Omat,2,sum))
# }	

thermal.sensitivity <- function(Pars){
  res <- wvlambda.projection(Kmat = K.matrix(Pars))
  Fmat <- Smat <- Omat <- Gmat <- matrix(NA, n+1, n+1)
  sx <- c(DEBsurvfun(e_m, Pars),DEBsurvfun(x, Pars))
  dsx <- dS.dT(c(e_m,x), Pars)
  fx <- c(DEBrepfun(0, Pars), DEBrepfun(x, Pars))
  dfx <- dR.dT(c(0,x), Pars)
  offs <- c(0,DEBoff.size(Pars))
  doff <- c(0,dO.dT(Pars))
  gx <- dgx <- matrix(0,n+1,n+1)
  gx[2:(n+1),2:(n+1)] <- DEBgrowthfun(x, Pars) # unlike vindenens 2014, the DEBgrowthfun returns all y-distributions
  for(j in 2:(n+1)){
    dgx[2:(n+1),j] <- dG.dT(x[j-1], Pars)
  }
  gx[,n+1] <- gx[,n]
  dgx[,n+1] <- dgx[,n]
  for(i in 1:(n+1)){			
    Fmat[i,] <- res$v[i]*res$w*(dfx*sx[i])*dx
    Omat[i,] <- res$v[i]*res$w*(doff[i]*el_surv)*dx
    Smat[i,] <- res$v[i]*res$w*(dsx*gx[i,])*dx
    Gmat[i,] <- res$v[i]*res$w*(dgx[i,]*sx[i])*dx
  }
  list("Fecundity contribution"= apply(Fmat,2,sum),"Survival contribution"=apply(Smat,2,sum),"Growth contribution"= apply(Gmat,2,sum),"Offspring length contribution"= apply(Omat,2,sum),"Total contribution"=apply(Fmat,2,sum)+apply(Smat,2,sum)+apply(Gmat,2,sum)+apply(Omat,2,sum))
}	

Tsensit <-thermal.sensitivity(GR_pars)
# 1=fec sens, 2=surv sens, 3=growth sens, 4=offs sens,

par(mfrow=c(2,2), las=1, bty="l")
plot(c(0,x), Tsensit[[1]], type="l", ylab="", ylim=c(-1e1,1e5),main=expression(paste("Fecudity contributions to d",  lambda, "/dT")), col=1, lwd=2, xlab="x")
plot(c(0,x), Tsensit[[2]], type="l", ylab="", ylim=c(-1e-5,1e-5), main=expression(paste("Survival contributions to d",  lambda, "/dT")), col=2, lwd=2, xlab="Size x")
plot(c(0,x), Tsensit[[3]], type="l", ylab="", ylim=c(-1e-6,1e-5),main=expression(paste("Growth Contributions to d",  lambda, "/dT")), col=3, lwd=2, xlab="Size x")
plot(c(0,x), Tsensit[[4]], type="l", ylab="", ylim=c(0,1e-9),main=expression(paste("Offspring size contributions to d",  lambda,"/dT")), col=2, lwd=2, xlab="Size x")
#plot(c(0,x), Tsensit[[5]], type="l", ylab="", ylim=c(0,1e-9),main=expression(paste("Contributions to d",  lambda, "/d",mu[z])), col=2, lwd=2, xlab="Size x")

Size.sensitivity <- function(Pars){
  res <- wvlambda.projection(Kmat = K.matrix(Pars))
  Fmat <- Smat <- Gmat <- matrix(NA, n+1, n+1)
  sx <- c(DEBsurvfun(e_m, Pars),DEBsurvfun(x, Pars))
  dsx <- dS.dm(c(e_m,x), Pars)
  fx <- c(DEBrepfun(0, Pars), DEBrepfun(x, Pars))
  dfx <- dR.dm(c(0,x), Pars)
  gx <- dgx <- matrix(0,n+1,n+1)
  gx[2:(n+1),2:(n+1)] <- DEBgrowthfun(x, Pars) # unlike vindenens 2014, the DEBgrowthfun returns all y-distributions
  for(j in 2:(n+1)){
    dgx[2:(n+1),j] <- dG.dm(x[j-1], Pars)
  }
  gx[,n+1] <- gx[,n]
  dgx[,n+1] <- dgx[,n]
  for(i in 1:(n+1)){			
    Fmat[i,] <- res$v[i]*res$w*(dfx*sx[i])*dx
    Smat[i,] <- res$v[i]*res$w*(dsx*gx[i,])*dx
    Gmat[i,] <- res$v[i]*res$w*(dgx[i,]*sx[i])*dx
  }
  list("Fecundity contribution"= apply(Fmat,2,sum),"Survival contribution"=apply(Smat,2,sum),"Growth contribution"= apply(Gmat,2,sum),"Total contribution"=apply(Fmat,2,sum)+apply(Smat,2,sum)+apply(Gmat,2,sum))
}	

Ssensit <- Size.sensitivity(GR_pars)
# 1=fec sens, 2=surv sens, 3=growth sens,4 = total

par(mfrow=c(2,2), las=1, bty="l")
plot(c(0,x), Ssensit[[1]], type="l", ylab="", ylim=c(-1e1,1e7),main=expression(paste("Fecudity contributions to d",  lambda, "/dW")), col=1, lwd=2, xlab="x")
plot(c(0,x), Ssensit[[2]], type="l", ylab="", ylim=c(-1e-8,1e-7), main=expression(paste("Survival contributions to d",  lambda, "/dW")), col=2, lwd=2, xlab="Size x")
plot(c(0,x), Ssensit[[3]], type="l", ylab="", ylim=c(0,1e-8),main=expression(paste("Growth Contributions to d",  lambda, "/dW")), col=3, lwd=2, xlab="Size x")

