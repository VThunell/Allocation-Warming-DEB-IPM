
############### Thunell et al. DEB-IPM script #########################
  
# This script contains the DEB model describing energy accumulation and allocation 
# to somatic growth and reproduction and functions for vital rates needed for 
# the IPM: variance function for somatic growth, a size and temperature dependent 
# survival function, function for temperature dependent mean offspring size distribution, 
# variance function for offspring size. Lastly, the script contains the 
# IPM kernel and projection matrix, calculation of lambda, stable population 
# structure and reproductive values

#setwd("C:/Users/vitl0001/VThunell_repos/Temperature-DEBIPM")
YVPike <- load("~/Manus2/R/Manus2R/PikeDataFiles.R")

### Packages ####
#install.packages("deSolve")
library(deSolve)   # for ode Solver
#install.packages("tidyverse")
library(tidyverse) # gglot etc.
#install.packages("grid")
library(grid)      # for grid.text
#install.packages("fields")
library(tictoc)
library(scales)
library(patchwork)

#VVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVV
### DEB Dynamic energy budget model ####
#VVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVV

#### Mass dependence (allometric) functions of vital rates #####
# a*(Mass^b), a = eps1 or rho1, b = eps2 or rho2

## Allometric parameters (mass dependence at the reference temperature)
# Parameters to fit growth and repro in Windermere pike with kappa=0.93
rho1 <- 0.02 # 0.011 for 287 # Maintenance allometric scalar, rescaled from 0.02 based on pike estimated by Lindmark from Armstrong 1992.
rho2 <- 0.76 # 0.83 for 287  # Maintenance allometric exponent, rescaled from 0.76 in Lindmark 2019 based on Diana 1982(?)
eps1 <- 0.69 #0.51 #45 #old 0.54 #614 #54 #64 # Intake allometric scalar, free parameter to fit grwth in Windermere pike.
eps2 <- 0.54 #0.574 #59 #old 0.58 #55 #58 # Intake allometric exponent, free parameter to fit grwth in Windermere pike.

## General DEB-parameters
alpha <- 0.4   # assimilation efficiency
#Y2 <- 1      # feeding level for mass < 1
sl <- 183      # Season length in days, i.e. number of growth time steps + 1
kap_fun <- function(m, kappa, ha=1.5){kappa*exp(-m/(ha*max(x)))} #for size dep. kappa

## Temperature scaling parameters
T0 <-  292            # Reference Temp
k  <-  8.617332e-05   # Boltzmann constant
cIa <- 0 #0.005    # linear interaction between size and temp for Maximum intake.
cM  <- 0 #.0026 # linear interaction between size and temp for Maintenance. 0.017 Lindmark unpub.
cId <- 0 # linear interaction effect (slope) between temp and mass for deactivation. -0.02 Guesstimate to fit assumptions on temperature dependent growth
EaI <- 0.58  # activation energy Intake, Lindmark bioRxiv 0.58 for sharpe and 0.69 for other
EaM <- 0.62 # activation energy Maintenance, Lindmark bioRxiv: 0.62
EaS <- 0.47  # activation energy Survival (mortality), Brown et al. 2004
EdI <- 2.64     # deactivation energy. 2 is guesstimate 
Td  <- T0+4  # deactivation temperature. 4 degrees in Lindmark unpub.?

DEBparams <- as.data.frame(rbind(rho1,rho2,eps1,eps2,alpha,sl,T0,
                                 k,cIa,cM,EaI,EaM,EaS,Td,cId,EdI))
### Temperature dependence functions of vital rates with interaction between Mass & Temp ####
# Arrhenius-Lindmark function (AL, Boltzmann-Arrhenius with interaction between Mass & Temp)
rM_T_AL2 <- function(T, m) { # AL Temp dependence of Maintenance 
  (m^(cM*(T-T0)))*exp(EaM*(T-T0)/(k*T*T0)) }

# Gardmark Unimodal function
rI_T_GU2 <- function(T, m) { # GU Temp dependence of intake 
  (m^(cIa*(T-T0))) * exp(EaI*(T-T0)/(k*T*T0)) *
    (m^(-(cId*(T-Td)))) * (1+exp(EdI*(T-Td)/(k*T*Td)))^(-1) * 
    (m^(cId*(T0-Td))) * (1+exp(EdI*(T0-Td)/(k*T0*Td))) }

# bTc=0.7
# rI_T_Pad <- function(T) { # GU Temp dependence of intake 
#   bTc * exp(EaI*(T-T0)/(k*T*T0)) *
#     (1+exp(EdI*(T-Td)/(k*T*Td)))^(-1)}
# rI_T_Pad(292)

###  TEMP & MASS DEPENDENT RATE FUNCTION ####
# ratefun() uses dwdt() through ode() and returns the mass and rate for intake and maintenance dependent on mass through
# integration over one growing season. ratefun() returns an 3d array where rows is days, cols is rates, dim is x
# if x is a single value, ratefun returns a matrix.
# Intake here refers to net energy available for growth, reproduction and maintenance

dwdt <- function(time, m, Pars) {
  maintenance <- (rho1*(m^rho2)*rM_T_AL2(Pars['T'], m))   # Maintenance rate at time, mass with Pars
  #ifelse(m >= 1,
  intake <- alpha*Pars['Y']*eps1*(m^eps2)*rI_T_GU2(Pars['T'], m)#,
  #intake <- (alpha*Y2*eps1*(m^eps2)*rI_T_GU2(Pars['T'], m))) # Intake energy rate at time, mass with Pars, NOTE intake rate multiplied with alpha and Y
  mass <- kap_fun(m,Pars["kappa"])*intake - maintenance    # Growth rate at time, mass with Pars
  #mass <- Pars['kappa']*intake - maintenance    # Growth rate at time, mass with Pars
  return(list(mass, maintenance, intake)) # return all rates
 }


ratefun <- function(m, Pars) { # mean function for y (dependent on mass) over one time step (one season). 
  yini <- m   # initial biomass in integration
  time <- seq(1, sl, by = 1) # time interval for an sl day long growing season, assumed for pike
  a <- sapply(yini, FUN=ode, time, func = dwdt, parms = Pars)#, simplify = "array")
  a <- array(as.numeric(unlist(a)), dim=c(sl, 4, length(yini)), 
             dimnames = list(c(1:(sl)), c("day","mass","maintenance","intake"), m))
  return(a) 
}

#VVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVV
### IPM Integral projection model ####
#VVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVV

### SIZES & PARAMETERS FOR IPM ####
mmin <- 1     # min weight
mmax <- 18000 # max weight
n <- 500    # row/col number of the discretization-matrix of the continuos rates
x <- seq(mmin,mmax, length=n)
dx <- x[2] - x[1] # step size (grams) in the continuos size spectra

e_m <- 0.00351388  # Egg mass calculated from windermere pike: mean(FData$Egg.weight))) 
# e_mvar <- 0.0009163516^2 # variance of egg weight from windermer pike sd(FData$Egg.weight)
# e_myd <- dlnorm(agg, meanlog = log(e_m) - .5*log(1 + e_mvar/e_m^2), sdlog = sqrt(log(1 + e_mvar/e_m^2))) #Distribution of y
el_surv = 1.9e-4#1e-5 # egg & larvae survival. losely based on Kipling and Frost 1970 1/50000 from laid egg to age 2. 1.9e-4 in Vindenes 2014
IPMparams <- as.data.frame(rbind(mmin,mmax,n,dx,e_m,el_surv))

# Length weight relationship from Windermere by YV
ltow_AYV  <-  function(l){exp(-6.49286)*l^3.4434} #length to weight function adults
wtol_AYV <- function(w){(w/exp(-6.49286))^(1/3.4434)}
#plot(wtol_AYV(x),x,type="l", lwd=2)
#lines(1:200,ltow_AYV(1:200),col="red", type="l",lty=2, lwd=2)
# Length weight relationship for juveniles from fishbase (Carlander 1969) 
# ltow_J <- function(l){0.01101*l^2.69} #length to weight function juveniles

### VITAL RATES FUNCTIONS ####

### GROWTH g(x,T) ####
# DEBgrowthfun returns a vector with a size distribution if given a single size x
# or if x is a vector, a matrix of all size distributions y.

DEBgrowthfun <- function(m, Pars, y=x) {
  xsd <- function(mu, nu.var = -0.0001,  sdres = 100){ # generates sd for mu
    sdres * exp(nu.var*mu) }
  ydf <- function(mu, y=x) { # generates m:s size distributions for next year
    var <- xsd(mu)^2 # variance
    yd <- dlnorm(y, meanlog = log(mu) - .5*log(1 + var/mu^2), sdlog = sqrt(log(1 + var/mu^2))) #Distribution of y
    if(sum(yd*dx) == 0){
      c(rep(0, n-1), 1/dx)
    } else { yd/sum(yd*dx) } } #scale it to one
  mu <- ratefun(m, Pars)[sl,2,] # mass on the last day of growth season for m
  sapply(mu, FUN = ydf)
}

### REPRODUCTION f(x,T) ####
# DEBrepfun uses the energy investment in reproduction dependent on mass and T 
# over t to produce eggs in t+1 for an individual that in the start of t is x. 
# Allocation to reproduction starts when an individuals reaches maturation size, 
# even if that occurs during t -> t+1 making reproductive reserve in t+1 
# dependent on number of days as mature in t.
# Mean egg weight 0.00351388 in Lake Windermere data

DEBrepfun <- function(m, Pars, e_m = 0.00351388) { # Note that intake refers to net energy available for growth, reproduction and maintenance
   Devm = 374 # *(1 - pars[["kappa"]])
   fecund <- function(m, Pars, e_m) {
    f <- ratefun(m, Pars)
    f <- ifelse(f[,"mass",] < Devm,  0, #is it mature? # see text for the following conditions...
           ifelse(kap_fun(m,Pars["kappa"])*f[,"intake",] >= f[,"maintenance",], f[,"intake",]*(1-kap_fun(m,Pars["kappa"])), #can available energy from intake cover maintenance, NOTE:intake is already scaled with Y*alpha
             ifelse(kap_fun(m,Pars["kappa"])*f[,"intake",] < f[,"maintenance",] & (f[,"maintenance",] <= f[,"intake",]),
                   (1-kap_fun(m,Pars["kappa"]))*f[,"intake",] + kap_fun(m,Pars["kappa"])*f[,"intake",] - f[,"maintenance",], 0)))
    sum(f)/e_m * 0.5 #- 10*y^0.6# Summed energy allocated to reproduction reserve in one season diveded by egg weight and 0.5 from sex ratio (half of the eggs are female)
   }
   sapply(m, fecund, Pars, e_m)
 }

### OFFSPRING SIZE DISTRIBUTION o(eggsize,T)####
# Probability of an egg to hatch and grow into size y in t+1?
# From Vindenes 2014: "Density distribution for offspring size y (we assume parental size x does not affect offspring size)."
# Temp dependent mean weight of offspring is calculated from YV 2014 linear relationship
# Variance from Windermere: sd(ltow_J((GData[GData$Age==1,]$Length))) (3.72 for length in Vindenes 2014)

DEBoff.size <- function(Pars, y=x, offvar = 24.05458^2){ #lognormal distribution of offspring weights
  mu <- ratefun(e_m, Pars)[round(sl/2),2,] 
  yd <- dlnorm(y, meanlog=log(mu) - .5*log(1 + offvar/(mu^2)), sdlog=sqrt(log(1 + offvar/(mu^2))))
      yd/sum(yd*dx) # scaled distributon to 1, dx is step size in the k-matrix
}

### SURVIVAL s(x,T) ####
DEBsurvfun <- function(m, Pars) { # Mass-Temp dependence of yearly Mortality 
  sx <- function(m, Pars){ # natural baseline survival
    exp(-3*m^-0.2*exp(EaS*(Pars[['T']]-T0)/(k*Pars[['T']]*T0))) #3*^-.288 is yearly mortality rate for 1 gram unit weight from Lorenzen 1996
    } #conf limits Lorenzon 1996: 0.315, 0.261
  Vsurvfun <- function(x, z=10.34){
    sxV <- function(m, z=10.34){ # sx from Vindenes 2014
      1/(1+exp(13.53316 - 0.50977*wtol_AYV(m) - (-0.00393)*wtol_AYV(m)^2 - 0.19312*z - (-0.00679)*wtol_AYV(m)*z))
    }
    xmax <- x[which(sxV(x) == max(sxV(x)))]
    ifelse(x < xmax, sxV(x), sxV(xmax))
  }
  sx.firstyear <- function(m, Pars){ # first year survival
    el_surv#*sx(ratefun(e_m,Pars)[round(sl/2),2,], Pars)
    }
  starvx <- function(m,Pars) { # starvation survival
    starvS <- ratefun(m, Pars)[sl,2,]
    ifelse(starvS < m | starvS > 18000, 0, 1)
    }
  #ifelse(m == e_m, sx.firstyear(m,Pars), sx(m,Pars)*starvx(m, Pars))# multiply with F (fishing mortality function)
  ifelse(m < 1, sx.firstyear(m,Pars), Vsurvfun(m)*starvx(m, Pars))# multiply with F (fishing mortality function)
  #ifelse(m == e_m, sx.firstyear(m,Pars), 0.67*starvx(m, Pars))# multiply with F (fishing mortality function)
} 

# starvx(17212, Pars=c(T = 285, kappa = 0.83, Y = 1))
# ratefun(17220, Pars=c(T = 285, kappa = 0.83, Y = 1))[sl,2,]
# # survi<-as.tibble(cbind(x,sx(x,GR_pars),Vsurvfun(x)))
# ggplot(survi)+
#        geom_line(aes(x,V2))+
#        geom_line(aes(x,V3), color = "red")
# plot(x,sx(x,GR_pars),type="l", col="blue", ylim=c(0,1))
# lines(x,Vsurvfun(x),type="l", col="red")
# GR_pars <- c(T = 287,     # parameters for Temperature, feeding, allocation and Mass dependence
#              kappa = 0.93, # allocation to respiration (Growth and maintenance)
#              Y = 1)         # Feeding level
# lines(x,sx(x,GR_pars),type="l", col="blue")
# legend("bottomright", c("Size dep. Lorenzon, exp=0.20" , "Vindenes et al. 2014","Size dep. Lorenzon, exp=0.288"), lty=c(1,1,1), col = c("green","red", "blue"),  cex=0.7)

### Projection Matrix ####
# The projection or Kernel (K) matrix maps the size distribution in time t to time t+1. 
# Smat column 1 describes the probality of an egg to hatch and grow into size y in t+1 or the density distribution of offspring size y, 
# the other columns of Smat describe the growth and survival probability density of y (from x)
# the first row is 0 to give way for an egg stage (the number of eggs produced by the sizes in x) through Fmat because at census for x, size x will give eggs at t+1 whihc will be offspring in t+2
# Fmat is the egg stage, number of eggs produced for each x which becomes y1.

K.matrix <- function(Pars) {
  Smat <- Fmat <- matrix(0,n+1,n+1)
  surv_x <- c(DEBsurvfun(e_m, Pars), DEBsurvfun(x, Pars))  # survival vector
  Smat[2:(n+1),1] <- surv_x[1]*DEBoff.size(Pars) * dx #dx is dy, First column (and 2:101 row) of Smat probability of transitions from egg to sizes at age 1.  
  Smat[2:(n+1),2:(n+1)] <- t(surv_x[2:(n+1)]*t(DEBgrowthfun(x, Pars) ) * dx)
# Smat[2:(n+1),2:(n+1)] <- t(diag(sx[2:(n+1)])%*%t(DEBgrowthfun(x, Pars) )* dx)
  Fmat[1,2:(n+1)] <- surv_x[2:(n+1)]*DEBrepfun(x, Pars) #Production of eggs
  Smat+Fmat 
}

# hg<- wvlambda.projection(K.matrix(GR_pars))
# sum(hg$w)*dx

# This function uses function eigen() to calculate lambda, a stable structure "w" and reproductive values "v", for a projection matrix "Kmat". The stable structure is scaled so that sum(w*dx)=1 and the reproductive values are scaled so that sum(v*w*dx)=1 (see details in the article).
# For IPMs, including threshold for setting value to zero
# wvlambda <- function(Kmat, tol=1e-20) {
#   ev <- eigen(Kmat)
#   tev <- eigen(t(Kmat))
#   lmax <- which.max(Re(ev$values))
#   W <- ev$vectors
#   V <- tev$vectors
#   w <- as.matrix(abs(Re(W[, lmax]))/sum(abs(Re(W[, lmax]))))
#   w <- ifelse(w<tol,0,w)
#   w <- w/(sum(w*dx))
#   v <- as.matrix(abs(Re(V[, lmax])))
#   # v <- v/sum(w*v*dx) # YV 2014 code producing inf with my Kernel for some T & Kappa
#   v <- ifelse(w*v*dx > tol, v/sum(w*v*dx), v/tol) # avoid producing inf
#   # v <- ifelse(is.infinite(v/(w*v*dx)), 0, v/sum(w*v*dx)) # alternative to avoid producing inf
#   v <- ifelse(w*v <= 0, 0, v)
#   return(list("lambda" = max(Re(ev$values)), "w"=w, "v"=v)) }

wvlambda.projection <- function(Kmat, N0=rep(10,length(Kmat[1,])), tol=1e-6) {
  Nt <- Nt2 <- N0
  lam <- 1
  prev.lam <- .5
  tKmat <- t(Kmat)
  while( abs(lam-prev.lam) > tol ){
    prev.lam <- lam
    Nt.new <- Kmat %*% Nt
    Nt.new2 <- tKmat %*% Nt2
    lam <- sum(Nt.new) / sum(Nt)
    Nt <- Nt.new
    Nt2 <- Nt.new2
  }
  w <- Nt/sum(Nt)
  w <- w/sum(w*dx)
  v <- Nt2/sum(Nt.new2)
  v <- v/sum(v*w*dx)
  return(list("lambda"=lam, "w"=w,"v"=v))
}

# tic()
# wvlambda.projection(K.matrix(Pars <- c(T = 287,     # parameters for Temperature, feeding, allocation and Mass dependence
#                                        kappa = 0.93, # allocation to respiration (Growth and maintenance)
#                                      Y = 1) ))[1]
# toc()

### COHORT PROJECTION MODEL ####
K.matrixCohort <- function(Pars) {
  Smat <- Fmat <- matrix(0,n+1,n+1)
  surv_x <- c(DEBsurvfun(e_m, Pars), DEBsurvfun(x, Pars))  # survival vector
  Smat[2:(n+1),1] <- surv_x[1]*DEBoff.size(Pars) * dx #dx is dy, First column (and 2:101 row) of Smat probability of transitions from egg to sizes at age 1.  
  Smat[2:(n+1),2:(n+1)] <- t(surv_x[2:(n+1)]*t(DEBgrowthfun(x, Pars) ) * dx)
  Fmat[1,2:(n+1)] <- 0
  Smat+Fmat 
}

ProjCoh <- function(Pars,t,I_pop=c(5000000,rep(0,n))) {
  Coh<-K.matrixCohort(Pars)
  for (i in 1:t){
    if(i==1){ y <- Coh%*%I_pop
    ay <- cbind(y,i) }
    else{y <- Coh%*%y
    ay <- rbind(ay,cbind(y,i)) }
  } 
  ay 
}

###CALCULATE LAMBDA, STABLE STRUCTURE AND REPRODUCTIVE VALUES FOR VARYING KAPPA AND TEMP VALUES ####

### Main RESULT - with temp x size interaction and size dep. kappa ####

# T <- seq(285,290,0.3) # temperature range
# kappa <- seq(0.6,1,0.04) # kappa range
# Y <- 1 # feeding levels

# Res_lam <- matrix(ncol = 3+1, nrow = length(kappa)*length(T)*length(Y)) # assuming 40 here from files
# Res_v  <-  matrix(ncol = 3+n+1, nrow = length(kappa)*length(T)*length(Y)) # assuming 40 here from files
# Res_w  <-  matrix(ncol = 3+n+1, nrow = length(kappa)*length(T)*length(Y)) # assuming 40 here from files

# # Res_lam <- matrix(ncol = 3+1, nrow = length(ha_Ad)*length(T)*length(Y)) # assuming 40 here from files
# # Res_v  <-  matrix(ncol = 3+n+1, nrow = length(ha_Ad)*length(T)*length(Y)) # assuming 40 here from files
# # Res_w  <-  matrix(ncol = 3+n+1, nrow = length(ha_Ad)*length(T)*length(Y)) # assuming 40 here from files

# parsK <- as.matrix(expand.grid(T,kappa,Y))
# colnames(parsK) <- c("T","kappa","Y")
# 
# tic()
# for (i in 1:nrow(parsK)){
#   res <- wvlambda.projection(K.matrix(parsK[i,]))
#   Res_lam[i,] <- c(parsK[i,],res$lam) #c(T[d],kappa[e],Y[f], res$lambda)
#   Res_v[i,]   <- c(parsK[i,],res$v) #REPRODUCTIVE VALUES
#   Res_w[i,]   <- c(parsK[i,],res$w) #STABLE STRUCTURE
#  }
# toc()
# 
# colnames(Res_lam) <- c("T","kappa","Y","Lambda")
# colnames(Res_v) <- c("T","kappa","Y","0",x)
# colnames(Res_w) <- c("T","kappa","Y","0",x)

# write.table(Res_lam, file="Res_lam0316_baseline_Onlyyellow_n500.txt",quote=TRUE, sep=",", row.names=TRUE)
# write.table(Res_v, file="Res_v0316__baseline_Onlyyellow_n500.txt",quote=TRUE, sep=",", row.names=TRUE)
# write.table(Res_w, file="Res_w0316__baseline_Onlyyellow_n500.txt",quote=TRUE, sep=",", row.names=TRUE)


# THERMAL or kappa SENSITIVITY OF LAMBDA ####
# Sensitivity analysis to study how change in Lambda can be
# attributed to change in different vital rates originating from 
# a change in kappa or temperature.
# The partial derivative of the output Y with respect to an input factor Xi

# Sensitivity based on Vindenes 2014:

dG.dm <- function(m, Pars){
  (DEBgrowthfun(m+0.1, Pars) - DEBgrowthfun(m, Pars))/0.1
}

#Sensitivity of the growth density function to temperature (approximate) 

dG.dT <- function(m, Pars){
  (DEBgrowthfun(m, c(T=Pars[["T"]]+0.001, 
                     kappa =Pars[["kappa"]], 
                     Y = Pars[["Y"]])) 
   - DEBgrowthfun(m, Pars))/0.001
}

##SURVIVAL SENSITIVITY
#Sensitivity of the survival probability function to temperature (approximate) 

dS.dT <- function(m, Pars){
  (DEBsurvfun(m, c(T=Pars[["T"]]+0.001, 
                   kappa =Pars[["kappa"]], 
                   Y = Pars[["Y"]])) 
   - DEBsurvfun(m, Pars))/0.001
}

#Sensitivity of the survival probability function to length (approximate) 
 
dS.dm <- function(m, Pars){
  (DEBsurvfun(m+0.1, Pars) - DEBsurvfun(m, Pars))/0.1
}

# OFFSPRING SIZE SENSITIVITY
#Sensitivity of the offspring length density function to temperature (approximate)

dO.dT <- function(Pars){
  (DEBoff.size(c(T=Pars[["T"]]+0.001, 
                 kappa =Pars[["kappa"]], 
                 Y = Pars[["Y"]]))
   - (DEBoff.size(Pars)))/0.001
}

# Sensitivity of the fecundity function to length (approximate) 
 
dR.dm <- function(m, Pars){
  (DEBrepfun(m+0.1, Pars) - DEBsurvfun(m, Pars))/0.1
}

# Sensitivity of the fecundity function to temperature (approximate) 

dR.dT <- function(m, Pars){
  (DEBrepfun(m, c(T=Pars[["T"]]+0.001, 
                  kappa =Pars[["kappa"]], 
                  Y = Pars[["Y"]])) 
   - DEBrepfun(m, Pars))/0.001
}

# Vindenes etal 2014 code
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
 
# Tsensit <-thermal.sensitivity(GR_pars)
# # 1=fec sens, 2=surv sens, 3=growth sens, 4=offs sens,
# 
# par(mfrow=c(2,2), las=1, bty="l")
# plot(c(0,x), Tsensit[[1]], type="l", ylab="", ylim=c(-1e1,1e5),main=expression(paste("Fecudity contributions to d",  lambda, "/dT")), col=1, lwd=2, xlab="x")
# plot(c(0,x), Tsensit[[2]], type="l", ylab="", ylim=c(-1e-5,1e-5), main=expression(paste("Survival contributions to d",  lambda, "/dT")), col=2, lwd=2, xlab="Size x")
# plot(c(0,x), Tsensit[[3]], type="l", ylab="", ylim=c(-1e-6,1e-5),main=expression(paste("Growth Contributions to d",  lambda, "/dT")), col=3, lwd=2, xlab="Size x")
# plot(c(0,x), Tsensit[[4]], type="l", ylab="", ylim=c(0,1e-9),main=expression(paste("Offspring size contributions to d",  lambda,"/dT")), col=2, lwd=2, xlab="Size x")
# #plot(c(0,x), Tsensit[[5]], type="l", ylab="", ylim=c(0,1e-9),main=expression(paste("Contributions to d",  lambda, "/d",mu[z])), col=2, lwd=2, xlab="Size x")

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

# Ssensit <- Size.sensitivity(GR_pars)
# # 1=fec sens, 2=surv sens, 3=growth sens,4 = total
# 
# par(mfrow=c(2,2), las=1, bty="l")
# plot(c(0,x), Ssensit[[1]], type="l", ylab="", ylim=c(-1e1,1e7),main=expression(paste("Fecudity contributions to d",  lambda, "/dT")), col=1, lwd=2, xlab="x")
# plot(c(0,x), Ssensit[[2]], type="l", ylab="", ylim=c(-1e-8,1e-7), main=expression(paste("Survival contributions to d",  lambda, "/dT")), col=2, lwd=2, xlab="Size x")
# plot(c(0,x), Ssensit[[3]], type="l", ylab="", ylim=c(0,1e-8),main=expression(paste("Growth Contributions to d",  lambda, "/dT")), col=3, lwd=2, xlab="Size x")

