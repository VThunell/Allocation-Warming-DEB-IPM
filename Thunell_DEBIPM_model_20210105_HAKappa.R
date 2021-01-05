
############### Thunell et al. DEB-IPM script #########################
  
# This script contains the DEB model describing mean somatic growth and reproduction
# and functions for vital rates needed for the IPM: variance function for somatic growth,
# a size and temperature dependent survival function, function for temperature 
# dependent mean offspring size distrubution, variance function for offspring size. 
# Lastly, the script contains the IPM kernel and projection matrix, calculation of lambda, 
# stable population structure and reproductive values

setwd("C:/Users/vitl0001/VThunell_repos/Temperature-DEBIPM")

## Windermere Pike data
YVPike <- load("~/Manus2/R/Manus2R/PikeDataFiles.R")

### Packages ####
#install.packages("deSolve")
library(deSolve)   # for ode Solver
#install.packages("tidyverse")
library(tidyverse) # gglot etc.
#install.packages("grid")
library(grid)      # for grid.text
#install.packages("fields")

#VVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVV
### DEB Dynamic energy budget model ####
#VVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVV

#### Mass dependence functions of vital rates #####
# a*(Mass^b), a = eps1 or rho1, b = eps2 or rho2

# Allometric parameters (mass dependence at the reference temperature)
rho1 <- 0.009  # Maintenance allometric scalar, rescaled from 0.02 based on pike estimated by Lindmark from Armstrong 1992.
rho2 <- 0.87   # Maintenance allometric exponent, rescaled from 0.76 in Lindmark 2019 based on Diana 1982(?)
eps1 <- 0.47   # Intake allometric scalar, free parameter to fit grwth in Windermere pike.
eps2 <- 0.58   # Intake allometric exponent, free parameter to fit grwth in Windermere pike.

# General DEB-parameters
alpha <- 0.4   # assimilation efficiency
Y2 <- 0.4      # feeding level for mass < 1
sl <- 183      # Season length in days
ha <- 4
#HR_exp = 0.999   # Hyper allometric relationship of kappa, e.g. Kappa = K*(Mass^HR_exp)/Mass
#plot(x, 0.8(x^HR_exp)/x, ylab = "Kappa", type="l", main = "Hyperallometric scaling of 1-Kappa")
#plot(x, 0.8*exp(-x/(3*max(x))), ylab = "Kappa", type="l", main = "Hyperallometric scaling of 1-Kappa")
#plot(x, -0.8*(x-)/x, ylab = "Kappa", type="l", main = "Hyperallometric scaling of 1-Kappa")

### Temperature dependence functions of vital rates with interaction between Mass & Temp ####
# Arrhenius-Lindmark function (AL, Boltzmnn-Arhhenius with interaction betwen Mass & Temp)
rM_T_AL2 <- function(T, m) { # AL Temp dependence of Maintenance 
  (m^(cM*(T-T0)))*exp(EaM*(T-T0)/(k*T*T0)) }

# Gardmark Unimodal function
rI_T_GU2 <- function(T, m) { # GU Temp dependence of intake 
  (m^(cIa*(T-T0))) * exp(EaI*(T-T0)/(k*T*T0)) *
    (m^(-(cId*(T-Td)))) * (1+exp(EdI*(T-Td)/(k*T*Td)))^(-1) * 
    (m^(cId*(T0-Td))) * (1+exp(EdI*(T0-Td)/(k*T0*Td))) }

# # Temperature parameters
T0 <-  283             # Reference Temp
k  <-  8.617332e-05    # Boltzmann constant
cIa <- 0     # linear interaction between size and temp for Max intake.
cM <-  0.017 # linear interaction between size and temp for Maintenance. 0.017 Lindmark unpub.
EaI <- 0.66  # activation energy Intake, Lindmark unpub. 0.66
EaM <- 0.61  # activation energy Maintenance, Lindmark unpub: 0.61
EaS <- 0.47  # activation energy Survival (mortality), Brown et al. 2004
Td  <- T0+6  # deactivation temperature. 6 degrees in Lindmark unpub.?
cId <- -0.02 # linear interaction effect (slope) between temp and mass for deactivation. -0.02 Guesstimate
EdI <- 2     # deactivation energy. 2 is guesstimate

###  TEMP & MASS DEPENDENT RATE FUNCTION ####
# ratefun() uses dwdt() through ode() and returns the mass and rate for intake and maintenance dependent on mass through
# integration over one growing season. ratefun() returns an 3d array where rows is days, cols is rates, dim is x
# if x is a single value, ratefun returns a matrix.
# Intake here refers to net energy available for growth, reproduction and maintenance

dwdt <- function(Time, Mass, Pars) { 
    Maintenance <- (rho1*(Mass^rho2)*rM_T_AL2(Pars['T'], Mass))   # Maintenance rate at Time, mass with Pars
    ifelse(Mass > 1, 
        Intake <- (alpha*Pars['Y']*eps1*(Mass^eps2)*rI_T_GU2(Pars['T'], Mass)),
        Intake <- (alpha*Y2*eps1*(Mass^eps2)*rI_T_GU2(Pars['T'], Mass))) # Intake energy rate at Time, mass with Pars, NOTE intake rate multiplied with alpha and Y
    G <- Pars['Kappa']*exp(-Mass/(ha*max(x)))*Intake - Maintenance    # Growth rate at Time, mass with Pars
    #G <- Pars['Kappa']*Intake - Maintenance    # Growth rate at Time, mass with Pars
     return(list(G, Maintenance, Intake))
    }# return all rates

ratefun <- function(x, Pars) { # mean function for y (dependent on mass) over one time step (one season). 
  yini <- x   # initial biomass in integration
  time <- seq(0, sl, by = 1) # time interval for an sl day long growing season, assumed for pike
  a <- sapply(yini, FUN=ode, time, func = dwdt, parms = Pars)#, simplify = "array")
  a <- array(as.numeric(unlist(a)), dim=c(sl+1, 4, length(yini)), 
             dimnames = list(c(1:(sl+1)), c("day","mass","maintenance","intake"), x))
  return(a) }

#view(ratefun(c(100,1000),GR_pars))
#VVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVV
### IPM Integral projection model ####
#VVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVV

### SIZES & PARAMETERS FOR IPM ####
mmin <- 1     # min weight
mmax <- 14000 # max weight
n <- 200     # row/col number of the discretization-matrix of the continuos rates
x <- seq(mmin,mmax, length=n)
dx <- x[2] - x[1] # step size (grams) in the continuos size spectra

o_e <- 0.00351388  # Egg weight calculated from windermere pike: mean(FData$Egg.weight))) 
# o_evar <- 0.0009163516^2 # variance of egg weight from windermer pike sd(FData$Egg.weight)
# o_eyd <- dlnorm(agg, meanlog = log(o_e) - .5*log(1 + o_evar/o_e^2), sdlog = sqrt(log(1 + o_evar/o_e^2))) #Distribution of y
el_surv = 2e-5 # egg & larvae survival. losely based on Kipling and Frost 1970 1/50000 from laid egg to age 2. 1.9e-4 in Vindenes 2014

ltow_A  <- function(l){0.00254*l^3.271} #length to weight function adults
# Length weight relationship from average in fishbase
ltow_J <- function(l){0.01101*l^2.69} #length to weight function juveniles
# Length weight relationship from fishbase (Carlander 1969) 

### VITAL RATES FUNCTIONS ####

### GROWTH g(x,T) ####
# DEBgrowthfun reurns a vector with a size distribution if given a single size x
# or if x is a vector, a matrix of all size distributions y.

DEBgrowthfun <- function(mass, Pars, y=x) {
  mean <- ratefun(mass, Pars)[sl+1,2,] #the last days mass for all x:s
  xsd <- function(mass, nu.var = -0.0001,  sdres = 100){
    sdres * exp(nu.var*mass) }
  ydf <- function(mu, y=x) {
    var <- xsd(mu)^2
    yd <- dlnorm(y, meanlog = log(mu) - .5*log(1 + var/mu^2), sdlog = sqrt(log(1 + var/mu^2))) #Distribution of y
    if(sum(yd*dx) == 0){
      return (c(rep(0, n-1), 1/dx))
    } else { return(yd/sum(yd*dx)) } #scale it to one
  } 
  return(sapply(mean, FUN = ydf))
}
#DEBgrowthfun(x,GR_pars)
### REPRODUCTION f(x,T) ####
# DEBrepfun uses the energy investment in reproduction dependent on mass and T 
# over t to produce eggs in t+1 for an individual that in the start of t is x. 
# Allocation to reproduction starts when an individuals reaches maturation size, 
# even if that occurs during t -> t+1 making reproductive reserve in t+1 
# dependent on number of days as mature in t.
# Mean egg weight 0.00351388 in Lake Windermere data

DEBrepfun <- function(x, Pars, o_E) { # Note that intake refers to net energy available for growth, reproduction and maintenance
  Devm = 374 # *(1 - pars[["Kappa"]])
  fecund <- function(m, Pars, o_E) { 
    f <- ratefun(m, Pars) 
    fe <- ifelse(f[,"mass",] < Devm,  0, #is it mature? # see text for the following conditions...
            ifelse(Pars[["Kappa"]]*exp(-x/(ha*max(x)))*f[,"intake",] >= f[,"maintenance",], f[,"intake",]*(1-Pars[["Kappa"]]*exp(-x/(ha*max(x)))), #can available energy from intake cover maintenance, NOTE:intake is already scaled with Y*alpha
              ifelse((Pars[["Kappa"]]*exp(-x/(ha*max(x)))*f[,"intake",] < f[,"maintenance",]) & (f[,"maintenance",] <= f[,"intake",]), 
                  f[,"intake",] - f[,"maintenance",], 0))) # can available energy from intake and reproduction cover maintenance, if not 0
         # expression for this reproduction energy in the text is simplified here to I(m,T)-M(m,T)
    sum(fe)/o_E * 0.5 #- 10*y^0.6# Summed energy allocated to reproduction reserve in one season diveded by egg weight and 0.5 from sex ratio (half of the eggs are female)
  }
  sapply(x, fecund, Pars, o_E)
}

### OFFSPRING SIZE DISTRIBUTION o(eggsize,T)####
# Probability of an egg to hatch and grow into size y in t+1?
# From Vindenes 2014: "Density distribution for offspring size y (we assume parental length x does not affect offspring length)."
# Temp dependent mean weight of offspring is calculated from YV 2014 linear relationship
# Variance from Windermere: sd(ltow_J((GData[GData$Age==1,]$Length))) (3.72 for length in Vindenes 2014)

DEBoff.size <- function(y, Pars, offvar = 24.05458^2){ #lognormal distribution of offspring weights
  mu <- ratefun(o_e, Pars)[sl/2,2,] 
  yd <- dlnorm(y, meanlog=log(mu) - .5*log(1 + offvar/(mu^2)), sdlog=sqrt(log(1 + offvar/(mu^2))))
      yd/sum(yd*dx) # scaled distributon to 1, dx is step size in the k-matrix
}

#DEBoff.size(y,GR_pars)
### SURVIVAL s(x,T) ####
DEBsurvfun <- function(Mass, Pars) { # Mass-Temp dependence of yearly Mortality 
    sx <- function(m, Pars){ # natural baseline survival  
      exp((-3*m^-.288)*exp(EaS*(Pars[['T']]-T0)/(k*Pars[['T']]*T0))) #3*^-.288 is yearly mortality rate for 1 gram unit weight from Lorenzen 1996 
    } 
    sx.firstyear <- function(m, Pars){ # first year survival
      el_surv*sx(ratefun(o_e,Pars)[sl/2,2,], Pars)
    }
    starvx <- function(m,Pars) { # starvation survival
      starvS <- ratefun(m, Pars)[sl+1,2,]
      ifelse(starvS < m, 0, 1)
    }
    ifelse(Mass == o_e, sx.firstyear(Mass,Pars), sx(Mass,Pars)*starvx(Mass, Pars))# multiply with F (fishing mortality function)
}

### Projection Matrix ####
# The projection or Kernel (K) matrix maps the size distribution in time t to time t+1. 

# Smat column 1 describes the probality of an egg to hatch and grow into size y in t+1 or the density distribution of offspring size y, 
# the other columns of Smat describe the growth and survival probability density of y (from x)
# the first row is 0 to give way for an egg stage (the number of eggs produced by the sizes in x) through Fmat because at census for x, size x will give eggs at t+1 whihc will be offspring in t+2

# Fmat is the egg stage, number of eggs produced for each x which becomes y1.

K.matrix <- function(Pars) {
  Smat <- Fmat <- matrix(0,n+1,n+1) 
  sx <- c(DEBsurvfun(o_e, Pars), DEBsurvfun(x, Pars)) # survival vector
  Smat[2:(n+1),1] <- sx[1]*DEBoff.size(x, Pars) * dx # First column (and 2:101 row) of Smat is sum of survival of laid eggs, larvae and fry to age 1 (offspring size),  
  Smat[2:(n+1),2:(n+1)] <- sx[2:(n+1)]*DEBgrowthfun(x, Pars, y=x) * dx
  Fmat[1,2:(n+1)] <- sx[2:(n+1)]*DEBrepfun(x, Pars, o_e) #Production of eggs from each size class, (the number of eggs that x in t (size in april) will release in april in t+1) * (probability to grow from egg to age 1 (in t+2-1)
  #Smat[,n+1] <- Smat[,n] # ????????
  Smat+Fmat }

# This function uses the R function "eigen" to calculate lambda, a stable structure "w" and reproductive values "v", for a projection matrix "Kmat". The stable structure is scaled so that sum(w*dx)=1 and the reproductive values are scaled so that sum(v*w*dx)=1 (see details in the article).
# For IPMs, including threshold for setting value to zero
wvlambda <- function(Kmat, tol=1e-20) {
  ev <- eigen(Kmat)
  tev <- eigen(t(Kmat))
  lmax <- which.max(Re(ev$values))
  W <- ev$vectors
  V <- tev$vectors
  w <- as.matrix(abs(Re(W[, lmax]))/sum(abs(Re(W[, lmax]))))
  w <- ifelse(w<tol,0,w)
  w <- w/(sum(w*dx))
  v <- as.matrix(abs(Re(V[, lmax])))
  v <- ifelse(w*v*dx > tol, v/sum(w*v*dx), v/tol) # avoid producing inf
  # v <- ifelse(is.infinite(v/(w*v*dx)), 0, v/sum(w*v*dx)) # alternative to avoid producing inf
  # v <- v/sum(w*v*dx) # YV 2014 original code
  v <- ifelse(w*v <= 0, 0, v)
  return(list("lambda" = max(Re(ev$values)), "w"=w, "v"=v)) }

###CALCULATE LAMBDA, STABLE STRUCTURE AND REPRODUCTIVE VALUES FOR VARYING KAPPA AND TEMP VALUES ####

# T <- seq(282,284,1) # temperature range
# Kappa <-seq(0.6,1,0.05) # Kappa range
# Y <- 1 # feeding levels
# 
# Res_lam <- matrix(ncol = 3+1, nrow = length(Kappa)*length(T)*length(Y)) # assuming 40 here from files
# Res_v  <-  matrix(ncol = 3+n+1, nrow = length(Kappa)*length(T)*length(Y)) # assuming 40 here from files
# Res_w  <-  matrix(ncol = 3+n+1, nrow = length(Kappa)*length(T)*length(Y)) # assuming 40 here from files
# 
# parsK <- as.matrix(expand.grid(T,Kappa,Y))
# colnames(parsK) <- c("T","Kappa","Y")
# 
# for (i in 1:nrow(parsK)){
#   res <- wvlambda(K.matrix(parsK[i,]))
#   Res_lam[i,] <- c(parsK[i,],res$lam) #c(T[d],Kappa[e],Y[f], res$lambda)
#   Res_v[i,]   <- c(parsK[i,],res$v) #REPRODUCTIVE VALUES
#   Res_w[i,]   <- c(parsK[i,],res$w) #STABLE STRUCTURE
#  }
# 
# colnames(Res_lam) <- c("T","Kappa","Y","Lambda")
# colnames(Res_v) <- c("T","Kappa","Y","0",x)
# colnames(Res_w) <- c("T","Kappa","Y","0",x)
#write.table(Res_lam, file="Res_lam1207.txt",quote=TRUE, sep=",", row.names=TRUE)

## Sensitivity analyses based on Merow et al. 2014 appendix (section 1.4.5)

#"The eigen-things can be combined to obtain the sensitivity and elasticity matrices."
# v.dot.w=sum(stable.dist*repro.val)*h
# sens=outer(repro.val,stable.dist)/v.dot.w
# elas=matrix(as.vector(sens)*as.vector(K)/lam,nrow=n)
