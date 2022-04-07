
############### THUNELL ET AL. 20XX DEBIPM R-script - MODEL #########################
  
# This script contains the DEB model describing energy accumulation and allocation 
# to somatic growth and reproduction and functions for demographic functions in
# the IPM:  Body size and temperature dependent mean and variance function for somatic
# growth g(m_s+1; m_s, T), fecundity g(m_s, T), survival function a(m_s,T) and
# temperature dependent size at age 1 o(m_s,T) distribution. The script also contains 
# the IPM kernel and projection matrix, calculation of lambda, stable population 
# structure and reproductive values and the numerical expression for the derivatives 
# for mass and temperature specific demographic functions with respect to $\kappa_0$

#setwd("C:/Users/vitl0001/VThunell_repos/Temperature-DEBIPM")
YVPike <- load("~/Manus2/R/Manus2R/PikeDataFiles.R")

### Packages ####
#install.packages("deSolve")
library(deSolve) # numerical integration of DEB
#install.packages("tidyverse")
library(tidyverse) # gglot etc.
#install.packages("grid")
library(tictoc) # clock model runs

#VVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVV
### DEB Dynamic energy budget model ####
#VVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVV

#### Mass dependence (allometric) functions DEB #####
# a*(Mass^b), a = eps1 or rho1, b = eps2 or rho2

## Allometric parameters (mass dependence at the reference temperature)
# Parameters to fit growth and repro in Windermere pike
rho1 <- 0.02 # Maintenance allometric scalar, based on pike estimated by Lindmark 2019 from Armstrong 1992.
rho2 <- 0.76 # Maintenance allometric exponent, rescaled from 0.76 in Lindmark 2019 based on Diana 1982
eps1 <- 1.55 # Consumption allometric scalar, free parameter to fit growth in Windermere pike.
eps2 <- 0.53 # Consumption allometric exponent, free parameter to fit growth in Windermere pike.

## General DEB-parameters 
alpha <- 0.4   # assimilation efficiency
sl <- 183      # Season length in days, i.e. number of growth time steps + 1
kap_fun <- function(m, kappa, ha=2){kappa*exp(-m/(ha*max(x)))} #for body size dep. kappa
#kap_fun <- function(m, kappa, ha=2){kappa*m/m} # for size independent kappa

## Temperature scaling parameter values from Lindmark et al. 2022, Global change biol.
T0 <-  292            # Reference Temp
k  <-  8.617333e-05   # Boltzmann constant
cM  <- 0.0026 # Linear interaction between size and temp for Maintenance
EaC <- 0.73 # activation energy Consumption for Sharpe-Schoolf
EaM <- 0.62 # activation energy Maintenance
EdI <- 1.89 # deactivation energy
Td  <- T0 + 0.75 # deactivation temperature
bTc <- 0.79 # rate at reference (a common) temperature Sharpe-Schoolf

# Parameters for Temp., allocation and feeding level (not used in analysis and excluded the main text) for testing demographic function
test_Pars <- c(T = 292, kappa = 0.8, Y = 1) 

### Temperature dependence functions of vital rates with interaction between Mass & Temp ####
# "Arrhenius-Lindmark" function (Boltzmann-Arrhenius with interaction (term cM) between Mass & Temp)
# gives Mass and temp dependence of Maintenance 
rM_T_AL <- function(T,m) { # 
  (m^(cM*(T-T0)))*
    exp(EaM*(T-T0)/(k*T*T0)) }

# Padfield unimodal function for unimodal consumption over temperature
rI_T_Pad <- function(T) { # GU Temp dependence of intake
  bTc * 
    exp(EaC*(T-T0)/(k*T*T0)) * (1+exp(EdI*(T-Td)/(k*T*Td)))^(-1)
  }

###  TEMP & MASS DEPENDENT RATE FUNCTION ####
# ratefun() uses dwdt() through ode() and returns the mass and rate for consumption 
# and maintenance through integration over one growing season. 
# ratefun() returns an 3d array where rows is days, cols is rates, dim is size in m_s
# if m_s is a single value, ratefun returns a matrix.
# Intake is not consumption but net energy available for growth, reproduction and maintenance

# dwdt <- function(time, m, Pars) {
#   with(as.list(c(m, Pars)), {
#   maintenance <- (rho1*(m^rho2)*rM_T_AL(T,m)) # Maintenance rate at time, mass with Pars
#   intake <- alpha*Y*eps1*(m^eps2)*rI_T_Pad(T) # Intake energy is not consumption
#   mass <- kap_fun(m,kappa)*intake - maintenance    # Growth rate at time, mass with Pars
#   return(list(mass, maintenance, intake)) # return all rates
#  })
# }

dwdt <- function(time, m, Pars) {
  with(as.list(c(m, Pars)), {
    maintenance <- (rho1*(m^rho2)*rM_T_AL(T,m)) # Maintenance rate at time, mass with Pars
    intake <- alpha*Y*eps1*(m^eps2)*rI_T_Pad(T) # Intake energy
    # Growth rate at time, mass with Pars:
    mass <- ifelse(kap_fun(m,kappa)*intake >= maintenance, kap_fun(m,kappa)*intake - maintenance,
                   ifelse(kap_fun(m,kappa)*intake < maintenance & maintenance <= intake, 0, intake - maintenance))
    return(list(mass, maintenance, intake)) # return all rates
  })
}

ratefun <- function(m, Pars) { # mean function for m_s+1 over one time step (one season). 
  yini <- m   # initial biomass in integration
  time <- seq(1, sl, by = 1) # time interval for an sl day long growing season, assumed for pike
  a <- sapply(yini, FUN=ode, time, func = dwdt, parms = Pars)#, simplify = "array")
  a <- array(as.numeric(unlist(a,F,F)), dim=c(sl, 4, length(yini)), 
             dimnames = list(c(1:(sl)), c("day","mass","maintenance","intake"), m))
  return(a) 
}

#VVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVV
### IPM Integral projection model ####
#VVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVV

### SIZES & PARAMETERS FOR IPM ####
mmin <- 1     # min weight
mmax <- 20000 # max weight
n <- 1400   # row/col number of the discretization-matrix of the continuous rates
x <- seq(mmin,mmax, length=n) # x is m_s in main text
dx <- x[2] - x[1] # step size (grams) in the continuous size spectra
dy <- dx # the IPM discretization is a square matrix, growth and age1 size is summed over y
x <- seq(mmin+dx/2,mmax+dx/2,dx) #for midpoint rule in IPM

e_m <- 0.00351388  # Egg mass calculated from Windermere pike: mean(FData$Egg.weight))) 
el_surv = 1.9e-4  # egg & larvae survival. Vindenes 2014 based on Kipling and Frost 1970

# Length-weight and weight-length relationship from Windermere
ltow  <-  function(l){exp(-6.49286)*l^3.4434} 
wtol <- function(w){(w/exp(-6.49286))^(1/3.4434)}

### DEMOGRAPHIC RATES FUNCTIONS ####

### GROWTH g(m_s+1; m_s, T) ####
# growthfun() returns a vector with a size distribution if given a single size m
# or if m is a vector, a matrix where columns are size distributions y from each entry in m.

growthfun <- function(m, Pars, y=x) {
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

### REPRODUCTION f(m_s,T) ####
# repfun uses the energy investment in reproduction dependent on mass and Temp 
# over s to produce eggs in s+1 for an individual that in the start of s is x. 
# Allocation to reproduction starts when an individual reaches maturation size, 
# even if that occurs during s -> s+1 making reproductive reserve in s+1 
# dependent on number of days as mature in s. Mean egg weight 0.00351388 in Lake Windermere data

repfun <- function(m, Pars, e_m = 0.00351388) { # Note that intake refers to net energy available for growth, reproduction and maintenance (i.e. alpha is accounted for)
   Devm = 417 # mean maturation size from data
   fecund <- function(m, Pars, e_m) {
     with(as.list(c(m, Pars)), {
      f <- ratefun(m, Pars)
      f <- ifelse(f[,"mass",] < Devm,  0, #is it mature? # see text for the following conditions...
             ifelse(kap_fun(f[,"mass",],kappa)*f[,"intake",] >= f[,"maintenance",], f[,"intake",]*(1-kap_fun(f[,"mass",],kappa)), #can available energy from intake cover maintenance, NOTE:intake is already scaled with Y*alpha
               ifelse(kap_fun(f[,"mass",],kappa)*f[,"intake",] < f[,"maintenance",] & (f[,"maintenance",] <= f[,"intake",]),
                     (1-kap_fun(f[,"mass",],kappa))*f[,"intake",] + kap_fun(f[,"mass",],kappa)*f[,"intake",] - f[,"maintenance",], 0)))
    sum(f)/e_m * 0.5 # Summed energy allocated to reproduction reserve in one season divided by egg weight and 0.5 from sex ratio (half of the eggs are female)
     })
   }
   sapply(X=m, FUN=fecund, Pars, e_m)
 }

### AGE 1 SIZE DISTRIBUTION o(e_m,T)####
# The age 1 size distribution is assumed to follow a lognormal density distribution.
# Mean size (mu) using ratefun() based on half a growing season.
# Variance from Windermere: sd(ltow((GData[GData$Age==1,]$Length)))

age1size <- function(Pars, y=x, offvar = 48.12555^2){ 
  mu <- ratefun(e_m, Pars = c(T = Pars[["T"]], kappa = Pars[["kappa"]], Y = Pars[["Y"]]))[round(sl/2),2,] 
  yd <- dlnorm(y, meanlog=log(mu) - .5*log(1 + offvar/(mu^2)), sdlog=sqrt(log(1 + offvar/(mu^2))))
      yd/sum(yd*dx) # scaled distributon to 1, dx is step size in the k-matrix
}

### SURVIVAL a(x,T) ####
survfun <- function(m, Pars) { # Mass-Temp dependence of yearly survival 
  Vsurvfun <- function(m, z=10.34){ # temperature independent survival function
    sxV <- function(m, z=10.34){ # sx from Vindenes 2014
      1/(1+exp(13.53316 - 0.50977*wtol(m) - (-0.00393)
               *wtol(m)^2 - 0.19312*z - (-0.00679)
               *wtol(m)*z))
    }
    xmax <- x[which(sxV(x) == max(sxV(x)))]
    ifelse(m < xmax, sxV(m), sxV(xmax))
  }
  VsurvfunT <- function(m, Pars){ # temperature dependent survival function
    zT=10.34+Pars[["T"]]-287 # make 287 equal to 283 survival in Vindenes 2014 since 283 is yearly mean Temp 287 is summer mean temp
    sxV <- function(m, z=zT){ # sx from Vindenes 2014
      1/(1+exp(13.53316 - 0.50977*wtol(m) - (-0.00393)
               *wtol(m)^2 - 0.19312*z - (-0.00679)
               *wtol(m)*z))
    }
    xmax <- x[which(sxV(x) == max(sxV(x)))]
    ifelse(m < xmax, sxV(m), sxV(xmax))
  }
  sx.firstyear <- function(m, Pars){ # first year survival
    el_surv }
  # Choose survival scenario below, compare fig. 4 main text.
  ifelse(m < mmin, sx.firstyear(m,Pars), VsurvfunT(m,Pars)) # main, i.e. size and temp depedent, survival
  #ifelse(m < mmin, sx.firstyear(m,Pars), Vsurvfun(m)) # temperature independent survival
  #ifelse(m < mmin, sx.firstyear(m,Pars), 0.68) # consant, i.e. size and temp INdepedent, survival
} 

### Projection Matrix ####
# The projection or Kernel (K) matrix maps the size distribution in time s to time s+1. 
# Smat column 1 describes the probability of an egg to hatch and grow into size y in s+1 or the density distribution of offspring size y, 
# the other columns of Smat describe the growth and survival probability density of y (from x)
# the first row is 0 to give way for an egg stage (the number of eggs produced by the sizes in x) through Fmat because at census for x, size x will give eggs at s+1 which will be offspring in s+2
# Fmat is the egg stage, number of eggs produced for each x which becomes y1.

K.matrix <- function(Pars) {
  Smat <- Fmat <- matrix(0,n+1,n+1)
  surv_x <- c(survfun(e_m, Pars), survfun(x, Pars))  # survival vector
  Smat[2:(n+1),1] <- surv_x[1]*age1size(Pars) * dx #dx is dy, First column (and 2:101 row) of Smat probability of transitions from egg to sizes at age 1.  
  Smat[2:(n+1),2:(n+1)] <- t(surv_x[2:(n+1)]*t(growthfun(x, Pars) ) * dx) # integrate over y as is rows, thus the transponate of DEBgrowth (???)
  Fmat[1,2:(n+1)] <- surv_x[2:(n+1)]*repfun(x, Pars) #Production of eggs
  Smat+Fmat 
}

# vwlambda.projection:function to estimate lambda, v and w
wvlambda.projection <- function(Kmat, N0=rep(10, length(Kmat[1,])), tol=1e-20) {
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
  w <- w/sum(w)
  v <- Nt2/sum(Nt.new2)
  v <- v/sum(v*w)
  return(list("lambda"=lam, "w"=w,"v"=v))
}

# Test the DEBIPM model for timing of each run and taht it sums to 1
# tic()
# onetest <- wvlambda.projection(K.matrix(test_Pars))$lambda
# toc()
# sum(onetest$w) 

### SENSITIVITY OF THE DEBIPM ####
## Numerical expression for the derivatives for mass and temperature specific 
## demographic functions with respect to $\kappa_0$...

# Growth function pert.
G_NumDer <- function(h ,m, Pars){
  (growthfun(m, Pars = c(T = Pars[["T"]],kappa = Pars[["kappa"]] + h,Y = Pars[["Y"]])) -
     growthfun(m, Pars = c(T = Pars[["T"]],kappa = Pars[["kappa"]],Y = Pars[["Y"]])))/h
} 
# Fecundity function pert.
F_NumDer <- function(h ,m, Pars){
  (repfun(m, Pars = c(T = Pars[["T"]],kappa = Pars[["kappa"]] + h,Y = Pars[["Y"]])) - 
     repfun(m, Pars = c(T = Pars[["T"]],kappa = Pars[["kappa"]],Y = Pars[["Y"]])))/h
}
# Age 1 size function pert.
A1_NumDer <- function(h ,m, Pars){
  (age1size(Pars = c(T = Pars[["T"]],kappa = Pars[["kappa"]] + h,Y = Pars[["Y"]])) - 
     age1size(Pars = c(T = Pars[["T"]],kappa = Pars[["kappa"]],Y = Pars[["Y"]])))/h
}
# the other calculations for the sensitivity analyses are found in the DEBIPM_plot script

# to perturb the kernel
K_NumDer <- function(h, Pars){
  (K.matrix(Pars = c(T = Pars[["T"]],kappa = Pars[["kappa"]] + h,Y = Pars[["Y"]])) -
     K.matrix(Pars = c(T = Pars[["T"]],kappa = Pars[["kappa"]],Y = Pars[["Y"]])))/h
}
