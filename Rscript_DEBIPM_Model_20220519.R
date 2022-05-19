
############### THUNELL ET AL. 20XX DEBIPM R-script - MODEL #########################
  
YVPike <- load("PikeDataFiles.R")

### Packages ####
#install.packages("deSolve")
library(deSolve) 
#install.packages("tidyverse")
library(tidyverse) 
#install.packages("grid")
library(tictoc) 

#VVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVV
### DEB Dynamic energy budget model ####
#VVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVV

## Allometric parameters 
rho1 <- 0.02 # Maintenance allometric scalar 
rho2 <- 0.76 # Maintenance allometric exponent
eps1 <- 1.55 # Consumption allometric scalar
eps2 <- 0.54 # Consumption allometric exponent

## DEB-parameters 
alpha <- 0.4   # assimilation efficiency
sl <- 183      # Season length in days, i.e. number of growth time steps
kap_fun <- function(m, kappa, ha=2){kappa*exp(-m/(ha*max(x)))} #for body size dep. kappa

## Temperature scaling parameter values, Lindmark et al. 2022, Global change biol.
T0 <-  292            # Reference Temp
k  <-  8.617333e-05   # Boltzmann constant
EaC <- 0.73 # activation energy Consumption
EaM <- 0.62 # activation energy Maintenance
EdI <- 1.89 # deactivation energy
Td  <- T0 + 0.75 # temperature at which half the rate is reduced due to temperature
bTc <- 0.79 # rate at reference (a common) temperature Sharpe-Schoolf

# Parameters for temperatur, allocation and feeding level (not used in analysis and excluded the main text) for testing demographic function
test_Pars <- c(T = 287, kappa = 0.8, Y = 1) 

### Temperature dependence ####

# Boltzmann-Arrhenius function
rM_T_A <- function(T) {  
  exp(EaM*(T-T0)/(k*T*T0)) }

# Sharpe-Schoolfield (as in Padfield) for unimodal consumption over temperature
rI_T_Pad <- function(T) { # GU Temp dependence of asim_energy
  bTc * exp(EaC*(T-T0)/(k*T*T0)) * (1+exp(EdI*(T-Td)/(k*T*Td)))^(-1) }

###  TEMP & MASS DEPENDENT RATE FUNCTION ####
# ratefun() uses dwdt() through ode() and returns the mass and rate for consumption 
# and maintenance through integration over one growing season. 
# ratefun() returns an 3d array where rows is days, cols is rates, dim is size in m_s
# if m_s is a single value, ratefun returns a matrix.
# asim_energy ís energy available for growth, reproduction and maintenance

dwdt <- function(time, m, Pars) {
  with(as.list(c(m, Pars)), {
    maintenance <- (rho1*(m^rho2)*rM_T_A(T)) # Maintenance rate at time, mass with Pars
    asim_energy <- alpha*Y*eps1*(m^eps2)*rI_T_Pad(T)
    mass <- ifelse(kap_fun(m,kappa)*asim_energy >= maintenance, kap_fun(m,kappa)*asim_energy - maintenance,
                   ifelse(kap_fun(m,kappa)*asim_energy < maintenance & maintenance <= asim_energy, 0, asim_energy - maintenance))
    return(list(mass, maintenance, asim_energy)) # return all rates ti ratefun
  })
}

ratefun <- function(m, Pars) { # mean function for m_s+1 over one time step (one season). 
  yini <- m   # initial biomass in integration
  time <- seq(1, sl, by = 1) # time interval for an sl day long growing season, assumed for pike
  a <- sapply(yini, FUN=ode, time, func = dwdt, parms = Pars)
  a <- array(as.numeric(unlist(a,F,F)), dim=c(sl, 4, length(yini)), 
             dimnames = list(c(1:(sl)), c("day","mass","maintenance","asim_energy"), m))
  return(a) 
}

#VVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVV
### IPM Integral projection model ####
#VVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVV

### SIZES & PARAMETERS FOR IPM ####
mmin <- 1     # min weight
mmax <- 20000 # max weight
n <- 1400   # the size class bin width, or matrix mesh size
x <- seq(mmin,mmax, length=n) # x is m_s in main text
dx <- x[2] - x[1] # step size (grams) in the continuous size spectra
dy <- dx # the IPM discretization is a square matrix, growth and age1 size is summed over y
x <- seq(mmin+dx/2,mmax+dx/2,dx) #for midpoint rule in IPM

e_m <- 0.00351388  # Egg mass calculated from Windermere pike: mean(FData$Egg.weight))) 
el_surv = 1.9e-4  # egg & larvae survival. Vindenes 2014 based on Kipling and Frost 1970

# Length-weight and weight-length relationship from Windermere
ltow  <-  function(l){exp(-6.49286)*l^3.4434} 
wtol <- function(w){(w/exp(-6.49286))^(1/3.4434)}

### DEMOGRAPHIC RATE FUNCTIONS ####

### GROWTH g(m_s+1; m_s, T) ####
# growthfun() returns a vector with a size distribution if given a single size m
# or if m is a vector, a matrix where columns are size distributions y from each entry in m.

growthfun <- function(m, Pars, y=x) {
  xsd <- function(mu, nu.var = -0.0001,  sdres = 300){ # sd for mu
    sdres * exp(nu.var*mu) }
  ydf <- function(mu, y=x) { # m:s size distributions for next year
    var <- xsd(mu)^2 # variance
    yd <- dlnorm(y, meanlog = log(mu) - .5*log(1 + var/mu^2), sdlog = sqrt(log(1 + var/mu^2))) #Distribution of y
    if(sum(yd*dx) == 0){
      c(rep(0, n-1), 1/dx)
    } else { yd/sum(yd*dx) } } #scale to one
  mu <- ratefun(m, Pars)[sl,2,] # mass on the last day of growth season for m
  sapply(mu, FUN = ydf)
}

### REPRODUCTION f(m_s,T) ####
# Allocation to reproduction starts when an individual reaches maturation size, 
# even if that occurs during s -> s+1 making reproductive reserve in s+1 
# dependent on number of days as mature in s. Mean egg weight 0.00351388 in Lake Windermere data

repfun <- function(m, Pars, e_m = 0.00351388) { #
   Devm = 417 # mean maturation size from data
   fecund <- function(m, Pars, e_m) {
     with(as.list(c(m, Pars)), {
      f <- ratefun(m, Pars)
      f <- ifelse(f[,"mass",] < Devm,  0, #is it mature? # see text for the following conditions...
             ifelse(kap_fun(f[,"mass",],kappa)*f[,"asim_energy",] >= f[,"maintenance",], f[,"asim_energy",]*(1-kap_fun(f[,"mass",],kappa)), #can available energy from asim_energy cover maintenance, NOTE:asim_energy is already scaled with Y*alpha
               ifelse(kap_fun(f[,"mass",],kappa)*f[,"asim_energy",] < f[,"maintenance",] & (f[,"maintenance",] <= f[,"asim_energy",]),
                     (1-kap_fun(f[,"mass",],kappa))*f[,"asim_energy",] + kap_fun(f[,"mass",],kappa)*f[,"asim_energy",] - f[,"maintenance",], 0)))
    sum(f)/e_m * 0.5 # Summed energy allocated to reproduction reserve in one season divided by egg weight and 0.5 from sex ratio (half of the eggs are female)
     })
   }
   sapply(X=m, FUN=fecund, Pars, e_m)
 }

### AGE 1 SIZE DISTRIBUTION o(e_m,T)####
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
    zT=10.34+Pars[["T"]]-287 # make 287 equal to 283 survival in Vindenes 2014 since 283 is yearly mean temp, 287 is summer mean temp whihc we use for the growth function
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

# Test the DEBIPM model for timing of each run and check that it sums to 1
# tic()
# onetest <- wvlambda.projection(K.matrix(test_Pars))
# toc()
sum(onetest$w) # Hold my breath as I wish for death

### SENSITIVITY OF THE DEBIPM ####
## Numerical expression for the derivatives for mass and temperature specific 
## demographic functions with respect to kappa_0

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
