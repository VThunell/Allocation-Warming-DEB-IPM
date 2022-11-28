
############### THUNELL ET AL. DEB-IPM R-script - MODEL #########################
  
data.g <-read.csv("PikeGrowthData1944_1995.csv", sep=",")
data.f <-read.csv("Windermere_Pike_Fecundity_and_Egg_Data_1963_to_2003.csv", sep=",")

# Packages ####

library(deSolve) 
library(tidyverse) 
#install.packages("tictoc")
library(tictoc) 

# DEB Dynamic energy budget model ####

# Paramters kappa,ha,eps1,eps2 from optimization
load("optim_Pars.RData")

# Allometric parameters 
rho1 <- 0.02 # Maintenance allometric scalar 
rho2 <- 0.76 # Maintenance allometric exponent
eps1 <- deb.optim$par[3] # Consumption allometric scalar
eps2 <- deb.optim$par[4] # Consumption allometric exponent

# DEB-parameters 
alpha <- 0.4   # assimilation efficiency
sl <- 183      # Season length in days, i.e. number of growth time steps
kap_fun <- function(m, kappa, ha=deb.optim$par[2]){kappa*exp(-m/(ha*max(x)))} #for body size dep. kappa
# kap_fun <- function(m, kappa, ha=2){kappa} #for body size indep. kappa

# Temperature scaling parameter values from Lindmark et al. 2022, Global change biol.
T0 <-  292            # Reference Temp
k  <-  8.617333e-05   # Boltzmann constant
EaC <- 0.73 # activation energy Consumption
EaM <- 0.62 # activation energy Maintenance
EdC <- 1.89 # deactivation energy
Td  <- T0 + 0.75 # temperature at which half the rate is reduced due to temperature
bTc <- 1.824953 #0.79/rC_T_Pad(292) rate at reference (a common) temperature Sharpe-Schoolf

# Parameters for temperature, allocation and feeding level  Y (which is not used in analysis and excluded the main text) for testing demographic function
test_Pars <- c(T = 287, kappa = deb.optim$par[1], Y = 1) 

# Temperature dependence ####
# Boltzmann-Arrhenius function
rM_T_A <- function(T) {  
  exp(EaM*(T-T0)/(k*T*T0)) }

# Sharpe-Schoolfield (as in Padfield) for unimodal consumption over temperature
rC_T_Pad <- function(T) { 
  bTc * exp(EaC*(T-T0)/(k*T*T0)) * (1+exp(EdC*(T-Td)/(k*T*Td)))^(-1) }

###  TEMP & MASS DEPENDENT RATE FUNCTION ####
# ratefun() uses dwdt() through ode() and returns the mass and rate for consumption 
# and maintenance through integration over one growing season. 
# ratefun() returns an 3d array where rows is days, cols is rates, dim is size in m_s
# if m_s is a single value, ratefun returns a matrix.
# asim_energy ís energy available for growth, reproduction and maintenance
dwdt <- function(time, m, Pars) {
    maintenance <- (rho1*(m^rho2)*rM_T_A(Pars[["T"]])) # Maintenance rate at time, mass with Pars
    asim_energy <- alpha*eps1*(m^eps2)*rC_T_Pad(Pars[["T"]])
    mass <- ifelse(kap_fun(m,Pars[["kappa"]])*asim_energy >= maintenance, kap_fun(m,Pars[["kappa"]])*asim_energy - maintenance,
                   ifelse(kap_fun(m,Pars[["kappa"]])*asim_energy < maintenance & maintenance <= asim_energy, 0, asim_energy - maintenance))
    return(list(mass, maintenance, asim_energy)) # return all rates to ratefun
}

ratefun <- function(m, Pars) { # mean function for m_s+1 over one time step (one season). 
  yini <- m   # initial biomass in integration
  time <- seq(1, sl, by = 1) # time interval for an sl day long growing season, assumed for pike
  a <- sapply(yini, FUN=ode, time, func = dwdt, parms = Pars)
  a <- array(as.numeric(unlist(a,F,F)), dim=c(sl, 4, length(yini)), 
             dimnames = list(c(1:(sl)), c("day","mass","maintenance","asim_energy"), m))
  return(a) 
}

# IPM Integral projection model ####

# SIZES & PARAMETERS FOR IPM ####
mmin <- 1     # min weight
mmax <- 20000 # max weight
n <- 1000  # the size class bin width, or matrix mesh size
x <- seq(mmin,mmax, length=n) # x is m_s in main text
dx <- x[2] - x[1] # step size (grams) in the continuous size spectra
dy <- dx # the IPM discretization is a square matrix, growth and age1 size is summed over y
x <- seq(mmin+dx/2,mmax+dx/2,dx) #for midpoint rule in IPM

e_m <- 0.00351388  # Egg mass calculated from Windermere pike: mean(data.f$Egg_weight))) 
el_surv = 1.9e-4  # egg & larvae survival. Vindenes 2014 based on Kipling and Frost 1970

# Length-weight and weight-length relationship from Windermere
ltow  <-  function(l){exp(-6.49286)*l^3.4434} 
wtol <- function(w){(w/exp(-6.49286))^(1/3.4434)}

# DEMOGRAPHIC RATE FUNCTIONS ####

## GROWTH g(m_s+1; m_s, T) ####
# growthfun() returns a vector with a size distribution if given a single size m
# or if m is a vector, a matrix where columns are size distributions y from each entry in m.
growthfun <- function(m, y_g, Pars, y=x) {
  xsd <- function(mu, nu.var = -0.0001,  sdres = 300){ # sd for mu
    sdres * exp(nu.var*mu) }
  ydf <- function(mu, y=x) { # m:s size distributions for next year
    var <- xsd(mu)^2 # variance
    yd <- dlnorm(y, meanlog = log(mu) - .5*log(1 + var/mu^2), sdlog = sqrt(log(1 + var/mu^2))) #Distribution of y
    if(sum(yd*dx) == 0){
      c(rep(0, n-1), 1/dx)
    } else { yd/sum(yd*dx) } } #scale to one
  mu <- y_g[sl,2,] # mass on the last day of growth season for m
  sapply(mu, FUN = ydf)
}

### REPRODUCTION f(m_s,T) ####
# Allocation to reproduction starts when an individual reaches maturation size, 
# even if that occurs during s -> s+1 making reproductive reserve in s+1 
# dependent on number of days as mature in s. Mean egg weight 0.00351388 in Windermere data
repfun <- function(m, y_g, Pars, e_m = 0.00351388) { #
  Devm = 417 # mean maturation size from data
  with(as.list(c(m, y_g, Pars)), {
    f=NULL
    for(i in 1:length(m)){
      ph = ifelse(y_g[,"mass",i] < Devm,  0, #is it mature? # see text for the following conditions...
                  ifelse(kap_fun(y_g[,"mass",i],kappa)*y_g[,"asim_energy",i] >= y_g[,"maintenance",i], y_g[,"asim_energy",i]*(1-kap_fun(y_g[,"mass",i],kappa)), #can available energy from asim_energy cover maintenance, NOTE:asim_energy is already scaled with Y*alpha
                         ifelse(kap_fun(y_g[,"mass",i],kappa)*y_g[,"asim_energy",i] < y_g[,"maintenance",i] & (y_g[,"maintenance",i] <= y_g[,"asim_energy",i]),
                                (1-kap_fun(y_g[,"mass",i],kappa))*y_g[,"asim_energy",i] + kap_fun(y_g[,"mass",i],kappa)*y_g[,"asim_energy",i] - y_g[,"maintenance",i], 0)))
      f=c(f,sum(ph)/e_m * 0.5) # Summed energy allocated to reproduction reserve in one season divided by egg weight and 0.5 from sex ratio (half of the eggs are female)
    } 
    return(f)
  }) 
}

### AGE 1 SIZE DISTRIBUTION o(e_m,T)####
# Mean size (mu) using ratefun() based on half a growing season.
# Variance from Windermere: sd(ltow((data.g[data.g$Age==1,]$Length)))
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
    zT=10.34+Pars[["T"]]-287 # make 287 equal to 283 survival in Vindenes 2014 since 283 is yearly mean temp, 287 is summer mean temp which we use for the growth function
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
  y_g = ratefun(x, Pars)
  Smat = Fmat = matrix(0,n+1,n+1)
  surv_x = c(survfun(e_m, Pars), survfun(x, Pars))  # survival vector
  Smat[2:(n+1),1] = surv_x[1]*age1size(Pars) * dx #dx is dy, First column (and 2:101 row) of Smat probability of transitions from egg to sizes at age 1.  
  Smat[2:(n+1),2:(n+1)] = t(surv_x[2:(n+1)]*t(growthfun(x, y_g, Pars) ) * dx) # integrate over y as is rows, thus the transponate of DEBgrowth (???)
  Fmat[1,2:(n+1)] = surv_x[2:(n+1)]*repfun(x, y_g, Pars) #Production of eggs
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
# sum(onetest$w) # Hold my breath as I wish for death
# onetest$lambda # mean fitness for test_Pars

### SENSITIVITY OF THE DEBIPM ####
## Numerical expression for the derivatives for mass and temperature specific 
## demographic functions with respect to kappa_0

# Growth function pert.
G_NumDer <- function(h ,m, Pars){
  (growthfun(m, y_g = ratefun(m, Pars = c(T = Pars[["T"]],kappa = Pars[["kappa"]] + h,Y = Pars[["Y"]])),
             Pars = c(T = Pars[["T"]],kappa = Pars[["kappa"]] + h,Y = Pars[["Y"]])) -
     growthfun(m, y_g = ratefun(m, Pars = c(T = Pars[["T"]],kappa = Pars[["kappa"]],Y = Pars[["Y"]])),
               Pars = c(T = Pars[["T"]],kappa = Pars[["kappa"]],Y = Pars[["Y"]])))/h
} 
# Fecundity function pert.
F_NumDer <- function(h ,m, Pars){
  (repfun(m, y_g = ratefun(m, Pars = c(T = Pars[["T"]],kappa = Pars[["kappa"]] + h,Y = Pars[["Y"]])),
          Pars = c(T = Pars[["T"]],kappa = Pars[["kappa"]] + h,Y = Pars[["Y"]])) - 
     repfun(m, y_g = ratefun(m, Pars = c(T = Pars[["T"]],kappa = Pars[["kappa"]],Y = Pars[["Y"]])),
            Pars = c(T = Pars[["T"]],kappa = Pars[["kappa"]],Y = Pars[["Y"]])))/h
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
