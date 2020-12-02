
############### Thunell et al. DEB-IPM script #########################

# This script contains the DEB model describing mean somatic growth and reproduction
# and functions for vital rates needed for the IPM: variance function for somatic growth,
# a size and temperature dependent survival function, function for temperature 
# dependent mean offspring size distrubution, variance function for offspring size. 
# Lastly, the script contains the IPM kernel and projection matrix, calculation of lambda, 
# stable population structure and reproductive values

setwd("C:/Users/vitl0001/VThunell_repos/Temperature-DEBIPM")

## Windermere Pike data
#YVPike <- load("~/Manus2/R/Manus2R/PikeDataFiles.R")
# str(YVPike)

### Packages ####
#install.packages("deSolve")
library(deSolve)   # for ode Solver
#install.packages("tidyverse")
library(tidyverse) # gglot etc.
#install.packages("grid")
library(grid)      # for grid.text
#install.packages("fields")
library(fields)      # for image.plot

#VVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVV
### DEB Dynamic energy budget model ####
#VVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVV

#### Mass dependence functions of vital rates #####
# a*(Mass^b), a = eps1 or rho1, b = eps2 or rho2

# Allometric parameters (mass dependence at the reference temperature)
rho1 = 0.009  # Maintenance allometric scalar, rescaled from 0.02 based on pike estimated by Lindmark from Armstrong 1992.
rho2 = 0.87   # Maintenance allometric exponent, rescaled from 0.76 in Lindmark 2019 based on Diana 1982(?)
eps1 = 0.64   # Intake allometric scalar, free parameter to fit grwth in Windermere pike.
eps2 = 0.54   # Intake allometric exponent, free parameter to fit grwth in Windermere pike.

# General DEB-parameters
Kap  = 0.8    # Kappa
alpha = 0.4   # assimilation efficiency
Y = 1         # feeding level
Y2 = 0.4      # feeding level for mass < 1
sl = 183      # Season length in days
#HR_exp = 0.99    # Hyper allometric relationship of kappa, e.g. Kappa = K*(Mass^HR_exp)/Mass
#plot(Mass, 0.8*(Mass^HR_exp)/Mass, ylab = "Kappa", type="l", main = "Hyperallometric scaling of 1-Kappa")

### Temperature dependence functions of vital rates with interaction between Mass & Temp ####

# Arrhenius-Lindmark function (AL, Boltzmnn-Arhhenius with interaction betwen Mass & Temp)
rM_T_AL2 <- function(T, m, T_p) { # AL Temp dependence of Maintenance 
  (m^(T_p[["cM"]]*(T-T0)))*exp(T_p[["EaM"]]*(T-T0)/(k*T*T0)) }

# Gardmark Unimodal function
rI_T_GU2 <- function(T, m, T_p) { # GU Temp dependence of intake 
  (m^(T_p[["cIa"]]*(T-T0))) * exp(T_p[["EaI"]]*(T-T0)/(k*T*T0)) *
    (m^(-(T_p[["cId"]]*(T-T_p[["Td"]])))) * (1+exp(T_p[["EdI"]]*(T-T_p[["Td"]])/(k*T*T_p[["Td"]])))^(-1) * 
    (m^(T_p[["cId"]]*(T0-T_p[["Td"]]))) * (1+exp(T_p[["EdI"]]*(T0-T_p[["Td"]])/(k*T0*T_p[["Td"]]))) }

# General temperature parameter values
T0 =  283             # Reference Temp
k  =  8.617332e-05    # Boltzmann constant

# Temperature parameters for intake and maintenance
T_par <- c(cIa = 0,     # linear interaction between size and temp for Max intake.
           cM =  0.017, # linear interaction between size and temp for Maintenance. 0.017 Lindmark unpub.
           EaI = 0.66,  # activation energy Intake, Lindmark unpub. 0.66
           EaM = 0.61,  # activation energy Maintenance, Lindmark unpub: 0.61
           Td  = T0+6,  # deactivation temperature. 6 degrees in Lindmark unpub.?
           cId = -0.02, # linear interaction effect (slope) between temp and mass for deactivation. -0.02 Guesstimate
           EdI = 2)     # deactivation energy. 2 is guesstimate

###  TEMP & MASS DEPENDENT RATE FUNCTION ####
# ratefun() uses dwdt() through ode() and returns the mass and rate for intake and maintenance dependent on mass through
# integration over one growing season.
# Intake here refers to net energy available for growth, reproduction and maintenance

dwdt <- function(Time, Mass, Pars) { # ode() requires that Time, mass and pars is passed to growthfunction
  with(as.list(c(Mass, Pars)),{
    Maintenance <- (rho1*(Mass^rho2)*rM_T_AL2(T, Mass, T_par))   # Maintenance rate at Time, mass with Pars
    ifelse(Mass > 1, Intake <- (alpha*Y*eps1*(Mass^eps2)*rI_T_GU2(T, Mass, T_par)),
           Intake <- (alpha*Y2*eps1*(Mass^eps2)*rI_T_GU2(T, Mass, T_par))) # Intake energy rate at Time, mass with Pars, NOTE intake rate multiplied with alpha and Y
    G <- Kappa*Intake - Maintenance    # Growth rate at Time, mass with Pars
     return(list(G, Maintenance, Intake)) # return all rates
  })
}

ratefun <- function(x, Pars) { # mean function for y (dependent on mass) over one time step (one season). 
  yini <- x    # initial biomass in integration
  time <- seq(0, sl, by = 1) # time interval for sl day long growing season, assumed for pike
  a <- ode(y = yini, time, func = dwdt, parms = Pars)
  b <- data.frame("day" = a[,1],
                  "mass" = a[,2],
                  "maintenance"= a[,3],
                  "intake" = a[,4])
  return(b)
}

#VVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVV
### IPM Integral projection model ####
#VVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVV

### SIZES & PARAMETERS FOR IPM ####
mmin = 1     # min weight
mmax = 14000 # max weight
n = 150      # row/col number of the discretization-matrix of the continuos rates
x <- seq(mmin,mmax, length=n)
dx <- x[2] - x[1] # step size (grams) in the continuos size spectra

o_e = 0.00351388  # Egg weight calculated from windermere pike: mean(FData$Egg.weight))) 
egg_surv = 1e-5 # egg survival. losely based on Kipling and Frost 1970 1/50000 from laid egg to age 2. 1.9e-4 in Vindenes 2014

ltow_A  <- function(l){0.00254*l^3.271} #length to weight function adults
# Lenght weight relationship from average in fishbase
ltow_J <- function(l){0.01101*l^2.69} #length to weight function juveniles
# Lenght weight relationship from fishbase (Carlander 1969) 

### VITAL RATES FUNCTIONS ####

### GROWTH g(x,y,T) ####
DEBgrowthfun <- function(x, y, Pars) {
  mu <- ratefun(x, Pars)[sl+1,2]
  xsd <- function(x, nu.var = -0.0001,  sdres = 100){
          sdres *exp(nu.var*x)
        }
  var <- xsd(mu)^2  # the variation in x becomes the variation in y-x
  yd <- dlnorm(y, meanlog = log(mu) - .5*log(1 + var/mu^2), sdlog = sqrt(log(1 + var/mu^2))) #Distribution of y
  if(sum(yd*dx) == 0){
    return (c(rep(0, n-1), 1/dx))
    } else { return(yd/sum(yd*dx)) } #scale it to one
  
}

### REPRODUCTION f(x,T) ####
# DEBrepfun uses the energy investment in reproduction dependent on mass and T 
# over t to produce eggs in t+1 for an individual that in the start of t is x. 
# Allocation to reproduction starts when an individuals reaches maturation size, 
# even if that occurs during t -> t+1 making reproductive reserve in t+1 
# dependent on number of days as mature in t.
# Mean egg weight 0.00351388 in Lake Windermere data

 DEBrepfun <- function(x, Pars, o_E) { # Note that intake refers to net energy available for growth, reproduction and maintenance
  R <- ratefun(x, Pars) 
  Devm = 374 # *(1 - pars[["Kappa"]])
  R <- ifelse(R$mass < Devm,  0, #is it mature? # see text for the following conditions...
              ifelse(Pars[["Kappa"]]*R$intake >= R$maintenance, R$intake*(1-Pars[["Kappa"]]), #can available energy from intake cover maintenance, NOTE:intake is already scaled with Y*alpha
                     ifelse((Pars[["Kappa"]]*R$intake < R$maintenance) & (R$maintenance <= R$intake), 
                            R$intake - R$maintenance, 0))) # can available energy from intake and reproduction cover maintenance, if not 0
                            # expression for this reproduction energy in the text is simplifiedhere to I(m,T)-M(m,T)
  ifelse(x <= 1, return(0),
         (sum(R)/o_E)*0.5) # Summed energy allocated to reproduction in one season diveded by egg weight and 0.5 from sex ratio (half of the eggs are female)
}

### OFFSPRING SIZE DISTRIBUTION o(eggsize,T)####
# Probability of an egg to hatch and grow into size y in t+1?
# From Vindenes 2014: "Density distribution for offspring size y (we assume parental length x does not affect offspring length)."
# Temp dependent mean weight of offspring is calculated from YV 2014 linear relationship
# Variance from Windermere: sd(ltow_J((GData[GData$Age==1,]$Length))) (3.72 for length in Vindenes 2014)

DEBoff.size <- function(y, Pars, offvar = 24.05458^2){ #lognormal distribution of offspring lengths
  mu <- ratefun(o_e, Pars)[sl/2,2] #mean offspring length is mass after 134 days of growth (i.e. 184-50 days)
  yd <- dlnorm(y, meanlog=log(mu) - .5*log(1 + offvar/(mu^2)), sdlog=sqrt(log(1 + offvar/(mu^2))))
      yd/sum(yd*dx) # scaled distributon to 1, dx is step size in the k-matrix
}

### SURVIVAL s(x,T) ####
DEBsurvfun <- function(Mass, Pars){ # Mass-Temp dependence of yearly Mortality 
       sx <- function(m, Pars){ # natural baseline survival  
         exp((-3*m^-.288)*exp(0.47*(Pars[['T']]-T0)/(k*Pars[['T']]*T0))) #3*^-.288 is yearly mortality rate for 1 gram unit weight from Lorenzen 1996 
       } 
       starvx <- function(m,Pars) { # starvation survival
         starvS <- NULL 
         for (i in m){
           starvS <- c(starvS, ifelse(ratefun(i, Pars)[sl+1,2] < i, 0, 1))
           } 
         return(starvS)
       }
      sx(Mass,Pars)*starvx(Mass, Pars)# multiply with F (fishing mortality function)
}

### Projection Matrix ####
# The projection or Kernel (K) matrix maps the size distribution in time t to time t+1. 

# Smat column 1 describes the probality of an egg to hatch and grow into size y in t+1 or the density distribution of offspring size y, 
# the other columns of Smat describe the growth and survival probability density of y (from x)
# the first row is 0 to give way for an egg stage (the number of eggs produced by the sizes in x) through Fmat because at census for x, size x will give eggs at t+1 whihc will be offspring in t+2

# Fmat is the egg stage, number of eggs produced for each x which becomes y1.

K.matrix <- function(Pars) {
  Smat <- Fmat <- matrix(0,n+1,n+1) 
  sx <- c(egg_surv, DEBsurvfun(x, Pars)) # survival vector
  Smat[2:(n+1),1] <- sx[1]*DEBoff.size(y=x, Pars) * dx # First column (and 2:101 row) of Smat is sum of survival of laid eggs, larvae and fry to age 1 (offspring size),  
  for (i in 2:(n+1)) {
    Smat[2:(n+1),i] <- sx[i]*DEBgrowthfun(x[i-1], y=x, Pars) * dx
    Fmat[1,i] <- sx[i]*DEBrepfun(x[i-1], Pars, o_e) #Production of eggs from each size class, (the number of eggs that x in t (size in april) will release in april in t+1) * (probability to grow from egg to age 1 (in t+2-1)
  }
   Smat[,n+1] <- Smat[,n] # TO MAKE UP FOR THE EGGSTAGE that Fmat adds?
   Kmat <- Smat+Fmat
}

# This function uses the R function "eigen" to calculate lambda, a stable structure "w" and reproductive values "v", for a projection matrix "Kmat". The stable structure is scaled so that sum(w*dx)=1 and the reproductive values are scaled so that sum(v*w*dx)=1 (see details in the article).

# For IPMs, including threshold for setting value to zero
wvlambda <- function(Kmat, tol=1e-20){
  ev <- eigen(Kmat)
  tev <- eigen(t(Kmat))
  lmax <- which.max(Re(ev$values))
  W <- ev$vectors
  V <- tev$vectors
  w <- as.matrix(abs(Re(W[, lmax]))/sum(abs(Re(W[, lmax]))))
  w <- ifelse(w<tol,0,w)
  w <- w/(sum(w*dx))
  v <- as.matrix(abs(Re(V[, lmax])))
  v <- v/sum(w*v*dx)
  v <- ifelse(w*v <= 0, 0, v)
  return(list("lambda" = max(Re(ev$values)), "w"=w, "v"=v))
}


###CALCULATE LAMBDA, STABLE STRUCTURE AND REPRODUCTIVE VALUES FOR VARYING KAPPA AND TEMP VALUES ####

Kappa <- seq(0.5,1,0.05)    # Kappa values to explore
T <- 283:286  # temperature range
Res_lamA = NULL
Res_vA = NULL
Res_wA = NULL

for (d in T) {
  for (e in Kappa) {
    
    parsK <- c(T = d,         # parameters for Temperature, feeding, allocation and Mass dependence
               Kappa = e,     # allocation to respiration (Growth and maintenance)
               Y = Y,         # Feeding level
               alpha=alpha,     # assimilation efficiency
               eps1=eps1,      # Intake allometric scalar, 0.248 based on roach by Lindmark assuming that intake scales equally with mass between species in his system
               eps2=eps2,      # Intake allometric exponent, 0.64 from Lindmark unpub. 0.725 to get an length asymtopte Lm for pike at 292 K. 0.767 based on roach by Lindmark assuming that intake scales equally with mass between species in his system
               rho1=rho1,      # Maintenance allometric scalar, Guesstimate 0.0204. 0.02 based on pike estimated by Lindmark
               rho2=rho2)
    
    K <- K.matrix(parsK)
    res <- wvlambda(K)
    Res_lamA <- rbind(Res_lamA, c(parsK, lam = res$lambda))
    Res_vA <- rbind(Res_vA, c(parsK, v = res$v)) #REPRODUCTIVE VALUES
    Res_wA <- rbind(Res_wA, c(parsK, w = res$w)) #STABLE STRUCTURE
    
    }
}
Res_lam <- as.data.frame(Res_lamA)  # Make a df instead of matrix
Res_v <- as.data.frame(Res_vA) # Make a df instead of matrix
Res_w <- as.data.frame(Res_wA) # Make a df instead of matrix
#colnames(Res_v) <- c("T","Y", "Kappa","alph", "eps1", "eps2", "rho1","rho2", as.character(round(x))) #minus recruits in x
#colnames(Res_w) <- c("T","Y", "Kappa","alph", "eps1", "eps2", "rho1","rho2", as.character(round(x))) #minus recruits in x

#rowSums(Res_w[,9:ncol(Res_w)])*dx

