# Code to fit DEB parameters

# Created: October 3, 2022 by EPD
# Last modified: October 3, 2022 by EPD

# Set working directory
#setwd("/Users/epdus/OneDrive/Breathless/Code/Packets/deb/")

# Load packages
library(deSolve)
library(optimParallel)
library(tidyverse)

# Contents (ctrl-f):
#	0a. Common values
#	0b. Common functions
#	I. Load data
#	II. Process data
#	III. Set parameters
#	IV. Optimization function
#	V. Optimize
#	VI. Bootstrap confidence intervals
#	VII. Plot with confidence intervals

########## 0a. Common values ##########

# Maintenance allometric scalar
rho1 = 0.02  

# Maintenance allometric exponent
rho2 = 0.76

# Assimilation efficiency
alpha = 0.4 

# Feeding level
Y = 1 

# Season length in days, i.e. number of growth time steps in a season
sl = 183

# Maximum body size 
mmax = 20000

# Reference temperature
T0 = 292         

# Boltzmann constant   
k = 8.617333e-05

# Temperature scaling parameter values, Lindmark et al. 2022, Global change biol.
EaC = 0.73 # Aactivation energy Consumption
EaM = 0.62 # Activation energy Maintenance
EdI = 1.89 # Deactivation energy

# temperature at which half the rate is reduced due to temperature
Td = T0 + 0.75 

# Rate at reference (a common) temperature Sharpe-Schoolfield
bTc <- 1.824953


########## 0b. Common functions ##########

# For body size dependent kappa
#	m: mass
#	kappa: baseline metabolism (?)
#	returns metabolism or maintenance
kap_fun = function(m, kappa, ha){
	kappa*exp(-m/(ha*mmax))
}

# Boltzmann-Arrhenius function
# 	T: temperature
#	returns scaling value
rM_T_A = function(T) {  
	exp(EaM*(T-T0)/(k*T*T0))
}

# Sharpe-Schoolfield (as in Padfield) for unimodal consumption over temperature
#	T: temperature
#	returns temperature
rC_T_Pad <- function(T) { 
		bTc * exp(EaC*(T-T0)/(k*T*T0)) * (1+exp(EdI*(T-Td)/(k*T*Td)))^(-1) 
}

# For body size dependent assimilation
#	m: mass
#	kappa: baseline metabolism (?)
#	returns assimilation
asim_fun = function(m, eps1, eps2){
	eps1*(m^eps2)
}

# Temperature & mass dependent rate function
#	time: continuous time value
#	m: mass
#	returns all rates
dwdt <- function(time, m, Pars) {
	with(as.list(c(m, Pars)), {
		maintenance = (rho1*(m^rho2)*rM_T_A(T))
		asim_energy = alpha*Y*asim_fun(m,eps1,eps2)*rC_T_Pad(T)
		mass = kap_fun(m,kappa,ha)*asim_energy - maintenance
		return(list(mass, maintenance, asim_energy))
	})
}

# Growth solver function
#	m: mass
#	Pars: parameters for ODE
#	returns growth trajectories
ratefun = function(m, Pars) { # mean function for m_s+1 over one time step (one season). 
	yini = m   # initial biomass in integration
	time = seq(1, sl, by = 1) # time interval for an sl day long growing season, assumed for pike
	a = sapply(yini, ode, time, func = dwdt, parms = Pars)
	a = array(as.numeric(unlist(a,F,F)), dim=c(sl, 4, length(yini)), 
		dimnames = list(c(1:(sl)), c("day","mass","maintenance","asim_energy"), m))
	return(a) 
}

# Modeled growth and fecundity
# 	parEst: parameter estimates
#	Temp: temperature
#	t: number of time steps
#	gdat: growth data
#	fdat: fecundity data
#	returns projected growth and fecundity
project = function(parEst, Temp, t, gdat, fdat) {

	# Assign parameters
	kappa = parEst[1]
	ha = parEst[2]
	eps1 = parEst[3]
	eps2 = parEst[4]

	# Initial growth and fecundity
	g = mean(gdat[gdat$Age==1, ]$mass)
	f = 0

	# Project growth and fecundity
	for(i in 2:t){
		ph = ratefun(g[i-1], Pars = c(T=Temp,kappa=kappa,ha=ha,eps1=eps1,eps2=eps2,Y=1))
		g[i] = ph[sl,2,]
		f[i] = sum(ph[,"asim_energy",]*(1-kap_fun(ph[,"mass",],kappa,ha)))/e_m
	}
	
	# Convert output to data.frames
	g = as.data.frame(cbind(Age=1:length(g),Mass=g))
	f = as.data.frame(cbind(mass=g$Mass, fecund=f))
	
	return(list(g = g, f = f))
}


########## I. Load data ##########

# Read in growth data
growth = read.csv("PikeGrowthData1944_1995.csv")

# Read in fecundity data
fecundity = read.csv("Windermere_Pike_Fecundity_and_Egg_Data_1963_to_2003.csv")

# Read in temepartrue data
temperature = read.csv("WindermereMonthlyTemp1946_2012.csv")

########## II. Process data ##########

# Mean growth season temperature matching growth data
temperature %>%
  #filter(Year >= 1963, Year <= 1995, Month %in% c("Apr","May","Jun","Jul","Aug","Sep")) %>%
  filter(Year <= 2003, Month %in% c("Apr","May","Jun","Jul","Aug","Sep")) %>%
  .$Temp %>%
  mean()

# Annual growth season temperature
mean_temp =
  temperature %>%
  filter(Month %in% c("Apr","May","Jun","Jul","Aug","Sep")) %>%
  #filter(Year >= 1963, Year <= 1995, Month %in% c("Apr","May","Jun","Jul","Aug","Sep")) %>%
  aggregate(Temp ~Year, data=., mean)
hist(mean_temp$Temp, 15)

# Merge mean growth season temperature and growth data
growthTemp = merge(growth,mean_temp,by="Year") # we lose 1944 & 1945 as they are not in temp data
# Calculate mean individual temperature experienced
growthTemp = merge(growthTemp, aggregate(Temp ~ Ind, data = growthTemp, mean), by="Ind", suffixes = c(".mean",".Indmean"))

# Select growth data +/- 1 degree from mean temperature or belonging to the second and third quartile
qT = quantile(growthTemp$Temp.Indmean, probs = c(.25, .5, .75))
growthTemp_select = 
 growthTemp %>%
  #filter(Temp.Indmean > 14.37871-1 & Temp.Indmean < 14.37871+1) 
  filter(Temp.Indmean > qT[1] & Temp.Indmean < qT[3]) 

# Merge mean growth season temperature and fecundity data
fecTemp = merge(fecundity,mean_temp,by="Year") # we lose 1944 & 1945 as they are not in temp data

# Select fecundity data +/- 1 degree from mean temperature or belonging to the second and third quartile
fecTemp_select = 
  fecTemp %>%
  filter(Temp > qT[1] & Temp < qT[3]) 
  #filter(Temp > 14.37871-1 & Temp < 14.37871+1) 

# Length-weight parameters (w = a*L^b)
a = exp(-6.49286)
b = 3.4434

# Convert length to weight
growthTemp_select$mass = a * growthTemp_select$Length ^ b
#growth$mass = a * growth$Length ^ b
fecTemp_select$mass = a * fecTemp_select$Length ^ b
#fecundity$mass = a * fecundity$Length ^ b

# Get age for each individual
birthyear = unlist(lapply(split(growthTemp_select, growthTemp_select$Ind), function(x) {min(x$Year)}))
growthTemp_select$Age = growthTemp_select$Year - birthyear[as.character(growthTemp_select$Ind)] + 1

# Plot data
par(mfrow = c(1,2))
plot(mass ~ Age, growthTemp_select, main = "Growth")
e_m = mean(fecTemp_select$Egg_weight) # mean egg mass
plot(Egg_number*Egg_weight/e_m ~ mass, fecTemp_select, main = "Fecundity")

########## III. Set parameters ##########

# Consumption allometric scalar starting value
eps1_start = 1

# Consumption allometric exponent starting value
eps2_start = 0.54 

# Baseline metabolism starting value
kappa_start = 0.8

# Metabolic scaling steepness starting value
ha_start = 2

# Choose temperature
Temp = 287

# Choose number of years
t = 20


########## IV. Optimization function ##########

# Function to optimize (find the minimum)
# 	par: a vector of four values, namely kappa, ha, eps1, and eps2 in that order
#	gw: the weighting for growth error fitting, a value in [0,1]; fecundity will be 1-gw
#	Temp: temperature in Kelvin
#	t: number of time steps
#	gdat: growth data
#	fdat: fecundity data
#	returns weighted sum of log error in growth and fecundity
fit_deb = function(par, gw, Temp, t, gdat, fdat) {
	
	# Set model parameters
	kappa = par[1]
	ha = par[2]
	eps1 = par[3]
	eps2 = par[4]
	
	# Set error weight
	fw = 1 - gw
	
	# Initial body size
	g = mean(gdat[gdat$Age==1, ]$mass)
	
	# Fecundity at age 1 is zero
	f = 0
	
	# Get modeled growth and fecundity
	for(i in 2:t){
		ph = ratefun(g[i-1], Pars = c(T=Temp,kappa=kappa,ha=ha,eps1=eps1,eps2=eps2,Y=1))
		g[i] = ph[sl,2,]
		f[i] = sum(ph[,"asim_energy",]*(1-kap_fun(ph[,"mass",],kappa,ha)))/e_m
	}
	g = as.data.frame(cbind(Age=1:length(g),Mass=g))
	f = as.data.frame(cbind(mass=g$Mass, fecund=f))
	
	# Get observed growth
	gmod = g[gdat$Age, ]$Mass
	
	# Growth error
	g_error = sum((gdat$mass - gmod)^2)
	
	# Get observed fecundity
	func = splinefun(x = f$mass, y = f$fecund, method = "fmm",  ties = mean)
	fmod = func(fdat$mass) #?
	
	# Fecundity error
	f_error = sum((fdat$Egg_number*fdat$Egg_weight/e_m - fmod)^2)
	
	# Final weighted error
	error = log(g_error)*gw + log(f_error)*fw
	
	return(error)
}


########## V. Optimize ##########

# Arrange starting values
par = c(kappa_start, ha_start, eps1_start, eps2_start)

# Set up parallel cores
noCores = detectCores() - 1
cl = makeCluster(noCores, setup_timeout = 0.5)
setDefaultCluster(cl = cl)
clusterExport(cl, as.list(ls()))
clusterEvalQ(cl, {
	library(deSolve)
})
deb.optim = optimParallel(par = par, fit_deb, gw = 0.5, Temp = Temp, t = t, gdat = growthTemp_select, fdat = fecTemp_select, lower = c(0.01,0.01,0.01,0), upper = c(10,100,10,2), method = "L-BFGS-B", parallel = list(cl=cl,loginfo=T), hessian = T)
stopCluster(cl)

# Get modeled growth and fecundity
out.optim = project(deb.optim$par, Temp, t, growthTemp_select, fecTemp_select)

save(deb.optim, file = "optim_Pars.RData")

########## VI. Bootstrap confidence intervals ##########

# Choose number of iterations
n = 10

# Output data frame
boot.deb = matrix(nrow = n, ncol = 4)
colnames(boot.deb) = c("kappa", "ha", "eps1", "eps2")

# Set seed
set.seed(5573)

# Resampling loop
for(i in 1:n) {

	# Sample growth and fecundity data
	gdat = growthTemp_select[sample(seq(nrow(growthTemp_select)), replace = T), ]
	fdat = fecTemp_select[sample(seq(nrow(fecTemp_select)), replace = T), ]
	
	# Run optimization procedure
	noCores = detectCores() - 1
	cl = makeCluster(noCores, setup_timeout = 0.5)
	setDefaultCluster(cl = cl)
	clusterExport(cl, as.list(ls()))
	clusterEvalQ(cl, {
		library(deSolve)
	})
	opt = optimParallel(par = par, fit_deb, gw = 0.5, Temp = Temp, t = t, gdat = gdat, fdat = fdat, lower = c(0.01,0.01,0.01,0), upper = c(10,100,10,2), method = "L-BFGS-B", parallel = list(cl=cl,loginfo=T))
	stopCluster(cl)
	
	# Save parameter estimates
	boot.deb[i,] = opt$par
}

# Predict and plot
preds = apply(unname(boot.deb), 1, project, Temp = Temp, t = t, gdat = growthTemp_select, fdat = fecTemp_select)

# Combine predictions
gboot = do.call(rbind, lapply(preds, function(x) x$g))
fboot = do.call(rbind, lapply(preds, function(x) x$f))

# Get 2.5% and 97.5% quantiles for growth
gci_list = lapply(split(gboot, gboot$Age), function(x) quantile(x$Mass, probs = c(0.025,0.975)))
gci = do.call(rbind, gci_list)

#  Get 2.5% and 97.5% quantiles for fecundity over ranges of mass from original optim
flist = split(fboot, cut(fboot$mass, c(0,out.optim$f$mass), include.lowest = FALSE))
fci_list = lapply(flist, function(x) quantile(x$fecund, probs = c(0.025,0.975)))
fci = do.call(rbind, fci_list)

########## VII. Plot with confidence intervals ##########

# Initialize plot
par(mfrow = c(1,2))

# Plot growth
plot(Mass ~ Age, out.optim$g, type = "l", col = "red", main = "Growth")
points(growthTemp_select$Age, growthTemp_select$mass)
polygon(x = c(out.optim$g$Age,rev(out.optim$g$Age)), y = c(gci[,1],rev(gci[,2])), col = rgb(1,0,0,0.2), border = NA)

# Plot fecundity
plot(fecund ~ mass, out.optim$f, type = "l", col = "red", main = "Fecundity")
points(fecTemp_select$mass, fecTemp_select$Egg_number*fecTemp_select$Egg_weight/e_m)
polygon(x = c(out.optim$f$mass,rev(out.optim$f$mass)), y = c(fci[,1],rev(fci[,2])), col = rgb(1,0,0,0.2), border = NA)

par(mfrow = c(1,1))
