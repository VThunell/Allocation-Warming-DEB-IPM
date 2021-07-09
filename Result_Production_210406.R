# Results production

# 1 - Baseline (Robust, Vindenes_surv) temp 282-294
# 2 - MortSizedep_Highexp (, Lorenzon_surv)
# 3 - MortSizedep_Lowexp (, Lorenzon_surv)
# 4 - Tempdepmortality (, Lorenzon_surv)
# 5 - SizeTempInt (cIa=0.01, Vindenes_surv)

# 7 - Allokering i ha och inte kappa
# 6 - SizeTempInt_weak (c:id=-0.01, Vindenes_surv)
# 8 - SizeTempInt_zero (c:id=0.0, Vindenes_surv)


### RESULT 1 ####

DEBparams <- as.data.frame(rbind(rho1,rho2,eps1,eps2,alpha,sl,T0,
                                 k,cIa,cM,EaI,EaM,EaS,Td,cId,EdI))
IPMparams <- as.data.frame(rbind(mmin,mmax,n,dx,e_m,el_surv))
meta <- rbind(DEBparams,IPMparams, kapfun=kap_fun(100))
write.table(meta, file="Res_meta_1_0406.txt",quote=TRUE, sep=",", row.names=TRUE)

T <- seq(284,294,0.25) # temperature range
kappa <- seq(0.8,1,0.025) # kappa range
Y <- 1 # feeding levels

Res_lam <- matrix(ncol = 3+1, nrow = length(kappa)*length(T)*length(Y)) # assuming 40 here from files
Res_v  <-  matrix(ncol = 3+n+1, nrow = length(kappa)*length(T)*length(Y)) # assuming 40 here from files
Res_w  <-  matrix(ncol = 3+n+1, nrow = length(kappa)*length(T)*length(Y)) # assuming 40 here from files

#length(Res_lam[,1])

parsK <- as.matrix(expand.grid(T,kappa,Y))
colnames(parsK) <- c("T","kappa","Y")

for (i in 1:nrow(parsK)){
  res <- wvlambda.projection(K.matrix(parsK[i,]))
  Res_lam[i,] <- c(parsK[i,],res$lam) #c(T[d],kappa[e],Y[f], res$lambda)
  Res_v[i,]   <- c(parsK[i,],res$v) #REPRODUCTIVE VALUES
  Res_w[i,]   <- c(parsK[i,],res$w) #STABLE STRUCTURE
}

colnames(Res_lam) <- c("T","kappa","Y","Lambda")
colnames(Res_v) <- c("T","kappa","Y","0",x)
colnames(Res_w) <- c("T","kappa","Y","0",x)

write.table(Res_lam, file="Res_lam_1_0406.txt",quote=TRUE, sep=",", row.names=TRUE)
write.table(Res_v, file="Res_v_1_0406.txt",quote=TRUE, sep=",", row.names=TRUE)
write.table(Res_w, file="Res_w_1_0406.txt",quote=TRUE, sep=",", row.names=TRUE)


### RESULT 2 ####

DEBsurvfun <- function(m, Pars) { # Mass-Temp dependence of yearly Mortality 
  sx <- function(m, Pars){ # natural baseline survival
    exp(-3*m^-0.315)#*exp(EaS*(Pars[['T']]-T0)/(k*Pars[['T']]*T0))) #3*^-.288 is yearly mortality rate for 1 gram unit weight from Lorenzen 1996
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
  ifelse(m == e_m, sx.firstyear(m,Pars), sx(m,Pars)*starvx(m, Pars))# multiply with F (fishing mortality function)
  #ifelse(m == e_m, sx.firstyear(m,Pars), Vsurvfun(m)*starvx(m, Pars))# multiply with F (fishing mortality function)
} 

DEBparams <- as.data.frame(rbind(rho1,rho2,eps1,eps2,alpha,sl,T0,
                                 k,cIa,cM,EaI,EaM,EaS,Td,cId,EdI))
IPMparams <- as.data.frame(rbind(mmin,mmax,n,dx,e_m,el_surv))
meta <- rbind(DEBparams,IPMparams, kapfun=kap_fun(100))
write.table(meta, file="Res_meta_2_0406.txt",quote=TRUE, sep=",", row.names=TRUE)

T <- seq(284,294,0.25) # temperature range
kappa <- seq(0.8,1,0.025) # kappa range
Y <- 1 # feeding levels

Res_lam <- matrix(ncol = 3+1, nrow = length(kappa)*length(T)*length(Y)) # assuming 40 here from files
Res_v  <-  matrix(ncol = 3+n+1, nrow = length(kappa)*length(T)*length(Y)) # assuming 40 here from files
Res_w  <-  matrix(ncol = 3+n+1, nrow = length(kappa)*length(T)*length(Y)) # assuming 40 here from files

parsK <- as.matrix(expand.grid(T,kappa,Y))
colnames(parsK) <- c("T","kappa","Y")

for (i in 1:nrow(parsK)){
  res <- wvlambda.projection(K.matrix(parsK[i,]))
  Res_lam[i,] <- c(parsK[i,],res$lam) #c(T[d],kappa[e],Y[f], res$lambda)
  Res_v[i,]   <- c(parsK[i,],res$v) #REPRODUCTIVE VALUES
  Res_w[i,]   <- c(parsK[i,],res$w) #STABLE STRUCTURE
}

colnames(Res_lam) <- c("T","kappa","Y","Lambda")
colnames(Res_v) <- c("T","kappa","Y","0",x)
colnames(Res_w) <- c("T","kappa","Y","0",x)

write.table(Res_lam, file="Res_lam_2_0406.txt",quote=TRUE, sep=",", row.names=TRUE)
write.table(Res_v, file="Res_v_2_0406.txt",quote=TRUE, sep=",", row.names=TRUE)
write.table(Res_w, file="Res_w_2_0406.txt",quote=TRUE, sep=",", row.names=TRUE)


### RESULT 3 ####

DEBsurvfun <- function(m, Pars) { # Mass-Temp dependence of yearly Mortality 
  sx <- function(m, Pars){ # natural baseline survival
    exp(-3*m^-0.261)#*exp(EaS*(Pars[['T']]-T0)/(k*Pars[['T']]*T0))) #3*^-.288 is yearly mortality rate for 1 gram unit weight from Lorenzen 1996
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
  ifelse(m == e_m, sx.firstyear(m,Pars), sx(m,Pars)*starvx(m, Pars))# multiply with F (fishing mortality function)
  #ifelse(m == e_m, sx.firstyear(m,Pars), Vsurvfun(m)*starvx(m, Pars))# multiply with F (fishing mortality function)
} 

DEBparams <- as.data.frame(rbind(rho1,rho2,eps1,eps2,alpha,sl,T0,
                                 k,cIa,cM,EaI,EaM,EaS,Td,cId,EdI))
IPMparams <- as.data.frame(rbind(mmin,mmax,n,dx,e_m,el_surv))
meta <- rbind(DEBparams,IPMparams, kapfun=kap_fun(100))
write.table(meta, file="Res_meta_3_0406.txt",quote=TRUE, sep=",", row.names=TRUE)

T <- seq(284,294,0.25) # temperature range
kappa <- seq(0.8,1,0.025) # kappa range
Y <- 1 # feeding levels

Res_lam <- matrix(ncol = 3+1, nrow = length(kappa)*length(T)*length(Y)) # assuming 40 here from files
Res_v  <-  matrix(ncol = 3+n+1, nrow = length(kappa)*length(T)*length(Y)) # assuming 40 here from files
Res_w  <-  matrix(ncol = 3+n+1, nrow = length(kappa)*length(T)*length(Y)) # assuming 40 here from files

#length(Res_lam[,1])

parsK <- as.matrix(expand.grid(T,kappa,Y))
colnames(parsK) <- c("T","kappa","Y")

for (i in 1:nrow(parsK)){
  res <- wvlambda.projection(K.matrix(parsK[i,]))
  Res_lam[i,] <- c(parsK[i,],res$lam) #c(T[d],kappa[e],Y[f], res$lambda)
  Res_v[i,]   <- c(parsK[i,],res$v) #REPRODUCTIVE VALUES
  Res_w[i,]   <- c(parsK[i,],res$w) #STABLE STRUCTURE
}

colnames(Res_lam) <- c("T","kappa","Y","Lambda")
colnames(Res_v) <- c("T","kappa","Y","0",x)
colnames(Res_w) <- c("T","kappa","Y","0",x)

write.table(Res_lam, file="Res_lam_3_0406.txt",quote=TRUE, sep=",", row.names=TRUE)
write.table(Res_v, file="Res_v_3_0406.txt",quote=TRUE, sep=",", row.names=TRUE)
write.table(Res_w, file="Res_w_3_0406.txt",quote=TRUE, sep=",", row.names=TRUE)


### RESULT 4 ####

DEBsurvfun <- function(m, Pars) { # Mass-Temp dependence of yearly Mortality 
  sx <- function(m, Pars){ # natural baseline survival
    exp(-3*m^-0.261)*exp(EaS*(Pars[['T']]-T0)/(k*Pars[['T']]*T0)) #3*^-.288 is yearly mortality rate for 1 gram unit weight from Lorenzen 1996
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
  ifelse(m == e_m, sx.firstyear(m,Pars), sx(m,Pars)*starvx(m, Pars))# multiply with F (fishing mortality function)
  #ifelse(m == e_m, sx.firstyear(m,Pars), Vsurvfun(m)*starvx(m, Pars))# multiply with F (fishing mortality function)
} 

DEBparams <- as.data.frame(rbind(rho1,rho2,eps1,eps2,alpha,sl,T0,
                                 k,cIa,cM,EaI,EaM,EaS,Td,cId,EdI))
IPMparams <- as.data.frame(rbind(mmin,mmax,n,dx,e_m,el_surv))
meta <- rbind(DEBparams,IPMparams, kapfun=kap_fun(100))
write.table(meta, file="Res_meta_4_0406.txt",quote=TRUE, sep=",", row.names=TRUE)

T <- seq(284,294,0.25) # temperature range
kappa <- seq(0.8,1,0.025) # kappa range
Y <- 1 # feeding levels

Res_lam <- matrix(ncol = 3+1, nrow = length(kappa)*length(T)*length(Y)) # assuming 40 here from files
Res_v  <-  matrix(ncol = 3+n+1, nrow = length(kappa)*length(T)*length(Y)) # assuming 40 here from files
Res_w  <-  matrix(ncol = 3+n+1, nrow = length(kappa)*length(T)*length(Y)) # assuming 40 here from files

parsK <- as.matrix(expand.grid(T,kappa,Y))
colnames(parsK) <- c("T","kappa","Y")

for (i in 1:nrow(parsK)){
  res <- wvlambda.projection(K.matrix(parsK[i,]))
  Res_lam[i,] <- c(parsK[i,],res$lam) #c(T[d],kappa[e],Y[f], res$lambda)
  Res_v[i,]   <- c(parsK[i,],res$v) #REPRODUCTIVE VALUES
  Res_w[i,]   <- c(parsK[i,],res$w) #STABLE STRUCTURE
}

colnames(Res_lam) <- c("T","kappa","Y","Lambda")
colnames(Res_v) <- c("T","kappa","Y","0",x)
colnames(Res_w) <- c("T","kappa","Y","0",x)

write.table(Res_lam, file="Res_lam_4_0406.txt",quote=TRUE, sep=",", row.names=TRUE)
write.table(Res_v, file="Res_v_4_0406.txt",quote=TRUE, sep=",", row.names=TRUE)
write.table(Res_w, file="Res_w_4_0406.txt",quote=TRUE, sep=",", row.names=TRUE)

### RESULT 5 ####

cIa <- 0.008

DEBsurvfun <- function(m, Pars) { # Mass-Temp dependence of yearly Mortality 
  sx <- function(m, Pars){ # natural baseline survival
    exp(-3*m^-0.261)*exp(EaS*(Pars[['T']]-T0)/(k*Pars[['T']]*T0)) #3*^-.288 is yearly mortality rate for 1 gram unit weight from Lorenzen 1996
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
  ifelse(m == e_m, sx.firstyear(m,Pars), Vsurvfun(m)*starvx(m, Pars))# multiply with F (fishing mortality function)
} 

DEBparams <- as.data.frame(rbind(rho1,rho2,eps1,eps2,alpha,sl,T0,
                                 k,cIa,cM,EaI,EaM,EaS,Td,cId,EdI))
IPMparams <- as.data.frame(rbind(mmin,mmax,n,dx,e_m,el_surv))
meta <- rbind(DEBparams,IPMparams, kapfun=kap_fun(100))
write.table(meta, file="Res_meta_5_0406.txt",quote=TRUE, sep=",", row.names=TRUE)

T <- seq(284,294,0.25) # temperature range
kappa <- seq(0.8,1,0.025) # kappa range
Y <- 1 # feeding levels

Res_lam <- matrix(ncol = 3+1, nrow = length(kappa)*length(T)*length(Y)) # assuming 40 here from files
Res_v  <-  matrix(ncol = 3+n+1, nrow = length(kappa)*length(T)*length(Y)) # assuming 40 here from files
Res_w  <-  matrix(ncol = 3+n+1, nrow = length(kappa)*length(T)*length(Y)) # assuming 40 here from files

parsK <- as.matrix(expand.grid(T,kappa,Y))
colnames(parsK) <- c("T","kappa","Y")

for (i in 1:nrow(parsK)){
  res <- wvlambda.projection(K.matrix(parsK[i,]))
  Res_lam[i,] <- c(parsK[i,],res$lam) #c(T[d],kappa[e],Y[f], res$lambda)
  Res_v[i,]   <- c(parsK[i,],res$v) #REPRODUCTIVE VALUES
  Res_w[i,]   <- c(parsK[i,],res$w) #STABLE STRUCTURE
}

colnames(Res_lam) <- c("T","kappa","Y","Lambda")
colnames(Res_v) <- c("T","kappa","Y","0",x)
colnames(Res_w) <- c("T","kappa","Y","0",x)

write.table(Res_lam, file="Res_lam_5_0406.txt",quote=TRUE, sep=",", row.names=TRUE)
write.table(Res_v, file="Res_v_5_0406.txt",quote=TRUE, sep=",", row.names=TRUE)
write.table(Res_w, file="Res_w_5_0406.txt",quote=TRUE, sep=",", row.names=TRUE)
