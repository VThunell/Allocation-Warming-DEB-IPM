
#### RUN DEBIPM TO PRODUCE RESLUTS ####

# Run main (1) and the two (2,3) contrasting survival rate models
# 1 - Baseline (Vindenes_Tsurv)
# 2 - Temp independent surv (Vsurvfun(m))
# 3 - Size independent surv (0.68)

# RUNNING THE MODEL FOR ONE SET OF PARAMETERS TAKES 20 SECONDS ON A GOOD COMPUTER, 
# RUNNING length(T) * length(kappa) NUMBER OF MODELS THUS TAKES TIME.

date <- Sys.Date()

T <- seq(283.5,294.5,0.25) # temperature range
kappa <- seq(0.59,1,0.01) # kappa range
Y <- 1 # feeding levels

### RESULT 1 ####
Res_lam <- matrix(ncol = 3+1, nrow = length(kappa)*length(T)*length(Y)) # assuming 40 here from files
Res_v  <-  matrix(ncol = 3+n+1, nrow = length(kappa)*length(T)*length(Y)) # assuming 40 here from files
Res_w  <-  matrix(ncol = 3+n+1, nrow = length(kappa)*length(T)*length(Y)) # assuming 40 here from files

parsK <- as.matrix(expand.grid(T,kappa,Y)) #1890 runs 
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

write.table(Res_lam, file=paste("Res_lam_1_",date,".txt", sep = ""),quote=TRUE,  sep=",", row.names=TRUE)
write.table(Res_v, file=paste("Res_v_1_",date,".txt", sep = ""),quote=TRUE, sep=",", row.names=TRUE)
write.table(Res_w, file=paste("Res_w_1_",date,".txt", sep = ""),quote=TRUE, sep=",", row.names=TRUE)

### RESULT 2 ####

### SURVIVAL a(x) ####
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
  #ifelse(m < mmin, sx.firstyear(m,Pars), VsurvfunT(m,Pars)) # main, i.e. size and temp depedent, survival
  ifelse(m < mmin, sx.firstyear(m,Pars), Vsurvfun(m)) # temperature independent survival
  #ifelse(m < mmin, sx.firstyear(m,Pars), 0.68) # consant, i.e. size and temp INdepedent, survival
} 

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

write.table(Res_lam, file=paste("Res_lam_2_",date,".txt", sep = ""),quote=TRUE,  sep=",", row.names=TRUE)
write.table(Res_v, file=paste("Res_v_2_",date,".txt", sep = ""),quote=TRUE, sep=",", row.names=TRUE)
write.table(Res_w, file=paste("Res_w_2_",date,".txt", sep = ""),quote=TRUE, sep=",", row.names=TRUE)


### RESULT 3 ####

### SURVIVAL a = 0.68 ####
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
  #ifelse(m < mmin, sx.firstyear(m,Pars), VsurvfunT(m,Pars)) # main, i.e. size and temp depedent, survival
  #ifelse(m < mmin, sx.firstyear(m,Pars), Vsurvfun(m)) # temperature independent survival
  ifelse(m < mmin, sx.firstyear(m,Pars), 0.68) # consant, i.e. size and temp INdepedent, survival
} 

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

write.table(Res_lam, file=paste("Res_lam_3_",date,".txt", sep = ""),quote=TRUE,  sep=",", row.names=TRUE)
write.table(Res_v, file=paste("Res_v_3_",date,".txt", sep = ""),quote=TRUE, sep=",", row.names=TRUE)
write.table(Res_w, file=paste("Res_w_3_",date,".txt", sep = ""),quote=TRUE, sep=",", row.names=TRUE)
