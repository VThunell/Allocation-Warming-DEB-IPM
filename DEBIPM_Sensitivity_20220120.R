
#### Sensitivity of lambda to kappa - decomposition into contributions 
#  from demographic functions ###

# Sensitivity (perturbation) analysis across different temperatures to study which
# demographic processes, and the size at which they occur, affect how optimal energy 
# allocation changes with warming. We analyze the partial contribution of each of 
# o(T) (size at age 1), f(x,T) (fecundity) and g(y:x,T) (somatic growth) to
# sensitivity of λ to changes in energy allocation (κ_0) at three different 
# temperatures.to understand how size dependent vital rates affects optimum 

library(tidyverse) # gglot etc.

# Temporary - create results for OptPars
Res_v287 <- c(T=OptPars[1,1], kappa=OptPars[1,2], Y=1,
            wvlambda.projection(K.matrix(OptPars[1,]))$v)
Res_w287 <- c(T=OptPars[1,1], kappa=OptPars[1,2], Y=1,
              wvlambda.projection(K.matrix(OptPars[1,]))$w)

Res_v289 <- c(T=OptPars[2,1], kappa=OptPars[2,2], Y=1,
              wvlambda.projection(K.matrix(OptPars[2,]))$v)
Res_w289 <- c(T=OptPars[2,1], kappa=OptPars[2,2], Y=1,
              wvlambda.projection(K.matrix(OptPars[2,]))$w)

Res_v291 <- c(T=OptPars[3,1], kappa=OptPars[3,2], Y=1,
              wvlambda.projection(K.matrix(OptPars[3,]))$v)
Res_w291 <- c(T=OptPars[3,1], kappa=OptPars[3,2], Y=1,
              wvlambda.projection(K.matrix(OptPars[3,]))$w)

Res_v <- as.data.frame(rbind(Res_v287, Res_v289, Res_v291))
colnames(Res_v)[4:ncol(Res_v)] <- c(0,x)
Res_w <- as.data.frame(rbind(Res_w287, Res_w289, Res_w291))
colnames(Res_w)[4:ncol(Res_w)] <- c(0,x)

# # Model output on lambda, w and v:
# Res_lam <- read.delim("//storage-og.slu.se/home$/vitl0001/My Documents/Manus2/Results/Results1108/Res_lam_1_1108.txt", sep = ",")
# Res_v <- read.delim("//storage-og.slu.se/home$/vitl0001/My Documents/Manus2/Results/Results1108/Res_v_1_1108.txt", sep = ",")
# colnames(Res_v)[4:ncol(Res_v)] <- sub("X", "", colnames(Res_v)[4:ncol(Res_v)])
# Res_w <- read.delim("//storage-og.slu.se/home$/vitl0001/My Documents/Manus2/Results/Results1108/Res_w_1_1108.txt", sep = ",")
# colnames(Res_w)[4:ncol(Res_w)] <- sub("X", "", colnames(Res_w)[4:ncol(Res_w)])

maxl_1 <- as.data.frame(Res_lam) %>% 
  group_by(T) %>%
  slice_max(Lambda)
maxl_1["Sur_f"] <- "Vin"

# optimum kappa dependent on temperature:
OptPars <- 
  as.data.frame(maxl_1[,1:3]) %>%
  filter(T %in% c(287, 289, 291))

#### Functions for numerical differentiation (difference forward differentiation ####
# i.e. ( f(x + h) - f(x) ) / h) ) of dGrow/dk_0, dFec/dk_0 and dAge1/k_0. 
# Kernel level pert.
K_NumDer <- function(h, Pars){
  (K.matrix(Pars = c(T = Pars[["T"]],kappa = Pars[["kappa"]] + h,Y = Pars[["Y"]])) -
     K.matrix(Pars = c(T = Pars[["T"]],kappa = Pars[["kappa"]],Y = Pars[["Y"]])))/h
}
# Growth function pert.
G_NumDer <- function(h ,m, Pars){
  (DEBgrowthfun(m, Pars = c(T = Pars[["T"]],kappa = Pars[["kappa"]] + h,Y = Pars[["Y"]])) -
     DEBgrowthfun(m, Pars = c(T = Pars[["T"]],kappa = Pars[["kappa"]],Y = Pars[["Y"]])))/h
} 
# Fecundity function pert.
F_NumDer <- function(h ,m, Pars){
  (DEBrepfun(m, Pars = c(T = Pars[["T"]],kappa = Pars[["kappa"]] + h,Y = Pars[["Y"]])) - 
   DEBrepfun(m, Pars = c(T = Pars[["T"]],kappa = Pars[["kappa"]],Y = Pars[["Y"]])))/h
}
# Age 1 size function pert.
A1_NumDer <- function(h ,m, Pars){
  (DEBage1.size(Pars = c(T = Pars[["T"]],kappa = Pars[["kappa"]] + h,Y = Pars[["Y"]])) - 
     DEBage1.size(Pars = c(T = Pars[["T"]],kappa = Pars[["kappa"]],Y = Pars[["Y"]])))/h
}

#### 1 - Calculate dlambda/dkappa_0 at 287 K ####
# "en perturbasjon av kappa_0 direkte og måle sensitiviteten i lambda 
# (da ser vi hva bidragene burde summeres til). "
(wvlambda.projection(K.matrix(Pars = c(T = OptPars[1,1],kappa = OptPars[1,2] 
                                       + 0.00001, Y = OptPars[1,3])))$lam
 - 
 # wvlambda.projection(K.matrix(Pars = c(T = OptPars[1,1],kappa = OptPars[1,2], 
 #                                        Y = OptPars[1,3])))$lam) / 0.00001
 #   
    Res_lam[c(Res_lam$T == OptPars[1,1] & Res_lam$kappa == OptPars[1,2]
           & Res_lam$Y == OptPars[1,3]),4]) / 0.00001
# gives -0.04407772


#### 2 - Calculate total size dependent sensitivity at 287 K by perturbing K ####
# and calculate: sum_i(sum_j(v*w*dK/dkapa_0)) at kappa_0* at 287 K
#res_v_t <- as.numeric(Res_v %>% filter(T == OptPars[1,1] & kappa==OptPars[1,2]))[4:ncol(Res_v)]
#res_w_t <- as.numeric(Res_w %>% filter(T == OptPars[1,1] & kappa==OptPars[1,2]))[4:ncol(Res_w)]

Sens_Kmat <- outer(res_v_t, res_w_t, "*") * K_NumDer(0.00001, OptPars[1,])
sum(apply(Sens_Kmat,1,sum)) # very close to equal to dlambda/dkappa_0
plot(c(0,x),apply(Sens_Kmat,1,sum), type = "l")
plot(x,apply(Sens_Kmat[2:(n+1),],1,sum), type = "l") # eggs excluded and this looks the same as in plot of 3.

#### 3 - Contributions to size dependent sensitivity from demographic functions ####
Res_Sens <- NULL
for (i in 1:nrow(OptPars)){
  res_v <- as.numeric(Res_v %>% filter(T == OptPars[i,1] & kappa==OptPars[i,2]))[4:ncol(Res_v)]
  res_w <- as.numeric(Res_w %>% filter(T == OptPars[i,1] & kappa==OptPars[i,2]))[4:ncol(Res_w)]

  Gr_cont <- bind_cols(Size=x, Temp=OptPars[i,1], kappa=OptPars[i,2], 
                       dfun="Gr", outer(res_v[2:(n+1)], res_w[2:(n+1)], "*") 
                       * t(DEBsurvfun(x, OptPars[i,]) * t(G_NumDer(h=0.00001, x, OptPars[i,])))*dx)
  Gr_cont["Sens"] <- apply(Gr_cont[5:ncol(Gr_cont)],1,FUN=sum)# sum the distribution y to get sensitivity for each size 
  Gr_cont<-Gr_cont[,c(1:4,ncol(Gr_cont))]
  
  F_cont  <- bind_cols(Size=x, Temp=OptPars[i,1], kappa=OptPars[i,2],
                       dfun="F",  Sens=res_v[1]*res_w[2:(n+1)] 
                       * DEBsurvfun(x, OptPars[i,]) * F_NumDer(h=0.00001, x, OptPars[i,]))
  
  A1_cont <- bind_cols(Size=x, Temp=OptPars[i,1], kappa=OptPars[i,2], 
                       dfun="A1", Sens=res_v[2:(n+1)]*res_w[1] 
                       * el_surv * A1_NumDer(h=0.00001, x, OptPars[i,])*dx)
  
  Res_Sens <- rbind(rbind(Gr_cont,F_cont,A1_cont), Res_Sens) 
}

# sum the contributions from demographic functions at each size
Sum_cont <- 
  Res_Sens %>%           
  group_by(Size, Temp) %>%
  summarise(sum.t = sum(Sens, na.rm = TRUE))

# compare with kernel perturbation
plot(Sum_cont[Sum_cont$Temp==287,]$Size,Sum_cont[Sum_cont$Temp==287,]$sum.t, type = "l", col="red")
lines(x,apply(sens_sum[2:(n+1),],1,sum), type = "l") # eggs excluded

sum(Sum_cont[Sum_cont$Temp==287,]$sum.t) # Same as in both in 1 and 2 ~ -0.04404539

#### 4 - Plot sensitivity contributions (figure  main text)
s_sum <- Res_Sens %>%
  group_by(dfun, Temp) %>%
  summarise(sum.c = sum(Sens, na.rm = TRUE)) #%>%

s_sum_text <- data.frame(Size = c(1500,1500,1500,3500,3500,3500,5500,5500,5500), Sens = 0.075,
                         dfun = s_sum$dfun, Temp = s_sum$Temp, label = round(s_sum$sum.c,4))

cont.labs <- c("Age 1 size", "Fecundity", "Growth")
names(cont.labs) <- c("A1", "F", "Gr") 

Fig_sens <-
Res_Sens %>%
  mutate(dfun=as.factor(dfun)) %>%
  ggplot(.) +
  geom_path(aes(Size, Sens, color = dfun)) +
  facet_grid(Temp~., scales="fixed") +
  xlim(0,12500) +
  geom_vline(xintercept=400, linetype="dashed") +
  ylab("Sensitivity") +
  xlab("Weight [g]") +
  scale_color_grey(labels= cont.labs, name= "",) +
  geom_text(data=s_sum_text, aes(Size, Sens, label = label, 
                              colour=dfun), size=4, show.legend = FALSE) +
  theme_bw() +
  theme(panel.background = element_rect(fill = "white", colour = "black"),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA))
Fig_sens

  
# Gr_AnDer <- function(m, Pars){
#   Gr_s <- unname(ratefun(m, Pars)[sl,2,]) #mass
#   Gr_s <- exp(-Gr_s/(max(x)*2))*alpha*(eps1*(Gr_s^eps2)*rI_T_Pad(Pars[["T"]], Gr_s)) #2 is slope k_m
#   return(Gr_s)
# }
# 
# F_AnDer <- function(m, Pars){
#   f_s <- unname(ratefun(m, Pars)[183,2,]) #mass
#   f_s <- (-0.5*e_m)*exp(-f_s/(max(x)*2))*alpha*(eps1*(f_s^eps2)*rI_T_Pad(Pars[["T"]], f_s)) #2 is slope k_m
#   return(f_s)
# }
# 
# 
# A1_AnDer <- function(m, Pars){
#   a1_s <- unname(ratefun(m, Pars)[round(sl/2),2,]) #mass
#   a1_s <- exp(-a1_s/(max(x)*2))*alpha*(eps1*(a1_s^eps2)*rI_T_Pad(Pars[["T"]], a1_s)) #2 is slope k_m
#   return(a1_s)
# }
