### PLOT DEB IPM 2021 04 26 ####  

#setwd("C:/Users/vitl0001/VThunell_repos/Temperature-DEBIPM")
#source("Thunell_DEBIPM_model_20210105_HAkappa.R")

library(scales)
library(patchwork)
library(gridExtra)
library(RColorBrewer)

#Plot-parameter s for IPM vital rate functions needed for some plots
GR_pars <- c(T = 287,     # parameters for Temperature, feeding, allocation and Mass dependence
             kappa = 0.83, # allocation to respiration (Growth and maintenance)
             Y = 1)         # Feeding level

## PLOT MASS AND TEMPERATURE FUNCTIONS ####

##### PLOT ENERGY BUDGET OVER TEMP AND MASS ####
Temp <- 283:294
T_rates = NULL

for(i in Temp) {
  maint <- rho1*(x^rho2)*rM_T_AL2(i, x)
  intake  <- eps1*(x^eps2)*rI_T_GU2(i, x) 
  T_rates  <- rbind(T_rates,
                    cbind(Temp = i, mass = x, 
                          maint, intake, 
                          growth_E = kap_fun(x,GR_pars["kappa"])*alpha*GR_pars["Y"]*intake - maint, 
                          repro_E = alpha*GR_pars["Y"]*(1-kap_fun(x,GR_pars["kappa"]))*intake))
}

T_rates_long <- gather(as.data.frame(T_rates), key = "E_type", value ="g_day", c(5,6)) # gather only column 6,7
## Over mass for three temperatures
T_rates_long %>%
  filter( Temp %in% c(285,287,289)) %>%
  ggplot(.) +
  geom_line(size=1 , aes(mass, g_day, color = as.factor(Temp), 
                         linetype = as.factor(E_type))) +
  #ylim(0,10) +
  theme_bw() +
  scale_colour_brewer(palette="Dark2")

## Over temperatures for three sizes
T_rates_long %>%
  filter( mass %in% c(x[2], x[20], x[200])) %>%
  ggplot(.) +
  geom_line(size=1 , aes(Temp, g_day, color = as.factor(mass), 
                         linetype = as.factor(E_type))) +
  geom_hline(yintercept = 0) +
  ylim(0,15) +
  theme_bw() +
  scale_colour_brewer(palette="Dark2")

### PLOT GROWTH / SIZE AT AGE DEPENDENT ON MASS, TEMPERATURE & kAPPA ####
#### PLOT RATE FUNCTION FOR IPM ####

### Plot growth fun ####
plot(x,  DEBgrowthfun(14000, y=x, GR_pars)*dx, type="l", ylab="DEBgrowth(y;x)", 
     xlab="y", main = "Prob. density of DEBgrowth(x,y)")
lines(x, DEBgrowthfun(10000, y=x, GR_pars)*dx, lty=2)
lines(x, DEBgrowthfun(1000, y=x, GR_pars)*dx, lty=3)
lines(x, DEBgrowthfun(100, y=x, GR_pars)*dx, lty=4) # values >~8100 are highly unlikley to grow into as y=x??
legend("topleft", c("14000","10000","1000","100"), lty=c(1,2,3,4), cex = 0.7, title="mass [g]")

# Plot egg size distribution in Windermere data ####
# agg <- seq(min(FData$Egg.weight),max(FData$Egg.weight), 0.0001)
# plot(agg, dnorm(agg,mean(FData$Egg.weight),sd(FData$Egg.weight)),type="l", xlab = "egg sizes", ylab= "density")
# plot(agg, o_eyd, type="l")

# Plot DEBage1.size ####
T_o <- NULL
xT_o <- 1:300
for(i in seq(285,293,2)){
  T_o_pars <- c(T = i,     # parameters for Temperature, feeding, allocation and Mass dependence
                 kappa = 0.83, # allocation to respiration (Growth and maintenance)
                 Y=1)         # Feeding level
  T_o <- as.data.frame(rbind(T_o, c(T_o_pars, DEBage1.size(T_o_pars, y=xT_o))))
}
colnames(T_o)[4:ncol(T_o)] <- xT_o
T_o_long <- gather(T_o, key = "Size", value ="biom", c(4:ncol(T_o))) # not use column 

oT <-
  T_o_long %>%
  #filter(Size > 500) %>%
  ggplot(., aes(as.numeric(Size), biom, color = as.factor(T))) +
  geom_line(size=0.85) +
  #scale_colour_discrete() +
  scale_x_continuous("Mass [g]", limits = c(0,300),breaks = seq(0,300,100)) +
  annotate(geom="text", -Inf, Inf, label="D", hjust = -5, vjust = 3, size= 4, fontface = "bold") +
  ylab(expression(paste(italic("o"),"(",italic("y:T"),")"))) +
  scale_colour_brewer(name="Temp [K]", palette="Dark2")+
  theme_bw()
oT

# Plot DEBsurv ####
DEBsurvfun <- function(m, Pars) { # Mass-Temp dependence of yearly Mortality 
  sx <- function(m, Pars){ # natural baseline survival
    exp(-3*m^-0.2*exp(EaS*(Pars[['T']]-T0)/(k*Pars[['T']]*T0))) #3*^-.288 is yearly mortality rate for 1 gram unit weight from Lorenzen 1996
  } #conf limits Lorenzon 1996: 0.315, 0.261
  Vsurvfun <- function(m, z=10.34){
    sxV <- function(m, z=10.34){ # sx from Vindenes 2014
      1/(1+exp(13.53316 - 0.50977*wtol_AYV(m) - (-0.00393)
               *wtol_AYV(m)^2 - 0.19312*z - (-0.00679)
               *wtol_AYV(m)*z))
    }
    xmax <- x[which(sxV(x) == max(sxV(x)))]
    ifelse(m < xmax, sxV(m), sxV(xmax))
  }
  VsurvfunT <- function(m, Pars){
    zT=10.34+Pars[["T"]]-287 # make 287 equal to 283 survival in Vindenes 2014 since 283 is yearly mean Temp 287 is summer mean temp
    sxV <- function(m, z=zT){ # sx from Vindenes 2014
      1/(1+exp(13.53316 - 0.50977*wtol_AYV(m) - (-0.00393)
               *wtol_AYV(m)^2 - 0.19312*z - (-0.00679)
               *wtol_AYV(m)*z))
    }
    xmax <- x[which(sxV(x) == max(sxV(x)))]
    ifelse(m < xmax, sxV(m), sxV(xmax))
  }
  sx.firstyear <- function(m, Pars){ # first year survival
    el_surv#*sx(ratefun(e_m,Pars)[round(sl/2),2,], Pars)
  }
  starvx <- function(m,Pars) { # starvation survival
    starvS <- ratefun(m, Pars)[sl,2,]
    ifelse(starvS < m | starvS > 18000, 0, 1)
  }
  #ifelse(m < mmin, sx.firstyear(m,Pars), sx(m,Pars))#*starvx(m, Pars))# multiply with F (fishing mortality function)
  #ifelse(m < mmin, sx.firstyear(m,Pars), Vsurvfun(m))#*starvx(m, Pars))# multiply with F (fishing mortality function)
  ifelse(m < mmin, sx.firstyear(m,Pars), VsurvfunT(m,Pars))#*starvx(m, Pars))# multiply with F (fishing mortality function)
  #ifelse(m < mmin, sx.firstyear(m,Pars), 0.68)*starvx(m, Pars))# multiply with F (fishing mortality function)
  #ifelse(m < mmin, sx.firstyear(m,Pars), 0.68)# multiply with F (fishing mortality function)
} 
Tsur <- NULL
for(i in seq(285,293,2)){
  Tsur_pars <- c(T = i,     # parameters for Temperature, feeding, allocation and Mass dependence
                 kappa = GR_pars[["kappa"]], # allocation to respiration (Growth and maintenance)
                 Y=GR_pars[["Y"]])
  Tsur <- as.data.frame(rbind(Tsur, c(Tsur_pars,DEBsurvfun(x, Tsur_pars))))
}
colnames(Tsur)[4:ncol(Tsur)] <- c(x)
Tsur_long <- pivot_longer(Tsur, cols = c(4:ncol(Tsur)), names_to = "Size", values_to = "biom") # not use column 10 & 11 (eggstage and recruits)
surT <-
  Tsur_long %>%
  ggplot(., aes(as.numeric(Size), biom, color = as.factor(T))) +
  #geom_line(size=0.85, linetype= "longdash") +
  geom_line(size=0.85, linetype= "solid") +
  #geom_line(size=0.85, linetype= "dotdash") +
  ylab(expression(paste(italic("a"),"(",italic("x,T"),")"))) +
  xlab("Mass [g]") +
  annotate(geom="text", -Inf, Inf, label="C", hjust = -2, vjust = 3, size= 4, fontface = "bold") +
  ylim(0,1) +
  scale_colour_brewer(name="Temp [K]", palette="Dark2") +
  theme_bw()+
  theme(legend.position="none")
surT
#vT<-surT
#vi<-surT
#co<-surT

# Plot individual growth trajectories or size at age of Windermere pike ####
kappp <- 0.83
Temp <- seq(285,293,2)
t <- 20

traj <- NULL
for(i in Temp){
  trajT <- as.data.frame(cbind(ProjCoh(c(T = i, kappa = kappp, Y = 1), t, I_pop=c(1000000,rep(0,n))),rep(c(0,round(x)),t),i, kappp))
  colnames(trajT) <- c("abund","t","mass","Temp", "Kappa")
  trajT <- 
    trajT %>% 
    group_by(t) %>%
    slice(which.max(abund)) 
  traj <- rbind(traj, trajT)
} 
trajG <- traj

groT <- ggplot(trajG) +
  #geom_point(data=GData, aes(Age, ltow_AYV(Length)), alpha=0.1) +
  geom_line(data = trajG, size=0.85, aes(t, mass, color = as.factor(Temp))) +
  ylab("Mass [g]") +
  xlab("Age [years]") +
  annotate(geom="text", -Inf, Inf, label="A", hjust = -2, vjust = 3, size= 4, fontface = "bold")+
  theme_bw() +
  scale_colour_brewer(palette="Dark2",name = "Temp [K]") +
  theme(legend.position="none")
groT

# Plot fecundity from DEB and in Lake Windermere
rep=NULL
for(i in 1:nrow(trajG)){
  if(trajG$t[i]==1){
    rep_tplus1 <- cbind(trajG$mass[i], DEBrepfun(e_m, c(T = trajG$Temp[i], kappa = kappp, Y = 1)),trajG$t[i],trajG$Temp[i], kappp)} 
  else{ rep_tplus1 <- cbind(trajG$mass[i], 
                            DEBrepfun(trajG$mass[i-1], c(T = trajG$Temp[i], kappa = kappp, Y = 1)),trajG$t[i],trajG$Temp[i], kappp)}
  
  colnames(rep_tplus1) <- c("mass","fec","t","Temp", "Kappa")
  rep <- rbind(rep, rep_tplus1)
}

# Plot DEBrepfun and windermere data ###
# plot(FData$Weight,FData$Eggs,ylab="#eggs in t+1", xlab = "weight in t+1",
#      main = "fecundity in t+1", ylim = c(0,600000),xlim = c(0,18000)) # Windermere Pike data
# abline(lm(FData$Eggs~FData$Weight), col="blue",lwd=1, lty=2)      
# as.tibble(rep) %>%
#   filter(Temp == 287) %>%
#   lines(fec~mass, data=., lwd = 2, type = "l", col = "red")
# legend("topleft", c("DEB", "Windermere Pike", "lm(Windermere Pike)"), 
#        lty = c(1, NA), pch = c(NA,1), col= c("red","black", "blue"),cex=.7)

repT <- ggplot()+
  #geom_point(data=FData, aes(Weight, Eggs), alpha=0.1)+
  geom_line(data=as.tibble(rep), size=0.85, aes(mass, fec, color = as.factor(Temp)))+
#  geom_smooth(data=FData, aes(Weight, Eggs),method=lm, lwd=0.5, linetype="dashed")+  
  scale_colour_brewer(palette="Dark2", name = "Temp [K]")+
  annotate(geom="text", -Inf, Inf, label="B", hjust = -2, vjust = 3, size= 4, fontface = "bold") +
  ylab(expression(paste(italic("f"),"(",italic("x,T"),")"))) +
  xlab("Mass [g]")+
  theme_bw()
repT

pdf("DEBIPM_vitalrates_20210629_nopoints.pdf", width = 8, height = 6)
(groT + repT) / (surT + oT)
dev.off()

groT
par(mfrow=c(2,2), las=1, bty="l")
par(mfrow=c(1,1), las=1, bty="l")
groT

#VVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVV
#### ### ## # PLOT RESULTS #### ### ## # ## ### ####
#VVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVV
## READ in results from high res Main and Contrast results ##

# Results Baseline scenario
Res_lam <- read.delim("//storage-og.slu.se/home$/vitl0001/My Documents/Manus2/Results/Results0624/Res_lam_1_0624.txt", sep = ",")
#Res_lam <- Res_lam[1:1149,]
Res_v <- read.delim("//storage-og.slu.se/home$/vitl0001/My Documents/Manus2/Results/Results0624/Res_v_1_0624.txt", sep = ",")
colnames(Res_v)[4:ncol(Res_v)] <- sub("X", "", colnames(Res_v)[4:ncol(Res_v)])
#Res_v <- Res_lam[1:1149,]
Res_w <- read.delim("//storage-og.slu.se/home$/vitl0001/My Documents/Manus2/Results/Results0624/Res_w_1_0624.txt", sep = ",")
colnames(Res_w)[4:ncol(Res_w)] <- sub("X", "", colnames(Res_w)[4:ncol(Res_w)])

### PLOT LAMBDA AS A FUNCTION OF KAPPA AND/OR TEMP ###
#as_tibble(Res_lam) %>%
#   filter(T %in% seq(min(Res_lam$T),max(Res_lam$T), length.out=13)) %>%
#   ggplot(., aes(kappa, Lambda, color = as.factor(T))) +
#   geom_line(size=0.5) +
#   ylab(expression(lambda~(Fitness))) +
#   xlab(expression(kappa~(Growth~Allocation))) +
#   labs(color = "Temp [K]") +
#   theme_bw()


### FIG 4A - PLOT A SURFACE OF FITNESS OVER TEMP (x) AND KAPPA (y) ####
Res_lam_1 <- read.delim("//storage-og.slu.se/home$/vitl0001/My Documents/Manus2/Results/Results0624/Res_lam_1_0624.txt", sep = ",")
maxl_1 <- as.data.frame(Res_lam_1) %>% 
  group_by(T) %>%
  slice_max(Lambda)
maxl_1["Sur_f"] <- "Vin"

Fig4a <- 
  as_tibble(Res_lam_1) %>%
  filter(kappa <= 1) %>%
  #mutate(kappa = kappa > 1) # set kappa > 1 to 1
  ggplot(., aes(T,kappa)) +
  geom_raster(aes(fill=round(Lambda, 4)))+#, interpolate = TRUE) +
  #geom_smooth(data=maxl_1 ,aes(T,kappa), size=0.85, se=FALSE) +
  #ggplot(., aes(T,kappa, z=round(Lambda, 4))) +#, color = as.factor(T))) +#, linetype = as.factor(kappa))) +
  # geom_contour_filled(breaks = 
  #                       c(0,seq(1,round(max(as_tibble(Res_lam)$Lambda),1),0.05),
  #                               max(as_tibble(Res_lam)$Lambda)))+
  geom_line(data=maxl_1 ,aes(T,kappa),size=0.85) +
  scale_fill_gradientn(
    limits = c(0,max(as.data.frame(Res_lam_1[,'Lambda'])+0.0001)),
    colours = c("black","white","yellow","red"),
    values = rescale(c(0,0.9999,1, round(max(as.data.frame(Res_lam_1[,'Lambda'])), 4)),
                     to=c(0,1)),
    breaks = c(0,1,round(max(as.data.frame(Res_lam_1[,'Lambda'])), 2)),
    labels=c("0",1,round(max(as.data.frame(Res_lam_1[,'Lambda'])),3)),
    name = expression(lambda)) +
  scale_x_continuous(expand = c(0,0), name = "Temperature [K]", limits = c(284,294),breaks = seq(285,293,2))+
  scale_y_continuous(expand = c(0,0), name=expression(paste(kappa," (Growth allocation)")), limits = c(0.6,1))+
  coord_cartesian(xlim = c(284,294)) +
  annotate(geom="text", -Inf, Inf, label="A", hjust = -26, vjust = 3, size= 4, fontface = "bold")+
  theme_bw()
Fig4a

### FIG 4B - OPTIMAL KAPPA FOR THREE SURVIVAL SCENARIOS ####
Res_lam_2 <- read.delim("//storage-og.slu.se/home$/vitl0001/My Documents/Manus2/Results/Results0624/Res_lam_2_0624.txt", sep = ",")
maxl_2 <- as.data.frame(Res_lam_2) %>% 
  group_by(T) %>%
  slice_max(Lambda)
maxl_2["Sur_f"] <- "NotempVin"

Res_lam_3 <- read.delim("//storage-og.slu.se/home$/vitl0001/My Documents/Manus2/Results/Results0624/Res_lam_3_0624.txt", sep = ",")
maxl_3 <- as.data.frame(Res_lam_3) %>% 
  group_by(T) %>%
  slice_max(Lambda)
maxl_3["Sur_f"] <- "flat"

Res_lam_4 <- read.delim("//storage-og.slu.se/home$/vitl0001/My Documents/Manus2/Results/Results0624/Res_lam_1_0624.txt", sep = ",")
maxl_4 <- as.data.frame(Res_lam_4) %>% 
  group_by(T) %>%
  slice_max(Lambda)
maxl_4["Sur_f"] <- "Nostarv"

maxl_123 <- as_tibble(rbind(maxl_1,maxl_2,maxl_3))

Fig4b <- 
maxl_123 %>%
  filter(kappa <= 1) %>%
  ggplot(., aes(T,kappa, linetype = Sur_f)) +
  #geom_raster(aes(fill=round(Lambda, 4)))+#, interpolate = TRUE) +
  #geom_smooth(data=maxl_123, aes(T,kappa), size=0.85, se=FALSE, color="black") +
  geom_line(data=maxl_123, aes(T,kappa),size=0.85) +
  scale_y_continuous(expand = c(0,0), name=expression(paste(kappa," (Growth allocation)")), limits = c(0.6,1))+
  scale_x_continuous(expand = c(0,0), name = "Temperature [K]", limits = c(284,294),breaks = seq(283,295,2))+
  scale_linetype_manual(values = c("dashed","dotdash","solid"),
                        name= "Survival",
                        #name = expression(paste(italic("a"),"(",italic("T,x"),")")), 
                         labels = c("a=0.68",
                                    expression(paste(italic("a"),"(",italic("x"),")")),
                                    expression(paste(italic("a"),"(",italic("T,x"),")"))),
                         guide = guide_legend(reverse = TRUE) ) +
  annotate(geom="text", -Inf, Inf, label="B", hjust = -26, vjust = 3, size= 4, fontface = "bold") +
  coord_cartesian(xlim = c(284,294)) +
  theme_bw() +
  theme( legend.key.width = unit(1, 'cm')) 
Fig4b

### PLOT FIG 4C-E, THREE SURVIVAL SCENARIOS  ####
# Baseline model survival 
DEBsurvfun <- function(m, Pars) { # Mass-Temp dependence of yearly Mortality 
  sx <- function(m, Pars){ # natural baseline survival
    exp(-3*m^-0.2*exp(EaS*(Pars[['T']]-T0)/(k*Pars[['T']]*T0))) #3*^-.288 is yearly mortality rate for 1 gram unit weight from Lorenzen 1996
  } #conf limits Lorenzon 1996: 0.315, 0.261
  Vsurvfun <- function(m, z=10.34){
    sxV <- function(m, z=10.34){ # sx from Vindenes 2014
      1/(1+exp(13.53316 - 0.50977*wtol_AYV(m) - (-0.00393)
               *wtol_AYV(m)^2 - 0.19312*z - (-0.00679)
               *wtol_AYV(m)*z))
    }
    xmax <- x[which(sxV(x) == max(sxV(x)))]
    ifelse(m < xmax, sxV(m), sxV(xmax))
  }
  VsurvfunT <- function(m, Pars){
    zT=10.34+Pars[["T"]]-287 # make 287 equal to 283 survival in Vindenes 2014 since 283 is yearly mean Temp 287 is summer mean temp
    sxV <- function(m, z=zT){ # sx from Vindenes 2014
      1/(1+exp(13.53316 - 0.50977*wtol_AYV(m) - (-0.00393)
               *wtol_AYV(m)^2 - 0.19312*z - (-0.00679)
               *wtol_AYV(m)*z))
    }
    xmax <- x[which(sxV(x) == max(sxV(x)))]
    ifelse(m < xmax, sxV(m), sxV(xmax))
  }
  sx.firstyear <- function(m, Pars){ # first year survival
    el_surv#*sx(ratefun(e_m,Pars)[round(sl/2),2,], Pars)
  }
  starvx <- function(m,Pars) { # starvation survival
    starvS <- ratefun(m, Pars)[sl,2,]
    ifelse(starvS < m | starvS > 18000, 0, 1)
  }
  ifelse(m < mmin, sx.firstyear(m,Pars), VsurvfunT(m,Pars))#*starvx(m, Pars))# multiply with F (fishing mortality function)
} 
Tsur <- NULL
for(i in seq(285,293,2)){
  Tsur_pars <- c(T = i,     # parameters for Temperature, feeding, allocation and Mass dependence
                 kappa = GR_pars[["kappa"]], # allocation to respiration (Growth and maintenance)
                 Y=GR_pars[["Y"]])
  Tsur <- as.data.frame(rbind(Tsur, c(Tsur_pars,DEBsurvfun(x, Tsur_pars))))
}
colnames(Tsur)[4:ncol(Tsur)] <- c(x)
Tsur_long <- pivot_longer(Tsur, cols = c(4:ncol(Tsur)), names_to = "Size", values_to = "biom") # not use column 10 & 11 (eggstage and recruits)
Fig4c <-
  Tsur_long %>%
  ggplot(., aes(as.numeric(Size), biom, color = as.factor(T))) +
  geom_line(size=0.85, linetype= "solid") +
  ylab(expression(paste(italic("a"),"(",italic("x,T"),")"))) +
  xlab("Mass [g]") +
  annotate(geom="text", -Inf, Inf, label="C", hjust = -1.5, vjust = 2.5, size= 4, fontface = "bold") +
  ylim(0,1) +
  scale_colour_brewer(name="Temp [K]", palette="Dark2") +
  theme_bw() +
  theme( legend.key.width = unit(1, 'cm') ) 
Fig4c

# Temeperature independent Vindenes et al. 2014
DEBsurvfun <- function(m, Pars) { # Mass-Temp dependence of yearly Mortality 
  sx <- function(m, Pars){ # natural baseline survival
    exp(-3*m^-0.2*exp(EaS*(Pars[['T']]-T0)/(k*Pars[['T']]*T0))) #3*^-.288 is yearly mortality rate for 1 gram unit weight from Lorenzen 1996
  } #conf limits Lorenzon 1996: 0.315, 0.261
  Vsurvfun <- function(m, z=10.34){
    sxV <- function(m, z=10.34){ # sx from Vindenes 2014
      1/(1+exp(13.53316 - 0.50977*wtol_AYV(m) - (-0.00393)
               *wtol_AYV(m)^2 - 0.19312*z - (-0.00679)
               *wtol_AYV(m)*z))
    }
    xmax <- x[which(sxV(x) == max(sxV(x)))]
    ifelse(m < xmax, sxV(m), sxV(xmax))
  }
  VsurvfunT <- function(m, Pars){
    zT=10.34+Pars[["T"]]-287 # make 287 equal to 283 survival in Vindenes 2014 since 283 is yearly mean Temp 287 is summer mean temp
    sxV <- function(m, z=zT){ # sx from Vindenes 2014
      1/(1+exp(13.53316 - 0.50977*wtol_AYV(m) - (-0.00393)
               *wtol_AYV(m)^2 - 0.19312*z - (-0.00679)
               *wtol_AYV(m)*z))
    }
    xmax <- x[which(sxV(x) == max(sxV(x)))]
    ifelse(m < xmax, sxV(m), sxV(xmax))
  }
  sx.firstyear <- function(m, Pars){ # first year survival
    el_surv#*sx(ratefun(e_m,Pars)[round(sl/2),2,], Pars)
  }
  starvx <- function(m,Pars) { # starvation survival
    starvS <- ratefun(m, Pars)[sl,2,]
    ifelse(starvS < m | starvS > 18000, 0, 1)
  }
  ifelse(m < mmin, sx.firstyear(m,Pars), Vsurvfun(m))#*starvx(m, Pars))# multiply with F (fishing mortality function)
}

Tsur <- NULL
for(i in seq(285,293,2)){
  Tsur_pars <- c(T = i,     # parameters for Temperature, feeding, allocation and Mass dependence
                 kappa = GR_pars[["kappa"]], # allocation to respiration (Growth and maintenance)
                 Y=GR_pars[["Y"]])
  Tsur <- as.data.frame(rbind(Tsur, c(Tsur_pars,DEBsurvfun(x, Tsur_pars))))
}
colnames(Tsur)[4:ncol(Tsur)] <- c(x)
Tsur_long <- pivot_longer(Tsur, cols = c(4:ncol(Tsur)), names_to = "Size", values_to = "biom") # not use column 10 & 11 (eggstage and recruits)
Fig4d <-
  Tsur_long %>%
  ggplot(., aes(as.numeric(Size), biom, color = as.factor(T))) +
  #geom_line(size=0.85, linetype= "longdash") +
  geom_line(size=0.85, linetype= "dotdash") +
  #geom_line(size=0.85, linetype= "dotdash") +
  ylab(expression(paste(italic("a"),"(",italic(x),")"))) +
  xlab("Mass [g]") +
  annotate(geom="text", -Inf, Inf, label="D", hjust = -1.5, vjust = 2.5, size= 4, fontface = "bold") +
  ylim(0,1) +
  scale_colour_manual(name="Temp [K]", values=c("black","black","black","black","black")) +
  theme_bw() +
  theme( legend.key.width = unit(1, 'cm') )
Fig4d


# Constant survival (0.68)
DEBsurvfun <- function(m, Pars) { # Mass-Temp dependence of yearly Mortality 
  sx <- function(m, Pars){ # natural baseline survival
    exp(-3*m^-0.2*exp(EaS*(Pars[['T']]-T0)/(k*Pars[['T']]*T0))) #3*^-.288 is yearly mortality rate for 1 gram unit weight from Lorenzen 1996
  } #conf limits Lorenzon 1996: 0.315, 0.261
  Vsurvfun <- function(m, z=10.34){
    sxV <- function(m, z=10.34){ # sx from Vindenes 2014
      1/(1+exp(13.53316 - 0.50977*wtol_AYV(m) - (-0.00393)
               *wtol_AYV(m)^2 - 0.19312*z - (-0.00679)
               *wtol_AYV(m)*z))
    }
    xmax <- x[which(sxV(x) == max(sxV(x)))]
    ifelse(m < xmax, sxV(m), sxV(xmax))
  }
  VsurvfunT <- function(m, Pars){
    zT=10.34+Pars[["T"]]-287 # make 287 equal to 283 survival in Vindenes 2014 since 283 is yearly mean Temp 287 is summer mean temp
    sxV <- function(m, z=zT){ # sx from Vindenes 2014
      1/(1+exp(13.53316 - 0.50977*wtol_AYV(m) - (-0.00393)
               *wtol_AYV(m)^2 - 0.19312*z - (-0.00679)
               *wtol_AYV(m)*z))
    }
    xmax <- x[which(sxV(x) == max(sxV(x)))]
    ifelse(m < xmax, sxV(m), sxV(xmax))
  }
  sx.firstyear <- function(m, Pars){ # first year survival
    el_surv#*sx(ratefun(e_m,Pars)[round(sl/2),2,], Pars)
  }
  starvx <- function(m,Pars) { # starvation survival
    starvS <- ratefun(m, Pars)[sl,2,]
    ifelse(starvS < m | starvS > 18000, 0, 1)
  }
  ifelse(m < mmin, sx.firstyear(m,Pars), 0.68)#*starvx(m, Pars))# multiply with F (fishing mortality function)
} 
Tsur <- NULL
for(i in seq(285,293,2)){
  Tsur_pars <- c(T = i,     # parameters for Temperature, feeding, allocation and Mass dependence
                 kappa = GR_pars[["kappa"]], # allocation to respiration (Growth and maintenance)
                 Y=GR_pars[["Y"]])
  Tsur <- as.data.frame(rbind(Tsur, c(Tsur_pars,DEBsurvfun(x, Tsur_pars))))
}
colnames(Tsur)[4:ncol(Tsur)] <- c(x)
Tsur_long <- pivot_longer(Tsur, cols = c(4:ncol(Tsur)), names_to = "Size", values_to = "biom") # not use column 10 & 11 (eggstage and recruits)

Fig4e <-
  Tsur_long %>%
  ggplot(., aes(as.numeric(Size), biom, color = as.factor(T))) +
  geom_line(size=0.85, linetype= "dashed") +
  ylab("a=0.68") +
  xlab("Mass [g]") +
  annotate(geom="text", -Inf, Inf, label="E", hjust = -1.5, vjust = 2.5, size= 4, fontface = "bold") +
  ylim(0,1) +
  scale_colour_manual(name="Temp [K]", values=c("black","black","black","black","black")) +
  theme_bw() +
  theme( legend.key.width = unit(1, 'cm')) 
Fig4e

#pdf("DEBIPM_fig4_20210630.pdf", width = 10, height = 6)
(Fig4a / Fig4b) | (Fig4c/Fig4d/Fig4e)
dev.off()

### PLOT STABLE STRUCTURE W ####
Kopt <- 
  as_tibble(maxl_1) %>%
  filter(T %in% c(287,289,291)) #%>%

Res_w_long <- pivot_longer(as.data.frame(Res_w), cols = c(4:ncol(Res_w)), 
                           names_to = "Size", values_to ="biom") # not use column 4 & 5 (eggstage and recruits)

Fig3a_1 <- Res_w_long %>%
  filter(T %in% c(285, 287, 289, 291)) %>%
  filter(as.character(kappa) %in% 0.83) %>% #Floating point issue when comparing vector, therefore the use of %in%, can also use near()
  ggplot(., aes(as.numeric(Size), biom, group = rev(as.factor(T)), color = as.factor(T))) +
  geom_line(size=0.8)+
  ylab("Stable structure w")+
  xlab("Weight [g]")+
  ylim(0,6e-7)+
  annotate(geom="text", -Inf, Inf, label="A", hjust = -24, vjust = 2, size= 4, fontface = "bold") +
  scale_colour_brewer(palette="Dark2",name = "Temp [K]") +
  theme_bw() +
  labs(color = "Temp [K]")
Fig3a_1

Fig3a_2 <- Res_w_long %>%
  filter(T %in% c(285, 287, 289, 291)) %>%
  filter(as.character(kappa) %in% 0.83) %>% #Floating point issue when comparing vector, therefore the use of %in%, can also use near()
  ggplot(., aes(as.numeric(Size), biom, group = rev(as.factor(T)), color = as.factor(T))) +
  geom_line(size=0.8)+
  ylab("Stable structure w")+
  xlab("Weight [g]")+
  ylim(0,5.5e-5)+
  scale_x_continuous(breaks = c(0,200,400), limits=c(0,400))+
  scale_colour_brewer(palette="Dark2",name = "Temp [K]") +
  theme_bw() +
  theme(legend.position="none",
        axis.title.x=element_blank(),
        plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
        axis.title.y=element_blank())
Fig3a_2

mycolors <- rev(colorRampPalette(brewer.pal(9, "RdPu"))(9))
Fig3b <- 
  as_tibble(Res_w_long) %>%
  filter(kappa %in% 0.83) %>%
  mutate(biom_log=log(biom)) %>%
  ggplot(., aes(T,as.double(Size), z=biom_log, colour=biom_log))+
  geom_contour_filled(breaks = c(min(log(Res_w_long$biom)), #-12
                                seq(-25,-15,length.out = 7),
                                max(log(Res_w_long$biom))))+ #-160
  annotate(geom="text", -Inf, Inf, label="B", hjust = -24, vjust = 2, size= 4, fontface = "bold")+
  scale_fill_manual(values = rev(mycolors), name="Relative densities", 
                     labels=c(paste("< ",signif(exp(-25),5)), 
                              signif(seq(exp(-25),exp(-15),length.out = 6),5),
                              #expression("">=signif(max(Res_w_long$biom),5))))+
                              paste(">",signif(max(Res_w_long$biom),5)))) + 
  scale_x_continuous(expand = c(0,0), name = "Temperature [K]", limits = c(284,294),breaks = seq(285,293,2))+
  scale_y_continuous(expand = c(0,0), name = "Weight [g]", limits = c(0,18000))+
  coord_cartesian(xlim = c(284,294)) +
  theme_bw()
Fig3b


pdf("DEBIPM_w_v_1_20210630.pdf", width = 10, height = 4)
subvp <- viewport(width = 0.19, height = 0.42, x = 0.25, y = 0.65)
(Fig3a_1+Fig3b) #/ (d_1+d_2)
print(Fig3a_2, vp = subvp)
dev.off()



  maxbio<- as.data.frame(Res_w_long) %>% 
  filter(kappa %in% 0.83) %>%
  filter(as.double(Size) > 400 & as.double(Size) < 1500) %>%
  group_by(T) %>%
  slice_max(biom)

### PLOT REPRODUCTIVE VALUES V ####
#colnames(Res_v) <- c("T","kappa","Y","0",x)
Res_v_long <- pivot_longer(as_tibble(Res_v), cols = c(4:ncol(Res_v)), 
                             names_to = "Size", values_to ="biom") # not use column 4 & 11 (eggstage and recruits)
  
#   OpCov_287 <- Res_v_long %>%
#     filter(T %in% Kopt$T[1]) %>%
#     filter(kappa %in% c(Kopt$kappa[1],
#                         Kopt$kappa[1]+0.01,0.94))
#   OpCov_289 <- Res_v_long %>%
#     filter(T %in% Kopt$T[2]) %>%
#     filter(kappa %in% c(Kopt$kappa[2],
#                         Kopt$kappa[2]+0.01,0.93))
#   OpCov_291 <- Res_v_long %>%
#     filter(T %in% Kopt$T[3]) %>%
#     filter(kappa %in% c(Kopt$kappa[3],
#                         Kopt$kappa[3]+0.01,91))
#   OptCov <- rbind(OpCov_287,OpCov_289,OpCov_291)
#   
#   
# OptCov %>%
#   ggplot(., aes(as.numeric(Size), biom, color = as.factor(kappa), linetype = as.factor(T))) +#, linetype = as.factor(kappa))) +
#   ggtitle("Repro. values over Size at 289") + 
#   geom_line(size=0.8) +
#   ylab("Reproductive value V") +
#   xlab("Size")+
#   theme_bw() +
#   labs(color = "Kappa (growth alloc.)")

d_1 <-  
Res_v_long %>%
  filter(T %in% c(285,287,289,291)) %>%
  #filter(T %in% c(Kopt$T[2],Kopt$T[2]+1,Kopt$T[2]-1)) %>%  
  #filter(kappa %in% c(Kopt$kappa[2])) %>%
  filter(kappa == 0.83) %>%
  #filter(T %in% c(Kopt$T[2],Kopt$T[2]+1,Kopt$T[2]-1)) %>%
  ggplot(., aes(as.numeric(Size), biom, color = as.factor(T))) +#, linetype = as.factor(kappa))) +
  geom_line(size=0.8)+
  ylab("Reproductive value v")+
  xlab("Weight [g]")+
  annotate(geom="text", -Inf, Inf, label="C", hjust = -18, vjust = 2, size= 4, fontface = "bold")+
  theme_bw()+
  scale_colour_brewer(palette="Dark2",name = "Temp [K]") #+
d_1 

d_2 <-
Res_v_long %>%
  filter(kappa == 0.83) %>%
  ggplot(., aes(T,as.numeric(Size), z=biom)) +#, color = as.factor(T))) +#, linetype = as.factor(kappa))) +
  geom_contour_filled()+
  scale_x_continuous(name = "Temperature [K]", limits = c(283,293),breaks = seq(283,293,2))+
  annotate(geom="text", -Inf, Inf, label="D", hjust = -20, vjust = 2, size= 4, fontface = "bold")+
  ylab("Weight [g]")+
  #scale_fill_manual(values = rev(mycolors), name="Reproductive value v)")+
  theme_bw()
d_2

# grid.arrange(
# tableGrob(DEBparams, rownames(DEBparams)),
# tableGrob(IPMparams, rownames(IPMparams)), nrow=1
#)

### PLOT RESULTING OPTIMAL W, V & growth/fecundity trajectories ####

Res_w_long %>%
  semi_join(maxl_1, Res_w_long, by=c('kappa', 'T')) %>%
  filter(T %in% c(285, 287, 289, 291, 293)) %>%
  #filter(kappa %in% maxl_1$kappa) #%>% #Floating point issue when comparing vector, therefore the use of %in%, can also use near()
  ggplot(., aes(as.numeric(Size), biom, group = rev(as.factor(T)), color = as.factor(T))) +
  geom_line(size=0.8)+
  ylab("Stable structure w")+
  xlab("Weight [g]")+
  ylim(0,2.5e-8)+
  #annotate(geom="text", -Inf, Inf, label="A", hjust = -21, vjust = 2, size= 4, fontface = "bold") +
  scale_colour_brewer(palette="Dark2",name = "Temp [K]") +
  theme_bw() +
  labs(color = "Temp [K]")

Res_v_long %>%
  semi_join(maxl_1, Res_v_long, by=c('kappa', 'T')) %>%
  filter(T %in% c(285, 287, 289, 291, 293)) %>%
  ggplot(., aes(as.numeric(Size), biom, color = as.factor(T))) +#, linetype = as.factor(kappa))) +
  geom_line(size=0.8) +
  ylab("Reproductive value v") +
  xlab("Weight [g]")+
  #annotate(geom="text", -Inf, Inf, label="C", hjust = -18, vjust = 2, size= 4, fontface = "bold")+
  theme_bw()+
  scale_colour_brewer(palette="Dark2",name = "Temp [K]") #+

# Optimal growth trajectories

t <- 20
OptPars <- 
as.data.frame(maxl_1[,1:3]) %>%
  filter(T %in% c(285, 287, 289, 291, 293))

trajOpt <- NULL
for(i in 1:nrow(OptPars)){
  trajT <- as.data.frame(cbind(ProjCoh(unlist(OptPars[i,]), t, I_pop=c(1000000,rep(0,n))),rep(c(0,round(x)),t),OptPars[i,1],OptPars[i,2]))
  colnames(trajT) <- c("abund","t","mass","Temp", "Kappa")
  trajT <- 
    trajT %>% 
    group_by(t) %>%
    slice(which.max(abund)) 
  trajOpt <- rbind(trajOpt, trajT)
} 

trajGOpt <- trajOpt
trajGOpt["strat"] <- "Opt"
trajG["strat"] <- "def"
trajBoth <- rbind(trajGOpt,trajG)

str(trajGOpt)
groTopt <- 
trajBoth %>%
  ggplot(.) +
  geom_line(size=0.85,  aes(t, mass, color = as.factor(Temp), alpha = strat)) +
  ylab("Mass [g]") +
  xlab("Age [years]") +
  scale_alpha_discrete(range=c(0.4,1), name = expression(kappa["0"]),labels = c(expression(kappa["0"]^"*"),"0.83")) +
  annotate(geom="text", -Inf, Inf, label="A", hjust = -2, vjust = 3, size= 4, fontface = "bold")+
  theme_bw() +
  scale_colour_brewer(palette="Dark2",name = "Temp [K]") +
  theme(legend.position="none")
groTopt

# Plot optimal fecundity from DEB and in Lake Windermere
repOpt=NULL
for(i in 1:nrow(trajGOpt)){
  if(trajGOpt$t[i]==1){
    rep_tplus1opt <- cbind(trajGOpt$mass[i], DEBrepfun(e_m, c(T = trajGOpt$Temp[i], kappa = trajGOpt$Kappa[i], Y = 1)), trajGOpt$t[i],trajGOpt$Temp[i], trajGOpt$Kappa[i])} 
  else{rep_tplus1opt <- cbind(trajGOpt$mass[i], 
                              DEBrepfun(trajGOpt$mass[i-1], c(T = trajGOpt$Temp[i], kappa = trajGOpt$Kappa[i], Y = 1)), trajGOpt$t[i], trajGOpt$Temp[i], trajGOpt$Kappa[i])}
  
  colnames(rep_tplus1opt) <- c("mass","fec","t","Temp", "Kappa")
  repOpt <- rbind(repOpt, rep_tplus1opt)
}

repOpt<- as.data.frame(repOpt)
rep <- as.data.frame(rep)
repOpt["strat"] <- "opt"
rep["strat"] <- "def"
repBoth <- rbind(repOpt,rep)


repTopt <- 
ggplot()+
  geom_line(data=as.tibble(repBoth), size=0.85, aes(mass, fec, color = as.factor(Temp), alpha= strat))+
  scale_alpha_discrete(range=c(0.4,1), name = expression(kappa["0"]),labels = c(expression(kappa["0"]^"*"),"0.83")) +
  scale_colour_brewer(palette="Dark2", name = "Temp [K]")+
  annotate(geom="text", -Inf, Inf, label="B", hjust = -2, vjust = 3, size= 4, fontface = "bold") +
  ylab(expression(paste(italic("f"),"(",italic("x,T"),")"))) +
  xlab("Mass [g]")+
  theme_bw()
repTopt

groTopt+repTopt

# Res_v_long %>%
  #   filter(kappa %in% c(0.7, 0.8)) %>%
  #   filter(T %in% c(283))  %>% 
  #   ggplot(., aes(as.numeric(Size), biom, color = as.factor(kappa))) +
  #   ggtitle("Repro. values over Size with three kappas") + 
  #   geom_line(size=0.5) +
  #   ylab("Reproductive value V") +
  #   xlab("Size")+
  #   xlim(0,14000) +
  #   theme_bw() +
  #   labs(color = expression(kappa))
  

### COHORT PROJECTPONS #### 

t=20
T285 <- cbind(ProjCoh(c(T = 288.2, kappa = 0.83, Y = 1),t),rep(c(0,round(x)),t))
colnames(T285) <- c("abund","t","size")
T287 <- cbind(ProjCoh(c(T = 287, kappa = 0.83, Y = 1),t),rep(c(0,round(x)),t))
colnames(T287) <- c("abund","t","size")
T289 <- cbind(ProjCoh(c(T = 289, kappa = 0.83, Y = 1),t),rep(c(0,round(x)),t))
colnames(T289) <- c("abund","t","size")
  
PT285 <- 
    as.tibble(T285) %>%
    filter(t != 1)%>%
    ggplot(., aes(x=size , y=abund)) +
    geom_line() +
    facet_wrap(.~t, nrow = max(t), strip.position= "right")+#, scales = "free_y") +#) +#, scales = "free")+
    geom_vline(xintercept = 1000, color="red") +
    ggtitle("285 K") + 
    #ylim(0,20) +
    theme_bw()
  PT287 <- 
    as.tibble(T287) %>%
    filter(t != 1)%>%
    ggplot(., aes(x=size , y=abund)) +
    geom_line() +
    facet_wrap(.~t, nrow = max(t), strip.position= "right")+#, scales = "free_y") +#, scales = "free")+
    geom_vline(xintercept = 1000, color="red") +
    ggtitle("287 K") + 
    #ylim(0,20) +
    theme_bw()
  PT289 <- 
    as.tibble(T289) %>%
    filter(t != 1)%>%
    ggplot(., aes(x=size , y=abund)) +
    geom_line() +
    facet_wrap(.~t, nrow = max(t), strip.position= "right")+#, scales = "free_y") +#, scales = "free")+
    geom_vline(xintercept = 1000, color="red") +
    ggtitle("289 K") + 
    theme_bw()
  
  #pdf("MS2_Vindenes2014survCoh.pdf", width = 8, height = 6)
  PT285+PT287+PT289
  #dev.off()
  
  
  # ### Size dependent cost #####
  # 
  # r0 <- 100
  # r1 <- 0.0001
  # eps1 <- 0.51  # Intake allometric scalar 
  # eps2 <- 0.58   # Intake allometric exponent
  # kap_fun <- function(m, kappa=0.89, ha=3){kappa*exp(-m/(ha*max(x)))} #for size dep. kappa
  # 
  # GR_pars <- c(T = 287,     # parameters for Temperature, feeding, allocation and Mass dependence
  #              kappa = 0.81, # allocation to respiration (Growth and maintenance)
  #              Y = 1)         # Feeding level
  # outR<-NULL
  # for(seas in 1:20) { 
  #   if(seas == 1) {
  #     offmean <- ratefun(e_m, GR_pars)[round(sl/2),2,]
  #     p <- data.frame(mass = max(ratefun(offmean, GR_pars)[,2,]),
  #                     re = DEBrepfun(offmean, GR_pars, e_m) - r0*max(ratefun(offmean, GR_pars)[,2,])^r1, 
  #                     age = seas+1,
  #                     Temp = GR_pars["T"],
  #                     Kappa = GR_pars["Kappa"],
  #                     Y = GR_pars["Y"], row.names = seas) #+1 as seas=1 is to age 2, offmean represents age 1
  #   } else p <- rbind(p,
  #                     c(max(ratefun(p[nrow(p),1], GR_pars)[,2,]),
  #                       DEBrepfun(p[nrow(p),1], GR_pars, e_m) - r0*max(ratefun(p[nrow(p),1], GR_pars)[,2,])^r1, 
  #                       seas+1, GR_pars[1], GR_pars[2], GR_pars[3]))
  # } 
  # p <- rbind(c(offmean, DEBrepfun(e_m, GR_pars, e_m), 
  #              1,GR_pars[1], GR_pars[2], GR_pars[3])
  #            , p ) # add Age 1
  # outR <- bind_rows(outR,p)
  # 
  # plot(FData$Weight,FData$Eggs,ylab="#eggs in t+1", xlab = "weight in t+1",
  #      main = "fecundity of y in t+1", ylim = c(0,400000)) # Windermere Pike data
  # outR %>%
  #   lines(re~mass, data=., lwd = 2, type = "l", col = "red")
  # outR %>%
  #   plot(mass~age, data=., lwd = 2, type = "l", col = "red", ylim = c(0,16000))
  # points(GData2$Age, ltow_AYV(GData2$Length))
  # 
  
  # Size distribution for Life cycle graph ####
Res_wLC <- bind_cols(wvlambda.projection(K.matrix(Pars <- c(T = 285,     # parameters for Temperature, feeding, allocation and Mass dependence
                                         kappa = 0.83, # allocation to respiration (Growth and maintenance)
                                         Y = 1) ))[2], c(0,x))

colnames(  Res_wLC)<-c("biom","Size")
  #c_1b <- 
Res_wLC %>%
    ggplot(., aes(as.numeric(Size), biom, color="red")) +
    geom_line(size=1.2)+
    ylab("Density")+
    xlab("Size")+
    #ylim(0,1.1e-7)+
    #annotate(geom="text", -Inf, Inf, label="A", hjust = -21, vjust = 2, size= 4, fontface = "bold") +
    #coord_cartesian(xlim = c(350,3000), expand = c(0,0)) +
    xlim(350,3000)+
    ylim(0,1.2e-8)+
    #xlim(0,4000)+
    #scale_x_continuous(breaks = c(0,500,1000), limits=c(0,1000))+
    #scale_x_continuous(limits = c(500,12000)) +
  #geom_vline(xintercept = 1000, color="red") +
    theme_bw() +
    theme(legend.position="none",
          #axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          #axis.ticks.x=element_blank(),
          axis.text.y=element_blank())
          #axis.ticks.y=element_blank(),
          #axis.title.y=element_blank())
  
  