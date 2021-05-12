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
### Plot mass functions
plot(x, kap_fun(x,kappa=0.83), xlab = "Weight", ylab = "kappa", type="l", ylim= c(0,1),main = "Size-dependent kappa")
lines(x, 1-kap_fun(x, kappa=0.83), ylab = "1-kappa", type="l", main = "Size-dependent 1-kappa", col="red")
legend("topright", c("Growth (kappa)", "Reproduction (1-kappa)"), lty=c(1,1), col = c("black" ,"red"),  cex=0.7)

par(mar = c(5,5,4,2) + 0.1) ## default is c(5,4,4,2) + 0.1
plot(x,kap_fun(x,kappa=0.83)*alpha*GR_pars[["Y"]]*eps1*x^eps2, #Imax, 
     type = "l", lty=2, xlab="Mass [g]", ylab=expression('gram day' ^-1),
     main= "Mass-scaling of rates and Growth energy at Tref", ylim = c(0,150))
lines(x,rho1*x^rho2, col= "red", type="l") #Maintenance
lines(x,eps1*x^eps2, col= "blue", type="l") #Imax, 
legend("topleft", c("Growth energy (Intake*kappa*alpha*Y)", "Maintenance rate","Max intake rate"), lty=c(2,1,1), col = c("black" ,"red", "blue"),  cex=0.7)

### Plot Temp functions 
plot(273:300, rI_T_GU2(273:300, 10), ylim=c(0,2), type = "l", xlab="Temperature [K]", ylab="Temp. effect", main= "Temperature effect on rates, Tref = 292 K")
lines(273:300, rM_T_AL2(273:300, 10), col="red")
  legend("topleft", c("Max Ingestion (Gardmark unimodal)", "Maintenance (Lindmark 2018)"),lty=c(1,1), col =c("black" ,"red"),  cex=0.7)

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

# Plot DEBoff.size ####
Toff <- NULL
xToff <- 1:300
for(i in seq(285,293,2)){
  Toff_pars <- c(T = i,     # parameters for Temperature, feeding, allocation and Mass dependence
                 kappa = 0.83, # allocation to respiration (Growth and maintenance)
                 Y=1)         # Feeding level
  Toff <- as.data.frame(rbind(Toff, c(Toff_pars, DEBoff.size(Toff_pars, y=xToff))))
}
colnames(Toff)[4:ncol(Toff)] <- xToff
Toff_long <- gather(Toff, key = "Size", value ="biom", c(4:ncol(Toff))) # not use column 

offT <-
  Toff_long %>%
  #filter(Size > 500) %>%
  ggplot(., aes(as.numeric(Size), biom, color = as.factor(T))) +
  geom_line(size=0.85) +
  #scale_colour_discrete() +
  scale_x_continuous("Weight [g]", limits = c(0,300),breaks = seq(0,300,50)) +
  annotate(geom="text", -Inf, Inf, label="D", hjust = -5, vjust = 2, size= 4, fontface = "bold") +
  ylab("Size dist. Age 1") +
  scale_colour_brewer(name="Temp [K]", palette="Dark2")+
  theme_bw()
offT

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
  #ifelse(m < mmin, sx.firstyear(m,Pars), sx(m,Pars)*starvx(m, Pars))# multiply with F (fishing mortality function)
  ifelse(m < mmin, sx.firstyear(m,Pars), Vsurvfun(m)*starvx(m, Pars))# multiply with F (fishing mortality function)
  #ifelse(m < mmin, sx.firstyear(m,Pars), VsurvfunT(m,Pars)*starvx(m, Pars))# multiply with F (fishing mortality function)
  #ifelse(m < mmin, sx.firstyear(m,Pars), 0.67*starvx(m, Pars))# multiply with F (fishing mortality function)
  #ifelse(m < mmin, sx.firstyear(m,Pars), 0.67)# multiply with F (fishing mortality function)
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
  ylab("Survival") +
  xlab("Weight [g]") +
  #annotate(geom="text", -Inf, Inf, label="C", hjust = -5, vjust = 3, size= 4, fontface = "bold") +
  ylim(0,1) +
  scale_colour_brewer(name="Temp [K]", palette="Dark2")+
  theme_bw()
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
  trajT <- as.data.frame(cbind(ProjCoh(c(T = i, kappa = kappp, Y = 1), t, I_pop=c(1000000,rep(0,n))),rep(c(0,round(x)),t),i,0))
  colnames(trajT) <- c("abund","t","mass","Temp","cIa")
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
  ylab("Weight [g]")+
  annotate(geom="text", -Inf, Inf, label="A", hjust = -5, vjust = 3, size= 4, fontface = "bold")+
  theme_bw()+ 
  scale_colour_brewer(palette="Dark2",name = "Temp [K]") #+
groT
#theme(legend.position="none") 

# Plot fecundity from DEB and in Lake Windermere
rep=NULL
for(i in 1:nrow(trajG)){
  if(trajG$t[i]==1){
    rep_tplus1 <- cbind(trajG$mass[i], DEBrepfun(e_m, c(T = trajG$Temp[i], kappa = kappp, Y = 1)),trajG$t[i],trajG$Temp[i])} 
  else{ rep_tplus1 <- cbind(trajG$mass[i], 
                            DEBrepfun(trajG$mass[i-1], c(T = trajG$Temp[i], kappa = kappp, Y = 1)),trajG$t[i],trajG$Temp[i])}
  
  colnames(rep_tplus1) <- c("mass","fec","t","Temp")
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
  geom_smooth(data=FData, aes(Weight, Eggs),method=lm, lwd=0.5, linetype="dashed")+  
  scale_colour_brewer(palette="Dark2", name = "Temp [K]")+
  annotate(geom="text", -Inf, Inf, label="B", hjust = -10, vjust = 3, size= 4, fontface = "bold") +
  ylab("Fecundity")+
  xlab("Weight [K]")+
  theme_bw()

pdf("DEBIPM_vitalrates_20210504_nopoints.pdf", width = 8, height = 6)
(groT + repT) / (surT + offT)
dev.off()


groT
par(mfrow=c(2,2), las=1, bty="l")
par(mfrow=c(1,1), las=1, bty="l")
groT
#VVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVV
#### ### ## # PLOT RESULTS #### ### ## # ## ### ####
#VVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVV
## READ in results from high res Main and Contrast results ##
# Results 1
Res_lam <- read.delim("//storage-og.slu.se/home$/vitl0001/My Documents/Manus2/Results/Results0427/Res_lam_1_0427.txt", sep = ",")
#Res_lam <- Res_lam[1:1149,]
Res_v <- read.delim("//storage-og.slu.se/home$/vitl0001/My Documents/Manus2/Results/Results0427/Res_v_1_0427.txt", sep = ",")
colnames(Res_v)[4:ncol(Res_v)] <- sub("X", "", colnames(Res_v)[4:ncol(Res_v)])
#Res_v <- Res_lam[1:1149,]
Res_w <- read.delim("//storage-og.slu.se/home$/vitl0001/My Documents/Manus2/Results/Results0427/Res_w_1_0427.txt", sep = ",")
colnames(Res_w)[4:ncol(Res_w)] <- sub("X", "", colnames(Res_w)[4:ncol(Res_w)])
#Res_w <- Res_lam[1:1149,]
# # Results 2 
# Res_lam <- read.delim("Res_lam_2_0406.txt", sep = ",")
# Res_v <- read.delim("Res_v_2_0406.txt", sep = ",")
# Res_w <- read.delim("Res_w_2_0406.txt", sep = ",")

### PLOT LAMBDA AS A FUNCTION OF KAPPA AND/OR TEMP ####
a_1 <- as_tibble(Res_lam) %>%
  filter(T %in% seq(min(Res_lam$T),max(Res_lam$T), length.out=13)) %>%
  ggplot(., aes(kappa, Lambda, color = as.factor(T))) +
  geom_line(size=0.5) +
  ylab(expression(lambda~(Fitness))) +
  xlab(expression(kappa~(Growth~Allocation))) +
  labs(color = "Temp [K]") +
  theme_bw()
a_1

### PLOT A SURFACE OF FITNESS OVER TEMP (x) AND KAPPA (y)
Res_lam_1 <- read.delim("//storage-og.slu.se/home$/vitl0001/My Documents/Manus2/Results/Results0427/Res_lam_1_0427.txt", sep = ",")
maxl_1 <- as.data.frame(Res_lam) %>% 
  group_by(T) %>%
  slice_max(Lambda)
maxl_1["Sur_f"] <- "Vin"
# as_tibble(maxl_1) %>%
#   ggplot(., aes(T, kappa))+
#   geom_line(size=0.85)+
#   scale_x_continuous(limits = c(283,293),breaks = seq(283,293,2))+
#   ylab(expression(paste(italic(Kappa[Optimum]))))+
#   xlab("Temp [K]")+
#   theme_bw()
b_1 <- 
  as_tibble(Res_lam) %>%
  filter(kappa <= 1) %>%
  ggplot(., aes(T,kappa)) +
  geom_raster(aes(fill=round(Lambda, 4)))+#, interpolate = TRUE) +
  #geom_smooth(data=maxl_1 ,aes(T,kappa), size=0.85, se=FALSE) +
  geom_line(data=maxl_1 ,aes(T,kappa),size=0.85) +
  scale_fill_gradientn(
    limits = c(0,max(as.data.frame(Res_lam[,'Lambda'])+0.0001)),
    colours = c("black","white","yellow","red"),
    values = rescale(c(0,0.9999,1, round(max(as.data.frame(Res_lam[,'Lambda'])), 4)),
                     to=c(0,1)),
    breaks = c(0,1,round(max(as.data.frame(Res_lam[,'Lambda'])), 2)),
    labels=c("0",1,round(max(as.data.frame(Res_lam[,'Lambda'])),3)),
    name = expression(lambda)) +
  scale_x_continuous(name = "Temperature [K]", limits = c(283,295),breaks = seq(283,295,2))+
  scale_y_continuous(expand = c(0,0), name=expression(paste(kappa," (Growth allocation)")))+
  theme_bw()
b_1

Res_lam_2 <- read.delim("//storage-og.slu.se/home$/vitl0001/My Documents/Manus2/Results/Results0427/Res_lam_2_0427.txt", sep = ",")
maxl_2 <- as.data.frame(Res_lam_2) %>% 
  group_by(T) %>%
  slice_max(Lambda)
maxl_2["Sur_f"] <- "tempVin"

Res_lam_3 <- read.delim("//storage-og.slu.se/home$/vitl0001/My Documents/Manus2/Results/Results0427/Res_lam_3_0427.txt", sep = ",")
maxl_3 <- as.data.frame(Res_lam_3) %>% 
  group_by(T) %>%
  slice_max(Lambda)
maxl_3["Sur_f"] <- "flat"

Res_lam_4 <- read.delim("//storage-og.slu.se/home$/vitl0001/My Documents/Manus2/Results/Results0507/Res_lam_1_0507.txt", sep = ",")
maxl_4 <- as.data.frame(Res_lam_4) %>% 
  group_by(T) %>%
  slice_max(Lambda)
maxl_4["Sur_f"] <- "Nostarv"

maxl_123 <- as_tibble(rbind(maxl_1,maxl_2,maxl_3))

b_123 <- 
maxl_123 %>%
  filter(kappa <= 1) %>%
  ggplot(., aes(T,kappa, linetype=Sur_f)) +
  #geom_raster(aes(fill=round(Lambda, 4)))+#, interpolate = TRUE) +
  #geom_smooth(data=maxl_123, aes(T,kappa), size=0.85, se=FALSE, color="black") +
  geom_line(data=maxl_123, aes(T,kappa),size=0.85) +
  scale_y_continuous(name = expression(paste(kappa," (Growth allocation)"))) +
  scale_x_continuous(name = "Temperature [K]", limits = c(283,295),breaks = seq(283,295,4))+
  scale_linetype_manual(values = c("longdash","dotdash","solid"),
                         name = "a(T,m)", 
                         labels = c("Constant, 0.67","s(T,m)","s(m)") 
                         )+
  ylim(0.65,1)+
  theme_bw()
b_123
b_1 + (b_123/b_poprep)

#pdf("DEBIPM_fig4_20210511.pdf", width = 10, height = 6)
(b_1 / b_123) | (vi/vT/co)
#dev.off()

Kopt <- 
  as_tibble(maxl_1) %>%
  filter(T %in% c(287,289,291)) #%>%

Res_w_long <- pivot_longer(as.data.frame(Res_w), cols = c(5:ncol(Res_w)), 
                           names_to = "Size", values_to ="biom") # not use column 4 & 5 (eggstage and recruits)

OpCow_287 <- Res_w_long %>%
  filter(T %in% Kopt$T[1]) %>%
  filter(kappa %in% c(Kopt$kappa[1],
                      Kopt$kappa[1]+0.01,0.94))
OpCow_289 <- Res_w_long %>%
  filter(T %in% Kopt$T[2]) %>%
  filter(kappa %in% c(Kopt$kappa[2],
                      Kopt$kappa[2]+0.01,0.93))
OpCow_291 <- Res_w_long %>%
  filter(T %in% Kopt$T[3]) %>%
  filter(kappa %in% c(Kopt$kappa[3],
                      Kopt$kappa[3]+0.01,91))
OptCow <- rbind(OpCow_287,OpCow_289,OpCow_291)

sum(Res_w[450,4:ncol(Res_w)])*dx

OptCow %>%
  ggplot(., aes(as.numeric(Size), biom, color = as.factor(kappa), linetype = as.factor(T))) +
  ggtitle("Stable structure for small Sizes, kappa=") + 
  geom_line(size=0.5) +
  ylab("Stable structure w") +
  xlab("Size") +
  xlim(0,1000) +
  theme_bw() +
  labs(color = "Temp [K]")

c_1a <- Res_w_long %>%
  filter(T %in% c(285, 287, 289, 291)) %>%
  filter(as.character(kappa) %in% 0.83) %>% #Floating point issue when comparing vector, therefore the use of %in%, can also use near()
  ggplot(., aes(as.numeric(Size), biom, color = as.factor(T))) +
  #ggtitle("Stable structure > Age 2, kappa = 0.83") +
  geom_line(size=0.8)+
  ylab("Stable structure w")+
  xlab("Weight [g]")+
  ylim(0,2e-8)+
  #xlim(300,NA)+
  annotate(geom="text", -Inf, Inf, label="A", hjust = -21, vjust = 2, size= 4, fontface = "bold") +
  #coord_cartesian(xlim = c(100,15000))+#, expand = c(50,), ) +
  #scale_x_continuous(breaks = c(100,seq(3000,15000,3000)))+
  scale_colour_brewer(palette="Dark2",name = "Temp [K]") +
  #geom_vline(xintercept = 1000, color="red") +
  theme_bw() +
  #theme(legend.position="none",
  #plot.title = element_text(size = 12))+
  labs(color = "Temp [K]")
c_1a

c_1b <- Res_w_long %>%
  filter(T %in% c(285, 287, 289, 291)) %>%
  filter(as.character(kappa) %in% 0.83) %>% #Floating point issue when comparing vector, therefore the use of %in%, can also use near()
  ggplot(., aes(as.numeric(Size), biom, color = as.factor(T))) +
  geom_line(size=0.8)+
  ylab("Stable structure w")+
  xlab("Weight [g]")+
  #ylim(0,1.1e-7)+
  #annotate(geom="text", -Inf, Inf, label="A", hjust = -21, vjust = 2, size= 4, fontface = "bold") +
  #coord_cartesian(xlim = c(0,1750))+#, expand = c(50,), ) +
  ylim(0,2.7e-6)+
  scale_x_continuous(breaks = c(0,500,1000), limits=c(0,1000))+
  scale_colour_brewer(palette="Dark2",name = "Temp [K]") +
  #geom_vline(xintercept = 1000, color="red") +
  theme_bw() +
  theme(legend.position="none",
        axis.title.x=element_blank(),
        axis.title.y=element_blank())
c_1b

mycolors <- rev(colorRampPalette(brewer.pal(9, "BuGn"))(10))
c_2 <- as_tibble(Res_w_long) %>%
  filter(kappa %in% 0.83) %>%
  mutate(biom_log=log10(biom)) %>%
  ggplot(., aes(T,as.double(Size), z=biom_log, colour=biom_log))+
  geom_contour_filled(breaks = c(max(log10(Res_w_long$biom)), #-12
                                 seq(-8,-15,length.out = 9),
                                 min(log10(Res_w_long$biom))))+ #-160
  # geom_contour_filled(breaks = c(max(log(Res_w_long$biom)), #-12
  #                                seq(-20,-50,length.out = 10),
  #                                min(log(Res_w_long$biom))))+ #-160
  ylab("Weight [g]")+
  annotate(geom="text", -Inf, Inf, label="B", hjust = -20, vjust = 2, size= 4, fontface = "bold")+
  scale_fill_manual(values = mycolors, name="Stable size distr.",
                    labels=c(max(Res_w_long$biom),
                             seq(1e-8,1e-15,length.out = 9),
                             min(Res_w_long$biom)))+
  scale_x_continuous(name="Temperature [K]",limits = c(283,293),breaks = seq(283,293,2))+
  theme_bw()
c_2

  maxbio<- as.data.frame(Res_w_long) %>% 
  filter(kappa %in% 0.83) %>%
  filter(as.double(Size) > 400 & as.double(Size) < 1500) %>%
  group_by(T) %>%
  slice_max(biom)

# OptCow %>%
#     filter(T %in% c(291, 292, 293, 294, 295)) %>%
#     #filter(kappa %in% 0.8) %>% #Floating point issue when comparing vector, therefore the use of %in%, can also use near()
#     filter(as.character(kappa) %in% 0.93) %>% #Floating point issue when comparing vector, therefore the use of %in%, can also use near()
#     ggplot(., aes(as.numeric(Size), biom, color = as.factor(kappa), linetype = as.factor(T))) +
#     ggtitle("Stable struct w/o Age 1, T=289") + 
#     geom_line(size=0.5) +
#     ylab("Stable structure w") +
#     xlab("Size") +
#     ylim(0,2.5e-8)+
#     xlim(0,3000)+
#     geom_vline(xintercept = 1000, color="red") +
#     theme_bw()+
#     theme(#legend.position="none", 
#       plot.title = element_text(size = 12)) +
#     labs(color = "Kappa")

  # Res_w_long %>%
  #   filter(kappa %in% c(0.75,0.8,0.85)) %>% #Floating point issue when comparing vector, therefore the use of near()
  #   filter(T %in% 283) %>%
  #   ggplot(., aes(as.numeric(Size), biom, color = as.factor(kappa))) +
  #   ggtitle("Stable structure over Size with three kappas at T0") +
  #   geom_line(size=0.8) +
  #   ylab("Stable structure w") +
  #   xlab("Size") +
  #   ylim(0,5e-9)+
  #   xlim(1000,NA)+
  # theme_bw()+
  #   labs(color = expression(kappa))
  
  
### PLOT REPRODUCTIVE VALUES V AS A FUNCTION OF KAPPA AND/OR TEMP ####
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
  scale_fill_manual(values = rev(mycolors), name="Reproductive value v)")+
  theme_bw()
d_2

b_poprep <-
as.tibble(Res_v_long) %>%
  group_by(T, kappa) %>%
  summarise(PopRep = sum(biom)) %>%#, groups = c(T, kappa)) 
  group_by(T) %>%
  slice_max(PopRep) %>%
  ggplot(., aes(T,kappa))+#, 
      geom_line()+
      scale_y_continuous(limits = c(0.7,1),name = expression(paste(kappa," (Growth allocation)")))+
      scale_x_continuous(name = "Temperature [K]", limits = c(283,295),breaks = seq(283,295,4))+
      theme_bw()


pdf("DEBIPM_w_v_1_20210512.pdf", width = 10, height = 6)
subvp <- viewport(width = 0.17, height = 0.20, x = 0.25, y = 0.82)
(c_1a+c_2) / (d_1+d_2)
print(c_1b, vp = subvp)
dev.off()
# grid.arrange(
# tableGrob(DEBparams, rownames(DEBparams)),
# tableGrob(IPMparams, rownames(IPMparams)), nrow=1
#)

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
    #ylim(0,20) +
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
  