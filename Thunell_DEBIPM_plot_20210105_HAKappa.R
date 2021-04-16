
### PLOT DEB IPM 2020 11 18 ####  

#setwd("C:/Users/vitl0001/VThunell_repos/Temperature-DEBIPM")
#source("Thunell_DEBIPM_model_20210105_HAkappa.R")

library(scales)
library(patchwork)

#Plot-parameter s for IPM vital rate functions needed for some plots
GR_pars <- c(T = 287,     # parameters for Temperature, feeding, allocation and Mass dependence
             kappa = 0.93, # allocation to respiration (Growth and maintenance)
             Y = 1)         # Feeding level

## PLOT MASS AND TEMPERATURE FUNCTIONS ####
### Plot mass functions
plot(x, kap_fun(x,kappa=0.93), yla = "kappa", type="l", ylim= c(0,1),main = "Size-dependent kappa")
lines(x, 1-kap_fun(x, kappa=0.93), ylab = "1-kappa", type="l", main = "Size-dependent 1-kappa", col="red")

par(mar = c(5,5,4,2) + 0.1) ## default is c(5,4,4,2) + 0.1
plot(x,kap_fun(x,kappa=0.93)*alpha*Y*eps1*x^eps2, #Imax, 
     type = "l", lty=2, xlab="Mass [g]", ylab=expression('gram day' ^-1),
     main= "Mass-scaling of rates and Growth energy at Tref", ylim = c(0,150))
lines(x,rho1*x^rho2, col="red") #Maintenance
lines(x,eps1*x^eps2, col= "blue") #Imax, 
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
                         growth_E = kap_fun(x,GR_pars["kappa"])*alpha*Y*intake - maint, 
                         repro_E = alpha*Y*(1-kap_fun(x,GR_pars["kappa"]))*intake))
  }

T_rates_long <- gather(as.data.frame(T_rates), key = "E_type", value ="g_day", c(5,6)) # gather only column 6,7

## Over mass for three temperatures
T_rates_long %>%
  filter( Temp %in% c(283,284,285)) %>%
  ggplot(.) +
  geom_line(size=1 , aes(mass, g_day, color = as.factor(Temp), 
                         linetype = as.factor(E_type))) +
  #ylim(0,10) +
  theme_bw() +
  scale_colour_brewer(palette="Dark2")

## Over temepratures for three sizes
T_rates_long %>%
  filter( mass %in% c(x[4], x[40], x[400])) %>%
  ggplot(.) +
  geom_line(size=1 , aes(Temp, g_day, color = as.factor(mass), 
                         linetype = as.factor(E_type))) +
  geom_hline(yintercept = 0) +
  ylim(-5,25) +
  theme_bw() +
  scale_colour_brewer(palette="Dark2")

### PLOT GROWTH / SIZE AT AGE DEPENDENT ON MASS, TEMPERATURE & kAPPA ####
# Build a dataframe of sizes, fecundity, temperatures and kappas using ratefun()

Temp = seq(285,289,2)   # temperatures
Kap  = 0.93 #c(0.75,0.8,0.85)  # kappas, allocation to growth (roughly 0.85 in Windermere pike)
Y = 1
parsP <- as.matrix(expand.grid(Temp,Kap,Y))
colnames(parsP) <- c("T","kappa","Y")
outR = NULL
for (i in 1:nrow(parsP)){
  plotpars <- parsP[i,]
  for(seas in 1:20) { # growth from age 1 happens in this for-loop, 25 years is from age 2 to 26
    if(seas == 1) {
        offmean <- ratefun(e_m, plotpars)[round(sl/2),2,]
        p <- data.frame(mass = ratefun(offmean, plotpars)[sl,2,],
                        re = DEBrepfun(offmean, plotpars), 
                        age = seas+1,
                        Temp = parsP[i,1],
                        kappa = parsP[i,2],
                        Y = parsP[i,3], row.names = seas) #+1 as seas=1 is to age 2, offmean represents age 1
        } else p <- rbind(p,
                      c(max(ratefun(p[nrow(p),1], plotpars)[,2,]),
                        DEBrepfun(p[nrow(p),1], plotpars, e_m), 
                        seas+1, parsP[i,1], parsP[i,2], parsP[i,3]))
  } 
  p <- rbind(c(offmean, DEBrepfun(e_m, plotpars), 
               1,parsP[i,1], parsP[i,2], parsP[i,3])
               , p ) # add Age 1
  outR <- bind_rows(outR,p)
}

# Plot mass over Age for the different Temps
outR %>%
  #filter(Temp %in% c(283, 285, 287, 289,291)) %>% #Floating point issue when comparing vector, therefore the use of %in%, can also use near()
  ggplot(.) +
  geom_point(data=GData, aes(Age, ltow_AYV(Length)), alpha=0.1) +
  geom_line(size=1,  aes(age, mass, color = as.factor(Temp))) +
  facet_wrap(vars(kappa), labeller = label_both, scales="free_y") +
  xlab("age") +
  ylab("mass [g]") +
 # xlim(0,15)+
  #coord_cartesian(ylim = c(0,16000)) +
  theme_bw() +  
  #theme(legend.position="none")+
  scale_colour_brewer(palette="Dark2")
  
plot(FData$Weight,FData$Eggs,ylab="#eggs in t+1", xlab = "weight in t+1",
     main = "fecundity in t+1", ylim = c(0,600000),xlim = c(0,18000)) # Windermere Pike data
abline(lm(FData$Eggs~FData$Weight), col="blue",lwd=1)      
#lines(ratefun(x,GR_pars)[183,2,],DEBrepfun(x,GR_pars),type = "l", col="blue")
outR %>%
  filter(Temp == 287) %>%
  lines(re~mass, data=., lwd = 2, type = "l", col = "red")
legend("topleft", c("DEB", "Windermere Pike", "lm(Windermere Pike)"), 
       lty=c(1, NA), pch =c(NA,1), col= c("red","black", "blue"),cex=.7)


#### PLOT RATE FUNCTION FOR IPM ####

### Plot growth fun ####
plot(x,  DEBgrowthfun(14000, y=x, GR_pars)*dx, type="l", ylab="DEBgrowth(y;x)", 
     xlab="y", main = "Prob. density of DEBgrowth(x,y)")
lines(x, DEBgrowthfun(10000, y=x, GR_pars)*dx, lty=2)
lines(x, DEBgrowthfun(1000, y=x, GR_pars)*dx, lty=3)
lines(x, DEBgrowthfun(100, y=x, GR_pars)*dx, lty=4) # values >~8100 are highly unlikley to grow into as y=x??
legend("topleft", c("14000","10000","1000","100"), lty=c(1,2,3,4), cex = 0.7, title="mass [g]")

# Plot growth against Windermere pikes ####
outR %>%
 filter(Temp == 287) %>%
 plot(mass~age, data=., lwd = 2, type = "l", col = "red", ylim = c(0,17000))
 points(GData2$Age, ltow_AYV(GData2$Length))
 legend("bottomright", c("DEB", "Windermere Pike"), col=c("red","black"),lty=c(1, NA), pch =c(NA,1,NA), cex = 0.75)

# Plot individual growth trajectories or size at age of Windermere pike
outR %>%
 filter(Temp == 287) %>%
 ggplot(., aes(age, mass)) +
 geom_line(data=GData, aes(Age, ltow_AYV(Length), group = as.factor(Ind))) +
 geom_line(size=1, color = "red") +
 theme_bw() +  
 scale_colour_brewer(palette="Dark2")

# Plot egg size distribution in Windermere data ####
# agg <- seq(min(FData$Egg.weight),max(FData$Egg.weight), 0.0001)
# plot(agg, dnorm(agg,mean(FData$Egg.weight),sd(FData$Egg.weight)),type="l", xlab = "egg sizes", ylab= "density")
# plot(agg, o_eyd, type="l")

# Plot DEBrepfun and windermere data ####
plot(FData$Weight,FData$Eggs,ylab="#eggs in t+1", xlab = "weight in t+1",
     main = "fecundity in t+1", ylim = c(0,600000),xlim = c(0,18000)) # Windermere Pike data
abline(lm(FData$Eggs~FData$Weight), col="blue",lwd=1)      
#lines(ratefun(x,GR_pars)[183,2,],DEBrepfun(x,GR_pars),type = "l", col="blue")
outR %>%
  filter(Temp == 287) %>%
  lines(re~mass, data=., lwd = 2, type = "l", col = "red")
  legend("topleft", c("DEB", "Windermere Pike", "lm(Windermere Pike)"), 
       lty=c(1, NA), pch =c(NA,1), col= c("red","black", "blue"),cex=.7)

outR %>%
  filter(Temp %in% c(283,285,287,289)) %>%
  ggplot(.) +
    geom_point(data=FData, aes(Weight, Eggs), alpha=0.1) +
    geom_line(size=1,  aes(mass, re, color = as.factor(Temp))) +
    xlim(0,20000)+
    theme_bw() +  
    scale_colour_brewer(palette="Dark2")

# Plot DEBoff.size ####
Toff <- NULL
xToff <- 1:300
for(i in seq(275,290,2)){
  Toff_pars <- c(T = i,     # parameters for Temperature, feeding, allocation and Mass dependence
               kappa = 0.8, # allocation to respiration (Growth and maintenance)
               Y=Y)         # Feeding level
  Toff <- as.data.frame(rbind(Toff, c(Toff_pars, DEBoff.size(Toff_pars, y = xToff))))
}
colnames(Toff)[4:ncol(Toff)] <- c(xToff)
Toff_long <- gather(Toff, key = "Size", value ="biom", c(10:ncol(Toff))) # not use column 10 & 11 (eggstage and recruits)
Toff_long %>%
#filter(near(kappa, 0.8)) %>%
ggplot(., aes(as.numeric(Size), biom, color = as.factor(T))) +
  geom_line(size=0.5) +
  ylab("Offspring size dist.") +
  xlab("Size") +
  xlim(0,300) +
  theme_bw()

# Plot DEBsurv ####
Tsur <- NULL
for(i in seq(282,295,2)){
  Tsur_pars <- c(T = i,     # parameters for Temperature, feeding, allocation and Mass dependence
                 kappa = Kap, # allocation to respiration (Growth and maintenance)
                 Y=Y)
  Tsur <- as.data.frame(rbind(Tsur, c(Tsur_pars,DEBsurvfun(x, Tsur_pars))))
}
colnames(Tsur)[4:ncol(Tsur)] <- c(x)
Tsur_long <- pivot_longer(Tsur, cols = c(4:ncol(Tsur)), names_to = "Size", values_to = "biom") # not use column 10 & 11 (eggstage and recruits)
Tsur_long %>%
  filter(near(kappa, 0.8)) %>%
  ggplot(., aes(as.numeric(Size), biom, color = as.factor(T))) +
  geom_line(size=0.5) +
  ylab("Survival") +
  xlab("Size") +
  xlim(0,16000) +
  theme_bw()

#VVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVV
#### ### ## # PLOT RESULTS #### ### ## # ## ### ####
#VVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVV
## READ in results from high res Main and Contrast results ##
# Results 1
Res_lam <- read.delim("//storage-og.slu.se/home$/vitl0001/My Documents/Manus2/Results/Results0414/Res_lam_flatsurv_0414.txt", sep = ",")
Res_lam<-Res_lam[1:642,]
Res_v <- read.delim("//storage-og.slu.se/home$/vitl0001/My Documents/Manus2/Results/Results0412/Res_v_YVmorttemp_0412.txt", sep = ",")
Res_v <- Res_v[1:220,]
colnames(Res_v)[4:ncol(Res_v)] <- sub("X", "", colnames(Res_v)[4:ncol(Res_v)])
Res_w <- read.delim("//storage-og.slu.se/home$/vitl0001/My Documents/Manus2/Results0408/Res_w_1_0408.txt", sep = ",")
colnames(Res_w)[4:ncol(Res_w)] <- sub("X", "", colnames(Res_w)[4:ncol(Res_w)])
# # Results 2 
# Res_lam <- read.delim("Res_lam_2_0406.txt", sep = ",")
# Res_v <- read.delim("Res_v_2_0406.txt", sep = ",")
# Res_w <- read.delim("Res_w_2_0406.txt", sep = ",")

### PLOT LAMBDA AS A FUNCTION OF KAPPA AND/OR TEMP ####
#pdf("\\\\storage-og.slu.se/home$/vitl0001/Desktop/Survival_Sizeexp0261N01_20210216.pdf", width = 8, height = 3.5)
#pdf("\\\\storage-og.slu.se/home$/vitl0001/Desktop/kappaTlamdba_20201120.pdf", width = 8, height = 3.5)
a <- as_tibble(Res_lam) %>%
  #filter(T %in% c(277, 279, 281, 283, 285, 287, 289)) %>%
  ggplot(., aes(kappa, Lambda, color = as.factor(T))) +
  geom_line(size=0.5) +
  ylab(expression(lambda~(Fitness))) +
  xlab(expression(kappa~(Growth~Allocation))) +
  labs(color = "Temp [K]") +
  theme_bw()
a

### PLOT A SURFACE OF FITNESS OVER TEMP (x) AND KAPPA (y)
maxl<- as.data.frame(Res_lam) %>% 
        group_by(T) %>%
        slice_max(Lambda)
# maxl <- maxl %>%
#   filter(T > 283) %>%
#   filter(kappa < 1)


b <- 
  as_tibble(Res_lam) %>%
#  filter(T > 283) %>%
 # filter(kappa < 1.025) %>%
   ggplot(., aes(T,kappa)) +
   geom_raster(aes(fill=round(Lambda, 4))) + #, interpolate = TRUE) +
   geom_smooth(data=maxl ,aes(T,kappa), size=1, se=FALSE) +
   geom_point(data=maxl ,aes(T,kappa),size=0.1) +
   # scale_fill_gradientn(
   #   limits = c(0,max(as.data.frame(Res_lam[,'Lambda'])+0.0001)),
   #   colours = c("black","white","yellow","red"),
   #   values = rescale(c(0,0.9999,1, round(max(as.data.frame(Res_lam[,'Lambda'])), 4)),
   #                    to=c(0,1)),
   #   breaks = c(0,1,round(max(as.data.frame(Res_lam[,'Lambda'])), 2)),
   #   labels=c("0",1,round(max(as.data.frame(Res_lam[,'Lambda'])),3)),
   #   name = expression(lambda)) +
  scale_x_continuous(breaks = seq(min(Res_lam[,"T"]), max(Res_lam[,"T"]), by = 1),expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  #ggtitle("kappa Size-dependent") +
  #ggtitle("Temperature effect independent of size ") +
  theme_bw()
b


  # breaks = c(0,1, round(max(as.data.frame(Res_lam[,'Lambda'])), 2)),
  # labels=c("0",1,round(max(as.data.frame(Res_lam[,'Lambda'])),3)),


## PLOT STABLE STRUCTURE W AS A FUNCTION OF KAPPA AND/OR TEMP ####
#colnames(Res_w) <- c("T","kappa","Y","0",x)
as.data.frame(Res_w) %>%
  filter(T %in% 283) %>%
  filter(kappa %in% 0.92) #%>%
#plot()  
Res_w_long <- pivot_longer(as.data.frame(Res_w), cols = c(5:ncol(Res_w)), 
                           names_to = "Size", values_to ="biom") # not use column 4 & 5 (eggstage and recruits)

c <- Res_w_long %>%
  #filter(T %in% c(277, 281, 283, 287, 290)) %>%
  #filter(as.character(T) %in% c(283,285,287,289, 291)) %>%
  filter(as.character(T) %in% c(285,287,289)) %>%
  #filter(kappa %in% 0.8) %>% #Floating point issue when comparing vector, therefore the use of %in%, can also use near()
  filter(as.character(kappa) %in% 0.95) #%>% #Floating point issue when comparing vector, therefore the use of %in%, can also use near()
  ggplot(., aes(as.numeric(Size), biom, color = as.factor(T))) +
  ggtitle("Stable structure for small Sizes, kappa=0.8") + 
  geom_line(size=0.5) +
  ylab("Stable structure w") +
  xlab("Size") +
  #xlim(0,2000) +
  theme_bw() +
  theme(legend.position="none",
  plot.title = element_text(size = 12))+
  labs(color = "Temp [K]")
c
d <- Res_w_long %>%
#  filter(T %in% c(283, 285, 287, 289, 290)) %>%
  #filter(kappa %in% 0.8) %>% #Floating point issue when comparing vector, therefore the use of %in%, can also use near()
  filter(as.character(kappa) %in% 0.95) %>% #Floating point issue when comparing vector, therefore the use of %in%, can also use near()
  filter(as.character(T) %in% c(285,287,289)) %>%
  #  filter(as.character(T) %in% c(283, 285, 287, 289, 290)) %>%
  #filter(kappa %in% 0.8) %>% #Floating point issue when comparing vector, therefore the use of %in%, can also use near()
  #filter(as.character(kappa) %in% "0.92") %>% #Floating point issue when comparing vector, therefore the use of %in%, can also use near()
  ggplot(., aes(as.numeric(Size), biom, color = as.factor(T))) +
  ggtitle("Stable struct w/o Age 1, kappa=0.8") + 
  geom_line(size=0.5) +
  ylab("Stable structure w") +
  xlab("Size") +
  ylim(0,5e-8)+
  theme_bw()+
  theme(legend.position="none", 
  plot.title = element_text(size = 12)) +
  labs(color = "Temp [K]")
d
# Res_w_long %>%
#   filter(kappa %in% c(0.75,0.8,0.85)) %>% #Floating point issue when comparing vector, therefore the use of near()
#   filter(T %in% 283) %>% 
#   ggplot(., aes(as.numeric(Size), biom, color = as.factor(kappa))) +
#   ggtitle("Stable structure over Size with three kappas at T0") + 
#   geom_line(size=0.8) +
#   ylab("Stable structure w") +
#   xlab("Size") +
#   ylim(0,2e-8)+
#   theme_bw()+
#   labs(color = expression(kappa))


### PLOT REPRODUCTIVE VALUES V AS A FUNCTION OF KAPPA AND/OR TEMP ####
colnames(Res_v) <- c("T","kappa","Y","0",x)
Res_v_long <- pivot_longer(as_tibble(Res_v), cols = c(5:ncol(Res_v)), 
                           names_to = "Size", values_to ="biom") # not use column 4 & 11 (eggstage and recruits)

e<-Res_v_long %>%
#  filter(T %in% c(277, 281, 283, 287, 290)) %>%
  #filter(kappa %in% 0.8) %>% #Floating point issue when comparing vector, therefore the use of %in%, can also use near()
  filter(as.character(T) %in% c(285, 287, 289)) %>%
  filter(as.character(kappa) %in% 0.85) %>% #Floating point issue when comparing vector, therefore the use of %in%, can also use near()
  ggplot(., aes(as.numeric(Size), biom, color = as.factor(T))) +#, linetype = as.factor(kappa))) +
  ggtitle("Repro. values over Size") + 
  geom_line(size=0.8) +
  ylab("Reproductive value V") +
  xlab("Size")+
  theme_bw() +
  labs(color = "Temp [K]")
e
library(gridExtra)
pdf("DEBIPM_SizeSurv_20210326_Robust.pdf", width = 8, height = 6)
(a+b) / (c+d+e)
# grid.arrange(
# tableGrob(DEBparams, rownames(DEBparams)),
# tableGrob(IPMparams, rownames(IPMparams)), nrow=1
#)
dev.off()


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
dev.off()


### COHORT PROJECTPONS #### 

T285 <- cbind(ProjCoh(c(T = 288.2, kappa = 0.93, Y = 1),t),rep(c(0,round(x)),t))
colnames(T285) <- c("abund","t","size")
T287 <- cbind(ProjCoh(c(T = 287, kappa = 0.93, Y = 1),t),rep(c(0,round(x)),t))
colnames(T287) <- c("abund","t","size")
T289 <- cbind(ProjCoh(c(T = 289, kappa = 0.93, Y = 1),t),rep(c(0,round(x)),t))
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



### Size dependent cost #####

r0 <- 100
r1 <- 0.0001
eps1 <- 0.51  # Intake allometric scalar 
eps2 <- 0.58   # Intake allometric exponent
kap_fun <- function(m, kappa=0.89, ha=3){kappa*exp(-m/(ha*max(x)))} #for size dep. kappa

GR_pars <- c(T = 287,     # parameters for Temperature, feeding, allocation and Mass dependence
             kappa = 0.81, # allocation to respiration (Growth and maintenance)
             Y = 1)         # Feeding level
outR<-NULL
for(seas in 1:20) { 
  if(seas == 1) {
    offmean <- ratefun(e_m, GR_pars)[round(sl/2),2,]
    p <- data.frame(mass = max(ratefun(offmean, GR_pars)[,2,]),
                    re = DEBrepfun(offmean, GR_pars, e_m) - r0*max(ratefun(offmean, GR_pars)[,2,])^r1, 
                    age = seas+1,
                    Temp = GR_pars["T"],
                    Kappa = GR_pars["Kappa"],
                    Y = GR_pars["Y"], row.names = seas) #+1 as seas=1 is to age 2, offmean represents age 1
  } else p <- rbind(p,
                    c(max(ratefun(p[nrow(p),1], GR_pars)[,2,]),
                      DEBrepfun(p[nrow(p),1], GR_pars, e_m) - r0*max(ratefun(p[nrow(p),1], GR_pars)[,2,])^r1, 
                      seas+1, GR_pars[1], GR_pars[2], GR_pars[3]))
} 
p <- rbind(c(offmean, DEBrepfun(e_m, GR_pars, e_m), 
             1,GR_pars[1], GR_pars[2], GR_pars[3])
           , p ) # add Age 1
outR <- bind_rows(outR,p)

plot(FData$Weight,FData$Eggs,ylab="#eggs in t+1", xlab = "weight in t+1",
     main = "fecundity of y in t+1", ylim = c(0,400000)) # Windermere Pike data
outR %>%
  lines(re~mass, data=., lwd = 2, type = "l", col = "red")
outR %>%
  plot(mass~age, data=., lwd = 2, type = "l", col = "red", ylim = c(0,16000))
points(GData2$Age, ltow_AYV(GData2$Length))

