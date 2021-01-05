
### PLOT DEB IPM 2020 11 18 ####

setwd("C:/Users/vitl0001/VThunell_repos/Temperature-DEBIPM")
#source("Thunell_DEBIPM_model_20210105_HAKappa.R")

library(scales)

plot(x, 0.85*exp(-x/(ha*max(x))), ylab = "Kappa", type="l", main = "Hyperallometric scaling")
plot(x, 1-0.85*exp(-x/(ha*max(x))), ylab = "Kappa", type="l", main = "Hyperallometric scaling")

## PLOT MASS AND TEMPERATURE FUNCTIONS ####
### Plot mass functions fo Kappa=0.8
Kappa <- 0.87
Y <- 1
par(mar = c(5,5,4,2) + 0.1) ## default is c(5,4,4,2) + 0.1
plot(x,Kappa*exp(-x/(ha*max(x)))*alpha*Y*eps1*x^eps2, #Imax, 
     type = "l", lty=2, xlab="Mass [g]", ylab=expression('gram day' ^-1), main= "Mass-scaling of rates and Growth energy at Tref", ylim = c(0,100))
lines(x,rho1*x^rho2, col="red") #Maintenance
lines(x,eps1*x^eps2, col= "blue") #Imax, 
legend("topleft", c("Growth energy (Intake*Kappa*alpha*Y)", "Maintenance rate","Max intake rate"), lty=c(2,1,1), col = c("black" ,"red", "blue"),  cex=0.7)

### Plot Temp functions 
plot(280:294, rI_T_GU2(280:294, 10), ylim=c(0,3), type = "l", xlab="Temperature [K]", ylab="Temp. effect", main= "Temperature effect on rates, Tref = 283 K")
lines(280:294, rM_T_AL2(280:294, 10), col="red")
legend("topleft", c("Max Ingestion (Gardmark unimodal)", "Maintenance (Lindmark 2018)"),lty=c(1,1), col =c("black" ,"red"),  cex=0.7)

##### PLOT ENERGY BUDGET OVER TEMP AND MASS ####
Temp <- 280:294
T_rates = NULL
for(i in Temp) {
  maint <- rho1*(x^rho2)*rM_T_AL2(i, x)
  intake  <- eps1*(x^eps2)*rI_T_GU2(i, x) 
  T_rates  <- rbind(T_rates,
                   cbind(Temp = i, mass = x, 
                         maint, intake, 
                         growth_E = Kappa*exp(-x/(ha*max(x)))*alpha*Y*intake - maint, 
                         repro_E = alpha*Y*(1-Kappa*exp(-x/(ha*max(x))))*intake))
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
  filter( mass %in% c(x[2], x[10], x[100])) %>%
  ggplot(.) +
  geom_line(size=1 , aes(Temp, g_day, color = as.factor(mass), 
                         linetype = as.factor(E_type))) +
  geom_hline(yintercept = 0) +
  ylim(-5,25) +
  theme_bw() +
  scale_colour_brewer(palette="Dark2")


#### PLOT GROWTH / SIZE AT AGE DEPENDENT ON MASS, TEMPERATURE & kAPPA ####

# Build a dataframe of sizes, fecundity, temperatures and kappas using ratefun()

Temp = 282:285   # temperatures
Kap  = Kappa #c(0.75,0.8,0.85)  # kappas, allocation to growth (roughly 0.85 in Windermere pike)
Y = 1
parsP <- as.matrix(expand.grid(Temp,Kap,Y))
colnames(parsP) <- c("T","Kappa","Y")
outR=NULL

for (i in 1:nrow(parsP)){
  plotpars <- parsP[i,]
  for(seas in 1:20) { # growth from age 1 happens in this for-loop, 25 years is from age 2 to 26
    if(seas == 1) {
        offmean <- ratefun(o_e, plotpars)[sl/2,2,]
        p <- data.frame(mass = max(ratefun(offmean, plotpars)[,2,]),
                        re = DEBrepfun(offmean, plotpars, o_e), 
                        age = seas+1,
                        Temp = parsP[i,1],
                        Kappa = parsP[i,2],
                        Y = parsP[i,3], row.names = seas) #+1 as seas=1 is to age 2, offmean represents age 1
        } else p <- rbind(p,
                      c(max(ratefun(p[nrow(p),1], plotpars)[,2,]),
                        DEBrepfun(p[nrow(p),1], plotpars, o_e), 
                        seas+1, parsP[i,1], parsP[i,2], parsP[i,3]))
  } 
  p <- rbind(c(offmean, DEBrepfun(o_e, plotpars, o_e), 
               1,parsP[i,1], parsP[i,2], parsP[i,3])
               , p ) # add Age 1
  outR <- bind_rows(outR,p)
}
  
# Plot mass over Age for the different Temps
ggplot(outR) +
  geom_point(data=GData, aes(Age, ltow_A(Length)), alpha=0.1) +
  geom_line(size=1.3,  aes(age, mass, color = as.factor(Temp))) +
  facet_wrap(vars(Kappa), labeller = label_both, scales="free_y") +
  xlab("age") +
  ylab("mass [g]") +
  theme_bw() +  
  scale_colour_brewer(palette="Dark2")

#### PLOT RATE FUNCTION FOR IPM ####

#Plot-parameters for IPM vital rate functions
GR_pars <- c(T = 283,     # parameters for Temperature, feeding, allocation and Mass dependence
             Kappa = 0.8, # allocation to respiration (Growth and maintenance)
             Y = 1)         # Feeding level

### Plot growth fun ####
plot(x,  DEBgrowthfun(14000, y=x, GR_pars)*dx, type="l", ylab="DEBgrowth(y;x)", 
     xlab="y", main = "Prob. density of DEBgrowth(x,y)")
lines(x, DEBgrowthfun(10000, y=x, GR_pars)*dx, lty=2)
lines(x, DEBgrowthfun(1000, y=x, GR_pars)*dx, lty=3)
lines(x, DEBgrowthfun(100, y=x, GR_pars)*dx, lty=4) # values >~8100 are highly unlikley to grow into as y=x??
legend("topleft", c("14000","10000","1000","100"), lty=c(1,2,3,4), cex = 0.7, title="mass [g]")

# Plot growth against Windermere pikes ####
outR %>%
 filter(Temp == 283) %>%
 plot(mass~age, data=., lwd = 2, type = "l", col = "red", ylim = c(0,13000))
 points(GData2$Age, ltow_A(GData2$Length))
 legend("bottomright", c("DEB", "Windermere Pike"), col=c("red","black"),lty=c(1, NA), pch =c(NA,1,NA), cex = 0.75)

# Plot individual growth trajectories or size at age of Windermere pike
outR %>%
 filter(Temp == 283) %>%
 ggplot(., aes(age, mass)) +
 geom_line(data=GData, aes(Age, ltow_A(Length), group = as.factor(Ind))) +
 geom_line(size=1, color = "red") +
 theme_bw() +  
 scale_colour_brewer(palette="Dark2")

# Plot egg size distribution in Windermere data ####
# agg <- seq(min(FData$Egg.weight),max(FData$Egg.weight), 0.0001)
# plot(agg, dnorm(agg,mean(FData$Egg.weight),sd(FData$Egg.weight)),type="l", xlab = "egg sizes", ylab= "density")
# plot(agg, o_eyd, type="l")

# Plot DEBrepfun and windermere data ####
plot(FData$Weight,FData$Eggs,ylab="#eggs in t+1", xlab = "weight in t+1",
     main = "fecundity of y in t+1", ylim = c(0,400000)) # Windermere Pike data
outR %>%
  filter(Temp == 283) %>%
  lines(re~mass, data=., lwd = 2, type = "l", col = "red")
  legend("topleft", c("DEB", "Windermere Pike"), lty=c(1, NA), pch =c(NA,1), col= c("red","black"),cex=.7)
  
# Plot DEBoff.size ####
Toff <- NULL
xToff <- 1:200
for(i in seq(282,290,2)){
  Toff_pars <- c(T = i,     # parameters for Temperature, feeding, allocation and Mass dependence
               Kappa = Kap, # allocation to respiration (Growth and maintenance)
               Y=Y)         # Feeding level
  Toff <- as.data.frame(rbind(Toff, c(Toff_pars, DEBoff.size(xToff, Toff_pars))))
}
colnames(Toff)[4:ncol(Toff)] <- c(xToff)
Toff_long <- gather(Toff, key = "Size", value ="biom", c(10:ncol(Toff))) # not use column 10 & 11 (eggstage and recruits)
Toff_long %>%
#filter(near(Kappa, 0.8)) %>%
ggplot(., aes(as.numeric(Size), biom, color = as.factor(T))) +
  geom_line(size=0.5) +
  ylab("Offspring size dist.") +
  xlab("Size") +
  xlim(0,200) +
  theme_bw()
  
# Plot DEBsurv ####
Tsur <- NULL
for(i in seq(282,295,2)){
  Tsur_pars <- c(T = i,     # parameters for Temperature, feeding, allocation and Mass dependence
                 Kappa = Kap, # allocation to respiration (Growth and maintenance)
                 Y=Y)
  Tsur <- as.data.frame(rbind(Tsur, c(Tsur_pars,DEBsurvfun(x, Tsur_pars))))
}
colnames(Tsur)[4:ncol(Tsur)] <- c(x)
Tsur_long <- pivot_longer(Tsur, cols = c(4:ncol(Tsur)), names_to = "Size", values_to = "biom") # not use column 10 & 11 (eggstage and recruits)
Tsur_long %>%
  filter(near(Kappa, 0.8)) %>%
  ggplot(., aes(as.numeric(Size), biom, color = as.factor(T))) +
  geom_line(size=0.5) +
  ylab("Survival") +
  xlab("Size") +
  xlim(0,14000) +
  theme_bw()

#VVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVV
#### ### ## # PLOT RESULTS #### ### ## # ## ### ####
#VVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVV

### PLOT LAMBDA AS A FUNCTION OF KAPPA AND/OR TEMP ####
#pdf("\\\\storage-og.slu.se/home$/vitl0001/Desktop/kappaTlamdba_20201120.pdf", width = 8, height = 3.5)
ggplot(as.data.frame(Res_lam), aes(Kappa, Lambda, color = as.factor(T))) +
  geom_line(size=0.5) +
  ylab(expression(lambda~(Fitness))) +
  xlab(expression(Kappa~(Growth~Allocation))) +
  labs(color = "Temp [K]") +
  theme_bw()

### PLOT KAPPA at max(LAMBDA) (y) AS A FUNCTION TEMP (x) ####
as.data.frame(Res_lam) %>% 
  mutate(lam.1 = ifelse(Lambda < 1, "extin", "grow")) %>% 
  group_by(T) %>%
  slice_max(Lambda) %>%
  ggplot(., aes(T, Kappa, color = lam.1)) +
  ggtitle("Kappa value that gives highest lambda over Temp") + 
  geom_line(size=0.5) +
  geom_vline(xintercept = 283, colour="blue") +
  geom_text(aes(x=283.5, label="T0", y=0.78), colour="blue") +
  ylab(expression(Kappa~at~max~lambda)) +
  scale_color_manual(values=c("black","red")) +
  xlab("Temperature") +
  labs(color = "") +
  theme_bw()

### PLOT A SURFACE OF FITNESS OVER TEMP (x) AND KAPPA (y)
ggplot(as.data.frame(Res_lam), aes(T,Kappa)) +
  geom_raster(aes(fill=round(Lambda, 4))) + #, interpolate = TRUE) +
  scale_fill_gradientn(
    colours = c("black","white","yellow","red"),
    values = rescale(c(0,0.9999,1, round(max(as.data.frame(Res_lam[,'Lambda'])), 4)), to=c(0,1)),
    breaks = c(0,1, round(max(as.data.frame(Res_lam[,'Lambda'])), 2)),
    labels=c("0",1,round(max(as.data.frame(Res_lam[,'Lambda'])),3)),
    name = expression(lambda)) +
  theme_bw()

## PLOT STABLE STRUCTURE W AS A FUNCTION OF KAPPA AND/OR TEMP ####
Res_w_long <- pivot_longer(as_tibble(Res_w), cols = c(6:ncol(Res_w)), 
                           names_to = "Size", values_to ="biom") # not use column 4 & 5 (eggstage and recruits)
Res_w_long %>%
  filter(T %in% c(282,283)) %>% #Floating point issue when comparing vector, therefore the use of %in%, can also use near()
  filter(Kappa %in% 0.8) %>% #Floating point issue when comparing vector, therefore the use of %in%, can also use near()
  ggplot(., aes(as.numeric(Size), biom, color = as.factor(T))) +
  ggtitle("Stable structure over Size with temperatures") + 
  geom_line(size=0.5) +
  ylab("Stable structure w") +
  xlab("Size") +
  xlim(0,14000) +
  theme_bw()+
  labs(color = "Temp [K]")

Res_w_long %>%
  filter(Kappa %in% c(0.7,0.8)) %>% #Floating point issue when comparing vector, therefore the use of near()
  filter(T %in% 283) %>% 
  ggplot(., aes(as.numeric(Size), biom, color = as.factor(Kappa))) +
  ggtitle("Stable structure over Size with three Kappas") + 
  geom_line(size=0.5) +
  ylab("Stable structure w") +
  xlab("Size") +
  xlim(0,14000) +
  theme_bw()+
  labs(color = expression(kappa))

### PLOT REPRODUCTIVE VALUES V AS A FUNCTION OF KAPPA AND/OR TEMP ####
Res_v_long <- pivot_longer(as_tibble(Res_v), cols = c(5:ncol(Res_v)), 
                           names_to = "Size", values_to ="biom") # not use column 4 & 11 (eggstage and recruits)
Res_v_long %>%
  filter(Kappa %in% 0.8) %>% #Floating point issue when comparing vector, therefore the use of %in%, can also use near()
  ggplot(., aes(as.numeric(Size), biom, color = as.factor(T), linetype = as.factor(Kappa))) +
  ggtitle("Repro. values over Size with temperatures") + 
  geom_line(size=0.5) +
  ylab("Reproductive value V") +
  xlab("Size")+
  xlim(0,14000) +
  theme_bw() +
  labs(color = "Temp [K]")

Res_v_long %>%
  filter(Kappa %in% c(0.7,0.8,0.9)) %>%
  filter(T %in% c(283))  %>% 
  ggplot(., aes(as.numeric(Size), biom, color = as.factor(Kappa))) +
  ggtitle("Repro. values over Size with three Kappas") + 
  geom_line(size=0.5) +
  ylab("Reproductive value V") +
  xlab("Size")+
  xlim(0,14000) +
  theme_bw() +
  labs(color = expression(kappa))
#dev.off()
