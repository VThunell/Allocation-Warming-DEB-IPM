############### THUNELL ET AL. DEBIPM R-script - PLOT #########################

library(scales)
library(patchwork) 
library(grid)   # for viewport()
library(tidyverse)

YVPike <- load("PikeDataFiles.R")
#Colorscheme from blue to red scale_color_manual(values = c('#4575b4','#91bfdb','#e0f3f8','#fee090','#fc8d59','#d73027'))+

## Load in results for Main and Contrast results ####
# Results Baseline scenario
Res_lam_1 <- read.delim("//storage-og.slu.se/home$/vitl0001/My Documents/Manus2/Results/Results220413/Res_lam_1_2022-04-21.txt", sep = ",")
Res_v_1 <- read.delim("//storage-og.slu.se/home$/vitl0001/My Documents/Manus2/Results/Results220413/Res_v_1_2022-04-21.txt", sep = ",")
colnames(Res_v_1)[4:ncol(Res_v_1)] <- sub("X", "", colnames(Res_v_1)[4:ncol(Res_v_1)])
Res_w_1 <- read.delim("//storage-og.slu.se/home$/vitl0001/My Documents/Manus2/Results/Results220413/Res_w_1_2022-04-21.txt", sep = ",")
colnames(Res_w_1)[4:ncol(Res_w_1)] <- sub("X", "", colnames(Res_w_1)[4:ncol(Res_w_1)])

maxl_1 <- as.data.frame(Res_lam_1) %>% # get max lambda for each temperature 
  group_by(T) %>%
  slice_max(Lambda)
maxl_1["Sur_f"] <- "main_s"

# Results for size-, but not temperature-dependent survival
Res_lam_2 <- read.delim("//storage-og.slu.se/home$/vitl0001/My Documents/Manus2/Results/Results220413/Res_lam_2_2022-04-21.txt", sep = ",")
maxl_2 <- as.data.frame(Res_lam_2) %>% # get max lambda for each temperature
  group_by(T) %>%
  slice_max(Lambda)
maxl_2["Sur_f"] <- "tind_s"

# Results for constant survival (i.e. independent of size and temperature)
Res_lam_3 <- read.delim("//storage-og.slu.se/home$/vitl0001/My Documents/Manus2/Results/Results220413/Res_lam_3_2022-04-21.txt", sep = ",")
maxl_3 <- as.data.frame(Res_lam_3) %>% # get max lambda fr each temperature
  group_by(T) %>%
  slice_max(Lambda)
maxl_3["Sur_f"] <- "const_s"

maxl_123 <- as_tibble(rbind(maxl_1,maxl_2,maxl_3))
maxl_123$Sur_f <- factor(maxl_123$Sur_f, levels = c("main_s", "tind_s", "const_s"))

### PLOT DEMOGRAPHIC FUNCTIONS FOR IPM (AND FIG. 3A-D)####

### Plot individual growth, g(m_s+1;m_s,T) via growthfun() ####
plot(x,  growthfun(14000, y=x, test_Pars)*dx, type="l", ylab="growthfun(y;x)",
     xlab="y", main = "Prob. density of DEBgrowth(x,y)")
lines(x, growthfun(10000, y=x, test_Pars)*dx, lty=2)
lines(x, growthfun(1000, y=x, test_Pars)*dx, lty=3)
lines(x, growthfun(100, y=x, test_Pars)*dx, lty=4) # values >~8100 are highly unlikley to grow into as y=x??
legend("topright", c("14000","10000","1000","100"), lty=c(1,2,3,4), cex = 0.7, title="Mass x [g]")

kappp = 0.8
Temp <- seq(285,293,2)
t <- 25

traj <- data.frame() 
 for(i in 1:length(Temp)){
   traj[1,i] <- x[which.max(age1size(Pars = c(T = Temp[i], kappa = 0.8, Y = 1) ))]
    for(j in 2:t){ 
      traj[j,i] <- x[which.max(growthfun(traj[j-1,i],Pars = c(T = Temp[i], kappa = 0.8, Y = 1)))]
   }
 }
colnames(traj) <- Temp[1:5]
traj<-cbind(Age=1:nrow(traj),traj)
trajG <- pivot_longer(traj, cols=c(2:6),names_to = "Temp", values_to = "mass" )
trajG$Temp <- as.double(trajG$Temp)

groT <- ggplot(trajG) +
  geom_line(data = trajG, size=0.85, aes(Age, mass, color = as.factor(Temp))) +
  ylab(expression(paste("Mass ", italic("\u03BC")["g"]," [g]"))) +
  xlab("Age [years]") +
  xlim(0,20)+
  ylim(0,max(trajG$mass))+
  annotate(geom="text", -Inf, Inf, label="A", hjust = -2, vjust = 3, size= 4, fontface = "bold")+
  theme_bw() +
  scale_color_manual(values = c('#4575b4','#91bfdb','#fee090','#fc8d59','#d73027'), name="Temp [K]")+
  theme(legend.position="none")

#### Plot size dependent fecundity, f(m_s,T) via repfun() ####
rep=data.frame()
for(i in 2:ncol(traj)){ # the temps, from 2nd last column in traj
  for(j in 1:nrow(traj)){ # the sizes in traj, i.e. Age traj
    if(traj$Age[j]==1){
      rep_tplus1 <- cbind(traj[j,i], repfun(e_m, c(T = Temp[i-1], kappa = kappp, Y = 1)), Temp[i-1], kappp, "Y"=1)} 
    else{ rep_tplus1 <- cbind(traj[j,i], repfun(traj[j-1,i], c(T = Temp[i-1], kappa = kappp, Y = 1)), Temp[i-1], kappp, "Y"=1)}
    colnames(rep_tplus1) <- c("mass","fec","Temp", "Kappa", "Y")
    rep <- rbind(rep, rep_tplus1)
  }
}

repT <- ggplot() +
  geom_line(data=as_tibble(rep), size=0.85, aes(mass, fec, color = as.factor(Temp)))+
  scale_color_manual(values = c('#4575b4','#91bfdb','#fee090','#fc8d59','#d73027'), name="Temp [K]")+
  annotate(geom="text", -Inf, Inf, label="B", hjust = -2, vjust = 3, size= 4, fontface = "bold") +
  ylab(expression(paste(italic("f"),"(",italic("m"["s"]),italic(",T"),")"))) +
  xlab(expression(paste("Mass ",italic("m")[s]," [g]"))) +
  theme_bw()

#### Plot Age 1 size distribution, o(e_m,T) via age1size() ####
T_o <- NULL
xT_o <- 1:500
for(i in seq(285,293,2)){
  T_o_pars <- c(T = i, kappa = 0.8, Y=1)
  T_o <- as.data.frame(rbind(T_o, c(T_o_pars, age1size(T_o_pars, y=xT_o))))
}
colnames(T_o)[4:ncol(T_o)] <- xT_o
T_o_long <- gather(T_o, key = "Size", value ="biom", c(4:ncol(T_o))) # not use column 

oT <-
  T_o_long %>%
  ggplot(., aes(as.numeric(Size), biom, color = as.factor(T), linetype = as.factor(T))) +
  geom_line(size=0.85) +
  scale_x_continuous(expression(paste("Mass ",italic("m")[s]," [g]")), limits = c(0,300),breaks = seq(0,300,100)) +
  annotate(geom="text", -Inf, Inf, label="D", hjust = -4, vjust = 3, size= 4, fontface = "bold") +
  ylab(expression(paste(italic("o"),"(",italic("m"["s+1"]),italic(":T"),")"))) +
  scale_color_manual(values = c('#4575b4','#91bfdb','#fee090','#fc8d59','#d73027'), name="Temp [K]")+
  scale_linetype_manual(values = c('solid','solid','solid','dashed','solid'), guide =NULL )+
  theme_bw()

## Plot survival probability, DEBsurv() ####
Tsur <- NULL
for(i in seq(285,293,2)){
  Tsur_pars <- c(T = i,     # parameters for Temperature, feeding, allocation and Mass dependence
                 kappa = test_Pars[["kappa"]], # allocation to respiration (Growth and maintenance)
                 Y=test_Pars[["Y"]])
  Tsur <- as.data.frame(rbind(Tsur, c(Tsur_pars,survfun(x, Tsur_pars))))
}
colnames(Tsur)[4:ncol(Tsur)] <- c(x)
Tsur_long <- pivot_longer(Tsur, cols = c(4:ncol(Tsur)), names_to = "Size", values_to = "biom") # not use column 10 & 11 (eggstage and recruits)

surT <-
  Tsur_long %>%
  ggplot(., aes(as.numeric(Size), biom, color = as.factor(T))) +
  geom_line(size=0.85, linetype= "solid") +
  ylab(expression(paste(italic("a"),"(",italic("m"["s"]),italic(",T"),")"))) +
  xlab(expression(paste("Mass ",italic("m")[s]," [g]"))) +
  annotate(geom="text", -Inf, Inf, label="C", hjust = -2, vjust = 3, size= 4, fontface = "bold") +
  ylim(0,1) +
  scale_color_manual(values = c('#4575b4','#91bfdb','#fee090','#fc8d59','#d73027'), name="Temp [K]")+
  theme_bw()+
  theme(legend.position="none")

### FIG 3E,F - STABLE STRUCTURE W ####
Res_w_1_long <- pivot_longer(as.data.frame(Res_w_1), cols = c(4:ncol(Res_w_1)), 
                           names_to = "Size", values_to ="biom") # not use column 4 & 5 (eggstage and recruits)

Fig3e_1 <- Res_w_1_long %>%
  filter(T %in% c(285, 287, 289, 291, 293)) %>%
  filter(as.character(kappa) %in% 0.8) %>% #Floating point issue when comparing vector, therefore the use of %in%, can also use near()
  ggplot(., aes(as.numeric(Size), biom, group = rev(as.factor(T)), color = as.factor(T))) +
  geom_line(size=0.8)+
  ylab("Stable structure w") +
  xlab(expression(paste("Mass ",italic("m"["s"])," [g]"))) +
  ylim(0,5e-7)+
  xlim(0,10000)+
  annotate(geom="text", -Inf, Inf, label="E", hjust = -21, vjust = 3, size= 4, fontface = "bold") +
  scale_color_manual(values = c('#4575b4','#91bfdb','#fee090','#fc8d59','#d73027'), name="Temp [K]")+
  theme_bw() +
  labs(color = "Temp [K]") +
  theme(legend.position="none")

Fig3e_2 <- Res_w_1_long %>%
  filter(T %in% c(285, 287, 289, 291, 293)) %>%
  filter(as.character(kappa) %in% 0.8) %>% #Floating point issue when comparing vector, therefore the use of %in%, can also use near()
  ggplot(., aes(as.numeric(Size), biom, group = rev(as.factor(T)), color = as.factor(T))) +
  geom_line(size=0.8)+
  ylim(0,5e-5)+
  scale_x_continuous(breaks = c(0,200,400), limits=c(0,400))+
  scale_color_manual(values = c('#4575b4','#91bfdb','#fee090','#fc8d59','#d73027'), name="Temp [K]")+
  theme_bw() +
  theme(legend.position="none",
        axis.title.x=element_blank(),
        plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
        axis.title.y=element_blank())

mycolors <- colorRampPalette(colors=c("#feebe2","#ce1256"))(9)
Fig3f <- 
  as_tibble(Res_w_1_long) %>%
  filter(kappa %in% 0.8) %>%
  mutate(biom_log=log(biom)) %>% 
    
  ggplot(., aes(as.double(Size),T, z=biom_log, colour=biom_log))+
  geom_contour_filled(breaks = c(min(log(Res_w_1_long$biom)), 
                                 log(1e-10),log(1e-9),log(1e-8),log(5e-8),log(1e-7),log(5e-7),log(1e-6),
                                 max(log(Res_w_1_long$biom)))) + 
  annotate(geom="text", -Inf, Inf, label="F", hjust = -23, vjust = 3, size= 4, fontface = "bold")+
  scale_fill_manual(values = mycolors, name="Relative densities", 
                     labels=c(paste("<",1e-10), 
                              1e-9,1e-8,5e-8,1e-7,5e-7,1e-6,
                              paste("<",1))) + 
  scale_y_continuous(name = "Temperature [K]", limits = c(284,294),breaks = seq(285,293,2))+
  scale_x_continuous( name = expression(paste("Mass ",italic("m"["s"])," [g]")), limits = c(0,20000))+
  coord_cartesian(ylim = c(284,294)) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

### FIG 2 - SURFACE OF FITNESS OVER TEMP AND KAPPA ####

Fig2 <- 
  as_tibble(Res_lam_1) %>%
  filter(kappa > 0.65) %>%
  ggplot(., aes(T,kappa)) +
  geom_raster(aes(fill=round(Lambda, 4))) +
  geom_line(data=maxl_1,aes(T,kappa),size=0.85) +
  scale_fill_gradientn( #here Im choosing values for the scale that is within the kappa-temp space that Im plotting
    limits = c(min(as.data.frame(Res_lam_1[Res_lam_1$kappa>0.65 & Res_lam_1$T<293.5, ][,'Lambda'])),max(as.data.frame(Res_lam_1[,'Lambda']))),
    colours = c("black","white","#f1eef6","#ce1256"),
    values = rescale(c(min(as.data.frame(Res_lam_1[Res_lam_1$kappa>0.65 & Res_lam_1$T<293.5, ][,'Lambda'])),0.9999,1, max(as.data.frame(Res_lam_1[,'Lambda']))),
                     to=c(0,1)),
    breaks = c(min(as.data.frame(Res_lam_1[Res_lam_1$kappa>0.65 & Res_lam_1$T<293.5, ][,'Lambda'])),1,max(as.data.frame(Res_lam_1[,'Lambda']))),
    labels=c(round(min(as.data.frame(Res_lam_1[Res_lam_1$kappa>0.65 & Res_lam_1$T<293.5, ][,'Lambda'])), 2),
             1, round(max(as.data.frame(Res_lam_1[,'Lambda'])),2)),
    name = expression(lambda)) +
  scale_x_continuous(expand = c(0,0), name = "Temperature [K]", breaks = seq(284,292,2), limits = c(283,293.5)) +
  scale_y_continuous(expand = c(0,0), name=expression(paste(kappa[0]," (Growth allocation)"))) +
  coord_cartesian(xlim = c(283.5,292.5)) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

### FIG 4A - OPTIMAL KAPPA FOR THREE SURVIVAL SCENARIOS ####
Fig4A <- 
maxl_123 %>%
  ggplot(., aes(T,kappa, linetype = Sur_f)) +
  geom_line(data=maxl_123, aes(T,kappa), size=0.85) +
  scale_x_continuous(expand = c(0,0), name = "Temperature [K]", breaks = seq(284,292,2), limits = c(283,293.5)) +
  scale_y_continuous(expand = c(0,0), name=expression(paste(kappa[0]," (Growth allocation)")), limits = c(0.695,1.005))+
  scale_linetype_manual(values = c("solid","dotdash","dashed"),
                        name= "Survival",
                        labels = c(expression(paste(italic("a"),"(",italic("T,"),italic("m"["s"]),")")),
                                   expression(paste(italic("a"),"(",italic("m"["s"]),")")),
                                   "a=0.68")) +
  annotate(geom="text", -Inf, Inf, label="A", hjust = -2, vjust = 3, size= 4, fontface = "bold") +
  coord_cartesian(xlim = c(283.5,292.5)) +
  theme_bw() +
  theme( legend.key.width = unit(1, 'cm')) 


### FIG 4D,B AND C, THREE SURVIVAL SCENARIOS  ####
# Defnie survfun for plots, first Baseline model survival 4D
survfun <- function(m, Pars) { 
  Vsurvfun <- function(m, z=10.34){ 
    sxV <- function(m, z=10.34){
      1/(1+exp(13.53316 - 0.50977*wtol(m) - (-0.00393)
               *wtol(m)^2 - 0.19312*z - (-0.00679)
               *wtol(m)*z))
    }
    xmax <- x[which(sxV(x) == max(sxV(x)))]
    ifelse(m < xmax, sxV(m), sxV(xmax))
  }
  VsurvfunT <- function(m, Pars){ 
    zT=10.34+Pars[["T"]]-287 
    sxV <- function(m, z=zT){ 
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

Tsur <- NULL
for(i in seq(285,293,2)){
  Tsur_pars <- c(T = i,     # parameters for Temperature, feeding, allocation and Mass dependence
                 kappa = test_Pars[["kappa"]], # allocation to respiration (Growth and maintenance)
                 Y=test_Pars[["Y"]])
  Tsur <- as.data.frame(rbind(Tsur, c(Tsur_pars,survfun(x, Tsur_pars))))
}
colnames(Tsur)[4:ncol(Tsur)] <- c(x)
Tsur_long <- pivot_longer(Tsur, cols = c(4:ncol(Tsur)), names_to = "Size", values_to = "biom") # not use column 10 & 11 (eggstage and recruits)

Fig4D <-
  Tsur_long %>%
  ggplot(., aes(as.numeric(Size), biom, color = as.factor(T))) +
  geom_line(size=0.85, linetype= "solid") +
  ylab(expression(paste(italic("a"),"(",italic("m"["s"]),italic(", T"),")"))) +
  xlab(expression(paste("Mass ",italic("m"["s"])," [g]"))) +
  annotate(geom="text", -Inf, Inf, label="D", hjust = -1.5, vjust = 2.5, size= 4, fontface = "bold") +
  ylim(0,1) +
  scale_color_manual(values = c('#4575b4','#91bfdb','#fee090','#fc8d59','#d73027'), name="Temp [K]")+
  theme_bw() +
  theme( legend.key.width = unit(1, 'cm') ) 

# Define survfun for Temperature independent version of Vindenes et al. 2014
survfun <- function(m, Pars) { # Mass-Temp dependence of yearly Mortality 
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
    zT=10.34+Pars[["T"]]-287 # make 287 equal to 283 survival in Vindenes 2014 since 283 is yearly mean Temp 287 is summer mean temp
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

Tsur <- NULL
for(i in seq(285,293,2)){
  Tsur_pars <- c(T = i,     # parameters for Temperature, feeding, allocation and Mass dependence
                 kappa = test_Pars[["kappa"]], # allocation to respiration (Growth and maintenance)
                 Y=test_Pars[["Y"]])
  Tsur <- as.data.frame(rbind(Tsur, c(Tsur_pars,survfun(x, Tsur_pars))))
}
colnames(Tsur)[4:ncol(Tsur)] <- c(x)

Tsur_long <- pivot_longer(Tsur, cols = c(4:ncol(Tsur)), names_to = "Size", values_to = "biom") # not use column 10 & 11 (eggstage and recruits)
Fig4B <-
  Tsur_long %>%
  ggplot(., aes(as.numeric(Size), biom, color = as.factor(T))) +
  geom_line(size=0.85, linetype= "dotdash") +
  ylab(expression(paste(italic("a"),"(",italic("m"["s"]),")"))) +
  xlab(expression(paste("Mass ",italic("m"["s"])," [g]"))) +
  annotate(geom="text", -Inf, Inf, label="B", hjust = -1.5, vjust = 2.5, size= 4, fontface = "bold") +
  ylim(0,1) +
  scale_colour_manual(name="Temp [K]", values=c("black","black","black","black","black")) +
  theme_bw() +
  theme( legend.position="none",
    legend.key.width = unit(1, 'cm') )

# Constant survival (0.68)
survfun <- function(m, Pars) { # Mass-Temp dependence of yearly Mortality 
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
    zT=10.34+Pars[["T"]]-287 # make 287 equal to 283 survival in Vindenes 2014 since 283 is yearly mean Temp 287 is summer mean temp
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

Tsur <- NULL
for(i in seq(285,293,2)){
  Tsur_pars <- c(T = i, kappa = test_Pars[["kappa"]], Y=test_Pars[["Y"]])
  Tsur <- as.data.frame(rbind(Tsur, c(Tsur_pars,survfun(x, Tsur_pars))))
}
colnames(Tsur)[4:ncol(Tsur)] <- c(x)
Tsur_long <- pivot_longer(Tsur, cols = c(4:ncol(Tsur)), names_to = "Size", values_to = "biom") # not use column 10 & 11 (eggstage and recruits)

Fig4C <-
  Tsur_long %>%
  ggplot(., aes(as.numeric(Size), biom, color = as.factor(T))) +
  geom_line(size=0.85, linetype= "dashed") +
  ylab("a=0.68") +
  xlab(expression(paste("Mass ",italic("m"["s"])," [g]"))) +
  annotate(geom="text", -Inf, Inf, label="C", hjust = -1.5, vjust = 2.5, size= 4, fontface = "bold") +
  ylim(0,1) +
  scale_colour_manual(name="Temp [K]", values=c("black","black","black","black","black")) +
  theme_bw() +
  theme( legend.position="none",
    legend.key.width = unit(1, 'cm')) 

# define main model survfun again for use in further calculations
survfun <- function(m, Pars) { # Mass-Temp dependence of yearly Mortality 
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
    zT=10.34+Pars[["T"]]-287 # make 287 equal to 283 survival in Vindenes 2014 since 283 is yearly mean Temp 287 is summer mean temp
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

### FIG. 5 - CACLULATE SENSITIVITES AND PLOT ####
# Optimum kappa dependent on temperature:
OptPars <- 
  as.data.frame(maxl_1[,1:3]) %>%
  filter(T %in% c(287, 289, 291))

#Calculate sensitivities based on methods in supplement S4
Res_Sens <-NULL
for (i in 1:nrow(OptPars)){
  res_v <- as.numeric(Res_v_1 %>% filter(T == OptPars[i,1] & kappa==OptPars[i,2]))[4:ncol(Res_v_1)]
  res_w <- as.numeric(Res_w_1 %>% filter(T == OptPars[i,1] & kappa==OptPars[i,2]))[4:ncol(Res_w_1)]
  
  Gr_cont <- bind_cols(Size=x, Temp=OptPars[i,1], kappa=OptPars[i,2], 
                       dfun="Gr", outer(res_v[2:(n+1)], res_w[2:(n+1)], "*") 
                       * t(survfun(x, OptPars[i,]) * t(G_NumDer(h=0.00001, x, OptPars[i,])))*dx)
  Gr_cont["Sens"] <- apply(Gr_cont[5:ncol(Gr_cont)],1,FUN=sum)# sum the distribution y to get sensitivity for each size 
  Gr_cont<-Gr_cont[,c(1:4,ncol(Gr_cont))]
  
  F_cont  <- bind_cols(Size=x, Temp=OptPars[i,1], kappa=OptPars[i,2],
                       dfun="F",  Sens=res_v[1]*res_w[2:(n+1)] 
                       * survfun(x, OptPars[i,]) * F_NumDer(h=0.00001, x, OptPars[i,]))
  
  A1_cont <- bind_cols(Size=x, Temp=OptPars[i,1], kappa=OptPars[i,2], 
                       dfun="A1", Sens=res_v[2:(n+1)]*res_w[1] 
                       * el_surv * A1_NumDer(h=0.00001, x, OptPars[i,])*dx)
  
  Res_Sens <- rbind(rbind(Gr_cont,F_cont,A1_cont), Res_Sens) 
}

# Labels for sensitivities in plot
cont.labs <- c("Age 1 size", "Fecundity", "Growth")
names(cont.labs) <- c("A1", "F", "Gr") 

# Plot sensitivity contributions (figure 5 main text)
# Summarize sensitivity contributions for each demographic function
s_sum <- Res_Sens %>%
  group_by(dfun, Temp) %>%
  summarise(sum.c = sum(Sens, na.rm = TRUE)) #%>%

s_sum_text <- data.frame(Size = c(2000,2000,2000,4500,4500,4500,7000,7000,7000), Sens = 0.075,
                      dfun = s_sum$dfun, Temp = s_sum$Temp, label = round(s_sum$sum.c,3))

# create Fig 5
Fig5 <-
  Res_Sens %>%
  mutate(dfun=as.factor(dfun)) %>%
  ggplot(.) +
  geom_path(aes(Size, Sens, color = dfun)) +
  facet_grid(Temp~., scales="fixed") +
  xlim(0,20000) +
  geom_vline(xintercept=400, linetype="dashed", size=0.3) +
  ylab("Sensitivity contribution") +
  xlab("Mass [g]") +
  scale_color_manual(values = c('#E69F00','#000000','#ce1256'), labels= cont.labs, name= "")+
  geom_text(data=s_sum_text, aes(Size, Sens, label = label, 
                              colour=dfun), size=4, show.legend = FALSE) +
  theme_bw() +
  theme(panel.background = element_rect(fill = "white", colour = "black"),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA))

#### Write plots to pdf:s ####
 date <- Sys.Date()
# pdf(paste("DEBIPM_fig2_",date,".pdf", sep=""), width = 6, height = 4)
# Fig2
# dev.off()
# 
# pdf(paste("DEBIPM_fig3_",date,".pdf", sep=""), width = 8, height = 8)
# subvp <- viewport(width = 0.19, height = 0.18, x = 0.26, y = 0.20)
# (groT + repT) / (surT + oT) / (Fig3e_1+Fig3f) #/ (d_1+d_2)
# print(Fig3e_2, vp = subvp)
# dev.off()
# 
# pdf(paste("DEBIPM_fig4_",date,".pdf", sep=""), width = 10, height = 6)
#  Fig4A / (Fig4B+Fig4C+Fig4D)
#  dev.off()
# 
# pdf(paste("DEBIPM_fig5_",date,".pdf", sep=""), width = 8, height = 6)
# Fig5
# dev.off()