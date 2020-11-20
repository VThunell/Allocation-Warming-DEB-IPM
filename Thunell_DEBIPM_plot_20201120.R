
### PLOT DEB IPM 2020 11 18 ####

setwd("~/Manus2/R/Manus2R")
source("Thunell_DEBIPM_model_20201116.R")

## PLOT MASS AND TEMPERATURE FUNCTIONS ####

### Plot mass functions 
par(mar = c(5,5,4,2) + 0.1) ## default is c(5,4,4,2) + 0.1
plot(x,Kap*alpha*Y*eps1*x^eps2, #Imax, 
     type = "l", lty=2, xlab="Mass [g]", ylab=expression('gram day' ^-1), main= "Mass-scaling of rates and Growth energy at Tref", ylim = c(0,100))
lines(x,rho1*x^rho2, col="red") #Maintenance
lines(x,eps1*x^eps2, col= "blue") #Imax, 
legend("topleft", c("Growth energy (Intake*IM_r*Kap*alpha*Y)", "Maintenance rate","Max intake rate"), lty=c(2,1,1), col = c("black" ,"red", "blue"),  cex=0.7)

### Plot Temp functions 
plot(280:290, rI_T_GU2(280:290, 10, T_par), ylim=c(0,2), type = "l", xlab="Temperature [K]", ylab="Temp. effect", main= "Temperature effect on rates, Tref = 283 K")
lines(280:290, rM_T_AL2(280:290, 10, T_par), col="red")
legend("bottomleft", c("Max Ingestion (Gardmark unimodal)", "Maintenance (Lindmark 2018)"),lty=c(1,1), col =c("black" ,"red"),  cex=0.75)

##### PLOT ENERGY BUDGET OVER TEMP AND MASS ####
Temp <- 280:290
T_rates = NULL
for(i in Temp) {
  maint <- rho1*(x^rho2)*rM_T_AL2(i, x, T_par)
  intake  <- IM_r*eps1*(x^eps2)*rI_T_GU2(i, x, T_par) 
  T_rates  <- rbind(T_rates,
                   cbind(Temp = i, mass = x, 
                         maint, intake, 
                         growth_E = Kap*alpha*Y*intake - maint, 
                         repro_E = alpha*Y*(1-Kap)*intake))
  }

T_rates_long <- gather(as.data.frame(T_rates), key = "E_type", value ="g_day", c(5,6)) # gather only column 6,7

## Over mass for three temperatures
T_rates_long %>%
  filter( Temp %in% c(283,284,285)) %>%
  ggplot(.) +
  geom_line(size=1 , aes(mass, g_day, color = as.factor(Temp), 
                         linetype = as.factor(E_type))) +
  ylim(0,10) +
  theme_bw() +
  scale_colour_brewer(palette="Dark2")

## Over temepratures for three sizes
T_rates_long %>%
  filter( mass %in% c(x[2], x[10], x[100])) %>%
  ggplot(.) +
  geom_line(size=1 , aes(Temp, g_day, color = as.factor(mass), 
                         linetype = as.factor(E_type))) +
  geom_hline(yintercept = 0) +
  ylim(-5,15) +
  theme_bw() +
  scale_colour_brewer(palette="Dark2")


#### PLOT GROWTH / SIZE AT AGE DEPENDENT ON MASS, TEMPERATURE & kAPPA ####

# Build a dataframe of sizes, fecundity, temperatures and kappas using ratefun()
Temp = 283:285   # temperatures
Kap  = Kap #c(0.75,0.8,0.85)  # kappas, allocation to growth (roughly 0.85 in Windermere pike)
outR = NULL
for (j in Kap) {
  for (l in Temp) {
    plotpars <- c(T = l,     # parameters for Temperature, feeding, allocation and Mass dependence
                  Kappa = j, # allocation to respiration (Growth and maintenance)
                  Y=Y,         # Feeding level
                  alpha=alpha,     # assimilation efficiency
                  IM_r=IM_r,      # ratio Intake maintenance
                  eps1=eps1,      # Intake allometric scalar, 0.248 based on roach by Lindmark assuming that intake scales equally with mass between species in his system
                  eps2=eps2,      # Intake allometric exponent, 0.64 from Lindmark unpub. 0.725 to get an length asymtopte Lm for pike at 292 K. 0.767 based on roach by Lindmark assuming that intake scales equally with mass between species in his system
                  rho1=rho1,      # Maintenance allometric scalar, Guesstimate 0.0204. 0.02 based on pike estimated by Lindmark
                  rho2=rho2)      # Maintenance allometric exponent,0.77 from Lindmark unpub.Guesstimate 0.8. 0.76 based on pike estimated by Lindmark,  0.807 to get an length asymtopte Lm for pike at 292 K
    
  for(seas in 1:20) { # growth from age 1 happens in this for-loop, 25 years is from age 2 to 26
      if(seas == 1){
        offmean <- ratefun(o_e, plotpars)[134,2]
        p <- data.frame(mass = ratefun(offmean, plotpars)[sl+1,2],
                        re = DEBrepfun(offmean, plotpars, o_e), 
                        age = seas+1) #age 2 as offmean represents age 1
      } else p <- rbind(p,
                      c(ratefun(p[nrow(p),1], plotpars)[sl+1,2],
                        DEBrepfun(p[nrow(p),1], plotpars, o_e), 
                        seas+1))
    }
    p <- rbind(c(offmean, DEBrepfun(o_e, plotpars, o_e),1),
               p) # add Age 1
    y <- data.frame(p,
                  Temp = l,
                  Kappa = j)
    outR <- rbind(outR,y)
      }
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

#Plot parameters for IPM vital rate functions
GR_pars <- c(T =283,     # parameters for Temperature, feeding, allocation and Mass dependence
             Kappa = Kap, # allocation to respiration (Growth and maintenance)
             Y=Y,         # Feeding level
             alpha=alpha,     # assimilation efficiency
             IM_r=IM_r,      # ratio Intake maintenance
             eps1=eps1,      # Intake allometric scalar, 0.248 based on roach by Lindmark assuming that intake scales equally with mass between species in his system
             eps2=eps2,      # Intake allometric exponent, 0.64 from Lindmark unpub. 0.725 to get an length asymtopte Lm for pike at 292 K. 0.767 based on roach by Lindmark assuming that intake scales equally with mass between species in his system
             rho1=rho1,      # Maintenance allometric scalar, Guesstimate 0.0204. 0.02 based on pike estimated by Lindmark
             rho2=rho2)

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
 #geom_line(data=GData, aes(Age, ltow_A(Length), group = as.factor(Ind))) +
 geom_line(size=1, color = "red") +
 theme_bw() +  
 scale_colour_brewer(palette="Dark2")

# Plot DEBrepfun and windermere data ####
plot(FData$Weight,FData$Eggs,ylab="fecundity", main = "fecundity of x in s+1", ylim = c(0,400000)) # Windermere Pike data
outR %>%
  filter(Temp == 283) %>%
  lines(re~mass, data=., lwd = 2, type = "l", col = "red")
  legend("topleft", c("DEB", "Windermere Pike"), lty=c(1, NA), pch =c(NA,1), col= c("red","black"),cex=.7)

# Plot DEBoff.size ####
# Temp effect on mean size not plotted yet
  plot(x,DEBoff.size(x, GR_pars), type = "l", xlab = "y", ylab = "f(x)",
     main= "Offspring size distribution") 

# Plot DEBsurv ####
# Temp effect on surv not plotted yet
plot(x, DEBsurvfun(x, GR_pars), type="l", col="red",ylim=c(0,1))

#VVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVV
#### ### ## # PLOT RESULTS #### ### ## # ## ### ####
#VVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVV

### PLOT LAMBDA AS A FUNCTION OF KAPPA AND/OR TEMP ####
#pdf("\\\\storage-og.slu.se/home$/vitl0001/Desktop/kappaTlamdba_20201120.pdf", width = 8, height = 3.5)
ggplot(as.data.frame(Res_lam), aes(Kappa, lam, color = as.factor(T))) +
  geom_line(size=0.5) +
  ylab(expression(lambda~(Fitness))) +
  xlab(expression(Kappa~(Growth~Allocation))) +
  labs(color = "Temp [K]") +
  theme_bw()

  
### PLOT STABLE STRUCTURE W AS A FUNCTION OF KAPPA AND/OR TEMP ####
colnames(Res_w)[10:ncol(Res_w)] <- c(0,x)
Res_w_long <- gather(Res_w, key = "Size", value ="biom", c(11:ncol(Res_w))) # not use column 10 & 11 (eggstage and recruits)
str(Res_w_long)
Res_w_long %>%
  filter(near(Kappa, 0.8)) %>% #Floating point issue when comparing vector, therefore the use of near()
  filter(T == 283) %>% 
  ggplot(., aes(as.numeric(Size), biom), color = as.factor(T)) +
    ggtitle("Only for T=283 & Kappa=0.8 at the moment") + # for the main title
    geom_line(size=0.5) +
    ylab("w") +
    xlab("Size") +
    xlim(0,14000) +
    theme_bw()+
    labs(color = "Temp [K]")

# YVs plot method from res for comparison
#plot(x,res$w[2:(n+1)],type="l")# ,ylim=c(0,.0001)) # the stable structure "w"

### PLOT REPRODUCTIVE VALUES V AS A FUNCTION OF KAPPA AND/OR TEMP ####
colnames(Res_v)[10:ncol(Res_v)] <- c(0,x)
Res_v_long <- gather(Res_v, key = "Size", value ="biom", c(11:ncol(Res_v))) # not use column 9 & 10 (eggstage and recruits)
str(Res_v_long)
Res_v_long %>%
  filter(Kappa == 0.8) %>% 
  filter(T == c(283,284,285))  %>% 
  ggplot(., aes(as.numeric(Size), biom, color = as.factor(T), linetype = as.factor(Kappa))) +
    geom_line(size=0.5) +
    ylab("Reproductive value V") +
    xlab("Size")+
    xlim(0,14000) +
    theme_bw() +
    labs(color = "Temp [K]")
#dev.off()

# YVs plot method from res for comparison
#plot(x,res$v[2:(n+1)],type="l" ) # reproductive values "v"

# ### PLOT the K matrix #### NOT WORKING

# L <- as.data.frame(K)
# rownames(L)= c(x,x[ncol(L)-1]+1)
# colnames(L)= c(0.1,x)

# library(reshape2)
# library(ggplot2)
# longData<-melt(K)
# ggplot(longData, aes(x = Var2, y = Var1)) + 
#   geom_raster(aes(fill=value)) + 
#   scale_fill_gradient(low="yellow", high="red") +
#   labs(x="x", y="y", title="Kernel") +
#   theme_bw() + theme(axis.text.x=element_text(size=9, angle=0, vjust=0.3),
#                      axis.text.y=element_text(size=9),
#                      plot.title=element_text(size=11))
# 
# Using countor
# contour(t(Mknsd[(n+1):1,1:(n+1)]),lty=1,lwd=.2,xaxt="n+1",xaxs="i",yaxs="i",yaxt="n",xlab="From Weight",ylab="To Weight",levels=seq(0,.1,.01))
# # contour(t(K[(n+1):1,1:(n+1)]),lty=1,lwd=.2,add=T,levels=seq(0,.5,.001), xaxt="n",xaxs="i",yaxs="i",yaxt="n",xlab="From Weight",ylab="To Weight",main="Mean projection kernel")
# # contour(t(K[(n+1):1,1:(n+1)]),lty=1,lwd=.2,xaxt="n",xaxs="i",yaxs="i",yaxt="n",xlab="From Weight",ylab="To Weight")
