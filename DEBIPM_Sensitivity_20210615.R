
#### Sensitivity of lambda to kappa decomposition into vital rate contributions ###

# Analysis to understand how size dependent vital rates affects optimum 
# kappa dependent on temperature. I.e. How does changes in k_0 affect lambda 
# through different vital rates and dependent on temperature

# Model output on lambda, w and v:
Res_lam <- read.delim("//storage-og.slu.se/home$/vitl0001/My Documents/Manus2/Results/Results0624/Res_lam_1_0624.txt", sep = ",")
#Res_lam <- read.delim("Res_lam_1_0427.txt", sep = ",")
Res_v <- read.delim("//storage-og.slu.se/home$/vitl0001/My Documents/Manus2/Results/Results0624/Res_v_1_0624.txt", sep = ",")
#Res_v <- read.delim("Res_v_1_0427.txt", sep = ",")
colnames(Res_v)[4:ncol(Res_v)] <- sub("X", "", colnames(Res_v)[4:ncol(Res_v)])
Res_w <- read.delim("//storage-og.slu.se/home$/vitl0001/My Documents/Manus2/Results/Results0624/Res_w_1_0624.txt", sep = ",")
#Res_w <- read.delim("Res_w_1_0427.txt", sep = ",")
colnames(Res_w)[4:ncol(Res_w)] <- sub("X", "", colnames(Res_w)[4:ncol(Res_w)])

# Numerical differentiation (difference forward differentiation, 
# i.e. ( f(x + h) - f(x) ) / h) ) of dGrow/dk_0, dFec/dk_0 and dAge1/k_0. 
G_NumDer <- function(h ,m, Pars){
  (DEBgrowthfun(m, Pars = c(T = Pars[["T"]],kappa = Pars[["kappa"]] + h,Y = Pars[["Y"]])) -
     DEBgrowthfun(m, Pars = c(T = Pars[["T"]],kappa = Pars[["kappa"]],Y = Pars[["Y"]])))/h
} # Note that DEBgrowthfun gives the n * n matrix of prob. density distributions, 
# the rowSums thus gives density of a size class in t+1
#plot(x,G_NumDer(h=0.00001,x,GR_pars), type = "l") 
#lines(x,G_NumDer(h=0.0001,x,GR_pars), type = "l", col="red") # converges at h=0.0001
F_NumDer <- function(h ,m, Pars){
  (DEBrepfun(m, Pars = c(T = Pars[["T"]],kappa = Pars[["kappa"]] + h,Y = Pars[["Y"]])) - 
     DEBrepfun(m, Pars = c(T = Pars[["T"]],kappa = Pars[["kappa"]],Y = Pars[["Y"]])))/h
}
# plot(x,F_NumDer(h=0.000001,x,GR_pars), type = "l") 
# lines(x,F_NumDer(h=0.00001,x,GR_pars), type = "l", col="red") # converges at h=0.00001

A1_NumDer <- function(h ,m, Pars){
  (DEBage1.size(Pars = c(T = Pars[["T"]],kappa = Pars[["kappa"]] + h,Y = Pars[["Y"]])) - 
     DEBage1.size(Pars = c(T = Pars[["T"]],kappa = Pars[["kappa"]],Y = Pars[["Y"]])))/h
}
# plot(x,A1_NumDer(h=0.0001,x,GR_pars), type = "l")
# lines(x,A1_NumDer(h=0.001,x,GR_pars), type = "l", col="red") # converges at h=0.001


whichMinMax <- function(sens){
  sens[which(is.na(sens))] <-0 # replace na:s with 0
  sens[which.max( abs(sens) )] # max (or min) of each row 
}


# optimum kappa dependent on temperature:
OptPars <- 
  as.data.frame(maxl_1[,1:3]) %>%
  filter(T %in% c(287, 289, 291))

Sens_Pars <- rbind(OptPars, cbind(T=c(287,289,291),kappa=c(0.83,0.83,0.83), Y=c(1,1,1)))
#Res_Sens <- matrix(ncol = 3+1+n, nrow = nrow(Sens_Pars)) # assuming 40 here from files
#parsK <- as.matrix(expand.grid(T,kappa,Y))
#colnames(parsK) <- c("T","kappa","Y")

Res_Sens <-NULL
for (i in 1:nrow(Sens_Pars)){
  res_v <- as.numeric(Res_v %>% filter(T == Sens_Pars[i,1] & kappa==Sens_Pars[i,2]))[4:ncol(Res_v)]
  res_w <- as.numeric(Res_w %>% filter(T == Sens_Pars[i,1] & kappa==Sens_Pars[i,2]))[4:ncol(Res_w)]
 
  Gr_cont <- bind_cols(Size=x, Temp=Sens_Pars[i,1], kappa=Sens_Pars[i,2],
                          vr="Gr", outer(res_v[2:(n+1)], res_w[2:(n+1)]) * DEBsurvfun(x, Sens_Pars[i,]) * G_NumDer(h=0.00001,x,Sens_Pars[i,])
                         )
  Gr_cont["Sens"] <-apply(Gr_cont[5:ncol(Gr_cont)],1,FUN=whichMinMax)
  Gr_cont<-Gr_cont[,c(1:4,ncol(Gr_cont))]
  
  F_cont  <- bind_cols(Size=x, Temp=Sens_Pars[i,1], kappa=Sens_Pars[i,2],
                          vr="F",  Sens=res_v[1]*res_w[2:(n+1)] * DEBsurvfun(x, Sens_Pars[i,]) * F_NumDer(h=0.00001,x,Sens_Pars[i,]))
  A1_cont <- bind_cols(Size=x, Temp=Sens_Pars[i,1], kappa=Sens_Pars[i,2], 
                          vr="A1", Sens=res_v[2:(n+1)]*res_w[1] * el_surv * A1_NumDer(h=0.00001,x,Sens_Pars[i,]))
  
  Res_Sens <- rbind(rbind(Gr_cont,F_cont,A1_cont), Res_Sens) #c(T[d],kappa[e],Y[f], res$lambda)
}

Res_Sens <-
Res_Sens %>%
  mutate(Scen = ifelse(kappa==0.83,"Baseline k_0","Optimal k_0"))
       # colnames(Res_R0) <- c("T","kappa","Y","R0")

mycolors<-RColorBrewer::brewer.pal(4,"Dark2")[2:4]
cont.labs <- c("Age 1 size", "Fecundity", "Growth")
names(cont.labs) <- c("A1", "F", "Gr")
Fig5_sens <-
Res_Sens %>%
  mutate(vr=as.factor(vr)) %>%
  ggplot(.) +
  geom_line(aes(Size, Sens, linetype = Scen, color = as.factor(Temp))) +
  #facet_wrap(Scen~vr, nrow=2, scales = "free", labeller = labeller(vr=cont_labs)) +
  facet_grid(.+vr~Temp, labeller = labeller(vr=cont.labs)) +
  xlim(0,2000) +
  ylab("Sensitivity") +
  xlab("Weight [g]") +
  scale_colour_manual(values = mycolors, name= "Temperature [K]") +
  theme_bw() +
  theme(panel.background = element_rect(fill = "white", colour = "black"),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA))

#pdf("DEBIPM_Fig5_Sens__small_20210709.pdf", width = 8, height = 6)
pdf("DEBIPM_Fig5_Sens__all_20210709.pdf", width = 8, height = 6)
Fig5_sens
dev.off()

# Is the sum of sensitivites 0 at optimum kappa? For T=287:
sum(Gr_cont289,F_cont289,A1_cont289)
sum(Gr_cont287[4:ncol(Gr_cont287)],na.rm=TRUE)
sum(F_cont289$Sens,na.rm=TRUE)
sum(A1_cont289$Sens,na.rm=TRUE)
# NO, what am I doing wrong

sum(res_w287[2:(n+1)])




### OLD OLD ###

# THERMAL or kappa SENSITIVITY OF LAMBDA ####
# Sensitivity analysis to study how change in Lambda can be
# attributed to change in different vital rates originating from 
# a change in kappa or temperature.
# The partial derivative of the output Y with respect to an input factor Xi
# 
# # Sensitivity based on Vindenes 2014:
# 
# dG.dm <- function(m, Pars){
#   (DEBgrowthfun(m+0.1, Pars) - DEBgrowthfun(m, Pars))/0.1
# }
# 
# #Sensitivity of the growth density function to temperature (approximate) 
# 
# dG.dT <- function(m, Pars){
#   (DEBgrowthfun(m, c(T=Pars[["T"]]+0.001, 
#                      kappa =Pars[["kappa"]], 
#                      Y = Pars[["Y"]])) 
#    - DEBgrowthfun(m, Pars))/0.001
# }
# 
# ##SURVIVAL SENSITIVITY
# #Sensitivity of the survival probability function to temperature (approximate) 
# 
# dS.dT <- function(m, Pars){
#   (DEBsurvfun(m, c(T=Pars[["T"]]+0.001, 
#                    kappa =Pars[["kappa"]], 
#                    Y = Pars[["Y"]])) 
#    - DEBsurvfun(m, Pars))/0.001
# }
# 
# #Sensitivity of the survival probability function to length (approximate) 
#  
# dS.dm <- function(m, Pars){
#   (DEBsurvfun(m+0.1, Pars) - DEBsurvfun(m, Pars))/0.1
# }
# 
# # OFFSPRING SIZE SENSITIVITY
# #Sensitivity of the offspring length density function to temperature (approximate)
# 
# dO.dT <- function(Pars){
#   (DEBage1.size(c(T=Pars[["T"]]+0.001, 
#                  kappa =Pars[["kappa"]], 
#                  Y = Pars[["Y"]]))
#    - (DEBage1.size(Pars)))/0.001
# }
# 
# # Sensitivity of the fecundity function to length (approximate) 
#  
# dR.dm <- function(m, Pars){
#   (DEBrepfun(m+0.1, Pars) - DEBrepfun(m, Pars))/0.1
# }
# 
# # Sensitivity of the fecundity function to temperature (approximate) 
# 
# dR.dT <- function(m, Pars){
#   (DEBrepfun(m, c(T=Pars[["T"]]+0.001, 
#                   kappa =Pars[["kappa"]], 
#                   Y = Pars[["Y"]])) 
#    - DEBrepfun(m, Pars))/0.001
# }
# 
# # Vindenes etal 2014 code
# thermal.sensitivity <- function(Pars){
#   res <- wvlambda.projection(Kmat = K.matrix(Pars))
#   Fmat <- Smat <- Omat <- Gmat <- matrix(NA, n+1, n+1)
#   sx <- c(DEBsurvfun(e_m, Pars),DEBsurvfun(x, Pars))
#   dsx <- dS.dT(c(e_m,x), Pars)
#   fx <- c(DEBrepfun(0, Pars), DEBrepfun(x, Pars))
#   dfx <- dR.dT(c(0,x), Pars)
#   offs <- c(0,DEBage1.size(Pars))
#   doff <- c(0,dO.dT(Pars))
#   gx <- dgx <- matrix(0,n+1,n+1)
#   gx[2:(n+1),2:(n+1)] <- DEBgrowthfun(x, Pars) # unlike vindenens 2014, the DEBgrowthfun returns all y-distributions
#   for(j in 2:(n+1)){
#     dgx[2:(n+1),j] <- dG.dT(x[j-1], Pars)
#   }
#   gx[,n+1] <- gx[,n]
#   dgx[,n+1] <- dgx[,n]
#   for(i in 1:(n+1)){			
#     Fmat[i,] <- res$v[i]*res$w*(dfx*sx[i])*dx
#     Omat[i,] <- res$v[i]*res$w*(doff[i]*el_surv)*dx
#     Smat[i,] <- res$v[i]*res$w*(dsx*gx[i,])*dx
#     Gmat[i,] <- res$v[i]*res$w*(dgx[i,]*sx[i])*dx
#   }
#   list("Fecundity contribution"= apply(Fmat,2,sum),"Survival contribution"=apply(Smat,2,sum),"Growth contribution"= apply(Gmat,2,sum),"Offspring length contribution"= apply(Omat,2,sum),"Total contribution"=apply(Fmat,2,sum)+apply(Smat,2,sum)+apply(Gmat,2,sum)+apply(Omat,2,sum))
# }	
#  
# # Tsensit <- thermal.sensitivity(GR_pars)
# # # 1=fec sens, 2=surv sens, 3=growth sens, 4=offs sens,
# # 
# # par(mfrow=c(2,2), las=1, bty="l")
# # plot(c(0,x), Tsensit[[1]], type="l", ylab="", ylim=c(-1e1,1e5),main=expression(paste("Fecudity contributions to d",  lambda, "/dT")), col=1, lwd=2, xlab="x")
# # plot(c(0,x), Tsensit[[2]], type="l", ylab="", ylim=c(-1e-5,1e-5), main=expression(paste("Survival contributions to d",  lambda, "/dT")), col=2, lwd=2, xlab="Size x")
# # plot(c(0,x), Tsensit[[3]], type="l", ylab="", ylim=c(-1e-6,1e-5),main=expression(paste("Growth Contributions to d",  lambda, "/dT")), col=3, lwd=2, xlab="Size x")
# # plot(c(0,x), Tsensit[[4]], type="l", ylab="", ylim=c(0,1e-9),main=expression(paste("Offspring size contributions to d",  lambda,"/dT")), col=2, lwd=2, xlab="Size x")
# # #plot(c(0,x), Tsensit[[5]], type="l", ylab="", ylim=c(0,1e-9),main=expression(paste("Contributions to d",  lambda, "/d",mu[z])), col=2, lwd=2, xlab="Size x")
# 
# Size.sensitivity <- function(Pars){
#   res <- wvlambda.projection(Kmat = K.matrix(Pars))
#   Fmat <- Smat <- Gmat <- matrix(NA, n+1, n+1)
#   sx <- c(DEBsurvfun(e_m, Pars),DEBsurvfun(x, Pars))
#   dsx <- dS.dm(c(e_m,x), Pars)
#   fx <- c(DEBrepfun(0, Pars), DEBrepfun(x, Pars))
#   dfx <- dR.dm(c(0,x), Pars)
#   gx <- dgx <- matrix(0,n+1,n+1)
#   gx[2:(n+1),2:(n+1)] <- DEBgrowthfun(x, Pars) # unlike vindenens 2014, the DEBgrowthfun returns all y-distributions
#   for(j in 2:(n+1)){
#     dgx[2:(n+1),j] <- dG.dm(x[j-1], Pars)
#   }
#   gx[,n+1] <- gx[,n]
#   dgx[,n+1] <- dgx[,n]
#   for(i in 1:(n+1)){			
#     Fmat[i,] <- res$v[i]*res$w*(dfx*sx[i])*dx
#     Smat[i,] <- res$v[i]*res$w*(dsx*gx[i,])*dx
#     Gmat[i,] <- res$v[i]*res$w*(dgx[i,]*sx[i])*dx
#   }
#   list("Fecundity contribution"= apply(Fmat,2,sum),"Survival contribution"=apply(Smat,2,sum),"Growth contribution"= apply(Gmat,2,sum),"Total contribution"=apply(Fmat,2,sum)+apply(Smat,2,sum)+apply(Gmat,2,sum))
# }	
# 
# # Ssensit <- Size.sensitivity(GR_pars)
# # # 1=fec sens, 2=surv sens, 3=growth sens,4 = total
# # 
# # par(mfrow=c(1,1), las=1, bty="l")
# # plot(c(0,x), Ssensit[[1]], type="l", ylab="", ylim=c(-1e1,1e7),main=expression(paste("Fecudity contributions to d",  lambda, "/dT")), col=1, lwd=2, xlab="x")
# # plot(c(0,x), Ssensit[[2]], type="l", ylab="", ylim=c(-1e-8,1e-7), main=expression(paste("Survival contributions to d",  lambda, "/dT")), col=2, lwd=2, xlab="Size x")
# # plot(c(0,x), Ssensit[[3]], type="l", ylab="", ylim=c(0,1e-8),main=expression(paste("Growth Contributions to d",  lambda, "/dT")), col=3, lwd=2, xlab="Size x")
# 
