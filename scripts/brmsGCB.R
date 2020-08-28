#load packages and data----
library(dplyr)
library(brms)
library(tidyr)
library(rstan)
library(bayestestR)

#load scaled (mean centered) data-output from brmssetup.R
load(file= "data/GCB_brms_data.RData")

#brms models---- 
#fungi----
m2ord.ITS.rich <- bf(ITS.otus ~  (1 | site) + symbiont*veg +  temp*veg  + precip*veg + vwc +  soilCN + ph )
CN.mod <- bf(soilCN ~ (1 | site) + leafPCA1+ leafPCA2)
ph.mod <- bf(ph ~ (1 | site) + leafPCA1+ leafPCA2)
vwc.mod<-bf(vwc ~ (1 | site) + veg)

sem2.ITS.rich <- brm(m2ord.ITS.rich + CN.mod + ph.mod + vwc.mod + set_rescor(FALSE),
data = ITS.dat, control = list(adapt_delta=0.99, max_treedepth = 12), cores=3, chains=3, iter=10000)

summary(sem2.ITS.rich)#this will give you mean estimates and confidence interals for all parameters except symbiont type

#calculate symbiont estimates----
postITSrich<- as_tibble(fixef(sem2.ITS.rich, summary=F))
postITSrich$ITSotus_AMF<-postITSrich$ITSotus_vegS
postITSrich$ITSotus_ECM_ERM<-postITSrich$ITSotus_vegS + postITSrich$`ITSotus_symbiontECM.ERM:vegS`
postITSrich$ITSotus_Nfix<-postITSrich$ITSotus_vegS + postITSrich$`ITSotus_symbiontNfix:vegS`

#Calculate pairwise differences between symbiont types from model posteriors
pairwise<-postITSrich
pairwise$NfixAMF<-pairwise$ITSotus_Nfix-pairwise$ITSotus_AMF
pairwise$NfixECMERM<-pairwise$ITSotus_Nfix-pairwise$ITSotus_ECM_ERM
pairwise$AMFECMERM<-pairwise$ITSotus_AMF-pairwise$ITSotus_ECM_ERM

#select parameters 
postITSrich<-dplyr::select(postITSrich, -ITSotus_Intercept, -soilCN_Intercept, -ph_Intercept,-vwc_Intercept, -ITSotus_symbiontECM.ERM, 
                    -ITSotus_symbiontNfix, -ITSotus_vegS)
postITSrich$`ITSotus_symbiontECM.ERM:vegS`<-NULL
postITSrich$`ITSotus_symbiontNfix:vegS`<-NULL

pairwise<-dplyr::select(pairwise,NfixAMF, NfixECMERM, AMFECMERM)

#calculate 85, 90, 95% confidence intervals---- 
temp <- summary(as.mcmc(postITSrich), quantiles = c(0.05, 0.1, 0.15, 0.85, 0.9, 0.95))
fit.stats <- as.data.frame(signif(cbind(temp$statistics[,1:2], temp$quantiles),3))
colnames(fit.stats) <- c("mean","sd","x5","x10","x15","x85","x90","x95")
fit.stats$param<-row.names(fit.stats)

#calculate 'significance' levels
fit.stats<-mutate(fit.stats, sig85=case_when(x15>0&x85>0~"Y", 
                                           x15<0&x85<0~"Y", 
                                           x15>0&x85<0~"N", 
                                           x15<0&x85>0~"N"))
fit.stats<-mutate(fit.stats, sig90=case_when(x10>0&x90>0~"Y", 
                                             x10<0&x90<0~"Y", 
                                             x10>0&x90<0~"N", 
                                             x10<0&x90>0~"N"))
fit.stats<-mutate(fit.stats, sig95=case_when(x5>0&x95>0~"Y", 
                                             x5<0&x95<0~"Y", 
                                             x5>0&x95<0~"N", 
                                             x5<0&x95>0~"N"))
#reorder for Table S2----
fit.stats$param = factor(fit.stats$param, levels = c("ITSotus_AMF",  "ITSotus_ECM_ERM", "ITSotus_Nfix", "ITSotus_temp",
   "ITSotus_vegS:temp", "ITSotus_precip", "ITSotus_vegS:precip", "ITSotus_vwc",  "ITSotus_ph", "ITSotus_soilCN",
   "ph_leafPCA1", "ph_leafPCA2", "soilCN_leafPCA1","soilCN_leafPCA2","vwc_vegS")) 

#calculate for pairwise in Table S2
tempx <- summary(as.mcmc(pairwise), quantiles = c(0.05, 0.1, 0.15, 0.85, 0.9, 0.95))
fit.statsx <- as.data.frame(signif(cbind(tempx$statistics[,1:2], tempx$quantiles),3))
colnames(fit.statsx) <- c("mean","sd","x5","x10","x15","x85","x90","x95")
fit.statsx$param<-row.names(fit.statsx)

#bacteria----
m2ord.16S.rich <- bf(log(X16S.otus+1) ~  (1 | site) + symbiont*veg +  temp*veg  + precip*veg +  vwc +  soilCN + ph )
CN.mod <- bf(soilCN ~ (1 | site) +leafPCA1 + leafPCA2)
ph.mod <- bf(ph ~ (1 | site) +leafPCA1 + leafPCA2)
vwc.mod<-bf(vwc ~ (1 | site) +veg)

sem2.16S.rich <- brm(m2ord.16S.rich + CN.mod + ph.mod + vwc.mod+ set_rescor(FALSE),
                     data = X16S.dat, control = list(adapt_delta=0.99, max_treedepth = 12), cores=3, chains=3, iter=10000)
summary(sem2.16S.rich)#this will give you mean estimates and confidence interals for all parameters except symbiont type

#calculate symbiont estimates----
post16Srich<- as_tibble(fixef(sem2.16S.rich, summary=F))
post16Srich$logX16Sotus1_AMF<-post16Srich$logX16Sotus1_vegS
post16Srich$logX16Sotus1_ECM_ERM<-post16Srich$logX16Sotus1_vegS + post16Srich$`logX16Sotus1_symbiontECM.ERM:vegS`
post16Srich$logX16Sotus1_Nfix<-post16Srich$logX16Sotus1_vegS + post16Srich$`logX16Sotus1_symbiontNfix:vegS`


#Calculate pairwise differences between symbiont types from model posteriors
pairwise<-post16Srich
pairwise$NfixAMF<-pairwise$logX16Sotus1_Nfix-pairwise$logX16Sotus1_AMF
pairwise$NfixECMERM<-pairwise$logX16Sotus1_Nfix-pairwise$logX16Sotus1_ECM_ERM
pairwise$AMFECMERM<-pairwise$logX16Sotus1_AMF-pairwise$logX16Sotus1_ECM_ERM

#select parameters 
post16Srich<-dplyr::select(post16Srich, -logX16Sotus1_Intercept, -soilCN_Intercept, -ph_Intercept,-vwc_Intercept,
 -logX16Sotus1_symbiontECM.ERM, -logX16Sotus1_symbiontNfix, -logX16Sotus1_vegS)
post16Srich$`logX16Sotus1_symbiontECM.ERM:vegS`<-NULL
post16Srich$`logX16Sotus1_symbiontNfix:vegS`<-NULL

pairwise<-dplyr::select(pairwise,NfixAMF, NfixECMERM, AMFECMERM)

#calculate 85, 90, 95% confidence intervals----
temp2 <- summary(as.mcmc(post16Srich), quantiles = c(0.05, 0.1, 0.15, 0.85, 0.9, 0.95))
fit.stats2 <- as.data.frame(signif(cbind(temp2$statistics[,1:2], temp2$quantiles),3))
colnames(fit.stats2) <- c("mean","sd","x5","x10","x15","x85","x90","x95")
fit.stats2$param<-row.names(fit.stats2)

#calculate 'significance' levels
fit.stats2<-mutate(fit.stats2, sig85=case_when(x15>0&x85>0~"Y", 
                                             x15<0&x85<0~"Y", 
                                             x15>0&x85<0~"N", 
                                             x15<0&x85>0~"N"))
fit.stats2<-mutate(fit.stats2, sig90=case_when(x10>0&x90>0~"Y", 
                                             x10<0&x90<0~"Y", 
                                             x10>0&x90<0~"N", 
                                             x10<0&x90>0~"N"))
fit.stats2<-mutate(fit.stats2, sig95=case_when(x5>0&x95>0~"Y", 
                                             x5<0&x95<0~"Y", 
                                             x5>0&x95<0~"N", 
                                             x5<0&x95>0~"N"))
#reorder for Table S2
fit.stats2$param = factor(fit.stats2$param, levels = c("logX16Sotus1_AMF",  "logX16Sotus1_ECM_ERM", "logX16Sotus1_Nfix", 
"logX16Sotus1_temp", "logX16Sotus1_vegS:temp", "logX16Sotus1_precip", "logX16Sotus1_vegS:precip", "logX16Sotus1_vwc", 
"logX16Sotus1_ph","logX16Sotus1_soilCN", "ph_leafPCA1", "ph_leafPCA2", "soilCN_leafPCA1","soilCN_leafPCA2", "vwc_vegS")) 

#calculate for pairwise in Table S2----
tempx <- summary(as.mcmc(pairwise), quantiles = c(0.05, 0.1, 0.15, 0.85, 0.9, 0.95))
fit.stats2x <- as.data.frame(signif(cbind(tempx$statistics[,1:2], tempx$quantiles),3))
colnames(fit.stats2x) <- c("mean","sd","x5","x10","x15","x85","x90","x95")
fit.stats2x$param<-row.names(fit.stats2x)

#Figures----
#Fig 3 a----
source("scripts/get_sem_params.R")
## fungi 
P0 <- get_sem_params(sem2.ITS.rich) 
level1.stats <- P0[[1]]
vars.order1 <- c("veg.AMF","veg.ECM.ERM","veg.Nfix"
                 ,"temp","vegS.temp","precip","vegS.precip","leafPCA1","leafPCA2","vwc","ph","soilCN")
level1.stats <- level1.stats[rownames(level1.stats) %in% vars.order1,]
level1.stats <- level1.stats[match(vars.order1, rownames(level1.stats)),]
level1.stats<-level1.stats[-c(8:9), ]

## bacteria 
source("scripts/get_sem_params_16S.R")
P0<-get_sem_params_16S(sem2.16S.rich) 
level1.stats.bac <- P0[[1]]
level1.stats.bac <- level1.stats.bac[rownames(level1.stats.bac) %in% vars.order1,]
level1.stats.bac <- level1.stats.bac[match(vars.order1, rownames(level1.stats.bac)),]
level1.stats.bac<-level1.stats.bac[-c(8:9), ]

param.names1 <- c(
  "AMF effect"
  ,"ECM.ERM effect"
  ,"N-Fixer effect"
  
  ,"Temp"
  ,"Temp x Woody"
  , "Precip"
  , "Precip x Woody"
  ,"VWC"
  ,"pH"
  ,"Soil C:N"
)


coul <- brewer.pal(7, "Dark2")[c(1,2,3,4,7,5,6)]
par(omi=c(0,.5,0,0), mar=c(4,7,2,3))
xl2 <- min(level1.stats$x2.5)
xr2 <- max(level1.stats$x97.5)
xl <- xl2 - (xr2-xl2)*.2  # add a bit of extra width
xr <- xr2 + (xr2-xl2)*.1  # add a bit of extra width
xx<- seq(xl,xr,by=.01)
N <- dim(level1.stats)[1]
treat.col='black'
plot(range(xx), range(c(0,N+1)), xlab='', ylab='', bty='l', type='n', axes=FALSE)
axis(1); axis(2, labels=FALSE, at=1:N)
mtext('Effect on OTU richness', 1, 3, cex=1)
mtext(param.names1, 2, 1, at=N:1, cex=1, las=1)

for(i in 1:N){
  zz <- N-i+1 + 0.15
  points(level1.stats$mean[i],zz, pch=16, col=coul[1])
  #arrows(x0= level1.stats$x2.5[i], x1=level1.stats$x97.5[i] , y0= zz, y1=zz, code=3, angle=90, col=coul[1], length=0,lwd=1)
  arrows(x0= level1.stats$x5[i], x1=level1.stats$x95[i] , y0= zz, y1=zz, code=3, angle=90, col=coul[1], length=0,lwd=2.5)
  signif.symbol <- ifelse(level1.stats$p[i] < .15 , "*"," ")
  signif.symbol <- ifelse(level1.stats$p[i] < .10 , "**", signif.symbol)
  signif.symbol <- ifelse(level1.stats$p[i] < .05 , "***", signif.symbol)
  #if(i > 3){ 
  text(min(xx), zz, signif.symbol, cex=1, adj=0, col=coul[1])# } #
  #if(i <= 3){ text(min(xx), (N+1-i), "A", cex=1, adj=0) } # not putting signif on symbiont intercepts
  
  zz <- N-i+1 - 0.15
  points(level1.stats.bac$mean[i],zz, pch=16, col=coul[2])
  #arrows(x0= level1.stats.bac$x2.5[i], x1=level1.stats.bac$x97.5[i] , y0= zz, y1=zz, code=3, angle=90, col=coul[2], length=0,lwd=1)
  arrows(x0= level1.stats.bac$x5[i], x1=level1.stats.bac$x95[i] , y0= zz, y1=zz, code=3, angle=90, col=coul[2], length=0,lwd=2.5)
  signif.symbol <- ifelse(level1.stats.bac$p[i] < .15 , "*"," ")
  signif.symbol <- ifelse(level1.stats.bac$p[i] < .10 , "**", signif.symbol)
  signif.symbol <- ifelse(level1.stats.bac$p[i] < .05 , "***", signif.symbol)
  #if(i > 3){ 
  text(min(xx), zz, signif.symbol, cex=1, adj=0, col=coul[2]) #} #
  #if(i <= 3){ text(min(xx), (N+1-i), "A", cex=1, adj=0) } # not putting signif on symbiont intercepts
}
abline(v=0, lty=2)
legend("right",legend=c("Fungi","Bacteria"), lwd=2, bty = 'n', col=coul[1:2],  cex=.8)


#Fig 3 c,d----
library(tidybayes)
library(tidyr)
library(modelr)
library(ggsci)

fit_precip_plot<-brm(ITS.otus ~  precip*veg, data=ITS.dat)
fit_temp_plot<-brm(X16S.otus ~temp*veg, data=X16S.dat)

plot2<-
  ITS.dat %>%
  group_by(veg) %>%
  data_grid(precip = seq_range(precip, n = 101)) %>%
  add_fitted_draws(fit_precip_plot, n = 300) %>%
  ggplot(aes(x = precip, y = ITS.otus, color = veg)) +
  theme_classic()+ ylim(-2, 4)+ 
  geom_line(aes(y = .value, group = paste(veg, .draw)), alpha = .1) +
  geom_point(data = ITS.dat)  +
  scale_color_jco(name="Vegetation", labels = c("herb", "woody"))+
  xlab("MAP") + ylab("ITS OTU richness")


plot3<-
  X16S.dat %>%
  group_by(veg) %>%
  data_grid(temp = seq_range(temp, n = 101)) %>%
  add_fitted_draws(fit_temp_plot, n = 300) %>%
  ggplot(aes(x = temp, y = X16S.otus, color = veg))  + ylim(-2,4)+
  geom_line(aes(y = .value, group = paste(veg, .draw)), alpha = .1) +
  geom_point(data = X16S.dat) + theme_classic() +
  scale_color_jco(name="Vegetation",labels = c("herb", "woody"))+
  xlab("MAT") + ylab("16S OTU richness")


grid.arrange(plot2, plot3, ncol=2)


#Fig S6
#a 
postITSrichl<-pivot_longer(postITSrich, cols=c("ITSotus_Nfix" ,"ITSotus_ECM_ERM" ,"ITSotus_AMF"), 
                           names_to = "symbiont")%>%dplyr::select(value, symbiont)
figa<-ggplot(postITSrichl, aes(x=value)) + geom_density(aes(fill=symbiont), alpha=0.8) +
  theme_classic()+ scale_fill_brewer(palette = 'Set2', type='qual')+
  xlab("Effect on ITS OTU richness")+ 
  geom_vline(xintercept = 0, lty=2)
#b  
post16Srichl<-pivot_longer(post16Srich, cols=c("logX16Sotus1_Nfix" ,"logX16Sotus1_ECM_ERM" ,"logX16Sotus1_AMF"), 
                           names_to = "symbiont")%>%dplyr::select(value, symbiont)
figb<-ggplot(post16Srichl, aes(x=value)) + geom_density(aes(fill=symbiont), alpha=0.8) +
  theme_classic()+ scale_fill_brewer(palette = 'Set2', type='qual')+
  xlab("Effect on 16S OTU richness")+ geom_vline(xintercept = 0, lty=2)+
   theme(legend.position = "none")

grid.arrange(figa, figb ,ncol=2)
