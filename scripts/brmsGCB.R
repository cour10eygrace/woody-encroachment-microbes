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

