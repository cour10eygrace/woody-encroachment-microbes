#load packages 
library(GGally)
library(formatR)
library(brms)
library(lme4)
library(quantreg)
library(ggplot2)
library(RColorBrewer) 
library(gridExtra) # for grid.arrange 
library(tidyverse)
library(prettydoc)
library(tidybayes)
library(dplyr)
library(reshape2)
library(maps)
library(mapdata)

#setup data for brms 
soil.data<-read.csv("soil.data.csv")
soil.data$soilCN<-soil.data$soilC2/soil.data$soilN2
adiv<-read.csv("ITS_alphadiv_bysample.csv")
adiv_bac<-read.csv("16S_alphadiv_bysample.csv")
adiv_all<-merge(adiv, adiv_bac, by="Tube", all.x = T, all.y=T)

# Confirm Juniperus AMF and Salix species ECM
soil.data$root.microb[grep("Salix",soil.data$plant.species)] <- "ECM.ERM"
soil.data$root.microb[grep("Juniperus",soil.data$plant.species)] <- "AMF"

dat<-merge(soil.data, adiv_all, by="Tube", all.x = T, all.y=T)

site.table <- group_by(dat[dat$veg=="S",], site) %>%
  summarize(lat = mean(coord.Y)
            ,long=mean(coord.X)
            ,elev = mean(elevation)
            ,root.microb= unique(root.microb)
            ,plant.species= unique(plant.species)[1]
  ) %>%
  as.data.frame()
clim.data<-read.csv("updated2.1_climate.data.csv") #has WorldClim v2.1 
site.table <- full_join(site.table, clim.data,by="site")
site.table<-site.table[-14,]#take out blank row

leaf.data<-read.csv("traits.csv") #updated 1/20/20
#new Mexico values-LDMC & SLA and corrected SLA-Spain and LDMC-Spain Ordesa 

#load leaf trait ordinations-output from leaftraitsGCB.R
PCA_traits<-read.csv("PCA_traits_euclidian.csv")

leaf.data <- merge(leaf.data, dplyr::select(PCA_traits, PC1,PC2,PC3,SampleID), by = "SampleID")
leaf.data <- rename(leaf.data, site = Location)
leaf.data$site<-trimws(leaf.data$site) #need this or dat2 does not merge correctly 

#average leaf traits at site level 
leaf.summary <- group_by(leaf.data, site) %>%
  summarise(ldmc.mean = mean(LDMC..g.g., na.rm=TRUE)
            ,sla.mean = mean(SLA..cm.2.g., na.rm=TRUE)
            ,leaf15N.mean = mean(X15N, na.rm=TRUE)
            ,leaf13C.mean = mean(X13C, na.rm=TRUE)
            ,leafC.mean = mean(Wt..C., na.rm=TRUE)
            ,leafN.mean = mean(Wt..N., na.rm=TRUE)
            ,ldmc.sd = sd(LDMC..g.g., na.rm=TRUE)
            ,sla.sd = sd(SLA..cm.2.g., na.rm=TRUE)
            ,leaf15N.sd = sd(X15N, na.rm=TRUE)
            ,leaf13C.sd = sd(X13C, na.rm=TRUE)
            ,leafC.sd = sd(Wt..C., na.rm=TRUE)
            ,leafN.sd = sd(Wt..N., na.rm=TRUE)
            ,PCA1.mean = mean(PC1, na.rm=TRUE)
            ,PCA1.sd = sd(PC1, na.rm=TRUE)
            ,PCA2.mean = mean(PC2, na.rm=TRUE)
            ,PCA2.sd = sd(PC2, na.rm=TRUE)
  ) %>% 
  as.data.frame()

dat2 <- merge(dat, leaf.summary, by="site",all.x=TRUE)
head(dat2)
dat3 <- merge(dat2, site.table, by=c("site"),all.x=TRUE) %>%
  dplyr::select(-root.microb.x) %>%
  dplyr::rename(symbiont = root.microb.y) #use the site.table root microbe, so it's the same for the whole site

#separate data for fungi and bacteria
ITS.dat <- with(dat3,
                data.frame(
                  site = site
                  ,veg = veg
                  ,symbiont = symbiont
                  ,ITS.otus =scale(ITS_otus)
                  ,ph = scale(pH)
                  ,vwc = scale(VWC)
                  ,soilN = scale(soilN)
                  ,soilCN = scale(soilCN)
                  
                  ,sla = scale(sla.mean)
                  ,ldmc = scale(ldmc.mean)
                  ,leafC = scale(leafC.mean)
                  ,leafN = scale(leafN.mean)
                  ,leaf15N  = scale(leaf15N.mean)
                  ,leaf13C = scale(leaf13C.mean)
                  ,leafPCA1=scale(PCA1.mean) #traits PCA1
                  ,leafPCA2=scale(PCA2.mean) #traits PCA2 
                  
                  ,temp = scale(MAT) # deg C
                  ,precip = scale(MAP) # cm precip
                ))

ITS.dat <- dplyr::filter(ITS.dat, is.na(ITS.otus)==FALSE)  ## take out NA's in response

summary(ITS.dat)
X16S.dat <- with(dat3 
                 ,data.frame(site = site
                             ,veg = veg
                             ,symbiont = symbiont
                             ,X16S.otus =scale(X16S.otus)
                             ,ph = scale(pH)
                             ,vwc = scale(VWC)
                             ,soilN = scale(soilN)
                             ,soilCN = scale(soilCN)
                             
                             ,sla = scale(sla.mean)
                             ,ldmc = scale(ldmc.mean)
                             ,leafC = scale(leafC.mean)
                             ,leafN = scale(leafN.mean)
                             ,leaf15N  = scale(leaf15N.mean)
                             ,leaf13C = scale(leaf13C.mean)
                             ,leafPCA1=scale(PCA1.mean) #traits PCoA1
                             ,leafPCA2=scale(PCA2.mean) #traits PCoA2 
                             
                             ,temp = scale(MAT) # deg C
                             ,precip = scale(MAP) # cm precip
                 ))
X16S.dat <- dplyr::filter(X16S.dat, is.na(X16S.otus)==FALSE)  ## take out NA's in response
X16S.dat <- dplyr::filter(X16S.dat, !is.na(veg))  ## take out NA's in veg (2 NTC samples)

summary(X16S.dat)

save(ITS.dat, X16S.dat, file="GCB_brms_datav2.RData") 


#Fig 2 a,b 
microbe.data <- dplyr::select(dat3, site,veg,symbiont,ITS_otus,X16S.otus) %>%
  reshape2::melt(id.vars=c("site","veg","symbiont")) #%>%

itsraw <- ggplot(filter(microbe.data,variable=="ITS_otus"), aes(x=fct_reorder(site,as.integer(symbiont)),y=log(value))) + 
  geom_boxplot(aes(fct_reorder(site,as.integer(symbiont)), fill=veg, col=symbiont)) +
  scale_fill_brewer(palette="Paired", name = "", labels = c("Herb","Woody")) + 
  scale_colour_brewer(palette = 'Set2', type='qual')+
  theme_minimal()+
  labs(y = "Fungal richness (Log OTUs)", x='Site') + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) 

x16raw <- ggplot(filter(microbe.data,variable=="X16S.otus"),aes(x=fct_reorder(site,as.integer(symbiont)),y=log(value))) + 
  geom_boxplot(aes(fct_reorder(site,as.integer(symbiont)), fill=veg, col=symbiont)) +
  scale_fill_brewer(palette="Paired", name = "", labels = c("Herb","Woody")) + 
  scale_colour_brewer(palette = 'Set2', type='qual')+
  theme_minimal()+ ylim(5, 8)+
  labs(y = "Bacterial richness (Log OTUs)", x='Site') + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) #+ 

grid.arrange(itsraw,x16raw, nrow = 1)

#Fig S1
a<-ggplot(data=site.table, aes(x=(MAP/10), y=MAT))+
  geom_point(aes(x=(MAP/10), y=MAT, col=root.microb.x), shape=15, size=2)+ 
  geom_text(aes(x=(MAP/10)+1, label=site), size=2.5, nudge_y = -0.2, nudge_x = -0.2) +
  theme_classic( )+ ylab("MAT(C)") + xlab("MAP(cm)")+ theme(legend.position = "none")
b<-ggplot(data=site.table, aes(x=elev.x, y=MAT))+
  geom_point(aes(x=elev.x, y=MAT, col=root.microb.x), shape=15, size=2)+ 
  geom_text(aes(x=elev.x+100, label=site),size=2.5, nudge_y = -0.2, nudge_x = -0.2)+
    theme_classic()+ ylab("MAT(C)") + xlab("Elevation(m)")+ theme(legend.position = "none")
c<-ggplot(data=site.table, aes(x=elev.x, y=MAP/10))+
  geom_point(aes(x=elev.x, y=MAP/10, col=root.microb.x), shape=15, size=2)+ 
 geom_text(aes(x=elev.x+100, label=site), size=2.5, nudge_y = 3, nudge_x = -0.5)+
  theme_classic()+ ylab("MAP(cm)") + xlab("Elevation(m)")
  
grid.arrange(a,b,c, nrow=1)
