#load libraries 
library(MASS)
library(ggplot2)
library(grid)
library(dplyr)
library(reshape2)
library(lme4)
library(lmerTest)
library(optparse)
library(vegan)
library(ape)
library(gridExtra)

#load and organize data 
alldata<-read.csv("Alldata.csv")
alldata$soilCN<-alldata$soilC2/alldata$soilN2 
alldata$VWC<-as.numeric(alldata$VWC....)
alldata$X.1<-NULL
alldata$X<-NULL
alldata$X.2<-NULL
alldata$Notes<-NULL

#fungi 
bray<-read.csv("bray_curtis_dm.csv")
rownames(bray)<-bray$X
bray$X<-NULL
map<-read.csv("bray_map.csv") #sample ID by tube #
map<-merge(map, alldata, by="Tube")
map<-dplyr::select(map,soilCN, VWC, pH, root.microb3, Vegetation, Location, ï..Sample.ID)
rownames(map)<-map$ï..Sample.ID
map$ï..Sample.ID<-NULL
pmap<-match(rownames(bray), rownames(map))
map2b<-map[pmap, ]
bray1<-as.dist(as(bray,"matrix"))

ITS.clim<-read.csv("ITS.data.CC.csv")
rownames(ITS.clim)<-ITS.clim$X
ITS.clim<-ITS.clim[,c(13:14)]
map2b<-merge(map2b, ITS.clim, by="row.names")
row.names(map2b)<-map2b$Row.names
map2b$Row.names<-NULL

#leaf traits 
traits<-read.csv("PCA_traits_euclidian.csv")
traits<-dplyr::select(traits, PC1, PC2, X)
row.names(traits)<-traits$X  
map2bS<-merge(map2b, traits, by="row.names")
row.names(map2bS)<-map2b$X
map2bS$Row.names<-NULL
map2bS$X<-NULL

#bacteria
wunif<-read.csv("weighted_unifrac_dm.csv")
rownames(wunif)<-(wunif$ï..)
wunif$ï..<-NULL
colnames(wunif)<-row.names(wunif)
wunif$NTC1<-NULL
wunif$NTC2<-NULL
wunif<-wunif[-c(216,217), ]

map2unif<-dplyr::select(alldata,soilCN, VWC, pH, root.microb3, Vegetation, Location, Tube)
map2unif<-filter(map2unif,map2unif$Tube%in% row.names(wunif), .preserve = TRUE) 
rownames(map2unif)<-map2unif$Tube
wunif1<-as.dist(as(wunif,"matrix"))

X16S.clim<-read.csv("X16S.data.CC.csv")
row.names(X16S.clim)<-X16S.clim$X
X16S.clim<-X16S.clim[,c(13:14, 17)]
map2unif<-merge(map2unif, X16S.clim, by="row.names")
row.names(map2unif)<-map2unif$Tube

#leaf traits 
traits<-read.csv("PCA_traits_euclidian.csv")
traits<-dplyr::select(traits, PC1, PC2, Tube)
row.names(traits)<-traits$Tube  
map2unifS<-left_join(map2unif, traits)
map2unifS<-subset(map2unifS, Vegetation=='S')
row.names(map2unifS)<-map2unifS$Tube
map2unifS$Row.names<-NULL
map2unifS$X<-NULL

#run NMDS
#fungi
NMDS_bray<-metaMDS(bray, k=4, trymax=200)#stress=0.137 
NMDS_bray_coords<-as.data.frame(NMDS_bray$points)
NMDS_bray_coords<-merge(NMDS_bray_coords, map2b, by="row.names")
NMDS_bray_shrub<-left_join(NMDS_bray_coords, map2bS)#pull out shrub only for leaf traits 
NMDS_bray_shrub<-subset(NMDS_bray_shrub, Vegetation=="S")
bray1S<-filter(bray, row.names(bray) %in% NMDS_bray_shrub$Row.names)
bray1S<-dplyr::select(bray1S, one_of(NMDS_bray_shrub$Row.names))
bray1S<-as.dist(as(bray1S,"matrix"))


#bacteria
NMDS_wunif<-metaMDS(wunif, k=3, trymax=100)#stress=0.107
NMDS_wunif_coords<-as.data.frame(NMDS_wunif$points)
NMDS_wunif_coords<-merge(NMDS_wunif_coords, map2unif, by="row.names")
NMDS_wunif_shrub<-left_join(NMDS_wunif_coords, map2unifS) #pull out shrub only for leaf traits 
NMDS_wunif_shrub<-subset(NMDS_wunif_shrub, Vegetation=="S")
wunif1S<-dplyr::select(wunif, one_of(NMDS_wunif_shrub$Row.names))
wunif1S<-filter(wunif1S, row.names(wunif1S) %in% NMDS_wunif_shrub$Row.names)
x<-colnames(wunif1S)
rownames(wunif1S)<-x
wunif1S<-as.dist(as(wunif1S,"matrix"))
map2unifS<-subset(map2unif, Vegetation=="S")


#run adonis in Vegan 
#fungi 
adonis(bray1~ Vegetation,data=map2b,permutations=999, strata = map2b$Location)
adonis(bray1~ temp + precip  + pH + VWC+ soilCN + root.microb3  ,data=map2b,permutations=999, strata = map2b$Vegetation)
adonis(bray1S~ PC1 + PC2 ,data=NMDS_bray_shrub, permutations=999) 
#bacteria
adonis(wunif1~Vegetation, data=map2unif,permutations=999, strata = map2unif$Location)
adonis(wunif1~ temp + precip + pH + soilCN + VWC+ root.microb3, data=map2unif,permutations=999, strata = map2unif$Vegetation)
adonis(wunif1S~ PC1 + PC2 ,data=NMDS_wunif_shrub, permutations=999) 

cols<-brewer.pal(n = 8, name = 'Greens')

#Figures 
library(RColorBrewer)
n <- 60
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
col_vector<-col_vector[c(26,25,24,22,21,20,19,18,17,16,13,11,10,9,15,12,14)]
cols<-col_vector[c(1,3,4,5,8,10,11,12,13,14,15,16,17)]
cols2<-cols[c(13,12,11,10,9,8,7,6,5,4,3,2,1)]


#Fig 2 c, d
c<-ggplot(data=NMDS_bray_coords,aes(x=MDS1,y=-MDS2, colour=Location)) + 
  stat_ellipse(level = 0.85)+
  scale_color_manual(values = cols2, name="Site") + theme_classic()+ 
  geom_point(aes(MDS1,y=-MDS2, colour=Location, shape=Vegetation))

d<-ggplot(data=NMDS_wunif_coords,aes(x=MDS1,y=-MDS2, colour=Location)) + 
  stat_ellipse(level = 0.85)+
  scale_color_manual(values = cols2, name="Site") + theme_classic()+ 
  geom_point(aes(MDS1,y=-MDS2, colour=Location, shape=Vegetation))

grid.arrange(c, d, ncol=2)

#Fig 5
#Fungi
p1<-qplot(data=NMDS_bray_coords,x=MDS2,y=MDS3,colour=Vegetation)+ scale_color_brewer(palette="Paired",  labels = c("Herb","Woody")) +stat_ellipse() +theme_classic()
p2<-qplot(data=NMDS_bray_shrub,x=MDS1,y=MDS2,colour=root.microb3) +stat_ellipse()+ scale_colour_brewer(palette = 'Set2', type='qual', name="Root symbiont")+ theme_classic()
p3<-qplot(data=NMDS_bray_coords,x=MDS1,y=MDS2,colour=precip) + theme_classic() +
  scale_color_continuous(type = 'viridis', name="MAP")

grid.arrange( p1,p2,p3,ncol=3)

#Bacteria
p5<-qplot(data=NMDS_wunif_coords,x=MDS2,y=MDS3,colour=Vegetation) + scale_color_brewer(palette="Paired",  labels = c("Herb","Woody")) +stat_ellipse() +theme_classic()
p6<-qplot(data=(NMDS_wunif_shrub),x=MDS1,y=MDS2,colour=root.microb3)+stat_ellipse()+ scale_colour_brewer(palette = 'Set2', type='qual', name="Root symbiont")+ theme_classic()
p7<-qplot(data=(NMDS_wunif_coords),x=MDS1,y=MDS2,colour=pH)+ theme_classic()+
  scale_color_continuous(type = 'viridis', name="pH")

grid.arrange(p5, p6, p7,ncol=3)


