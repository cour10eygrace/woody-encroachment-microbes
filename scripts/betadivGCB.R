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
soil.data<-read.csv("data/soil.data.csv")
soil.data$soilCN<-soil.data$soilC2/soil.data$soilN2

#fungi 
bray<-read.csv("data/bray_curtis_dm.csv")
rownames(bray)<-bray$X
bray$X<-NULL
map<-read.csv("data/bray_map.csv") #sample ID by tube #
map<-merge(map, soil.data, by="Tube")
map<-dplyr::select(map,soilCN, VWC, pH, root.microb, veg, site, Sample.ID)
rownames(map)<-map$Sample.ID
map$Sample.ID<-NULL
pmap<-match(rownames(bray), rownames(map))
map2b<-map[pmap, ]
bray1<-as.dist(as(bray,"matrix"))

ITS.clim<-read.csv("data/ITS.data.CC.csv")
rownames(ITS.clim)<-ITS.clim$X
map2b<-merge(map2b, ITS.clim, by="row.names")
row.names(map2b)<-map2b$Row.names
map2b$Row.names<-NULL

#leaf traits 
traits<-read.csv("data/PCA_traits_euclidian.csv")
traits<-dplyr::select(traits, PC1, PC2, X, SampleID)
row.names(traits)<-traits$SampleID
map2bS<-merge(map2b, traits, by="row.names")
row.names(map2bS)<-map2b$Row.names
map2bS$Row.names<-NULL

#bacteria
wunif<-read.csv("data/weighted_unifrac_dm.csv")
rownames(wunif)<-(wunif$Tube)
wunif$Tube<-NULL
colnames(wunif)<-row.names(wunif)
#remove no template control samples
wunif$NTC1<-NULL
wunif$NTC2<-NULL
wunif<-wunif[-c(216,217), ]

map2unif<-dplyr::select(soil.data,soilCN, VWC, pH, root.microb, veg, site, Tube)
map2unif<-filter(map2unif,map2unif$Tube%in% row.names(wunif), .preserve = TRUE) 
rownames(map2unif)<-map2unif$Tube
wunif1<-as.dist(as(wunif,"matrix"))

X16S.clim<-read.csv("data/X16S.data.CC.csv")
row.names(X16S.clim)<-X16S.clim$X
map2unif<-merge(map2unif, X16S.clim, by="row.names")
row.names(map2unif)<-map2unif$Tube

#leaf traits 
traits<-read.csv("data/PCA_traits_euclidian.csv")
traits<-dplyr::select(traits, PC1, PC2, Tube)
row.names(traits)<-traits$Tube  
map2unifS<-left_join(map2unif, traits)
map2unifS<-subset(map2unifS, veg=='S')
row.names(map2unifS)<-map2unifS$Tube
map2unifS$Row.names<-NULL
map2unifS$X<-NULL

#run NMDS
#fungi
NMDS_bray<-metaMDS(bray, k=4, trymax=200)#stress=0.137 
NMDS_bray_coords<-as.data.frame(NMDS_bray$points)
NMDS_bray_coords<-merge(NMDS_bray_coords, map2b, by="row.names")
NMDS_bray_shrub<-left_join(NMDS_bray_coords, map2bS)#pull out shrub only for leaf traits 
NMDS_bray_shrub<-subset(NMDS_bray_shrub, veg=="S")
bray1S<-filter(bray, row.names(bray) %in% NMDS_bray_shrub$Row.names)
bray1S<-dplyr::select(bray1S, one_of(NMDS_bray_shrub$Row.names))
bray1S<-as.dist(as(bray1S,"matrix"))


#bacteria
NMDS_wunif<-metaMDS(wunif, k=3, trymax=100)#stress=0.107
NMDS_wunif_coords<-as.data.frame(NMDS_wunif$points)
NMDS_wunif_coords<-merge(NMDS_wunif_coords, map2unif, by="row.names")
NMDS_wunif_shrub<-left_join(NMDS_wunif_coords, map2unifS) #pull out shrub only for leaf traits 
NMDS_wunif_shrub<-subset(NMDS_wunif_shrub, veg=="S")
wunif1S<-dplyr::select(wunif, one_of(NMDS_wunif_shrub$Row.names))
wunif1S<-filter(wunif1S, row.names(wunif1S) %in% NMDS_wunif_shrub$Row.names)
x<-colnames(wunif1S)
rownames(wunif1S)<-x
wunif1S<-as.dist(as(wunif1S,"matrix"))
map2unifS<-subset(map2unif, veg=="S")


#run adonis in Vegan 
#fungi 
adonis(bray1~  veg,data=map2b,permutations=999, strata = map2b$Location)
adonis(bray1~ temp + precip  + pH + VWC+ soilCN + root.microb  ,data=map2b,permutations=999, strata = map2b$Vegetation)
adonis(bray1S~ PC1 + PC2 ,data=NMDS_bray_shrub, permutations=999) 
#bacteria
adonis(wunif1~veg, data=map2unif,permutations=999, strata = map2unif$Location)
adonis(wunif1~ temp + precip + pH + soilCN + VWC+ root.microb, data=map2unif,permutations=999, strata = map2unif$Vegetation)
adonis(wunif1S~ PC1 + PC2 ,data=NMDS_wunif_shrub, permutations=999) 

