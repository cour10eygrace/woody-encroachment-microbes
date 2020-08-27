#load libraries
library(MASS)
library(ggplot2)
library(dplyr)
library(reshape2)
library(lme4)
library(lmerTest)
library(devtools)
library(vegan)
library(ape)
library(stats)

#load data & organize 
data<-read.csv("Alldata.csv")
data$VWC<-as.numeric(data$VWC....)
adiv<-read.csv("ITS_alphadiv_bysample.csv")
adiv_bac<-read.csv("16S_alphadiv_bysample.csv")
adiv_all<-merge(adiv, adiv_bac, by="Tube", all.x = T, all.y=T)
alldata<-merge(data, adiv_all, by="Tube")
alldata$soilCN<-alldata$soilC2/alldata$soilN2

str(alldata)
alldata$X.1<-NULL
alldata$X<-NULL
alldata$X.2<-NULL
alldata$Notes<-NULL
alldata$root.microb<-NULL
alldata$root.microb2<-NULL
alldata$VWC....<-NULL
alldata$VWCplot<-NULL
alldata$pHplot<-NULL
alldata$soilCplot<-NULL
alldata$Biome2<-NULL
alldata$Biome3<-NULL

str(alldata)
shrubs<-subset(alldata,Vegetation=="S" )
herbs<-subset(alldata,Vegetation=="NS" )
shrubs<-subset(shrubs,soilCN!="Inf")
shrubs<-subset(shrubs,root.microb3!=" ")


#ORDINATE LEAF TRAITS----
#infill missing data-China and some NAs

library(mice)
shrubsfill<-mice(shrubs[, c(18:31)], m=100)
shrubs[42, 19]<-mean(as.numeric(shrubsfill$imp$LDMC..g.g.[1,]))
shrubs[45, 19]<-mean(as.numeric(shrubsfill$imp$LDMC..g.g.[2,]))
shrubs[107, 19]<-mean(as.numeric(shrubsfill$imp$LDMC..g.g.[3,]))
shrubs[108, 19]<-mean(as.numeric(shrubsfill$imp$LDMC..g.g.[4,]))
shrubs[109, 19]<-mean(as.numeric(shrubsfill$imp$LDMC..g.g.[5,]))
shrubs[110, 19]<-mean(as.numeric(shrubsfill$imp$LDMC..g.g.[6,]))
shrubs[111, 19]<-mean(as.numeric(shrubsfill$imp$LDMC..g.g.[7,]))
shrubs[112, 19]<-mean(as.numeric(shrubsfill$imp$LDMC..g.g.[8,]))
shrubs[113, 19]<-mean(as.numeric(shrubsfill$imp$LDMC..g.g.[9,]))
shrubs[114, 19]<-mean(as.numeric(shrubsfill$imp$LDMC..g.g.[10,]))
shrubs[115, 19]<-mean(as.numeric(shrubsfill$imp$LDMC..g.g.[11,]))
shrubs[116, 19]<-mean(as.numeric(shrubsfill$imp$LDMC..g.g.[12,]))

shrubs[42, 22]<-mean(as.numeric(shrubsfill$imp$SLA..cm.2.g.[1,]))
shrubs[45, 22]<-mean(as.numeric(shrubsfill$imp$SLA..cm.2.g.[2,]))

shrubs[10, 24]<-mean(as.numeric(shrubsfill$imp$X15N[1,]))
shrubs[42, 24]<-mean(as.numeric(shrubsfill$imp$X15N[2,]))
shrubs[45, 24]<-mean(as.numeric(shrubsfill$imp$X15N[3,]))
shrubs[77, 24]<-mean(as.numeric(shrubsfill$imp$X15N[4,]))
shrubs[78, 24]<-mean(as.numeric(shrubsfill$imp$X15N[5,]))
shrubs[107, 24]<-mean(as.numeric(shrubsfill$imp$X15N[6,]))
shrubs[108, 24]<-mean(as.numeric(shrubsfill$imp$X15N[7,]))
shrubs[109, 24]<-mean(as.numeric(shrubsfill$imp$X15N[8,]))
shrubs[110, 24]<-mean(as.numeric(shrubsfill$imp$X15N[9,]))
shrubs[111, 24]<-mean(as.numeric(shrubsfill$imp$X15N[10,]))
shrubs[112, 24]<-mean(as.numeric(shrubsfill$imp$X15N[11,]))
shrubs[113, 24]<-mean(as.numeric(shrubsfill$imp$X15N[12,]))
shrubs[114, 24]<-mean(as.numeric(shrubsfill$imp$X15N[13,]))
shrubs[115, 24]<-mean(as.numeric(shrubsfill$imp$X15N[14,]))
shrubs[116, 24]<-mean(as.numeric(shrubsfill$imp$X15N[15,]))

shrubs[10, 26]<-mean(as.numeric(shrubsfill$imp$X13C[1,]))
shrubs[42, 26]<-mean(as.numeric(shrubsfill$imp$X13C[2,]))
shrubs[45, 26]<-mean(as.numeric(shrubsfill$imp$X13C[3,]))
shrubs[77, 26]<-mean(as.numeric(shrubsfill$imp$X13C[4,]))
shrubs[78, 26]<-mean(as.numeric(shrubsfill$imp$X13C[5,]))
shrubs[107, 26]<-mean(as.numeric(shrubsfill$imp$X13C[6,]))
shrubs[108, 26]<-mean(as.numeric(shrubsfill$imp$X13C[7,]))
shrubs[109, 26]<-mean(as.numeric(shrubsfill$imp$X13C[8,]))
shrubs[110, 26]<-mean(as.numeric(shrubsfill$imp$X13C[9,]))
shrubs[111, 26]<-mean(as.numeric(shrubsfill$imp$X13C[10,]))
shrubs[112, 26]<-mean(as.numeric(shrubsfill$imp$X13C[11,]))
shrubs[113, 26]<-mean(as.numeric(shrubsfill$imp$X13C[12,]))
shrubs[114, 26]<-mean(as.numeric(shrubsfill$imp$X13C[13,]))
shrubs[115, 26]<-mean(as.numeric(shrubsfill$imp$X13C[14,]))
shrubs[116, 26]<-mean(as.numeric(shrubsfill$imp$X13C[15,]))

shrubs[10, 28]<-mean(as.numeric(shrubsfill$imp$Wt..N.[1,]))
shrubs[42, 28]<-mean(as.numeric(shrubsfill$imp$Wt..N.[2,]))
shrubs[45, 28]<-mean(as.numeric(shrubsfill$imp$Wt..N.[3,]))
shrubs[77, 28]<-mean(as.numeric(shrubsfill$imp$Wt..N.[4,]))
shrubs[78, 28]<-mean(as.numeric(shrubsfill$imp$Wt..N.[5,]))
shrubs[107, 28]<-mean(as.numeric(shrubsfill$imp$Wt..N.[6,]))
shrubs[108, 28]<-mean(as.numeric(shrubsfill$imp$Wt..N.[7,]))
shrubs[109, 28]<-mean(as.numeric(shrubsfill$imp$Wt..N.[8,]))
shrubs[110, 28]<-mean(as.numeric(shrubsfill$imp$Wt..N.[9,]))
shrubs[111, 28]<-mean(as.numeric(shrubsfill$imp$Wt..N.[10,]))
shrubs[112, 28]<-mean(as.numeric(shrubsfill$imp$Wt..N.[11,]))
shrubs[113, 28]<-mean(as.numeric(shrubsfill$imp$Wt..N.[12,]))
shrubs[114, 28]<-mean(as.numeric(shrubsfill$imp$Wt..N.[13,]))
shrubs[115, 28]<-mean(as.numeric(shrubsfill$imp$Wt..N.[14,]))
shrubs[116, 28]<-mean(as.numeric(shrubsfill$imp$Wt..N.[15,]))

shrubs[10, 30]<-mean(as.numeric(shrubsfill$imp$Wt..C.[1,]))
shrubs[42, 30]<-mean(as.numeric(shrubsfill$imp$Wt..C.[2,]))
shrubs[45, 30]<-mean(as.numeric(shrubsfill$imp$Wt..C.[3,]))
shrubs[77, 30]<-mean(as.numeric(shrubsfill$imp$Wt..C.[4,]))
shrubs[78, 30]<-mean(as.numeric(shrubsfill$imp$Wt..C.[5,]))
shrubs[107, 30]<-mean(as.numeric(shrubsfill$imp$Wt..C.[6,]))
shrubs[108, 30]<-mean(as.numeric(shrubsfill$imp$Wt..C.[7,]))
shrubs[109, 30]<-mean(as.numeric(shrubsfill$imp$Wt..C.[8,]))
shrubs[110, 30]<-mean(as.numeric(shrubsfill$imp$Wt..C.[9,]))
shrubs[111, 30]<-mean(as.numeric(shrubsfill$imp$Wt..C.[10,]))
shrubs[112, 30]<-mean(as.numeric(shrubsfill$imp$Wt..C.[11,]))
shrubs[113, 30]<-mean(as.numeric(shrubsfill$imp$Wt..C.[12,]))
shrubs[114, 30]<-mean(as.numeric(shrubsfill$imp$Wt..C.[13,]))
shrubs[115, 30]<-mean(as.numeric(shrubsfill$imp$Wt..C.[14,]))
shrubs[116, 30]<-mean(as.numeric(shrubsfill$imp$Wt..C.[15,]))

traits<-shrubs[, c(1, 7, 5, 15, 37, 19,22,24,26,28,30)]
row.names(traits)<-traits$SampleID
#write.csv(traits, "traits.csv")


#remove negatives and log for vegdist
traits$X15N<-log(traits$X15N+9)
traits$X13C<-log(traits$X13C+32)
traits$LDMC..g.g.<-log(traits$LDMC..g.g.)+4
traits$SLA..cm.2.g.<-log(traits$SLA..cm.2.g.)+1
traits$Wt..N.<-log(traits$Wt..N.)
traits$Wt..C.<-log(traits$Wt..C.)
traits$Wt..C.<-log(traits$Wt..C.)


#vegan
#normalize data (above)
#prcomp(scale=T)
head(traits)
traits2 = traits[ ,c(6:11)]
head(traits2)

##markos code update
PCA_trt = prcomp(traits2, scale.=TRUE)
biplot(PCA_trt)

##this gives loadings
PCA_trt
## PC1 loads with LDMC, SLA 15N and %N
## PC2 loads with  %C (and 13C if only 2 axes)
## PC3 loads with 13C 

##this give amount explained by each axis
summary(PCA_trt)
PCA_traitsdf<-as.data.frame(PCA_trt$x)
PCA_traitsdf<-PCA_traitsdf[,c(1:3)]

#merge back with metadata
PCA_traits_full<-merge(PCA_traitsdf, traits, by="row.names")
var_explained<-PCA_trt$sdev^2/sum(PCA_trt$sdev^2)
var_explained<-var_explained[1:5]
trait.scrs <- as.data.frame(PCA_trt$rotation)
trait.scrs <- cbind(trait.scrs, traits = c("LDMC", "SLA", "15N", "13C", "TotalN", "TotalC"))

#Fig S3
library(RColorBrewer)
n <- 60
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
col_vector<-col_vector[c(26,25,24,22,21,20,19,18,17,16,13,11,10,9,15,12,14)]
cols<-col_vector[c(1,3,4,5,8,10,11,12,13,14,15,16,17)]


ggplot(PCA_traits_full)+
  geom_point(aes(x=PC1,y=PC2, colour=Location, shape=root.microb3), size=2) + 
 scale_color_manual(values = cols[1:13])+ theme_classic()+
  labs(x=paste0("PC1: ",round(var_explained[1]*100,1),"%"),
       y=paste0("PC2: ",round(var_explained[2]*100,1),"%")) +
  geom_segment(data = trait.scrs, aes(x = 0, xend =PC1 , y = 0, yend = PC2),
               arrow = arrow(length = unit(0.25, "cm")), colour = "grey") +
  geom_text(data = trait.scrs, aes(x =PC1, y = PC2, label = traits),size = 4)

#write.csv(PCA_traits_full, "PCA_traits_euclidian.csv")

#Fig S5
#test for differences in leaf N among symbiont types
traits<-shrubs[, c(1, 7, 5, 15, 37, 19,22,24,26,28,30)]

nmod<-aov(log(Wt..N.)~root.microb3, traits)
TukeyHSD(nmod)

#plot 
plot(log(Wt..N.)~root.microb3, traits, ylab="leaf N %", xlab=" ")

