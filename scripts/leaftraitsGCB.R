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

#load leaf trait data & organize 
traits<-read.csv("data/traits.csv")

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

#write.csv(PCA_traits_full, "data/PCA_traits_euclidian.csv")

#Fig S5
#test for differences in leaf N among symbiont types
#leaf N is already logged above prior to PCA 
nmod<-aov(Wt..N.~root.microb3, traits)
TukeyHSD(nmod)

#plot 
plot(Wt..N.~root.microb3, traits, ylab="leaf N %", xlab=" ")

