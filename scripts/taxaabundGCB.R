#TAXA---- 
library(reshape2)
library(lme4)
#FUNGI 
#read in otu table 
otus <- read.csv("CCGlobalProject.otu_table.taxonomy.csv")
otus$OTUID<-otus$ï..OTU.ID
otus$ï..OTU.ID<-NULL
otus$Taxonomy<-NULL
otus$ITS.cc.NTC<-NULL

#read in updated taxonomy
taxmat <- read.csv("CCGlobalProject.taxonomy.fix.csv")
taxmat<-subset(taxmat, Kingdom=="Fungi")
taxmat$OTUID<-taxmat$ï..
taxmat$ï..<-NULL

#merge
otu_tax<-merge(taxmat, otus, by="OTUID")

#read in metadata 
map = read.csv("Alldata.csv" , header=T, row.names=43)
map = map[ ,c(42,4,18)]
map$Sample.ID<-row.names(map)

#make otu table long 
otul <- melt(otu_tax)
names(otul)[names(otul)=="variable"] <-"Sample.ID"
names(otul)[names(otul)=="value"] <- "itsreads"
#merge with metadata 
otu_all<-merge(otul,map, by="Sample.ID")

#filter out zero reads 
itsl_reads<-subset(otu_all, itsreads>0)

#subset by phyla
asco<-subset(itsl_reads, Phylum=="Ascomycota")
basidio<-subset(itsl_reads, Phylum=="Basidiomycota")
blast<-subset(itsl_reads, Phylum=="Blastocladiomycota")#not in both 
chyt<-subset(itsl_reads, Phylum=="Chytridiomycota")
glom<-subset(itsl_reads, Phylum=="Glomeromycota")
zyg<-subset(itsl_reads, Phylum=="Zygomycota")#not in both
mort<-subset(itsl_reads, Phylum=="Mortierellomycota")

##run mixed effects models 
hist(asco$itsreads)# not normal 
hist(log(asco$itsreads))#still not normal 
asco_veg<-glmer((log(itsreads)+1)~Vegetation + (1|Location), asco, family = Gamma) 
summary(asco_veg)#S<NS (p=0.08)

hist(basidio$itsreads)#not normal 
hist(log(basidio$itsreads))#still not normal 
basidio_veg<-glmer(log(itsreads)+1~Vegetation + (1|Location), basidio, family = Gamma) 
summary(basidio_veg)#S<NS
plot(log(itsreads)~Vegetation, basidio)

hist(chyt$itsreads)#not normal
hist(log(chyt$itsreads))#still not normal 
chyt_veg<-glmer((log(itsreads)+1)~Vegetation + (1|Location), chyt, family=Gamma)
summary(chyt_veg)#NS

hist(glom$itsreads)#not normal
hist(log(glom$itsreads))#still not normal 
glom_veg<-glmer((log(itsreads)+1)~Vegetation + (1|Location), glom, family=Gamma)
summary(glom_veg)#NS

hist(mort$itsreads)#not normal
hist(log(mort$itsreads))#still not normal 
mort_veg<-glmer((log(itsreads)+1)~Vegetation + (1|Location), mort, family=Gamma)
summary(mort_veg)#NS

#TAXA---- 
#read in otu table
otus_b <- read.csv("feature-table_Cfilt.csv")
otus_b$ITS.cc.NTC1<-NULL
otus_b$ITS.cc.NTC2<-NULL
#read in taxonomy
taxmat_b<- read.csv("taxonomy.csv")
taxmat_b<-subset(taxmat_b, Kingdom=="Bacteria")
taxmat_b$Feature.ID<-taxmat_b$ï..
taxmat_b$ï..<-NULL
otu_tax_b<-merge(otus_b, taxmat_b, by="Feature.ID")
#make otu table long 
otulb <- melt(otu_tax_b)
names(otulb)[names(otulb)=="variable"] <- "Sample.ID"
names(otulb)[names(otulb)=="value"] <- "X16Sreads"
#merge with metadata
otu_all_b<-merge(otulb,map, by="Sample.ID")
#filter out zero reads 
otulb_reads<-subset(otu_all_b, X16Sreads>0)

#subset by phylum 
acido<-subset(otulb_reads, Phylum=="Acidobacteria")#NS>S
actino<-subset(otulb_reads, Phylum=="Actinobacteria")#NS>S
proteo<-subset(otulb_reads, Phylum=="Proteobacteria")#NS>S
firm<-subset(otulb_reads, Phylum=="Firmicutes")#S>NS
very<-subset(otulb_reads, Phylum=="Verrucomicrobia")#NS>S
bact<-subset(otulb_reads, Phylum=="Bacteroidetes")#S>NS
planct<-subset(otulb_reads, Phylum=="Planctomycetes")#NS>S
chlor<-subset(otulb_reads, Phylum=="Chloroflexi")#not sig


hist(log(acido$X16Sreads))#normal 
acido_veg<-lmer(log(X16Sreads)~Vegetation + (1|Location), acido) 
summary(acido_veg)#sig
plot(log(X16Sreads)~Vegetation, acido)

hist(log(actino$X16Sreads))#normal 
actino_veg<-lmer(log(X16Sreads)~Vegetation + (1|Location), actino) 
summary(actino_veg)#sig
plot(log(X16Sreads)~Vegetation, actino)

hist(log(proteo$X16Sreads))#normal 
proteo_veg<-lmer(log(X16Sreads)~Vegetation + (1|Location), proteo) 
summary(proteo_veg)#sig
plot(log(X16Sreads)~Vegetation, proteo)

hist(log(firm$X16Sreads))#normal 
firm_veg<-lmer(log(X16Sreads)~Vegetation + (1|Location), firm) 
summary(firm_veg)#sig (p=0.075)
plot(log(X16Sreads)~Vegetation, firm)

hist(log(very$X16Sreads))#normal 
ver_veg<-lmer(log(X16Sreads)~Vegetation + (1|Location), very) 
summary(ver_veg)#sig
plot(log(X16Sreads)~Vegetation, very)

hist(log(chlor$X16Sreads))#normal 
chlor_veg<-lmer(log(X16Sreads)~Vegetation + (1|Location), chlor) 
summary(chlor_veg)#NS

hist(log(bact$X16Sreads))#normal 
bact_veg<-lmer(log(X16Sreads)~Vegetation + (1|Location), bact) 
summary(bact_veg) #sig (p=0.056)
plot(log(X16Sreads)~Vegetation, bact)

hist(log(planct$X16Sreads))#normal 
planct_veg<-lmer(log(X16Sreads)~Vegetation + (1|Location), planct) 
summary(planct_veg)#sig
plot(log(X16Sreads)~Vegetation, planct)




