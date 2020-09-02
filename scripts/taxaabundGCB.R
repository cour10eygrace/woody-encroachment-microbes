#TAXA---- 
library(reshape2)
library(lme4)
#FUNGI 
#read in otu table 
otus <- read.csv("data/CCGlobalProject.otu_table.taxonomy.csv")
otus$OTUID<-otus$OTU_ID
otus$OTU_ID<-NULL
otus$Taxonomy<-NULL
otus$ITS.cc.NTC<-NULL

#read in updated taxonomy
taxmat <- read.csv("data/CCGlobalProject.taxonomy.fix.csv")
taxmat<-subset(taxmat, Kingdom=="Fungi")

#merge
otu_tax<-merge(taxmat, otus, by="OTUID")

#read in metadata 
map = read.csv("data/soil.data.csv" , header=T)

#make otu table long 
otul <- melt(otu_tax)
names(otul)[names(otul)=="variable"] <-"SampleID"
names(otul)[names(otul)=="value"] <- "itsreads"
#merge with metadata 
otu_all<-merge(otul,map, by="SampleID")

#filter out zero reads 
itsl_reads<-subset(otu_all, itsreads>0)

#subset by phyla
asco<-subset(itsl_reads, Phylum=="Ascomycota")
basidio<-subset(itsl_reads, Phylum=="Basidiomycota")
chyt<-subset(itsl_reads, Phylum=="Chytridiomycota")
glom<-subset(itsl_reads, Phylum=="Glomeromycota")
mort<-subset(itsl_reads, Phylum=="Mortierellomycota")

##run mixed effects models 
hist(asco$itsreads)# not normal 
hist(log(asco$itsreads))#not normal 
#add 1 to all response vars for gamma requirement positive values 
asco_veg<-glmer((log(itsreads)+1)~veg + (1|site), asco, family = Gamma) 
summary(asco_veg)

hist(basidio$itsreads)#not normal 
hist(log(basidio$itsreads))#not normal 
basidio_veg<-glmer(log(itsreads)+1~veg + (1|site), basidio, family = Gamma) 
summary(basidio_veg)

hist(chyt$itsreads)#not normal
hist(log(chyt$itsreads))#not normal 
chyt_veg<-glmer((log(itsreads)+1)~veg + (1|site), chyt, family=Gamma)
summary(chyt_veg)

hist(glom$itsreads)#not normal
hist(log(glom$itsreads))#not normal 
glom_veg<-glmer((log(itsreads)+1)~veg + (1|site), glom, family=Gamma)
summary(glom_veg)

hist(mort$itsreads)#not normal
hist(log(mort$itsreads))#not normal 
mort_veg<-glmer((log(itsreads)+1)~veg + (1|site), mort, family=Gamma)
summary(mort_veg)

#TAXA---- 
#read in otu table
otus_b <- read.csv("data/feature-table_Cfilt.csv")
otus_b$ITS.cc.NTC1<-NULL
otus_b$ITS.cc.NTC2<-NULL
#read in taxonomy
taxmat_b<- read.csv("data/taxonomy.csv")
taxmat_b<-subset(taxmat_b, Kingdom=="Bacteria")
taxmat_b$Feature.ID<-taxmat_b$ID
taxmat_b$ID<-NULL
otu_tax_b<-merge(otus_b, taxmat_b, by="Feature.ID")
#make otu table long 
otulb <- melt(otu_tax_b)
names(otulb)[names(otulb)=="variable"] <- "SampleID"
names(otulb)[names(otulb)=="value"] <- "X16Sreads"
#merge with metadata
otu_all_b<-merge(otulb,map, by="SampleID")
#filter out zero reads 
otulb_reads<-subset(otu_all_b, X16Sreads>0)

#subset by phylum 
acido<-subset(otulb_reads, Phylum=="Acidobacteria")
actino<-subset(otulb_reads, Phylum=="Actinobacteria")
bact<-subset(otulb_reads, Phylum=="Bacteroidetes")
chlor<-subset(otulb_reads, Phylum=="Chloroflexi")
firm<-subset(otulb_reads, Phylum=="Firmicutes")
planct<-subset(otulb_reads, Phylum=="Planctomycetes")
proteo<-subset(otulb_reads, Phylum=="Proteobacteria")
very<-subset(otulb_reads, Phylum=="Verrucomicrobia")

hist(log(acido$X16Sreads))#normal 
acido_veg<-lmer(log(X16Sreads)~veg + (1|site), acido) 
summary(acido_veg)

hist(log(actino$X16Sreads))#normal 
actino_veg<-lmer(log(X16Sreads)~veg + (1|site), actino) 
summary(actino_veg)

hist(log(proteo$X16Sreads))#normal 
proteo_veg<-lmer(log(X16Sreads)~veg + (1|site), proteo) 
summary(proteo_veg)

hist(log(firm$X16Sreads))#normal 
firm_veg<-lmer(log(X16Sreads)~veg + (1|site), firm) 
summary(firm_veg)

hist(log(very$X16Sreads))#normal 
ver_veg<-lmer(log(X16Sreads)~veg + (1|site), very) 
summary(ver_veg)

hist(log(chlor$X16Sreads))#normal 
chlor_veg<-lmer(log(X16Sreads)~veg + (1|site), chlor) 
summary(chlor_veg)

hist(log(bact$X16Sreads))#normal 
bact_veg<-lmer(log(X16Sreads)~veg + (1|site), bact) 
summary(bact_veg) 

hist(log(planct$X16Sreads))#normal 
planct_veg<-lmer(log(X16Sreads)~veg + (1|site), planct) 
summary(planct_veg)




