## load libraries
library(permute)
library(indicspecies)

#FUNGI
#load otu table, .csv, remove #created from biom file first line, and remove #otu_id cell
otu<- read.csv("CCGlobalProject.otu_table.taxonomy.csv")
rownames(otu)<-otu$ï..OTU.ID
otu$ï..OTU.ID<-NULL
otu$Taxonomy<-NULL 
otu$ITS.cc.NTC<-NULL ##no template control

#load mapping file, remove 'sample ID' cell before
map = read.csv("Alldata.csv" , header=T, row.names=43)
map = map[ ,c(4,8,10,18)]

#make dataframe for taxa
taxa<- read.csv("CCGlobalProject.taxonomy.fix.csv", row.names = 1)
taxa<-subset(taxa, Kingdom=="Fungi")
taxa<-as.data.frame(taxa)

#transpose OTU table so that samples are rows and OTU IDs are columns
OTU.transp<-t(otu)
OTU.transp<-as.data.frame(OTU.transp)
str(OTU.transp)
head(OTU.transp)

#match order of mapping file with order of samples in OTU table
p<-match(row.names(OTU.transp),row.names(map))
map2<-map[p , ]
map2<-na.omit(map2)

#check
nrow(map2)==nrow(OTU.transp)

#create factor name to call for Uniqueclass, i.e., unique tree by site combination
Veg<-as.factor(map2$Vegetation)

#indicator species analysis by unique identifier
indval_veg<-multipatt(x=as.data.frame(OTU.transp),
                      cluster=Veg,
                      func="IndVal.g",
                      duleg=TRUE,
                      control=how(nperm=999))

#To get summary information of OTUs that are indicators 
summary.multipatt(indval_veg, minstat = 0.3)#stat=0.3 means taxa is found in 30% or more of the samples in that group  

itsveg<-as.data.frame(indval_veg$sign)%>%
  #subset(., p.value<0.05&stat>0.1)%>%
  merge(.,taxa,by='row.names')%>%
  subset(., Species!="")


#BACTERIA
#load otu table, .csv, remove #created from biom file first line, and remove #otu_id cell
otu_b<-read.csv("feature-table_Cfilt.csv")  
rownames(otu_b)<-otu_b$Feature.ID
otu_b$Feature.ID<-NULL
otu_b$ITS.cc.NTC1<-NULL ##no template control
otu_b$ITS.cc.NTC2<-NULL ##no template control

#make dataframe for taxa
taxa_b<-read.csv("taxonomy.csv")  
taxa_b<-subset(taxa_b, Kingdom=="Bacteria")
taxa_b<-as.data.frame(taxa_b)
rownames(taxa_b)<-taxa_b$ï..
taxa_b$ï..<-NULL

#transpose OTU table so that samples are rows and OTU IDs are columns
OTU.transp_b<-t(otu_b)
OTU.transp_b<-as.data.frame(OTU.transp_b)
str(OTU.transp_b)
head(OTU.transp_b)

#match order of mapping file with order of samples in OTU table
pb<-match(row.names(OTU.transp_b),row.names(map))
map2b<-map[pb , ]
map2b<-na.omit(map2b)

#check
nrow(map2b)==nrow(OTU.transp_b)

#create factor name to call for Uniqueclass, i.e., unique tree by site combination
Veg<-as.factor(map2b$Vegetation)

#indicator species analysis by unique identifier
indval_veg_b<-multipatt(x=as.data.frame(OTU.transp_b),
                      cluster=Veg,
                      func="IndVal.g",
                      duleg=TRUE,
                      control=how(nperm=999))
#To get summary information of OTUs that are indicators 
summary.multipatt(indval_veg_b, minstat = 0.3)#stat=0.3 means taxa is found in 30% or more of the samples in that group  

X16sveg<-as.data.frame(indval_veg_b$sign)%>%
  subset(., p.value<0.05&stat>0.2)%>%
  merge(.,taxa_b,by='row.names')%>%
  subset(., Genus!="")


