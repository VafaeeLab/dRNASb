# Cluster of Gene ---------------------------------------------------------

All.ides<-read.csv("All.ides.csv")
Mean<-read.csv("Mean of Normalise.csv")
All.mean<-merge(Mean,All.ides, by="Attributes")
rownames(All.mean)<-All.mean$Attributes
All.mean<-as.matrix((All.mean[,-c(1)]))
colnames(All.mean)<-c("0","2","4","8","16","24")
timepoint<-c(0,2,4,8,16,24)
All.mean.time<-rbind(timepoint,All.mean)
row.names(All.mean.time)[1]<-"time"
All.mean.time<-data.frame(All.mean.time)
tmp<-tempfile()
write.table(All.mean.time,file=tmp, sep='\t', quote=FALSE, col.names= NA)
data <- table2eset(tmp)
data.s<-standardise(data)
m1<-mestimate(data.s)
m1
cluster<-Dmin(data.s,m=m1,crange = seq(2,22,1), repeats = 3, visu = TRUE)
cluster<-data.frame(cluster)
clust=8  # m.ship= 1/8= 0.125
c<-mfuzz(data.s,c=clust, m=m1)
mfuzz.plot(data.s,cl=c,mfrow=c(4,4),min.mem=0.5,time.labels=c(0,2,4,8,16,24),new.window=FALSE)
fd<-data.frame(cor(t(c[[1]])))
acore<-acore(data.s,c,min.acore = 0.5)
acore_list<-do.call(rbind,lapply(seq_along(acore), function(i){data.frame(CLUSTER=i, acore[[i]])}))
colnames(acore_list)[2]<-"Attributes"
genelist<- acore(data.s,cl=c,min.acore=0.7)
temp <- do.call("rbind", lapply(genelist, FUN = function(x){
  return(paste0(as.character(x$NAME), collapse = ","))
}))

Cluster_list<-as.data.frame(temp)
Cluster_list<-str_split_fixed(Cluster_list$V1,",", n=Inf)
Cluster_list<-t(Cluster_list)
colnames(Cluster_list)<- c("Cluster1", "Cluster2","Cluster3","Cluster4","Cluster5", "Cluster6", "Cluster7","Cluster8")
rownames(Cluster_list)= NULL
Cluster_list<-as.data.frame(Cluster_list)


# Make list ---------------------------------------------------------------

GO.UNI<-read.csv("Bio_uni_MD.csv")
GO.UNI<-GO.UNI[,c(4,5)]
go<-gsub(" ", "", GO.UNI$Gene.Ontology.ID,fixed = TRUE)
go<-as.data.frame(go)
UNI<-cbind(go, GO.UNI)
GO.UNI<-UNI[,c(1,2)]
colnames(GO.UNI)[1]<-"Gene.Ontology.ID"
GO <- unique(GO.UNI$Gene.Ontology.ID)
Uniprot.ID <- sapply(1:length(GO), function(i) paste(gsub("[[:space:]]", "", GO.UNI[which(GO.UNI$Gene.Ontology.ID==GO[i]),]$Uniprot.ID),collapse=" "))
Uniprot.ID<-as.data.frame(Uniprot.ID)
GO_Pro_ID<-data.frame(GO.ID=unique(GO.UNI$Gene.Ontology.ID),
                      Uniprot.ID=Uniprot.ID)
write.csv(GO_Pro_ID,file ="GO_Pro_ID.csv", row.names = FALSE)

#  Make list of list  -------------------------------------------------------------------
Anno <- list()
#groupSize <- 436
GO_IDs <- as.vector(GO_Pro_ID[,1])

for (i in GO_IDs) {
  myindex <- which(GO_Pro_ID == i)
  Anno[i] <- strsplit(as.character(GO_Pro_ID[myindex, 2]), " ")
}


# Enrichment --------------------------------------------------------------


ce <- clustEnrichment(c, annotation=Anno, effectiveSize=c(2,100), pvalueCutoff=0.01)
