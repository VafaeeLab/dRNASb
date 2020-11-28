

# skip all inner join step for sig genes by just reading the end  --------

x <-fread("Sig_oneTimepoint_oneZerocount_logFC1.csv")
x <- as.data.frame(x)
#rownames(x)<- x$Attributes
#x<- x[,-1]


# jump Fuzzy c-means clustering  ------------------------------------------


# -------------------------------------------------------------------------


# find the gene ID if they are significant atleast in one time point --------
sig_2h<- read.csv("newdt2h.csv", header = TRUE)
sig_4h<- read.csv("newdt4h.csv", header = TRUE)
sig_8h<- read.csv("newdt8h.csv", header = TRUE)
sig_16h<- read.csv("newdt16h.csv", header = TRUE)
sig_24h<- read.csv("newdt24h.csv", header = TRUE)


# Full significant data counts with LogFC>1 to reduce count number --------

d2h_sub <- subset(sig_2h, abs(sig_2h$LogFC2h) > 1)
d4h_sub <- subset(sig_4h, abs(sig_4h$LogFC4h) > 1)
d8h_sub <- subset(sig_8h, abs(sig_8h$LogFC8h) > 1)
d16h_sub <- subset(sig_16h, abs(sig_16h$LogFC16h) > 1)
d24h_sub <- subset(sig_24h, abs(sig_24h$LogFC24h) > 1)



all.ides <-unique(
  union(union(union(union(d2h_sub$Attributes, d4h_sub$Attributes),d8h_sub$Attributes), d16h_sub$Attributes), d24h_sub$Attributes))

df_all.ides <- as.data.frame(all.ides)
names(df_all.ides)[1]<-paste("Attributes")

# Full significant data counts --------------------------------------------

#join<-full_join(sig_2h,sig_4h, by="Attributes")
#join2<-full_join(join,sig_8h, by="Attributes")
#join3<-full_join(join2,sig_16h, by="Attributes")
#join4<-full_join(join3,sig_24h, by="Attributes")

#uniq_genes <- unique(join4)

#sig_genes <- as.data.frame(join4[,1])
#names(sig_genes)[1]<-paste("Attributes")

N.M.D <- read.csv("Mean normalized data.csv", header = TRUE)
names(N.M.D)[1]<-paste("Attributes")

#inner_join sig genes by total genes to get normalised read count for all 6 time-points
x <- inner_join(N.M.D,df_all.ides, by="Attributes")


write.csv(x, file= "All.ides.csv", row.names = FALSE)
 



# Fuzzy c-means clustering  -----------------------------------------------
library(Mfuzz)
library(edgeR)
library(SummarizedExperiment)

y.dat<- x
rownames(y.dat)<- y.dat$Attributes
y.dat<- y.dat [,c(-1)]

y.dat <- as.matrix(y.dat)
#z_var <- apply(y.dat, 1, var)
#z_mean <- apply(y.dat,1,mean)
#y.dat <- y.dat[which(z_var>1 & z_mean>1), 1:6]

timepoint <- c(0,2,4,8,16,24)
y.dat <- rbind(timepoint, y.dat)
rownames(y.dat)[1]<- "time"

tmp<- tempfile()
write.table(y.dat, file = tmp, sep = '\t',quote = FALSE, col.names = NA)
# A expression matrix stored as a table (in a defined format) is r 
eset.data <- table2eset(tmp)

# standardise  ------------------------------------------------------------
data.z <- standardise(eset.data)
m1<-mestimate(data.z)
m1



# c-selection of NA for optimal cluster number ----------------------------
#cselection(data.z, m=m1,crange = seq(2,32,1), repeats = 5, visu = TRUE )

# Apply mfuzz function after optimal range dimensions ---------------------
Dmin(data.z, m=m1, crange=seq(2,22,1), repeats = 3, visu = TRUE)

clust=8
c<-mfuzz(data.z,c=clust, m=m1)
# mfuzz.plot(data.z,cl=c,mfrow=c(3,3),min.mem=0.3,time.labels=c(0,2,4,8,16,24),new.window=FALSE)

mfuzz.plot2(data.z,cl=c,mfrow=c(3,3),min.mem=0.3,
            ylim.set=c(0,0), xlab="Time",ylab="Expression changes",x11=TRUE,
            ax.col="black",bg = "white",col.axis="black",col.lab="white",
            col.main="black",col.sub="blue",col="blue",cex.main=1.3,cex.lab=1.1, 
            centre=TRUE,
            centre.col="black",centre.lwd=2,
            Xwidth=9,Xheight=9,single=FALSE)


# Cluster analysis(correlation between cluster centroids) -----------------
# (no two clusters should exhbit a correlation greater than 0.85)
fd<-data.frame(cor(t(c[[1]])))
colnames(fd)<- c("Cluster1", "Cluster2","Cluster3","Cluster4","Cluster5", "Cluster6", "Cluster7","Cluster8")
  rownames(fd)<- c("Cluster1", "Cluster2","Cluster3","Cluster4","Cluster5", "Cluster6", "Cluster7","Cluster8")
  write.csv(fd, file = "Correlation between Cluster centroid at 0.3 (post hoc test).csv", row.names = TRUE)
  

# Extract list of genes ---------------------------------------------------
acore<-acore(data.z,c,min.acore = 0.3)


acore_list_0.3<-do.call(rbind,lapply(seq_along(acore), function(i){data.frame(CLUSTER=i, acore[[i]])}))

#colnames(acore_list_0.3)[2]<-"Attributes"
write.csv(acore_list_0.3, file = "New Acore_list_0.3.csv", row.names = FALSE)


# OR ----------------------------------------------------------------------

acore_0.3<-acore(data.z,c,min.acore = 0.3)
temp <- do.call("rbind", lapply(acore_0.3, FUN = function(x){
  return(paste0(as.character(x$NAME), collapse = ","))
}))



v<- as.data.frame(temp)
write.csv(v, file = "Cluster acore FINAL for GENE Name Cluster innerjoin PPI Attributes.csv",row.names = FALSE)
#test<-str_split(test,",", n=Inf)

#need to convert text to colimn in excel 
test<- read.csv("Cluster acore.csv", header=TRUE)

#test<-t(test)
test<- as.data.frame(test)
colnames(test)<- c("Name")
rownames(test)= c("Cluster1", "Cluster2","Cluster3","Cluster4","Cluster5", "Cluster6", "Cluster7","Cluster8")

write.csv(test, file = "Temporal cluster genelist 0.3.csv",row.names = TRUE)


remove(x,sig_2h,sig_4h,sig_8h,sig_16h,sig_24h,d2h_sub,d4h_sub,d8h_sub,d16h_sub,d24h_sub,all.ides,df_all.ides,N.M.D,y.dat,tmp,eset.data,data.z,m1,c,fd,v,test,timepoint)


# Ad Go term --------------------------------------------------------------


join<-full_join(dat,acore_list_0.3, by="Uniprot.ID")

# Cluster enrichment (using 'ClueR') --------------------------------------
install.packages("ClueR")
library(ClueR)

# Read clusters and Preprocess



#make list of lists for cluster data
    newClusters.dat <- as.data.frame(fread("Temporal cluster genelist 0.3.csv", header = T))
    newClusters <- sapply(newClusters.dat[,2], function(x){
      return(str_split(x, ","))
    })
    names(newClusters) <- as.character(str_replace(newClusters.dat[,1], "Cluster",""))



#make list of lists for annotation
   dat <- as.data.frame(fread("Annotation.enrichment.csv"))
    GO.dat <- dat[,c(2,4,6)] %>% group_by(Gene.Ontology.Terms) %>% 
      summarise_each(funs(paste0(., collapse = "\t"))) %>% 
      as.data.frame()
    fwrite(GO.dat[,c(1,2)], file = "GO_uniprots_SL1344.csv", sep="\t", col.names = F,  quote = F) 
    
    
    GO <- sapply(GO.dat[,2], function(x){
      return(str_split(x, "\t"))
    })
    names(GO) <- as.character(str_replace(GO.dat[,1], "Gene.Ontology.terms",""))


clustEnrichment(clustObj, annotation= GO, effectiveSize=c(5, 100), pvalueCutoff=0.05)

cluster_enrich <- clustEnrichment(newClusters, GO, effectiveSize=8, pvalueCutoff = 0.05)




# OR this -----------------------------------------------------------------

c=fd

clustObj <- cmeans(c, centers=8, iter.max=10, m=1.77)
ce <-clustEnrichment(clustObj, annotation=GO, effectiveSize=c(5, 100), pvalueCutoff=0.05)





# Merge 8 clusters into annotation ----------------------------------------
colnames(acore_list_0.3)[2]<- c("Uniprot.ID")
join<-full_join(dat,acore_list_0.3, by="Uniprot.ID")
colnames(join)[7]<- c("Cluster.number")
colnames(join)[8]<- c("Cluster.Mem.ship")
write.csv(join, file = "Merge 8 clusters in annotation.csv", row.names = FALSE)





# Innerjoin wih PPI attributes --------------------------------------------



PPI_attributes <- read.csv (file= "PPI_Attributes.csv", header= TRUE)
colnames(PPI_attributes)[2] <-c("NAME")

GENE_NAME_CLUSTERS<- inner_join(acore_list_0.3,PPI_attributes, by="NAME")




# Innerjoin wih GO term annotation ----------------------------------------

dat <- as.data.frame(fread("Annotation.enrichment.csv"))

colnames(dat)[2] <-c("NAME")

GO_TERM_CLUSTERS<- inner_join(acore_list_0.3,dat, by="NAME")



