
# find the gene ID if they are significant atleast in one time poi --------
sig_2h<- read.csv("newdt2h.csv", header = TRUE)
sig_4h<- read.csv("newdt4h.csv", header = TRUE)
sig_8h<- read.csv("newdt8h.csv", header = TRUE)
sig_16h<- read.csv("newdt16h.csv", header = TRUE)
sig_24h<- read.csv("newdt24h.csv", header = TRUE)


#sig_4h<-data.frame(sig_4h[,c(1)])
#Attributes <-colnames(sig_4h)
#x<- merge(t1,t2, by= "Attributes")
#names(sig_4h)[1]<-paste("Attributes")
#yx<- rbind(sig_24h,sig_2h,sig_8h,sig_16h,sig_4h)
#rownames(sig_4h)<- sig_4h$Attributes
#uniq_yx <- unique(yx)



join<-full_join(sig_2h,sig_4h, by="Attributes")
join2<-full_join(join,sig_8h, by="Attributes")
join3<-full_join(join2,sig_16h, by="Attributes")
join4<-full_join(join3,sig_24h, by="Attributes")

 uniq_genes <- unique(join4)
 
sig_genes <- as.data.frame(join4[,1])
names(sig_genes)[1]<-paste("Attributes")

N.M.D <- read.csv("Mean normalized data.csv", header = TRUE)
names(N.M.D)[1]<-paste("Attributes")

#inner_join sig genes by total genes to get normalised read count for all 6 time-points
x <- inner_join(N.M.D,sig_genes, by="Attributes")

#inner_join x by z(genes extracted from cytoscape)
z <- read.csv("Significant genes from Cytoscape.csv", header = TRUE)
genes_431 <- inner_join(x,z, by="Attributes")

write.csv(x, file= "Sig_oneTimepoint_oneZerocount.csv", row.names = FALSE)


# Fuzzy c-means clustering  -----------------------------------------------
library(Mfuzz)
library(edgeR)
library(SummarizedExperiment)

y.dat<- genes_431
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

z.data <- table2eset(tmp)

# standardise  ------------------------------------------------------------
data.z <- standardise(z.data)
class(data.z)
m1 <-mestimate(data.z)
m1 ### 1.730301, 1.95498

Dmin(data.z, m=m1, crange=seq(2,22,1), repeats = 3, visu = TRUE)

clust=10
c<- mfuzz(data.z, c=clust, m=m1)

mfuzz.plot(data.z, cl=c, mfrow = c(1,1), time.labels = c(0,2,4,8,16,24), new.window = FALSE)



# Extract list of genes ---------------------------------------------------


genelist<- acore(data.z,cl=c,min.acore=0.7)

temp <- do.call("rbind", lapply(acore, FUN = function(x){
  return(paste0(as.character(x$NAME), collapse = ","))
}))


v<- as.data.frame(temp)
write.csv(v, file = "Cluster acore.csv",row.names = FALSE)
#test<-str_split(test,",", n=Inf)

#need to convert text to colimn in excel 
test<- read.csv("Cluster acore.csv", header=TRUE)


test<-t(test)
test<- as.data.frame(test)
colnames(test)<- c("Cluster1", "Cluster2","Cluster3","Cluster4","Cluster5", "Cluster6", "Cluster7","Cluster8","Cluster9","Cluster10","Cluster11", "Cluster12")
rownames(test)= NULL

write.csv(test, file = "Temporal cluster genelist 0.3.csv",row.names = FALSE)




 
