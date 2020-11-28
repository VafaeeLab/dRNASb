
# Load library ------------------------------------------------------------


library(ppclust)
library(factoextra)
library(cluster)
library(fclust)

library(edgeR)
library(SummarizedExperiment)
library(Mfuzz)




#Matrix to df from rownames to columns:
d <- dt
names <- rownames(d)
rownames(d) <- NULL
data <- cbind(names,d)

#df to matrix from columns to rownames:
rownames(data)<- data$Attributes
data <- data [,-1]





# Average of 3  replicates ------------------------------------------------
#before Normalization

data_mean=data.frame(apply(array(as.matrix(data[,-1]), c(nrow(data),3, ncol(data)/3)),3, rowMeans))
data_mean=cbind(data$names, data_mean)
colnames(data_mean)<- c("Attributes", "Mean_0h","Mean_2h","Mean_4h","Mean_8h","Mean_16h", "Mean_24h")
write.csv(data_mean, file = "Normalized data mean.csv", row.names = FALSE)

write.csv(data_mean, file = "data mean.csv", row.names = FALSE)


# Read data ---------------------------------------------------------------


mean.dat <- read.csv(file = "data mean.csv", header = TRUE)
rownames(mean.dat)<- mean.dat$Attributes
mean.dat.p <- mean.dat [,-1]



# Preprocessing, filtering  -----------------------------------------------

y.dat <- as.matrix((mean.dat.p))
y.dat <- DGEList(counts = y.dat, group=c(1,2,3,4,5,6))
y.dat <- calcNormFactors(y.dat,method="TMM")
z.dat <- cpm(y.dat, normalized.lib.size=TRUE)

z_var <- apply(z.dat, 1, var)
z_mean <- apply(z.dat, 1, mean)

z.test_data <- z.dat[which(z_var > 50 & z_mean > 50), 1:6]


#first get the time point data together:
timepoint <- c(0,2,4,8,16,24)
# bind that to the dataframe
z.test_data <- rbind(timepoint, z.test_data)
row.names(z.test_data)[1]<-"time"

#save it to a temp file so ti doesnt clutter up my blog directory
tmp <- tempfile()
write.table(z.test_data,file=tmp, sep='\t', quote = F, col.names=NA)

#read it back in as an expression set
z.data <- table2eset(file=tmp)


# -------------------------------------------------------------------------


#standardise
data.z <- standardise(z.data)


#estimate the fuzzifier( a fuzzifier of 1 is essentially hard clustering)
m1 <- mestimate(data.z)
m1


#cluster number determination
Dmin(data.z, m=m1, crange=seq(2,22,1), repeats=3, visu=TRUE)


clust=10

c <- mfuzz(data.z,c=clust,m=m1)

mfuzz.plot(data.z,cl=c,mfrow=c(1,1),time.labels=c(0,2,4,8,16,24),new.window=FALSE)


# To extract list of genes
cluster_grp_genelist<- acore(data.z,cl=c,min.acore=0.7)
cluster_grp_genelist<-data.frame(cluster_grp_genelist)
cluster_grp_genelist<-as.data.frame(cluster_grp_genelist)


write.table(cluster_grp_genelist, file = "cluster_grp_genelist.txt", sep = "\t", col.names = TRUE, row.names = false)
class(cluster_grp_genelist)

r<-ClusterApply(cluster_grp_genelist, write, "test3.txt", append=T, ncolumns=1000 )
)
install.packages(clusterApply)


library(clusterApply)
library(BiocParallel)
library(Biobase)
library(BiocGenerics)
library(BiocManager)
library(parallel)
library(future.apply)


















