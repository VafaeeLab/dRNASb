#Code developed by Mojdeh Dinarvand

# Load library  -----------------------------------------------------------

library(limma)
library(DESeq2)
library(dplyr)
library(tidyr)
library(igraph)
library(data.table)
library(Mfuzz)
library(edgeR)
library(stringr)
library(stringi)
library(ClueR)
library(VennDiagram)
library(igraph)
library(corrplot)
library(psych)

# Load data and statistics analysis---------------------------------------------------------------

dat <- read.csv("Data/Alldat.csv",row.names = 1)
pheno <- read.csv("Data/pheno.csv", row.names = 1)
fun<-read.csv("Data/GO_Function.csv")
ppi<-read.csv("Data/PPI.csv")
ann<-read.csv("Data/Annotation.csv")

# Normalise ---------------------------------------------------------------

c <- pheno[colnames(dat), "groups"]
y <- DGEList(counts=dat, group=c, genes=rownames(dat))
y <- cpm(calcNormFactors(y, method="TMM"), log = TRUE)


# Filter ------------------------------------------------------------------

keep <- filterByExpr(y, group = c, min.count = log2(10))
y <- y[keep,]
normlise.count.dat<-data.frame(y)


# Differential gene expression analysis using "limma" ---------------------

group = as.factor(c)
design <- model.matrix(~ 0 + group)
colnames(design) <- levels(group)
rownames(design) <- colnames(y)
fit <- lmFit(y, design = design)
cont.matrix <- makeContrasts(WT.02_h - WT.00_h,
                             WT.04_h - WT.00_h,
                             WT.08_h - WT.00_h,
                             WT.16_h - WT.00_h,
                             WT.24_h - WT.00_h,levels=design)
fit.cont <- eBayes(contrasts.fit(fit, cont.matrix))
summa.fit <- decideTests(fit.cont)


de.ppi <- function(fit.cont, coef=1, lfc = 1, adjP =0.05){
  wtdt <- topTable(fit.cont, n = Inf, coef = coef)
  updt = with(wtdt, logFC > lfc & adj.P.Val < adjP)
  downdt = with(wtdt, logFC < lfc & adj.P.Val < adjP)
  wtdt$col="Not sig"
  wtdt$col[updt] = "Up"
  wtdt$col[downdt] = "Down"
  return(wtdt)}

DE <- list()
for (n in 1:5){
  DE[[n]] <- de.ppi(fit.cont, coef = n)
}

D<-DE[[1]] # 1=2h, 2=4h, 3=8h, 4=16h, 5=24h
Gene.name<-as.data.frame(row.names(D))    # add row to data frame in r
D<-cbind(Gene.name,D)
colnames(D)[3]<-"Gene.name"
write.csv(D,file = "Differential.gene.expression.for.2h.csv",row.names = FALSE)


# Average replecates across each time -------------------------------------

d <- cbind(rowMeans(dat[,c(1:3)], na.rm = T),
           rowMeans(dat[,c(4:6)], na.rm = T),
           rowMeans(dat[,c(7:9)], na.rm = T),
           rowMeans(dat[,c(10:12)], na.rm = T),
           rowMeans(dat[,c(13:15)], na.rm = T),
           rowMeans(dat[,c(16:18)], na.rm = T))
colnames(d) <- c("Mean.0h","Mean.20h","Mean.4h","Mean.8h","Mean.16h","Mean.24h")

Gene.name<-as.data.frame(row.names(d))    # add row to data frame in r
d <- cbind(Gene.name,d)
colnames(d)[1]<-"Gene.name"
write.csv(d,file = "pathogen.mean.csv",row.names = FALSE)

# Clustering using "Mfuzz" ------------------------------------------------

d <- read.csv("pathogen.mean.csv", row.names = 1)

y.dat<- as.matrix(d)
y.dat <- y.dat[which(apply(y.dat, 1, var)>2 & apply(y.dat,1,mean)>2), 1:6]
timepoint <- c(0,2,4,8,16,24)
y.dat <- rbind(timepoint, y.dat)
rownames(y.dat)[1]<- "time"
tmp<- tempfile()
write.table(y.dat,file=tmp, sep='\t',quote=FALSE, col.names=NA)
z.data <- table2eset(tmp)
data.z <-standardise(z.data)
class(data.z)
m1 <-mestimate(data.z)
Dmin(data.z, m=m1, crange=seq(2,22,1), repeats = 3, visu = TRUE)
clust=8
# set.seed(123456)
c<- mfuzz(data.z, c=clust, m=m1)
mfuzz.plot(data.z,cl=c,mfrow=c(4,4),min.mem=0.5,time.labels=c(0,2,4,8,16,24),new.window=FALSE)
membership<-c$membership
membership<-data.frame(membership)
fd<-data.frame(cor(t(c[[1]])))
acore<-acore(data.z,c,min.acore = 0.5)
acore_list<-do.call(rbind,lapply(seq_along(acore), function(i){data.frame(CLUSTER=i, acore[[i]])}))
colnames(acore_list)[2]<-"gene_name"
genelist<- acore(data.z,cl=c,min.acore=0.7)
temp <- do.call("rbind", lapply(genelist, FUN = function(x){
  return(paste0(as.character(x$NAME), collapse = ","))
}))
Cluster_list<-as.data.frame(temp)
colnames(Cluster_list) <-"gene_name"
Cluster_list<-str_split_fixed(Cluster_list$gene_name,",", n=Inf)
Cluster_list<-t(Cluster_list)
colnames(Cluster_list)<- c("Cluster1", "Cluster2","Cluster3","Cluster4","Cluster5","Cluster6","Cluster7","Cluster8")

write.csv(acore_list, file = "acore_list.CSV", quote = F, row.names = F)

# Make list ---------------------------------------------------------------

anno<-read.csv("Data/Annotation.csv")
GO<-unique(anno$Gene.Ontology.ID)
Uniprot.ID <- sapply(1:length(GO), function(i) paste(gsub("[[:space:]]", "", anno[which(anno$Gene.Ontology.ID==GO[i]),]$Uniprot.ID),collapse=" "))
Uniprot.ID<-as.data.frame(Uniprot.ID)
GO_Pro_ID<-data.frame(GO.ID=unique(anno$Gene.Ontology.ID),
                      Uniprot.ID=Uniprot.ID)


# Make list of list -------------------------------------------------------

Anno <- list()
groupSize <- 422
GO_IDs <- as.vector(GO_Pro_ID[,1])

for (i in GO_IDs) {
  myindex <- which(GO_Pro_ID == i)
  Anno[i] <- strsplit(as.character(GO_Pro_ID[myindex, 2]), " ")
}


# Enrichment using "ClueR" ------------------------------------------------

ce <- clustEnrichment(c, annotation=Anno, effectiveSize=c(2,100), pvalueCutoff=0.01)

out <- c()
i <- 1
for (clus in ce$enrich.list) {
  clus<- cbind(rep(paste0("Cluster_",i), nrow(clus)), clus)
  out <- rbind(out,clus)
  i = i+1
}

write.csv(out, file = "./Enrich.csv", quote = F, row.names = F)


# Venn diagram ------------------------------------------------------------------


D<-DE[[1]] # 1=2h, 2=4h, 3=8h, 4=16h, 5=24h
Gene.name<-as.data.frame(row.names(D))    # add row to data frame in r
D<-cbind(Gene.name,D)
colnames(D)[3]<-"Gene.name"
d2h<-D%>% filter(logFC<(-1)) # downregulation
u2h<-D%>% filter(logFC>(1)) # upregulation
Q2<-d2h


D<-DE[[2]] # 1=2h, 2=4h, 3=8h, 4=16h, 5=24h
Gene.name<-as.data.frame(row.names(D))    # add row to data frame in r
D<-cbind(Gene.name,D)
colnames(D)[3]<-"Gene.name"
d4h<-D%>% filter(logFC<(-1)) # downregulation
u4h<-D%>% filter(logFC>(1)) # upregulation
Q4<-d4h

D<-DE[[3]] # 1=2h, 2=4h, 3=8h, 4=16h, 5=24h
Gene.name<-as.data.frame(row.names(D))    # add row to data frame in r
D<-cbind(Gene.name,D)
colnames(D)[3]<-"Gene.name"
d8h<-D%>% filter(logFC<(-1)) # downregulation
u8h<-D%>% filter(logFC>(1)) # upregulation
Q8<-d8h

D<-DE[[4]] # 1=2h, 2=4h, 3=8h, 4=16h, 5=24h
Gene.name<-as.data.frame(row.names(D))    # add row to data frame in r
D<-cbind(Gene.name,D)
colnames(D)[3]<-"Gene.name"
d16h<-D%>% filter(logFC<(-1)) # downregulation
u16h<-D%>% filter(logFC>(1)) # upregulation
Q16<-d16h

D<-DE[[5]] # 1=2h, 2=4h, 3=8h, 4=16h, 5=24h
Gene.name<-as.data.frame(row.names(D))    # add row to data frame in r
D<-cbind(Gene.name,D)
colnames(D)[3]<-"Gene.name"
d24h<-D%>% filter(logFC<(-1)) # downregulation
u24h<-D%>% filter(logFC>(1)) # upregulation
Q24<-d24h



A=data.frame(intersect(d2h$Gene.name, d4h$Gene.name))
B=data.frame(intersect(d2h$Gene.name,d8h$Gene.name))
C=data.frame(intersect(d2h$Gene.name,d16h$Gene.name))
D=data.frame(intersect(d2h$Gene.name,d24h$Gene.name))
E=data.frame(intersect(d4h$Gene.name,d8h$Gene.name))
FF=data.frame(intersect(d4h$Gene.name,d16h$Gene.name))
K=data.frame(intersect(d4h$Gene.name,d24h$Gene.name))
G=data.frame(intersect(d8h$Gene.name,d16h$Gene.name))
M=data.frame(intersect(d8h$Gene.name,d24h$Gene.name))
H=data.frame(intersect(d16h$Gene.name,d24h$Gene.name))
colnames(A)<-"Attributes"
colnames(B)<-"Attributes"
colnames(C)<-"Attributes"
colnames(D)<-"Attributes"
colnames(E)<-"Attributes"
colnames(FF)<-"Attributes"
colnames(K)<-"Attributes"
colnames(G)<-"Attributes"
colnames(M)<-"Attributes"
colnames(H)<-"Attributes"

A1=data.frame(intersect(A$Attributes,E$Attributes))
B1=data.frame(intersect(A$Attributes,FF$Attributes))
C1=data.frame(intersect(A$Attributes,K$Attributes))
D1=data.frame(intersect(B$Attributes,G$Attributes))
E1=data.frame(intersect(B$Attributes,M$Attributes))
FF1=data.frame(intersect(C$Attributes,H$Attributes))
K1=data.frame(intersect(E$Attributes,G$Attributes))
G1=data.frame(intersect(E$Attributes,M$Attributes))
M1=data.frame(intersect(FF$Attributes,H$Attributes))
H1=data.frame(intersect(G$Attributes,H$Attributes))
colnames(A1)<-"Attributes"
colnames(B1)<-"Attributes"
colnames(C1)<-"Attributes"
colnames(D1)<-"Attributes"
colnames(E1)<-"Attributes"
colnames(FF1)<-"Attributes"
colnames(K1)<-"Attributes"
colnames(G1)<-"Attributes"
colnames(M1)<-"Attributes"
colnames(H1)<-"Attributes"

A2=data.frame(intersect(A1$Attributes,G$Attributes))
B2=data.frame(intersect(A1$Attributes,M$Attributes))
C2=data.frame(intersect(B1$Attributes,H$Attributes))
D2=data.frame(intersect(D1$Attributes,H$Attributes))
E2=data.frame(intersect(K1$Attributes,H$Attributes))
FF2=data.frame(intersect(A1$Attributes,H$Attributes))
colnames(A2)<-"Attributes"
colnames(B2)<-"Attributes"
colnames(C2)<-"Attributes"
colnames(D2)<-"Attributes"
colnames(E2)<-"Attributes"
colnames(FF2)<-"Attributes"

grid.newpage()


#below line not running
venn.plot <- draw.quintuple.venn(
  area1 = Q2,
  area2 =Q4,
  area3 = Q8,
  area4 = Q16,
  area5 = Q24,
  n12 = A,
  n13 = B,
  n14 = C,
  n15 = D,
  n23 = E,
  n24 = FF,
  n25 = K,
  n34 = G,
  n35 =M,
  n45 = H,
  n123 = A1,
  n124 = B1,
  n125 = C1,
  n134 = D1,
  n135 = E1,
  n145 = FF1,
  n234 = K1,
  n235 = G1,
  n245 = M1,
  n345 = H2,
  n1234 = A2,
  n1235 = B2,
  n1245 = C2,
  n1345 = D2,
  n2345 = E2,
  n12345 = FF2,
  category = c("2h", "4h", "8h", "16h", "24h"),
  fill = c("#e1bebe", "darkgoldenrod1", "#c70000", "#ff99cc", "#9ecbff"),
  cat.col = c("#e1bebe", "darkgoldenrod1", "#c70000", "#ff99cc", "#9ecbff"),
  cat.cex = 2,
  margin = 0.05,
  cex = c(1.5, 1.5, 1.5, 1.5, 1.5, 1, 0.8, 1, 0.8, 1, 0.8, 1, 0.8, 1, 0.8,
          1, 0.55, 1, 0.55, 1, 0.55, 1, 0.55, 1, 0.55, 1, 1, 1, 1, 1, 1.5),
  ind = TRUE
);
# Writing to file
tiff(filename = "Quintuple_Venn_diagram.tiff", compression = "lzw");
grid.draw(venn.plot);
dev.off()


# Network analysis using igraph --------------------------------

net<-graph.data.frame(unique(ppi[,c(2,3)]),directed = FALSE)
V(net)$label<-V(net)$name

# V(net)$igraph::degre<-degree(net)

deg<-igraph::degree(net,v=V(net), mode = c("total"),
                    loops = TRUE,normalized = FALSE)
# Hub gene significance ---------------------------------------------------

Hub<-as.data.frame(deg[which(deg>=1)])
gene.name<-as.data.frame(row.names(Hub))    # add row to data frame in r
Hub<-cbind(gene.name,Hub)
colnames(Hub)[1]<-"Gene.name"
colnames(Hub)[2]<-"edge.number"
remove(gene.name)
hist(Hub$edge.number)
write.csv(Hub,file = "Hub.correlated.90.posetive.csv",row.names = FALSE)

# Network betweenness --------------------------------------------------

b<-igraph::betweenness(net, v = V(net), directed = TRUE, weights = NULL,
                       nobigint = TRUE, normalized = FALSE)

bt<-as.data.frame(b)
gene.name<-as.data.frame(row.names(bt))    # add row to data frame in r
v<-cbind(gene.name,bt)
colnames(v)[1]<-"Gene.name"
colnames(v)[2]<-"Betweenness"
write.csv(v,file = "Betweenness.correlated.90.posetive.csv",row.names = FALSE)

# Shortest (directed or undirected) paths between vertices
distance_table(net, directed = TRUE)
mean_distance(net, directed = TRUE, unconnected = TRUE)

d<-distances(
  net,
  v = V(net),
  to = V(net),
  mode = c("all", "out", "in"),
  weights = NULL,
  algorithm = c("automatic", "unweighted", "dijkstra", "bellman-ford", "johnson")
)
write.csv(d,file = "Shortest.paths.correlated.90.posetive.csv",row.names = TRUE)

# Network Modules ---------------------------------------------------------

cl<-cluster_louvain(net, weights = NULL)
t<-as.data.frame(cl$membership)
t1<-as.data.frame(cl$names)
t2<-cbind(t1,t)
colnames(t2)[1]<-"Gene.name"
colnames(t2)[2]<-"membership"
write.csv(t2,file = "Modules.in.Negative.R.corr.all.GC026-26695.network.csv",row.names = FALSE)

# Plot Modules ------------------------------------------------------------

g_grouped = net

for(i in unique(V(net)$community)){
  groupV = which(V(net)$community == i)
  g_grouped = add_edges(g_grouped, combn(groupV, 2), attr=list(weight = 2))
}

l <- layout_nicely(g_grouped)

plot(cl,net, layout = layout_with_fr,
     vertex.size =10,
     edge.width = 1,
     vertex.label.dist=0.001,
     vertex.color ='gold',
     vertex.frame.color="#555555",
     edge.label=net$v,
     vertex.size=1,
     edge.color="gray",
     vertex.label.font=1,
     edge.label.font =1,
     edge.label.cex = 1,
     edge.arrow.size=0.2,
     edge.curved=0,
     vertex.label=V(net)$v,
     vertex.label.color="black",
     vertex.label.cex=0.5,
     vertex.label.cex = 0.5 )


# cor.analysis ------------------------------------------------------------
d <- read.csv("pathogen.mean.csv", row.names = 1)
dt <- t(d)

# t<-read.csv("dt", header = T, row.names = 1, check.names = F)
t <- dt


All.cor<-corr.test(t[1:54],t[1:54], method="pearson",adjust="holm", ci=FALSE) # kendall , # pearson, # spearman
r.df <- as.data.frame(All.cor$r)
R.corr <- r.df %>%
  mutate(gene1 = row.names(r.df)) %>%
  tidyr::pivot_longer(-gene1,
               names_to = "gene2", names_ptypes = list(gene2=character()),
               values_to = "corr") %>%
  mutate(gene1 = str_replace(gene1, "\\.\\.", " ("),
         gene1 = str_replace(gene1, "\\.$", ")"),
         gene2 = str_replace(gene2, "\\.\\.", " ("),
         gene2 = str_replace(gene2, "\\.$", ")")) %>%
  mutate(comb = paste(gene1, "-", gene2))


#below line not running
corrplot(as.matrix(r.df), is.corr = FALSE, method = "circle", order = "hclust", type = "upper")


write.csv(R.corr,file = "R.corr.all.intersect.genes.csv",row.names = FALSE)

posetive<-R.corr%>% filter(corr>0.7)
negative<-R.corr %>% filter(corr<(-0.7))











