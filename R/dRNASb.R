

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
library(gdata)
library(venn)
library(UpSetR)
library(plyr)
library(ggplot2)

# importFrom(magrittr,"%>%")

library(magrittr)

dir_path <- "/home/abhivij/UNSW/VafaeeLab/HostPathogenInteraction/dRNASb/R"
setwd(dir_path)
# Load data and statistics analysis---------------------------------------------------------------

pheno <- read.csv("Inputs/Pheno.csv", row.names = 1)

### Host
dat<- read.csv("Inputs/Host.data.csv",row.names = 1)
fun<-read.csv("Inputs/Host.annotation.function.csv")
ppi<-read.csv("Inputs/Host.ppi.csv")

### Pathogen
dat.p<- read.csv("Inputs/Pathogen.data.csv",row.names = 1)
fun.p<-read.csv("Inputs/Pathogen.annotation.function.csv")
ppi.p<-read.csv("Inputs/Pathogen.ppi.csv")


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
  return(wtdt)}

DE <- list()
for (n in 1:5){
  DE[[n]] <- de.ppi(fit.cont, coef = n)
}

################### 2h
D<-DE[[1]]
Gene.name<-as.data.frame(row.names(D))
colnames(Gene.name)<-"Gene.name"
D<-cbind(Gene.name,D)
colnames(D)[2]<-"logFC.2h"
Q2<-subset(D,D$logFC.2h<(-1)|D$logFC.2h>1)
Q2<-Q2[,c(1,2)]
dQ2<-D%>% filter(logFC.2h<(-1))
dQ2<-dQ2[,c(1,2)]
uQ2<-D%>% filter(logFC.2h>(1))
uQ2<-uQ2[,c(1,2)]
write.csv(D,file = paste0("./Results/Differential_gene_expression_analysis/","Host.DE-2h.csv"), row.names = FALSE)

################### 4h
D<-DE[[2]]
Gene.name<-as.data.frame(row.names(D))
colnames(Gene.name)<-"Gene.name"
D<-cbind(Gene.name,D)
colnames(D)[2]<-"logFC.4h"
Q4<-subset(D,D$logFC.4h<(-1)|D$logFC.4h>1)
Q4<-Q4[,c(1,2)]
dQ4<-D%>% filter(logFC.4h<(-1))
dQ4<-dQ4[,c(1,2)]
uQ4<-D%>% filter(logFC.4h>(1))
uQ4<-uQ4[,c(1,2)]
write.csv(D,file = paste0("./Results/","./Differential_gene_expression_analysis/","Host.DE-4h.csv"), row.names = FALSE)

################### 8h
D<-DE[[3]]
Gene.name<-as.data.frame(row.names(D))
colnames(Gene.name)<-"Gene.name"
D<-cbind(Gene.name,D)
colnames(D)[2]<-"logFC.8h"
Q8<-subset(D,D$logFC.8h<(-1)|D$logFC.8h>1)
Q8<-Q8[,c(1,2)]
dQ8<-D%>% filter(logFC.8h<(-1))
dQ8<-dQ8[,c(1,2)]
uQ8<-D%>% filter(logFC.8h>(1))
uQ8<-uQ8[,c(1,2)]
write.csv(D,file = paste0("./Results/","./Differential_gene_expression_analysis/","Host.DE-8h.csv"), row.names = FALSE)

################### 16h
D<-DE[[4]]
Gene.name<-as.data.frame(row.names(D))
colnames(Gene.name)<-"Gene.name"
D<-cbind(Gene.name,D)
colnames(D)[2]<-"logFC.16h"
Q16<-subset(D,D$logFC.16h<(-1)|D$logFC.16h>1)
Q16<-Q16[,c(1,2)]
dQ16<-D%>% filter(logFC.16h<(-1))
dQ16<-dQ16[,c(1,2)]
uQ16<-D%>% filter(logFC.16h>(1))
uQ16<-uQ16[,c(1,2)]
write.csv(D,file = paste0("./Results/","./Differential_gene_expression_analysis/","Host.DE-16h.csv"), row.names = FALSE)

################### 24h
D<-DE[[5]]
Gene.name<-as.data.frame(row.names(D))
colnames(Gene.name)<-"Gene.name"
D<-cbind(Gene.name,D)
colnames(D)[2]<-"logFC.24h"
Q24<-subset(D,D$logFC.24h<(-1)|D$logFC.24h>1)
Q24<-Q24[,c(1,2)]
dQ24<-D%>% filter(logFC.24h<(-1))
dQ24<-dQ24[,c(1,2)]
uQ24<-D%>% filter(logFC.24h>(1))
uQ24<-uQ24[,c(1,2)]
write.csv(D,file = paste0("./Results/","./Differential_gene_expression_analysis/","Host.DE-24h.csv"), row.names = FALSE)


# Average replecates across each time -------------------------------------
d <- cbind(rowMeans(dat[,c(1:3)], na.rm = T),
           rowMeans(dat[,c(4:6)], na.rm = T),
           rowMeans(dat[,c(7:9)], na.rm = T),
           rowMeans(dat[,c(10:12)], na.rm = T),
           rowMeans(dat[,c(13:15)], na.rm = T),
           rowMeans(dat[,c(16:18)], na.rm = T))
colnames(d) <- c("Mean.0h","Mean.2h","Mean.4h","Mean.8h","Mean.16h","Mean.24h")

Gene.name<-as.data.frame(row.names(d))
da <- cbind(Gene.name,d)
colnames(da)[1]<-"Gene.name"
Mean.host<-da
write.csv(Mean.host,file = paste0("./Results/","./Average_data/","All.host.mean.data.csv"), row.names = FALSE)


# DE Matrix  -----------------------------------------------------

### Select gene that shown DE at less in one time point
Q2$logFC.2h [Q2$logFC.2h<0]<-1
Q4$logFC.4h [Q4$logFC.4h<0]<-1
Q8$logFC.8h [Q8$logFC.8h<0]<-1
Q16$logFC.16h [Q16$logFC.16h<0]<-1
Q24$logFC.24h [Q24$logFC.24h<0]<-1
Q<-rbind.fill(Q2,Q4,Q8,Q16,Q24)
Q[is.na(Q[1:6])]<-0
Q <- Q %>% pivot_longer(logFC.2h:logFC.24h, names_to = "timepoint", values_to = "sign") %>%
  dplyr::filter(sign == 1) %>%
  pivot_wider(names_from = timepoint, values_from = sign, values_fill = 0)

### Get mean value for DE genes
m<-merge(Q,da, by="Gene.name")
m<-m[,-c(2:7)]
ED.host<-m
write.csv(m,file = paste0("./Results/","./Average_data/","All.host.ED.gene.csv"), row.names = FALSE)

### Get transpose matrix
mt<-data.frame(t(m))
g<-as.data.frame(row.names(mt))
mt <- cbind(g,mt)
colnames(mt)[1]<-"Gene.name"
write.table(mt,file = paste0("./Results/","./Matrix_Transpose/","Transpose.host.gene.csv"),  sep=",",col.names= FALSE,row.names = FALSE)

# Mfuzz Clustering   ------------------------------------------------
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
clust=10
c<- mfuzz(data.z, c=clust, m=m1)

tiff(filename ="./Results/Mfuzz_Clustering/Host.mfuzz.plot.tiff", compression = "lzw")
mfuzz.plot(data.z,cl=c,mfrow=c(4,4),min.mem=0.5,time.labels=c(0,2,4,8,16,24),new.window=FALSE)
dev.off()

pdf(file = "./Results/Mfuzz_Clustering/Host.mfuzz.plot.pdf",width = 10, height = 10)
mfuzz.plot(data.z,cl=c,mfrow=c(4,4),min.mem=0.5,time.labels=c(0,2,4,8,16,24),new.window=FALSE)
dev.off()

membership<-c$membership
membership<-data.frame(membership)
fd<-data.frame(cor(t(c[[1]])))
acore<-acore(data.z,c,min.acore = 0.5)
acore_list<-do.call(rbind,lapply(seq_along(acore), function(i){data.frame(CLUSTER=i, acore[[i]])}))
colnames(acore_list)[2]<-"Gene.name"
genelist<- acore(data.z,cl=c,min.acore=0.7)
temp <- do.call("rbind", lapply(genelist, FUN = function(x){
  return(paste0(as.character(x$NAME), collapse = ","))
}))
Cluster_list<-as.data.frame(temp)
colnames(Cluster_list) <-"Gene.name"
Cluster_list<-str_split_fixed(Cluster_list$Gene.name,",", n=Inf)
Cluster_list<-t(Cluster_list)
colnames(Cluster_list)<- c("Cluster1", "Cluster2","Cluster3","Cluster4","Cluster5","Cluster6","Cluster7","Cluster8","Cluster9","Cluster10")
write.csv(acore_list,file = paste0("./Results/","./Mfuzz_Clustering/","Host.cluster.acore_list.csv"), quote = F, row.names = FALSE)



# Make list ---------------------------------------------------------------
anno<-fun
GO<-unique(anno$Gene.ontology.ID)
Uniprot.ID <- sapply(1:length(GO), function(i) paste(gsub("[[:space:]]", "", anno[which(anno$Gene.ontology.ID==GO[i]),]$Gene.name),collapse=" "))
Uniprot.ID<-as.data.frame(Uniprot.ID)
GO_Pro_ID<-data.frame(GO.ID=unique(anno$Gene.ontology.ID),
                      Uniprot.ID=Uniprot.ID)


# Make list of list -------------------------------------------------------
Anno <- list()
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

colnames(out) [1] <-"Cluster.number"
colnames(out) [2] <-"Gene.ontology.ID"
colnames(out) [5] <-"Overlap.Gene.name"
write.csv(out, file= paste0 ("./Results/","./Mfuzz_Clustering/","Host.cluster.enrichment.csv"),quote = F, row.names = FALSE)


# Function using GO.term --------------------------------------------------
fu<-fun[,c(1:4)]
u<-merge(out, fu, by="Gene.ontology.ID")
write.csv(u, file= paste0 ("./Results/","./Mfuzz_Clustering/","Host.cluster.enrichment.function.csv"),quote = F, row.names = FALSE)

#  Variable and  frequency of Function
cl<-as.data.frame(out[,c(5)])
colnames(cl)<-"Gene.name"
gene.list<-cl %>% separate_rows(Gene.name, sep = "\\|") %>% group_by(Gene.name)
freq.list.host.gene<-cl %>% separate_rows(Gene.name, sep = "\\|") %>% group_by(Gene.name)%>% dplyr::summarize(n = n()) %>% arrange(desc(n))
write.csv(freq.list.host.gene, file= paste0 ("./Results/","./Mfuzz_Clustering/","Variable.enriched.genes.and.their.frequency.in.host.csv"),quote = F, row.names = FALSE)



# Venn diagram ------------------------------------------------------------------

################################ Host.downrulated
A=data.frame(intersect(dQ2$Gene.name,dQ4$Gene.name))
B=data.frame(intersect(dQ2$Gene.name,dQ8$Gene.name))
C=data.frame(intersect(dQ2$Gene.name,dQ16$Gene.name))
D=data.frame(intersect(dQ2$Gene.name,dQ24$Gene.name))
E=data.frame(intersect(dQ4$Gene.name,dQ8$Gene.name))
FF=data.frame(intersect(dQ4$Gene.name,dQ16$Gene.name))
K=data.frame(intersect(dQ4$Gene.name,dQ24$Gene.name))
G=data.frame(intersect(dQ8$Gene.name,dQ16$Gene.name))
M=data.frame(intersect(dQ8$Gene.name,dQ24$Gene.name))
H=data.frame(intersect(dQ16$Gene.name,dQ24$Gene.name))
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
venn.plot <- draw.quintuple.venn(
  area1=1331,  #dQ2
  area2=1082,  #dQ4
  area3=1661,  #dQ8
  area4=2774,  #dQ16
  area5=3315,  #dQ24
  n12=970,    #A
  n13=933,    #B
  n14=1097,   #C
  n15=1162,   #D
  n23=838,    #E
  n24=969,     #FF
  n25=991,     #K
  n34=1502,    #G
  n35=1543,     #M
  n45=2653,    #H
  n123=771,   #A1
  n124=883,    #B1
  n125=903,    #C1
  n134=892,    #D1
  n135=912,    #E1
  n145=1080,    #FF1
  n234=819,    #K1
  n235=822,    #G1
  n245=957,    #M1
  n345=1481,    #H1
  n1234=756,   #A2
  n1235=759,   #B2
  n1245=873,  #C2
  n1345=889,  #D2
  n2345=815,  #E2
  n12345=753, #FF2
  category = c("2h", "4h", "8h", "16h", "24h"),
  fill = c("#84b3e7", "#317456", "#abcdef", "#ff99cc", "#bd0000"),
  cat.col = c("#84b3e7", "#317456", "#abcdef", "#ff99cc", "#bd0000"),
  cat.cex = 2,
  margin = 0.05,
  cex = c(1.5, 1.5, 1.5, 1.5, 1.5, 1, 0.8, 1, 0.8, 1, 0.8, 1, 0.8, 1, 0.8,
          1, 0.55, 1, 0.55, 1, 0.55, 1, 0.55, 1, 0.55, 1, 1, 1, 1, 1, 1.5),
  ind = TRUE
)

# Writing to file
tiff(filename ="./Results/Venn_diagram/Host.downregulate.Venn.diagram.tiff", compression = "lzw")
grid.draw(venn.plot)
dev.off()

pdf(file ="./Results/Venn_diagram/Host.downregulate.Venn.diagram.plot.pdf",width = 5, height = 5)
grid.draw(venn.plot)
dev.off()


################################ Host.uprulated
A=data.frame(intersect(uQ2$Gene.name,uQ4$Gene.name))
B=data.frame(intersect(uQ2$Gene.name,uQ8$Gene.name))
C=data.frame(intersect(uQ2$Gene.name,uQ16$Gene.name))
D=data.frame(intersect(uQ2$Gene.name,uQ24$Gene.name))
E=data.frame(intersect(uQ4$Gene.name,uQ8$Gene.name))
FF=data.frame(intersect(uQ4$Gene.name,uQ16$Gene.name))
K=data.frame(intersect(uQ4$Gene.name,uQ24$Gene.name))
G=data.frame(intersect(uQ8$Gene.name,uQ16$Gene.name))
M=data.frame(intersect(uQ8$Gene.name,uQ24$Gene.name))
H=data.frame(intersect(uQ16$Gene.name,uQ24$Gene.name))
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
venn.plot <- draw.quintuple.venn(
  area1 = 648,
  area2 =633,
  area3 = 1091,
  area4 = 1577,
  area5 = 1949,
  n12 = 427,
  n13 = 247,
  n14 = 406,
  n15 = 447,
  n23 = 287,
  n24 = 417,
  n25 = 440,
  n34 = 650,
  n35 = 710,
  n45 = 1443,
  n123 = 195,
  n124 = 300,
  n125 = 317,
  n134 = 204,
  n135 = 207,
  n145 = 395,
  n234 = 242,
  n235 = 243,
  n245 = 398,
  n345 = 597,
  n1234 = 167,
  n1235 = 167,
  n1245 = 292,
  n1345 = 198,
  n2345 = 232,
  n12345 = 162,
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
tiff(filename ="./Results//Venn_diagram/Host.upregulate.Venn.diagram.tiff", compression = "lzw")
grid.draw(venn.plot)
dev.off()

pdf(file ="./Results//Venn_diagram/Host.upregulate.Venn.diagram.plot.pdf",width = 5, height = 5)
grid.draw(venn.plot)
dev.off()




# Upset Plot --------------------------------------------------------------
### Select downregulated gene
dQ2$logFC.2h [dQ2$logFC.2h<0]<-1
dQ4$logFC.4h [dQ4$logFC.4h<0]<-1
dQ8$logFC.8h [dQ8$logFC.8h<0]<-1
dQ16$logFC.16h [dQ16$logFC.16h<0]<-1
dQ24$logFC.24h [dQ24$logFC.24h<0]<-1
dQ<-rbind.fill(dQ2,dQ4,dQ8,dQ16,dQ24)
dQ[is.na(dQ[1:6])]<-0
dQ <- dQ %>% pivot_longer(logFC.2h:logFC.24h, names_to = "timepoint", values_to = "sign") %>%
  dplyr::filter(sign == 1) %>%
  pivot_wider(names_from = timepoint, values_from = sign, values_fill = 0)%>% t()

write.table(dQ,file = paste0("./Inputs/","Host.all.downregulated.ED.csv"),  sep=",",col.names= FALSE,row.names = TRUE)

### Select upregulated gene
uQ2$logFC.2h [uQ2$logFC.2h>0]<-1
uQ4$logFC.4h [uQ4$logFC.4h>0]<-1
uQ8$logFC.8h [uQ8$logFC.8h>0]<-1
uQ16$logFC.16h [uQ16$logFC.16h>0]<-1
uQ24$logFC.24h [uQ24$logFC.24h>0]<-1
uQ<-rbind.fill(uQ2,uQ4,uQ8,uQ16,uQ24)
uQ[is.na(uQ[1:6])]<-0
uQ <- uQ %>% pivot_longer(logFC.2h:logFC.24h, names_to = "timepoint", values_to = "sign") %>%
  dplyr::filter(sign == 1) %>%
  pivot_wider(names_from = timepoint, values_from = sign, values_fill = 0)%>% t()

write.table(uQ,file = paste0("./Inputs/","Host.all.upregulated.ED.csv"),  sep=",",col.names= FALSE,row.names = TRUE)



### Downregulated plot
y <-t(read.csv("Inputs/Host.all.downregulated.ED.csv")) %>% as.data.frame()
colnames(y)<-y[1,]
y <- y[-1,] %>% mutate_if(is.character, as.numeric)
upset(y)

# Setting colors
main_bar_col <- c("blue4")
sets_bar_col <- c("coral1")
matrix_col <- c("forestgreen")
shade_col <- c("wheat4")


# Setting Set Variables
mb_ratio1 <- c(0.55,0.45)

tiff(filename ="./Results//Upset_plot/Host.downregulate.upset.plot.tiff", compression = "lzw")
upset(y,
      mb.ratio = mb_ratio1,
      mainbar.y.label = "Interaction of downregulated genes",
      sets.x.label = "Number of Genes",
      order.by = "freq",
      # show.numbers = TRUE,
      point.size = 2,
      line.size = 1,
      main.bar.color = main_bar_col,
      sets.bar.color = sets_bar_col,
      matrix.color = matrix_col,
      shade.color = shade_col )


dev.off()

pdf(file ="./Results//Upset_plot/Host.downregulate.upset.plot.pdf",width = 5, height = 5)
upset(y,
      mb.ratio = mb_ratio1,
      mainbar.y.label = "Interaction of downregulated genes",
      sets.x.label = "Number of Genes",
      order.by = "freq",
      # show.numbers = TRUE,
      point.size = 2,
      line.size = 1,
      main.bar.color = main_bar_col,
      sets.bar.color = sets_bar_col,
      matrix.color = matrix_col,
      shade.color = shade_col )
dev.off()



### Upregulate plot

y <-t(read.csv("Inputs/Host.all.upregulated.ED.csv")) %>% as.data.frame()
colnames(y)<-y[1,]
y <- y[-1,] %>% mutate_if(is.character, as.numeric)
upset(y)

# Setting colors
main_bar_col <- c("violetred4")
sets_bar_col <- c("turquoise4")
matrix_col <- c("slateblue4")
shade_col <- c("wheat4")


# Setting Set Variables
mb_ratio1 <- c(0.55,0.45)

tiff(filename ="./Results//Upset_plot/Host.upregulate.upset.plot.tiff", compression = "lzw")
upset(y,
      mb.ratio = mb_ratio1,
      mainbar.y.label = "Interaction of upregulated genes",
      sets.x.label = "Number of Genes",
      order.by = "freq",
      # show.numbers = TRUE,
      point.size = 2,
      line.size = 1,
      main.bar.color = main_bar_col,
      sets.bar.color = sets_bar_col,
      matrix.color = matrix_col,
      shade.color = shade_col )


dev.off()

pdf(file ="./Results//Upset_plot/Host.upregulate.upset.plot.pdf",width = 5, height = 5)
upset(y,
      mb.ratio = mb_ratio1,
      mainbar.y.label = "Interaction of upregulated genes",
      sets.x.label = "Number of Genes",
      order.by = "freq",
      # show.numbers = TRUE,
      point.size = 2,
      line.size = 1,
      main.bar.color = main_bar_col,
      sets.bar.color = sets_bar_col,
      matrix.color = matrix_col,
      shade.color = shade_col )
dev.off()


# Network analysis using igraph --------------------------------
net<-graph.data.frame(unique(ppi[,c(2,4)]),directed = FALSE)
V(net)$label<-V(net)$name
deg<-igraph::degree(net,v=V(net), mode = c("total"),
                    loops = TRUE,normalized = FALSE)

# Hub gene  ---------------------------------------------------
Hub<-as.data.frame(deg[which(deg>=100)])
gene.name<-as.data.frame(row.names(Hub))
Hub<-cbind(gene.name,Hub)
colnames(Hub)[1]<-"Gene.name"
colnames(Hub)[2]<-"edge.number"
Hub.host<-Hub
write.csv(Hub.host, file= paste0 ("./Results/","./Network_analysis/","Host.hub.more.100.csv"),quote = F, row.names = FALSE)

# Network betweennes --------------------------------------------------
b<-igraph::betweenness(net, v = V(net), directed = TRUE, weights = NULL,
                       nobigint = TRUE, normalized = FALSE)

b<-as.data.frame(b)
gene.name<-as.data.frame(row.names(b))
b<-cbind(gene.name,b)
colnames(b)[1]<-"Gene.name"
colnames(b)[2]<-"Betweennes"
Betweennes.host<-b%>% filter(Betweennes>1000)
write.csv(Betweennes.host,file = paste0 ("./Results/","./Network_analysis/","Host.betweennes.csv"),quote = F, row.names = FALSE)

# Network short distance --------------------------------------------------
d<-distances(
  net,
  v = V(net),
  to = V(net),
  mode = c("all", "out", "in"),
  weights = NULL,
  algorithm = c("automatic", "unweighted", "dijkstra", "bellman-ford", "johnson")
)
g<-as.data.frame(row.names(d))
colnames(g)[1]<-"Gene.name"
d<-cbind(g,d)
write.csv(d,file = paste0 ("./Results/","./Network_analysis/","Host.shortest.path.csv"),quote = F, row.names = FALSE)


# Network Modules ---------------------------------------------------------
cl<-cluster_louvain(net, weights = NULL)
t<-as.data.frame(cl$membership)
t1<-as.data.frame(cl$names)
t2<-cbind(t1,t)
colnames(t2)[1]<-"Gene.name"
colnames(t2)[2]<-"membership"
write.csv(t2,file = paste0 ("./Results/","./Network_analysis/","Host.modules.in.ppi.network.csv"),quote = F, row.names = FALSE)


# Plot Modules ------------------------------------------------------------
g_grouped = net

for(i in unique(V(net)$community)){
  groupV = which(V(net)$community == i)
  g_grouped = add_edges(g_grouped, combn(groupV, 2), attr=list(weight = 2))
}

l <- layout_nicely(g_grouped)

# Writing to file
pdf(file ="./Results/Network_analysis/Host.modules.in.ppi.network.plot.pdf",width = 5, height = 5)
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
dev.off()

tiff(filename ="./Results/Network_analysis/Host.modules.in.ppi.network.plot.tiff", compression = "lzw")
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

dev.off()


################################################### Pathogen #######################################################################



# Normalise ---------------------------------------------------------------
c.p <- pheno[colnames(dat.p), "groups"]
y.p <- DGEList(counts=dat.p, group=c.p, genes=rownames(dat.p))
y.p <- cpm(calcNormFactors(y.p, method="TMM"), log = TRUE)


# Filter ------------------------------------------------------------------
keep.p <- filterByExpr(y.p, group = c.p, min.count = log2(10))
y.p <- y.p[keep.p,]
normlise.count.dat.p<-data.frame(y.p)

# Differential gene expression analysis using "limma" ---------------------

group.p = as.factor(c.p)
design.p <- model.matrix(~ 0 + group.p)
colnames(design.p) <- levels(group.p)
rownames(design.p) <- colnames(y.p)
fit.p <- lmFit(y.p, design = design.p)
cont.matrix.p <- makeContrasts(WT.02_h - WT.00_h,
                               WT.04_h - WT.00_h,
                               WT.08_h - WT.00_h,
                               WT.16_h - WT.00_h,
                               WT.24_h - WT.00_h,levels=design.p)
fit.cont.p <- eBayes(contrasts.fit(fit.p, cont.matrix.p))
summa.fit.p <- decideTests(fit.cont.p)


de.ppi.p <- function(fit.cont, coef=1, lfc = 1, adjP =0.05){
  wtdt <- topTable(fit.cont, n = Inf, coef = coef)
  updt = with(wtdt, logFC > lfc & adj.P.Val < adjP)
  downdt = with(wtdt, logFC < lfc & adj.P.Val < adjP)
  return(wtdt)}

DE.p <- list()
for (n in 1:5){
  DE.p[[n]] <- de.ppi.p(fit.cont.p, coef = n)
}


################### 2h
D<-DE.p[[1]]
Gene.name<-as.data.frame(row.names(D))
colnames(Gene.name)<-"Gene.name"
D<-cbind(Gene.name,D)
colnames(D)[2]<-"logFC.2h"
Qp2<-subset(D,D$logFC.2h<(-1)|D$logFC.2h>1)
Qp2<-Qp2[,c(1,2)]
dpQ2<-D%>% filter(logFC.2h<(-1))
dpQ2<-dpQ2[,c(1,2)]
upQ2<-D%>% filter(logFC.2h>(1))
upQ2<-upQ2[,c(1,2)]
write.csv(D,file = paste0("./Results/","./Differential_gene_expression_analysis/","Pathogen.DE-2h.csv"), row.names = FALSE)

################### 4h
D<-DE.p[[2]]
Gene.name<-as.data.frame(row.names(D))
colnames(Gene.name)<-"Gene.name"
D<-cbind(Gene.name,D)
colnames(D)[2]<-"logFC.4h"
Qp4<-subset(D,D$logFC.4h<(-1)|D$logFC.4h>1)
Qp4<-Qp4[,c(1,2)]
dpQ4<-D%>% filter(logFC.4h<(-1))
dpQ4<-dpQ4[,c(1,2)]
upQ4<-D%>% filter(logFC.4h>(1))
upQ4<-upQ4[,c(1,2)]
write.csv(D,file = paste0("./Results/","./Differential_gene_expression_analysis/","Pathogen.DE-4h.csv"), row.names = FALSE)

################### 8h
D<-DE.p[[3]]
Gene.name<-as.data.frame(row.names(D))
colnames(Gene.name)<-"Gene.name"
D<-cbind(Gene.name,D)
colnames(D)[2]<-"logFC.8h"
Qp8<-subset(D,D$logFC.8h<(-1)|D$logFC.8h>1)
Qp8<-Qp8[,c(1,2)]
dpQ8<-D%>% filter(logFC.8h<(-1))
dpQ8<-dpQ8[,c(1,2)]
upQ8<-D%>% filter(logFC.8h>(1))
upQ8<-upQ8[,c(1,2)]
write.csv(D,file = paste0("./Results/","./Differential_gene_expression_analysis/","Pathogen.DE-8h.csv"), row.names = FALSE)

################### 16h
D<-DE.p[[4]]
Gene.name<-as.data.frame(row.names(D))
colnames(Gene.name)<-"Gene.name"
D<-cbind(Gene.name,D)
colnames(D)[2]<-"logFC.16h"
Qp16<-subset(D,D$logFC.16h<(-1)|D$logFC.16h>1)
Qp16<-Qp16[,c(1,2)]
dpQ16<-D%>% filter(logFC.16h<(-1))
dpQ16<-dpQ16[,c(1,2)]
upQ16<-D%>% filter(logFC.16h>(1))
upQ16<-upQ16[,c(1,2)]
write.csv(D,file = paste0("./Results/","./Differential_gene_expression_analysis/","Pathogen.DE-16h.csv"), row.names = FALSE)

################### 24h
D<-DE.p[[5]]
Gene.name<-as.data.frame(row.names(D))
colnames(Gene.name)<-"Gene.name"
D<-cbind(Gene.name,D)
colnames(D)[2]<-"logFC.24h"
Qp24<-subset(D,D$logFC.24h<(-1)|D$logFC.24h>1)
Qp24<-Qp24[,c(1,2)]
dpQ24<-D%>% filter(logFC.24h<(-1))
dpQ24<-dpQ24[,c(1,2)]
upQ24<-D%>% filter(logFC.24h>(1))
upQ24<-upQ24[,c(1,2)]
write.csv(D,file = paste0("./Results/","./Differential_gene_expression_analysis/","Pathogen.DE-24h.csv"), row.names = FALSE)


# Average replecates across each time -------------------------------------
d.p <- cbind(rowMeans(dat.p[,c(1:3)], na.rm = T),
             rowMeans(dat.p[,c(4:6)], na.rm = T),
             rowMeans(dat.p[,c(7:9)], na.rm = T),
             rowMeans(dat.p[,c(10:12)], na.rm = T),
             rowMeans(dat.p[,c(13:15)], na.rm = T),
             rowMeans(dat.p[,c(16:18)], na.rm = T))
colnames(d.p) <- c("Mean.0h","Mean.2h","Mean.4h","Mean.8h","Mean.16h","Mean.24h")

Gene.name<-as.data.frame(row.names(d.p))
da.p <- cbind(Gene.name,d.p)
colnames(da.p)[1]<-"Gene.name"
Mean.pathogen<-da.p
write.csv(Mean.pathogen,file = paste0("./Results/","./Average_data/","All.pathogen.mean.data.csv"), row.names = FALSE)



# ED Matrix -----------------------------------------------------

### Select gene that shown DE at less in one time point
Qp2$logFC.2h[Qp2$logFC.2h<0]<-1
Qp4$logFC.4h[Qp4$logFC.4h<0]<-1
Qp8$logFC.8h[Qp8$logFC.8h<0]<-1
Qp16$logFC.16h[Qp16$logFC.16h<0]<-1
Qp24$logFC.24h[Qp24$logFC.24h<0]<-1
Qp<-rbind.fill(Qp2,Qp4,Qp8,Qp16,Qp24)
Qp[is.na(Qp[1:6])]<-0
Qp <- Qp %>% pivot_longer(logFC.2h:logFC.24h, names_to = "timepoint", values_to = "sign") %>%
  dplyr::filter(sign == 1) %>%
  pivot_wider(names_from = timepoint, values_from = sign, values_fill = 0)

### Get mean value for ED genes
m.p<-merge(Qp,da.p, by="Gene.name")
m.p<-m.p[,-c(2:7)]
ED.pathogen<-m.p
write.csv(m.p,file = paste0("./Results/","./Average_data/","All.pathogen.ED.gene.csv"), row.names = FALSE)


### Get transpose matrix
mt<-data.frame(t(m.p))
g<-as.data.frame(row.names(mt))
mt <- cbind(g,mt)
colnames(mt)[1]<-"Gene.name"
write.table(mt,file = paste0("./Results/","./Matrix_Transpose/","Transpose.pathogen.gene.csv"),  sep=",",col.names= FALSE,row.names = FALSE)



# Mfuzz Clustering   ------------------------------------------------
y.dat<- as.matrix(d.p)
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
clust=10
c<- mfuzz(data.z, c=clust, m=m1)

tiff(filename ="./Results/Mfuzz_Clustering/Pathogen.mfuzz.plot.tiff", compression = "lzw")
mfuzz.plot(data.z,cl=c,mfrow=c(4,4),min.mem=0.5,time.labels=c(0,2,4,8,16,24),new.window=FALSE)
dev.off()

pdf(file = "./Results/Mfuzz_Clustering/Pathogen.mfuzz.plot.pdf",width = 10, height = 10)
mfuzz.plot(data.z,cl=c,mfrow=c(4,4),min.mem=0.5,time.labels=c(0,2,4,8,16,24),new.window=FALSE)
dev.off()

membership<-c$membership
membership<-data.frame(membership)
fd<-data.frame(cor(t(c[[1]])))
acore<-acore(data.z,c,min.acore = 0.5)
acore_list<-do.call(rbind,lapply(seq_along(acore), function(i){data.frame(CLUSTER=i, acore[[i]])}))
colnames(acore_list)[2]<-"Gene.name"
genelist<- acore(data.z,cl=c,min.acore=0.7)
temp <- do.call("rbind", lapply(genelist, FUN = function(x){
  return(paste0(as.character(x$NAME), collapse = ","))
}))
Cluster_list<-as.data.frame(temp)
colnames(Cluster_list) <-"Gene.name"
Cluster_list<-str_split_fixed(Cluster_list$Gene.name,",", n=Inf)
Cluster_list<-t(Cluster_list)
colnames(Cluster_list)<- c("Cluster1", "Cluster2","Cluster3","Cluster4","Cluster5","Cluster6","Cluster7","Cluster8","Cluster9","Cluster10")  ### ? how we can make it depanded to cluster number???

write.csv(acore_list,file = paste0("./Results/","./Mfuzz_Clustering/","Pathogen.cluster.acore_list.csv"), quote = F, row.names = FALSE)



# Make list ---------------------------------------------------------------
anno<-fun.p
GO<-unique(anno$Gene.ontology.ID)
Uniprot.ID <- sapply(1:length(GO), function(i) paste(gsub("[[:space:]]", "", anno[which(anno$Gene.ontology.ID==GO[i]),]$Uniprot.ID),collapse=" "))
Uniprot.ID<-as.data.frame(Uniprot.ID)
GO_Pro_ID<-data.frame(GO.ID=unique(anno$Gene.ontology.ID),
                      Uniprot.ID=Uniprot.ID)


# Make list of list -------------------------------------------------------
Anno <- list()
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

colnames(out) [1] <-"Cluster.number"
colnames(out) [2] <-"Gene.ontology.ID"
colnames(out) [5] <-"Overlap.Gene.name"
write.csv(out, file= paste0 ("./Results/","./Mfuzz_Clustering/","Pathogen.cluster.enrichment.csv"),quote = F, row.names = FALSE)


# Function using GO.term --------------------------------------------------
fu<-fun.p[,c(1,2)]
u<-merge(out, fu, by="Gene.ontology.ID")
write.csv(u, file= paste0 ("./Results/","./Mfuzz_Clustering/","Pathogen.cluster.enrichment.function.csv"),quote = F, row.names = FALSE)



#  Variable and  frequency of Function
cl<-as.data.frame(out[,c(5)])
colnames(cl)<-"Gene.name"
gene.list<-cl %>% separate_rows(Gene.name, sep = "\\|") %>% group_by(Gene.name)
freq.list.pathogen.gene<-cl %>% separate_rows(Gene.name, sep = "\\|") %>% group_by(Gene.name)%>% dplyr::summarize(n = n()) %>% arrange(desc(n))
write.csv(freq.list.pathogen.gene, file= paste0 ("./Results/","./Mfuzz_Clustering/","Variable.enriched.genes.and.their.frequency.in.pathogen.csv"),quote = F, row.names = FALSE)


# Venn diagram ------------------------------------------------------------------

################################ Pathogen.downrulated
A=data.frame(intersect(dpQ2$Gene.name,dpQ4$Gene.name))
B=data.frame(intersect(dpQ2$Gene.name,dpQ8$Gene.name))
C=data.frame(intersect(dpQ2$Gene.name,dpQ16$Gene.name))
D=data.frame(intersect(dpQ2$Gene.name,dpQ24$Gene.name))
E=data.frame(intersect(dpQ4$Gene.name,dpQ8$Gene.name))
FF=data.frame(intersect(dpQ4$Gene.name,dpQ16$Gene.name))
K=data.frame(intersect(dpQ4$Gene.name,dpQ24$Gene.name))
G=data.frame(intersect(dpQ8$Gene.name,dpQ16$Gene.name))
M=data.frame(intersect(dpQ8$Gene.name,dpQ24$Gene.name))
H=data.frame(intersect(dpQ16$Gene.name,dpQ24$Gene.name))
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
venn.plot <- draw.quintuple.venn(
  area1=303,  #dpQ2
  area2=208,  #dpQ4
  area3=206,  #dpQ8
  area4=134,  #dpQ16
  area5=152,  #dpQ24
  n12=141,    #A
  n13=136,    #B
  n14=80,     #C
  n15=68,     #D
  n23=137,    #E
  n24=95,     #FF
  n25=79,     #K
  n34=105,    #G
  n35=85,     #M
  n45=102,    #H
  n123=106,   #A1
  n124=71,    #B1
  n125=58,    #C1
  n134=69,    #D1
  n135=53,    #E1
  n145=57,    #FF1
  n234=86,    #K1
  n235=68,    #G1
  n245=72,    #M1
  n345=80,    #H1
  n1234=66,   #A2
  n1235=52,   #B2
  n1245=52,  #C2
  n1345=51,  #D2
  n2345=66,  #E2
  n12345=50, #FF2
  category = c("2h", "4h", "8h", "16h", "24h"),
  fill = c("#84b3e7", "#317456", "#abcdef", "#ff99cc", "#bd0000"),
  cat.col = c("#84b3e7", "#317456", "#abcdef", "#ff99cc", "#bd0000"),
  cat.cex = 2,
  margin = 0.05,
  cex = c(1.5, 1.5, 1.5, 1.5, 1.5, 1, 0.8, 1, 0.8, 1, 0.8, 1, 0.8, 1, 0.8,
          1, 0.55, 1, 0.55, 1, 0.55, 1, 0.55, 1, 0.55, 1, 1, 1, 1, 1, 1.5),
  ind = TRUE
)

# Writing to file
tiff(filename ="./Results/Venn_diagram/Pathogn.downregulate.Venn.diagram.tiff", compression = "lzw")
grid.draw(venn.plot)
dev.off()

pdf(file ="./Results/Venn_diagram/Pathogen.downregulate.Venn.diagram.plot.pdf",width = 5, height = 5)
grid.draw(venn.plot)
dev.off()


################################ Pathogen.uprulated
A=data.frame(intersect(upQ2$Gene.name,upQ4$Gene.name))
B=data.frame(intersect(upQ2$Gene.name,upQ8$Gene.name))
C=data.frame(intersect(upQ2$Gene.name,upQ16$Gene.name))
D=data.frame(intersect(upQ2$Gene.name,upQ24$Gene.name))
E=data.frame(intersect(upQ4$Gene.name,upQ8$Gene.name))
FF=data.frame(intersect(upQ4$Gene.name,upQ16$Gene.name))
K=data.frame(intersect(upQ4$Gene.name,upQ24$Gene.name))
G=data.frame(intersect(upQ8$Gene.name,upQ16$Gene.name))
M=data.frame(intersect(upQ8$Gene.name,upQ24$Gene.name))
H=data.frame(intersect(upQ16$Gene.name,upQ24$Gene.name))
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
venn.plot <- draw.quintuple.venn(
  area1=157,  #upQ2
  area2=180,  #upQ4
  area3=202,  #upQ8
  area4=225,  #upQ16
  area5=310,  #upQ24
  n12=122,    #A
  n13=119,    #B
  n14=124,    #C
  n15=96,     #D
  n23=145,    #E
  n24=151,    #FF
  n25=120,    #K
  n34=161,    #G
  n35=121,    #M
  n45=169,    #H
  n123=108,   #A1
  n124=110,   #B1
  n125=85,    #C1
  n134=111,   #D1
  n135=81,    #E1
  n145=91,    #FF1
  n234=134,   #K1
  n235=99,    #G1
  n245=110,   #M1
  n345=114,   #H1
  n1234=103,  #A2
  n1235=77,   #B2
  n1245=81,   #C2
  n1345=80,   #D2
  n2345=98,   #E2
  n12345=76,  #FF2
  category = c("2h", "4h", "8h", "16h", "24h"),
  fill = c("#e1bebe", "darkgoldenrod1", "#c70000", "#ff99cc", "#9ecbff"),
  cat.col = c("#e1bebe", "darkgoldenrod1", "#c70000", "#ff99cc", "#9ecbff"),
  cat.cex = 2,
  margin = 0.05,
  cex = c(1.5, 1.5, 1.5, 1.5, 1.5, 1, 0.8, 1, 0.8, 1, 0.8, 1, 0.8, 1, 0.8,
          1, 0.55, 1, 0.55, 1, 0.55, 1, 0.55, 1, 0.55, 1, 1, 1, 1, 1, 1.5),
  ind = TRUE
)

# Writing to file
tiff(filename ="./Results//Venn_diagram/Pathogen.upregulate.Venn.diagram.tiff", compression = "lzw")
grid.draw(venn.plot)
dev.off()

pdf(file ="./Results//Venn_diagram/Pathogen.upregulate.Venn.diagram.plot.pdf",width = 5, height = 5)
grid.draw(venn.plot)
dev.off()


# Upset Plot --------------------------------------------------------------
### Select downregulated gene
dpQ2$logFC.2h [dpQ2$logFC.2h<0]<-1
dpQ4$logFC.4h [dpQ4$logFC.4h<0]<-1
dpQ8$logFC.8h [dpQ8$logFC.8h<0]<-1
dpQ16$logFC.16h [dpQ16$logFC.16h<0]<-1
dpQ24$logFC.24h [dpQ24$logFC.24h<0]<-1
dpQ<-rbind.fill(dpQ2,dpQ4,dpQ8,dpQ16,dpQ24)
dpQ[is.na(dpQ[1:6])]<-0
dpQ <- dpQ %>% pivot_longer(logFC.2h:logFC.24h, names_to = "timepoint", values_to = "sign") %>%
  dplyr::filter(sign == 1) %>%
  pivot_wider(names_from = timepoint, values_from = sign, values_fill = 0)%>% t()

write.table(dpQ,file = paste0("./Inputs/","Pathogen.all.downregulated.ED.csv"),  sep=",",col.names= FALSE,row.names = TRUE)


### Select upregulated gene
upQ2$logFC.2h [upQ2$logFC.2h>0]<-1
upQ4$logFC.4h [upQ4$logFC.4h>0]<-1
upQ8$logFC.8h [upQ8$logFC.8h>0]<-1
upQ16$logFC.16h [upQ16$logFC.16h>0]<-1
upQ24$logFC.24h [upQ24$logFC.24h>0]<-1
upQ<-rbind.fill(upQ2,upQ4,upQ8,upQ16,upQ24)
upQ[is.na(upQ[1:6])]<-0
upQ <- upQ %>% pivot_longer(logFC.2h:logFC.24h, names_to = "timepoint", values_to = "sign") %>%
  dplyr::filter(sign == 1) %>%
  pivot_wider(names_from = timepoint, values_from = sign, values_fill = 0)%>% t()

write.table(upQ,file = paste0("./Inputs/","Pathogen.all.upregulated.ED.csv"),  sep=",",col.names= FALSE,row.names = TRUE)



### Downregulated
y <-t(read.csv("Inputs/Pathogen.all.downregulated.ED.csv")) %>% as.data.frame()
colnames(y)<-y[1,]
y <- y[-1,] %>% mutate_if(is.character, as.numeric)
upset(y)

# Setting colors
main_bar_col <- c("blue4")
sets_bar_col <- c("coral1")
matrix_col <- c("forestgreen")
shade_col <- c("wheat4")


# Setting Set Variables
mb_ratio1 <- c(0.55,0.45)

tiff(filename ="./Results//Upset_plot/Pathogen.downregulate.upset.plot.tiff", compression = "lzw")
upset(y,
      mb.ratio = mb_ratio1,
      mainbar.y.label = "Interaction of downregulated genes",
      sets.x.label = "Number of Genes",
      order.by = "freq",
      # show.numbers = TRUE,
      point.size = 2,
      line.size = 1,
      main.bar.color = main_bar_col,
      sets.bar.color = sets_bar_col,
      matrix.color = matrix_col,
      shade.color = shade_col )


dev.off()

pdf(file ="./Results//Upset_plot/Pathogen.downregulate.upset.plot.pdf",width = 5, height = 5)
upset(y,
      mb.ratio = mb_ratio1,
      mainbar.y.label = "Interaction of downregulated genes",
      sets.x.label = "Number of Genes",
      order.by = "freq",
      # show.numbers = TRUE,
      point.size = 2,
      line.size = 1,
      main.bar.color = main_bar_col,
      sets.bar.color = sets_bar_col,
      matrix.color = matrix_col,
      shade.color = shade_col )
dev.off()



### Upregulate

y <-t(read.csv("Inputs/Pathogen.all.upregulated.ED.csv")) %>% as.data.frame()
colnames(y)<-y[1,]
y <- y[-1,] %>% mutate_if(is.character, as.numeric)
upset(y)

# Setting colors
main_bar_col <- c("violetred4")
sets_bar_col <- c("turquoise4")
matrix_col <- c("slateblue4")
shade_col <- c("wheat4")


# Setting Set Variables
mb_ratio1 <- c(0.55,0.45)

tiff(filename ="./Results//Upset_plot/Pathogen.upregulate.upset.plot.tiff", compression = "lzw")
upset(y,
      mb.ratio = mb_ratio1,
      mainbar.y.label = "Interaction of upregulated genes",
      sets.x.label = "Number of Genes",
      order.by = "freq",
      # show.numbers = TRUE,
      point.size = 2,
      line.size = 1,
      main.bar.color = main_bar_col,
      sets.bar.color = sets_bar_col,
      matrix.color = matrix_col,
      shade.color = shade_col )


dev.off()

pdf(file ="./Results//Upset_plot/Pathogen.upregulate.upset.plot.pdf",width = 5, height = 5)
upset(y,
      mb.ratio = mb_ratio1,
      mainbar.y.label = "Interaction of upregulated genes",
      sets.x.label = "Number of Genes",
      order.by = "freq",
      # show.numbers = TRUE,
      point.size = 2,
      line.size = 1,
      main.bar.color = main_bar_col,
      sets.bar.color = sets_bar_col,
      matrix.color = matrix_col,
      shade.color = shade_col )
dev.off()


# Network analysis using igraph --------------------------------
net<-graph.data.frame(unique(ppi.p[,c(1,4)]),directed = FALSE)
V(net)$label<-V(net)$name
deg<-igraph::degree(net,v=V(net), mode = c("total"),
                    loops = TRUE,normalized = FALSE)

# Hub gene  ---------------------------------------------------
Hub.p<-as.data.frame(deg[which(deg>=10)])
gene.name<-as.data.frame(row.names(Hub.p))
Hub.p<-cbind(gene.name,Hub.p)
colnames(Hub.p)[1]<-"Gene.name"
colnames(Hub.p)[2]<-"edge.number"
Hub.pathogen<-Hub.p
write.csv(Hub.pathogen, file= paste0 ("./Results/","./Network_analysis/","Pathogen.hub.more.10.csv"),quote = F, row.names = FALSE)


# Network betweennes --------------------------------------------------
b<-igraph::betweenness(net, v = V(net), directed = TRUE, weights = NULL,
                       nobigint = TRUE, normalized = FALSE)

b.p<-as.data.frame(b)
gene.name<-as.data.frame(row.names(b.p))
b.p<-cbind(gene.name,b.p)
colnames(b.p)[1]<-"Gene.name"
colnames(b.p)[2]<-"Betweennes"
Betweennes.pathogen<-b.p%>% filter(Betweennes>100)
write.csv(Betweennes.pathogen,file = paste0 ("./Results/","./Network_analysis/","Pathogen.betweennes.csv"),quote = F, row.names = FALSE)

# Network short distance --------------------------------------------------
d<-distances(
  net,
  v = V(net),
  to = V(net),
  mode = c("all", "out", "in"),
  weights = NULL,
  algorithm = c("automatic", "unweighted", "dijkstra", "bellman-ford", "johnson")
)
g<-as.data.frame(row.names(d))
colnames(g)[1]<-"Gene.name"
d<-cbind(g,d)
write.csv(d,file = paste0 ("./Results/","./Network_analysis/","Pathogen.shortest.path.csv"),quote = F, row.names = FALSE)


# Network Modules ---------------------------------------------------------
cl<-cluster_louvain(net, weights = NULL)
t<-as.data.frame(cl$membership)
t1<-as.data.frame(cl$names)
t2<-cbind(t1,t)
colnames(t2)[1]<-"Gene.name"
colnames(t2)[2]<-"membership"
write.csv(t2,file = paste0 ("./Results/","./Network_analysis/","Pathogen.modules.in.ppi.network.csv"),quote = F, row.names = FALSE)


# Plot Modules ------------------------------------------------------------
g_grouped = net

for(i in unique(V(net)$community)){
  groupV = which(V(net)$community == i)
  g_grouped = add_edges(g_grouped, combn(groupV, 2), attr=list(weight = 2))
}

l <- layout_nicely(g_grouped)

# Writing to file
pdf(file ="./Results/Network_analysis/Pathogen.modules.in.ppi.network.plot.pdf",width = 5, height = 5)
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
dev.off()

tiff(filename ="./Results/Network_analysis/Pathogen.modules.in.ppi.network.plot.tiff", compression = "lzw")
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

dev.off()



# t: Matrix Transpose -----------------------------------------------------

#### Select host gene set for correlation analysis based on the network analysis
m<-read.csv("Results/Average_data/All.host.mean.data.csv")
m<-Mean.host

d<-read.csv("Results/Average_data/All.host.ED.gene.csv")
d<-ED.host
d<-as.data.frame(d[,c(1)])
colnames(d)[1]<-"Gene.name"


h<-read.csv("Results/Network_analysis/Host.hub.more.100.csv")
h<-Hub.host

b<-read.csv("Results/Network_analysis/Host.betweennes.csv")
b<-Betweennes.host

e<-read.csv("Results/Mfuzz_Clustering/Variable.enriched.genes.and.their.frequency.in.host.csv")
e<-freq.list.host.gene

m.p<-read.csv("Results/Average_data/All.pathogen.mean.data.csv")
m.p<-Mean.pathogen


d.p<-read.csv("Results/Average_data/All.pathogen.ED.gene.csv")
d.p<-ED.pathogen
d.p<-as.data.frame(d.p[,c(1)])
colnames(d.p)[1]<-"Gene.name"


h.p<-read.csv("Results/Network_analysis/Pathogen.hub.more.10.csv")
h.p<-Hub.pathogen

b.p<-read.csv("Results/Network_analysis/Pathogen.betweennes.csv")
b.p<-Betweennes.pathogen

e.p<-read.csv("Results/Mfuzz_Clustering/Variable.enriched.genes.and.their.frequency.in.pathogen.csv")
e.p<-freq.list.pathogen.gene


### Relation between hub and mean
s<-merge(d,m,by="Gene.name")
s<-merge(s,h,by="Gene.name")
s<-merge(s,b,by="Gene.name")
s<-merge(s,e,by="Gene.name")

s.p<-merge(d.p,m.p,by="Gene.name")
s.p<-merge(s.p,h.p,by="Gene.name")
s.p<-merge(s.p,b.p,by="Gene.name")
s.p<-merge(s.p,e.p,by="Gene.name")

l<-plyr::rbind.fill(s,s.p)
l<-l[,c(1:7)]

### Get transpose matrix
l<-data.frame(t(l))
g<-as.data.frame(row.names(l))
l <- cbind(g,l)
colnames(l)[1]<-"Gene.name"
write.table(l,file = paste0("./Results/","./Matrix_Transpose/","Transpose.host.pathogen.gene.csv"),  sep=",",col.names= FALSE,row.names = FALSE)


# cor.analysis with all gene in study  ------------------------------------------------------------



c<-read.csv("Inputs/Select.gene.set.for.correlation.study/All.mean.transpose.for.cor.csv", header = T, row.names = 1, check.names = F)
All.cor<-psych::corr.test(c[1:293],c[294:422], method="pearson",adjust="holm", ci=FALSE)

r.df <- as.data.frame(All.cor$r)
R.corr <- r.df %>%
  dplyr::mutate(gene1 = row.names(r.df)) %>%
  tidyr::pivot_longer(-gene1,
               names_to = "gene2", names_ptypes = list(gene2=character()),
               values_to = "corr") %>%
  dplyr::mutate(gene1 = stringr::str_replace(gene1, "\\.\\.", " ("),
         gene1 = stringr::str_replace(gene1, "\\.$", ")"),
         gene2 = stringr::str_replace(gene2, "\\.\\.", " ("),
         gene2 = stringr::str_replace(gene2, "\\.$", ")")) %>%
  dplyr::mutate(comb = paste(gene1, "-", gene2))

#######################################################################  Collect negative correlation
negative<-R.corr %>% dplyr::filter(corr<(-0.7))
net<-igraph::graph.data.frame(unique(negative[,c(1,2)]),directed = FALSE)

cl<-igraph::cluster_louvain(net, weights = NULL)
t<-as.data.frame(cl$membership)
t1<-as.data.frame(cl$names)
t2<-cbind(t1,t)
colnames(t2)[1]<-"Gene.name"
colnames(t2)[2]<-"membership"

g_grouped = net

for(i in unique(igraph::V(net)$community)){
  groupV = which(igraph::V(net)$community == i)
  g_grouped = igraph::add_edges(g_grouped, combn(groupV, 2), attr=list(weight = 2))
}

l <- igraph::layout_nicely(g_grouped)



tiff(filename ="./Results/Correlation_analysis/Negative.highly.corrplot.plot.tiff", compression = "lzw")
plot(cl,net, layout = igraph::layout_with_fr,
     vertex.size =10,
     edge.width = 1,
     vertex.label.dist=0.001,
     vertex.color ='gold',
     vertex.frame.color="#555555",
     edge.label=net$v,
     vertex.size=1,
     edge.color="gray",
     vertex.label.font=0.5,
     edge.label.font =0.5,
     edge.label.cex = 0.5,
     edge.arrow.size=0.2,
     edge.curved=0,
     vertex.label=igraph::V(net)$v,
     vertex.label.color="black",
     vertex.label.cex=0.5,
     vertex.label.cex = 0.5 )
dev.off()

pdf(file ="./Results//Correlation_analysis/Negative.highly.corrplot.plot.tiff.pdf",width = 5, height = 5)
plot(cl,net, layout = igraph::layout_with_fr,
     vertex.size =10,
     edge.width = 1,
     vertex.label.dist=0.001,
     vertex.color ='gold',
     vertex.frame.color="#555555",
     edge.label=net$v,
     vertex.size=1,
     edge.color="gray",
     vertex.label.font=0.5,
     edge.label.font =0.5,
     edge.label.cex = 0.5,
     edge.arrow.size=0.2,
     edge.curved=0,
     vertex.label=igraph::V(net)$v,
     vertex.label.color="black",
     vertex.label.cex=0.5,
     vertex.label.cex = 0.5 )
dev.off()

#########################################################################################  Collect positive correlation
posetive<-R.corr%>% dplyr::filter(corr>0.7, corr<1)

net<-igraph::graph.data.frame(unique(posetive[,c(1,2)]),directed = FALSE)

# cor.analysis Network Modules ---------------------------------------------------------

cl<-igraph::cluster_louvain(net, weights = NULL)
t<-as.data.frame(cl$membership)
t1<-as.data.frame(cl$names)
t2<-cbind(t1,t)
colnames(t2)[1]<-"Gene.name"
colnames(t2)[2]<-"membership"

# cor.analysis Plot Cluster ------------------------------------------------------------
g_grouped = net

for(i in unique(igraph::V(net)$community)){
  groupV = which(igraph::V(net)$community == i)
  g_grouped = igraph::add_edges(g_grouped, combn(groupV, 2), attr=list(weight = 2))
}

l <- igraph::layout_nicely(g_grouped)

tiff(filename ="./Results/Correlation_analysis/Posetive.highly.corrplot.plot.tiff", compression = "lzw")
plot(cl,net, layout = igraph::layout_with_fr,
     vertex.size =10,
     edge.width = 1,
     vertex.label.dist=0.001,
     vertex.color ='gold',
     vertex.frame.color="#555555",
     edge.label=net$v,
     vertex.size=1,
     edge.color="gray",
     vertex.label.font=0.5,
     edge.label.font =0.5,
     edge.label.cex = 0.5,
     edge.arrow.size=0.2,
     edge.curved=0,
     vertex.label=igraph::V(net)$v,
     vertex.label.color="black",
     vertex.label.cex=0.5,
     vertex.label.cex = 0.5 )
dev.off()

pdf(file ="./Results//Correlation_analysis/Posetive.highly.corrplot.plot.tiff.pdf",width = 5, height = 5)
plot(cl,net, layout = igraph::layout_with_fr,
     vertex.size =10,
     edge.width = 1,
     vertex.label.dist=0.001,
     vertex.color ='gold',
     vertex.frame.color="#555555",
     edge.label=net$v,
     vertex.size=1,
     edge.color="gray",
     vertex.label.font=0.5,
     edge.label.font =0.5,
     edge.label.cex = 0.5,
     edge.arrow.size=0.2,
     edge.curved=0,
     vertex.label=igraph::V(net)$v,
     vertex.label.color="black",
     vertex.label.cex=0.5,
     vertex.label.cex = 0.5 )
dev.off()



# cor.analysis with selected gene in study  ------------------------------------------------------------

#between host and pathogen

c<-read.csv("Inputs/Select.gene.set.for.correlation.study/Transpose.both gene set.csv", header = T, row.names = 1, check.names = F)
All.cor<-psych::corr.test(c[1:54],c[1:54], method="pearson",adjust="holm", ci=FALSE)

r.df <- as.data.frame(All.cor$r)
R.corr <- r.df %>%
  dplyr::mutate(gene1 = row.names(r.df)) %>%
  tidyr::pivot_longer(-gene1,
               names_to = "gene2", names_ptypes = list(gene2=character()),
               values_to = "corr") %>%
  dplyr::mutate(gene1 = stringr::str_replace(gene1, "\\.\\.", " ("),
         gene1 = stringr::str_replace(gene1, "\\.$", ")"),
         gene2 = stringr::str_replace(gene2, "\\.\\.", " ("),
         gene2 = stringr::str_replace(gene2, "\\.$", ")")) %>%
  dplyr::mutate(comb = paste(gene1, "-", gene2))


### Writing to file
tiff(filename ="./Results/Correlation_analysis/Corrplot.plot.tiff", compression = "lzw")
corrplot::corrplot(as.matrix(r.df), is.corr = FALSE, method = "circle", order = "hclust")
dev.off()

pdf(file ="./Results//Correlation_analysis/Corrplot.plot.pdf",width = 5, height = 5);
corrplot::corrplot(as.matrix(r.df), is.corr = FALSE, method = "circle", order = "hclust")
dev.off()

write.csv(R.corr,file = paste0 ("./Results/","/Correlation_analysis/","R.corr.host.pathogen.csv"),quote = F, row.names = FALSE)

#until here added 'package::' for correlation analysis.
#rest seems same with different correlation cutoff values


############################################################################# Negative correlated 70
negative<-R.corr %>% filter(corr<(-0.7))

net<-graph.data.frame(unique(negative[,c(1,2)]),directed = FALSE)

cl<-cluster_louvain(net, weights = NULL)
t<-as.data.frame(cl$membership)
t1<-as.data.frame(cl$names)
t2<-cbind(t1,t)
colnames(t2)[1]<-"Gene.name"
colnames(t2)[2]<-"membership"

g_grouped = net

for(i in unique(V(net)$community)){
  groupV = which(V(net)$community == i)
  g_grouped = add_edges(g_grouped, combn(groupV, 2), attr=list(weight = 2))
}

l <- layout_nicely(g_grouped)


tiff(filename ="./Results/Correlation_analysis/Negative.70.highly.corrplot.plot.tiff", compression = "lzw")
plot(cl,net, layout = layout_with_fr,
     vertex.size =10,
     edge.width = 1,
     vertex.label.dist=0.001,
     vertex.color ='gold',
     vertex.frame.color="#555555",
     edge.label=net$v,
     vertex.size=1,
     edge.color="gray",
     vertex.label.font=0.5,
     edge.label.font =0.5,
     edge.label.cex = 0.5,
     edge.arrow.size=0.2,
     edge.curved=0,
     vertex.label=V(net)$v,
     vertex.label.color="black",
     vertex.label.cex=0.5,
     vertex.label.cex = 0.5 )
dev.off()

pdf(file ="./Results//Correlation_analysis/Negative.70.highly.corrplot.plot.tiff.pdf",width = 5, height = 5)
plot(cl,net, layout = layout_with_fr,
     vertex.size =10,
     edge.width = 1,
     vertex.label.dist=0.001,
     vertex.color ='gold',
     vertex.frame.color="#555555",
     edge.label=net$v,
     vertex.size=1,
     edge.color="gray",
     vertex.label.font=0.5,
     edge.label.font =0.5,
     edge.label.cex = 0.5,
     edge.arrow.size=0.2,
     edge.curved=0,
     vertex.label=V(net)$v,
     vertex.label.color="black",
     vertex.label.cex=0.5,
     vertex.label.cex = 0.5 )
dev.off()

############################################################################# Negative correlation 90
negative<-R.corr %>% filter(corr<(-0.9))

net<-graph.data.frame(unique(negative[,c(1,2)]),directed = FALSE)

cl<-cluster_louvain(net, weights = NULL)
t<-as.data.frame(cl$membership)
t1<-as.data.frame(cl$names)
t2<-cbind(t1,t)
colnames(t2)[1]<-"Gene.name"
colnames(t2)[2]<-"membership"

g_grouped = net

for(i in unique(V(net)$community)){
  groupV = which(V(net)$community == i)
  g_grouped = add_edges(g_grouped, combn(groupV, 2), attr=list(weight = 2))
}

l <- layout_nicely(g_grouped)



tiff(filename ="./Results/Correlation_analysis/Negative.90.highly.corrplot.plot.tiff", compression = "lzw")
plot(cl,net, layout = layout_with_fr,
     vertex.size =10,
     edge.width = 1,
     vertex.label.dist=0.001,
     vertex.color ='gold',
     vertex.frame.color="#555555",
     edge.label=net$v,
     vertex.size=1,
     edge.color="gray",
     vertex.label.font=0.5,
     edge.label.font =0.5,
     edge.label.cex = 0.5,
     edge.arrow.size=0.2,
     edge.curved=0,
     vertex.label=V(net)$v,
     vertex.label.color="black",
     vertex.label.cex=0.5,
     vertex.label.cex = 0.5 )
dev.off()

pdf(file ="./Results//Correlation_analysis/Negative.90.highly.corrplot.plot.tiff.pdf",width = 5, height = 5)
plot(cl,net, layout = layout_with_fr,
     vertex.size =10,
     edge.width = 1,
     vertex.label.dist=0.001,
     vertex.color ='gold',
     vertex.frame.color="#555555",
     edge.label=net$v,
     vertex.size=1,
     edge.color="gray",
     vertex.label.font=0.5,
     edge.label.font =0.5,
     edge.label.cex = 0.5,
     edge.arrow.size=0.2,
     edge.curved=0,
     vertex.label=V(net)$v,
     vertex.label.color="black",
     vertex.label.cex=0.5,
     vertex.label.cex = 0.5 )
dev.off()


############################################################################## Positive correlated 70
posetive<-R.corr%>% filter(corr>0.7, corr<1)

net<-graph.data.frame(unique(posetive[,c(1,2)]),directed = FALSE)

# cor.analysis Network Modules ---------------------------------------------------------

cl<-cluster_louvain(net, weights = NULL)
t<-as.data.frame(cl$membership)
t1<-as.data.frame(cl$names)
t2<-cbind(t1,t)
colnames(t2)[1]<-"Gene.name"
colnames(t2)[2]<-"membership"

# cor.analysis Plot Cluster ------------------------------------------------------------
g_grouped = net

for(i in unique(V(net)$community)){
  groupV = which(V(net)$community == i)
  g_grouped = add_edges(g_grouped, combn(groupV, 2), attr=list(weight = 2))
}

l <- layout_nicely(g_grouped)

tiff(filename ="./Results/Correlation_analysis/Posetive.70.highly.corrplot.plot.tiff", compression = "lzw")
plot(cl,net, layout = layout_with_fr,
     vertex.size =10,
     edge.width = 1,
     vertex.label.dist=0.001,
     vertex.color ='gold',
     vertex.frame.color="#555555",
     edge.label=net$v,
     vertex.size=1,
     edge.color="gray",
     vertex.label.font=0.5,
     edge.label.font =0.5,
     edge.label.cex = 0.5,
     edge.arrow.size=0.2,
     edge.curved=0,
     vertex.label=V(net)$v,
     vertex.label.color="black",
     vertex.label.cex=0.5,
     vertex.label.cex = 0.5 )
dev.off()

pdf(file ="./Results//Correlation_analysis/Posetive.70.highly.corrplot.plot.tiff.pdf",width = 5, height = 5)
plot(cl,net, layout = layout_with_fr,
     vertex.size =10,
     edge.width = 1,
     vertex.label.dist=0.001,
     vertex.color ='gold',
     vertex.frame.color="#555555",
     edge.label=net$v,
     vertex.size=1,
     edge.color="gray",
     vertex.label.font=0.5,
     edge.label.font =0.5,
     edge.label.cex = 0.5,
     edge.arrow.size=0.2,
     edge.curved=0,
     vertex.label=V(net)$v,
     vertex.label.color="black",
     vertex.label.cex=0.5,
     vertex.label.cex = 0.5 )
dev.off()

############################################################################ Positive correlated 90

posetive<-R.corr%>% filter(corr>0.9, corr<1)

net<-graph.data.frame(unique(posetive[,c(1,2)]),directed = FALSE)

# cor.analysis Network Modules ---------------------------------------------------------

cl<-cluster_louvain(net, weights = NULL)
t<-as.data.frame(cl$membership)
t1<-as.data.frame(cl$names)
t2<-cbind(t1,t)
colnames(t2)[1]<-"Gene.name"
colnames(t2)[2]<-"membership"

# cor.analysis Plot Cluster ------------------------------------------------------------
g_grouped = net

for(i in unique(V(net)$community)){
  groupV = which(V(net)$community == i)
  g_grouped = add_edges(g_grouped, combn(groupV, 2), attr=list(weight = 2))
}

l <- layout_nicely(g_grouped)

tiff(filename ="./Results/Correlation_analysis/Posetive.90.highly.corrplot.plot.tiff", compression = "lzw")
plot(cl,net, layout = layout_with_fr,
     vertex.size =10,
     edge.width = 1,
     vertex.label.dist=0.001,
     vertex.color ='gold',
     vertex.frame.color="#555555",
     edge.label=net$v,
     vertex.size=1,
     edge.color="gray",
     vertex.label.font=0.5,
     edge.label.font =0.5,
     edge.label.cex = 0.5,
     edge.arrow.size=0.2,
     edge.curved=0,
     vertex.label=V(net)$v,
     vertex.label.color="black",
     vertex.label.cex=0.5,
     vertex.label.cex = 0.5 )
dev.off()

pdf(file ="./Results//Correlation_analysis/Posetive.90.highly.corrplot.plot.tiff.pdf",width = 5, height = 5)
plot(cl,net, layout = layout_with_fr,
     vertex.size =10,
     edge.width = 1,
     vertex.label.dist=0.001,
     vertex.color ='gold',
     vertex.frame.color="#555555",
     edge.label=net$v,
     vertex.size=1,
     edge.color="gray",
     vertex.label.font=0.5,
     edge.label.font =0.5,
     edge.label.cex = 0.5,
     edge.arrow.size=0.2,
     edge.curved=0,
     vertex.label=V(net)$v,
     vertex.label.color="black",
     vertex.label.cex=0.5,
     vertex.label.cex = 0.5 )
dev.off()



