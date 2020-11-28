dat<-read.csv(file = "GSE67758_Raw_gene_wise_quantifications_Human_and_Salmonella.csv", header = TRUE,sep = "\t")
dat.host=dat[grep("chr",dat$Sequence.name,fixed = TRUE),]

dat.host = dat.host[which(dat.host[,1]=="sense"),]
dat.host.gene = dat.host[grep("gene_status=KNOWN", dat.host$Attributes, fixed = T),]

dat.host.attributes<-dat.host.gene[,-c(1:9,11:28,47:61)]
Attributes <-dat.host.attributes[,c(1)]
Attributes<-str_split_fixed(Attributes,";", n=Inf)
Attributes<-Attributes[ , c(2)]
Attributes<-str_split_fixed(Attributes,"=", n=Inf)
Attributes<-Attributes[,-c(1)]
dat.host.attributes<-cbind(dat.host.attributes[,1],Attributes,dat.host.attributes[,2:19])
dat.host.attributes<-dat.host.attributes[,-c(1)]


write.csv(dat.host.attributes, file = "dat.H.csv", row.names = FALSE)




#mean of same uniprotID
raw.agg<-aggregate(.~ Attributes,dat.host.attributes, FUN="mean")
#making rownames 1,2,3 as uniprotID name
rownames(raw.agg) <- raw.agg[,1]
raw.agg <- raw.agg[,-1]
dat.host.uniq<-raw.agg
dt=dat.host.uniq[rowSums(is.na(dat.host.uniq[,-1])) !=ncol(dat.host.uniq[,-1]),]

data <-gsub('HeLa.S3.','',colnames(dt))

##Normdt <- read.csv("y.norm.csv", header = TRUE)


library(limma)
library(DESeq2)
library(dplyr)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(igraph)
library(network)
library(clusterProfiler)
library(data.table)
library(AnnotationDbi)
library(edgeR)
library(stringr)
library(igraph)
library(stringi)

library(EnhancedVolcano)
library(ggplot2)
library(BiocManager)


# Limma -------------------------------------------------------------------

#data.filter<-dat.host.na
#dt<-t(data.filter)


conditions<-gsub('HeLa.S3.','',colnames(dt))

c<-gsub('.R[0-9]','',conditions)


#dt <- dt + 10e-10 # option: try to add epsilon only to zeros
y <- DGEList(counts=dt,group=c,genes=rownames(dt))
y <- calcNormFactors(y, method="TMM")
y.norm <-log(cpm(y))
#replace -Inf with zero
keep = which(y.norm == -Inf)
y.norm[keep] = 0
dat.H.norm <- as.data.frame(y.norm)

write.csv(dat.host.attributes, file = "dat.H.norm.csv", row.names = FALSE)

# Limma DE analysis -------------------------------------------------------



require(limma)
group = as.factor(c)
design <- model.matrix(~ 0 + group)
colnames(design) <- levels(group)
rownames(design) <- colnames(y.norm)
fit <- lmFit(y.norm, design = design)
cont.matrix <- makeContrasts(WT.02_h - WT.00_h,
                             WT.04_h - WT.00_h,
                             WT.08_h - WT.00_h,
                             WT.16_h - WT.00_h,
                             WT.24_h - WT.00_h,levels=design)
fit.cont <- eBayes(contrasts.fit(fit, cont.matrix))
summa.fit <- decideTests(fit.cont)
summary(summa.fit)

#Host_DE_summary <- as.data.frame(summary(summa.fit))
#write.csv(Host_DE_summary, file = "Host_DE_summary.csv", row.names = FALSE)

# topTable(fit.cont, coef = 2)

out.Wt2h_0h <- topTable(fit.cont, n = Inf, coef = 1)
# add zero counts to each table
out.Wt4h_0h <- topTable(fit.cont, n = Inf, coef = 2)
# repeat for the rest
out.Wt8h_0h <- topTable(fit.cont, n = Inf, coef = 3)
out.Wt16h_0h <- topTable(fit.cont, n = Inf, coef = 4)
out.Wt24h_0h <- topTable(fit.cont, n = Inf, coef = 5)

#all.out = cbind(out.Wt2h,out.Wt4h)
write.csv(out.Wt2h_0h, file = "Host_2h_vs_0h.csv")
write.csv(out.Wt4h_0h, file = "Host_4h_vs_0h.csv")
write.csv(out.Wt8h_0h, file = "Host_8h_vs_0h.csv")
write.csv(out.Wt16h_0h, file = "Host_16h_vs_0h.csv")
write.csv(out.Wt24h_0h, file = "Host_24h_vs_0h.csv") 


#UPREGULATED ONES ONLY
Wt2h.up = subset(out.Wt2h, logFC > 0 & adj.P.Val < 0.05)
Wt4h.up = subset(out.Wt4h, logFC > 0 & adj.P.Val < 0.05)
Wt8h.up = subset(out.Wt8h, logFC > 0 & adj.P.Val < 0.05)
Wt16h.up = subset(out.Wt16h, logFC > 0 & adj.P.Val < 0.05)
Wt24h.up = subset(out.Wt24h, logFC > 0 & adj.P.Val < 0.05)




#Addind zero count in the table 
count0<-function(y.norm,col){
  x=y.norm[,col]
  x0=x==0
  n_0=apply(x0,1,sum)
  out=data.frame(n_0=n_0)
  out
}


y=count0(y.norm, col = 1:6)
y.norm=cbind(y.norm,zeros02=y[,1])

y=count0(y.norm, c(1:3,7:9))
y.norm=cbind(y.norm,zeros04=y[,1])

y=count0(y.norm, c(1:3,10:12))
y.norm=cbind(y.norm,zeros08=y[,1])

y=count0(y.norm, c(1:3,13:15))
y.norm=cbind(y.norm,zeros016=y[,1])

y=count0(y.norm, c(1:3,16:18))
y.norm=cbind(y.norm,zeros024=y[,1])

write.csv(y.norm, file = "Host_zerocount_y.norm.csv")



#GETTING MEAN OF REPLICATES
y.norm_ <-y.norm[, -c(19:24)]
y.norm_ <-as.data.frame(y.norm_)

names <- rownames(y.norm_)
rownames(y.norm_) <- NULL
data <- cbind(names,y.norm_)

data_mean=data.frame(apply(array(as.matrix(data[,-1]), c(nrow(data),3, ncol(data)/3)),3, rowMeans))
data_mean=cbind(data$names, data_mean)
colnames(data_mean)<- c("Attributes", "Mean_0h","Mean_2h","Mean_4h","Mean_8h","Mean_16h", "Mean_24h")
write.csv(data_mean, file = "Host Normalized data mean.csv", row.names = FALSE)

write.csv(data_mean, file = "data mean.csv", row.names = FALSE)




# End ---------------------------------------------------------------------

#Matrix to df

y.norm_zc <- tibble::rownames_to_column(y.norm, "VALUE")


#Rownames to Column
out.Wt2h_zc <- tibble::rownames_to_column(out.Wt2h, "Attributes")



#Combine by Value


both_dat=inner_join(out.Wt2h,y.norm,by= rownamnes
                    )
names <- rownames(y.norm)
rownames(y.norm) <- NULL
y.norm <- cbind(names,y.norm)

y.norm<- read.csv(file = "y.norm.csv", header = TRUE)
y.norm=y.norm [,20:23]
zero_count=y.norm[grep(y.norm,fixed = TRUE),]

# -------------------------------------------------------------------------




unq<-unique(c)
a<-unq[1]
res <- c()
for(b in unq[-1]){
  t.res <- exactTest(y,c(a,b), dispersion = "common")
  res <- cbind(res,t.res$table['PValue']$PValue)
}
rownames(res) <- names(dt)[-1]
colnames(res) <- unq[-1]
write.csv(res,file = "res.csv")

#cpm.mat<-log(cpm(exprs(data.filter)))
#keep<-rowSums(cpm(pathogen.zer.na.nor[,2:34])>5)>=3
#data_filt<-pathogen.zer.na.nor[keep,]
#data_filt1<-cbind(data_filt,Attributes,data_filt[,2:34])
#filter <- apply(dat.host.na, 1, function(x){length(x[x>10])>=3}) (based on the paper information # select rows with expression at least 10 in at least 3 replicates)
# var_rows = apply(data_filt[,-1], 1, var)                             (removal of rows with little variance)    
# data_filt= data_filt[-which(var_rows <= 2), ]
# data_filt <- apply(data_filt[,-1],1,sum)==0                           (remove rows whose sum is 0)

out.Wt2h_zc <- tibble::rownames_to_column(out.Wt2h, "VALUE")





# -------------------------------------------------------------------------


alldt2h <- read.csv(file = "Alldt2h.csv", header = TRUE)
alldt4h <- read.csv(file = "Alldt4h.csv", header = TRUE)
alldt8h <- read.csv(file = "Alldt8h.csv", header = TRUE)
alldt16h <- read.csv(file = "Alldt16h.csv", header = TRUE)
alldt24h <- read.csv(file = "Alldt24h.csv", header = TRUE)

mean.data <- read.csv(file = "Mean of sample replicates.csv", header = TRUE)
rownames(mean.data)<- mean.data$Attributes
mean.data<- mean

# Heatmap -----------------------------------------------------------------

mean_expr <- cbind(alldt2h[,c(1,3)], alldt4h$AveExpr4h, alldt8h$AveExpr8h, alldt16h$AveExpr16h, alldt24h$AveExpr24h)
colnames(mean_expr)<- c("Attributes", "2h","4h","8h","16h", "24h")
class(mean_expr)
#[1] "data.frame"

mean_expr_matrix <-data.matrix(mean_expr)
class(mean_expr_matrix)
#[1] "matrix"

mean_heatmap <- heatmap(mean.data, Rowv = NULL, Colv= NULL, col=cm.colors(256), scale= "column", margins=c(5,10))
#or, mean_heatmap <- heatmap(mean_expr_matrix, Rowv = NULL, Colv= NULL, col=heat.colors(256), scale= "column", margins=c(5,10))

heatmap3(mean_expr_matrix,showRowDendro=FALSE,colorCell=colorCell,
         highlightCell=highlightCell)





# Boxplot ~ Gordana -----------------------------------------------------------------


alldt2h$col="C"
alldt2h$col[alldt2h$LogFC2h>2&alldt2h$Adj.P.Val2h>0.05]="A"
alldt2h$col[alldt2h$LogFC2h<=2&alldt2h$Adj.P.Val2h>0.05]="B"

alldt2h$col=factor(alldt2h$col,levels=c("A","B","C"),labels = c("up","down","not sig"))

ggplot(alldt2h,aes(x=LogFC2h,y=-Adj.P.Val2h,color=col))+geom_point()+
  scale_color_manual(values=c("green","red","grey"))+
  theme_classic()

ggplot(alldt2h,aes(x=col,y=Adj.P.Val2h,fill=col, xlab= "logFC"))+geom_boxplot(title= "Boxplot 2h")


# 4h ------------------------

alldt4h$col="C"
alldt4h$col[alldt4h$LogFC4h>2&alldt4h$Adj.P.Val4h>0.05]="A"
alldt4h$col[alldt4h$LogFC4h<=2&alldt4h$Adj.P.Val4h>0.05]="B"

alldt4h$col=factor(alldt4h$col,levels=c("A","B","C"),labels = c("up","down","not sig"))

ggplot(alldt4h,aes(x=LogFC4h,y=-Adj.P.Val4h,color=col))+geom_point()+
  scale_color_manual(values=c("green","red","grey"))+
  theme_classic()

ggplot(alldt4h,aes(x=col,y=Adj.P.Val4h,fill=col, xlab= "logFC"))+geom_boxplot(title= "Boxplot 4h")


# 8h -----------------------


alldt8h$col="C"
alldt8h$col[alldt8h$LogFC8h>2&alldt8h$Adj.P.Val8h>0.05]="A"
alldt8h$col[alldt8h$LogFC8h<=2&alldt8h$Adj.P.Val8h>0.05]="B"

alldt8h$col=factor(alldt8h$col,levels=c("A","B","C"),labels = c("up","down","not sig"))

ggplot(alldt8h,aes(x=LogFC8h,y=-Adj.P.Val8h,color=col))+geom_point()+
  scale_color_manual(values=c("green","red","grey"))+
  theme_classic()

ggplot(alldt8h,aes(x=col,y=Adj.P.Val8h,fill=col, xlab= "logFC"))+geom_boxplot(title= "Boxplot 8h")

# 16h ----------------------

alldt16h$col="C"
alldt16h$col[alldt16h$LogFC16h>2&alldt16h$Adj.P.Val16h>0.05]="A"
alldt16h$col[alldt16h$LogFC16h<=2&alldt16h$Adj.P.Val16h>0.05]="B"

alldt16h$col=factor(alldt16h$col,levels=c("A","B","C"),labels = c("up","down","not sig"))

ggplot(alldt16h,aes(x=LogFC16h,y=-Adj.P.Val16h,color=col))+geom_point()+
  scale_color_manual(values=c("green","red","grey"))+
  theme_classic()

ggplot(alldt16h,aes(x=col,y=Adj.P.Val16h,fill=col, xlab= "logFC"))+geom_boxplot(title= "Boxplot 16h")

# 24h ---------------------


alldt24h$col="C"
alldt24h$col[alldt24h$LogFC24h>2&alldt24h$Adj.P.Val24h>0.05]="A"
alldt24h$col[alldt24h$LogFC24h<=2&alldt24h$Adj.P.Val24h>0.05]="B"

alldt24h$col=factor(alldt24h$col,levels=c("A","B","C"),labels = c("up","down","not sig"))

ggplot(alldt24h,aes(x=LogFC24h,y=-Adj.P.Val24h,color=col))+geom_point()+
  scale_color_manual(values=c("green","red","grey"))+
  theme_classic()

ggplot(alldt24h,aes(x=col,y=Adj.P.Val24h,fill=col, xlab= "logFC"))+geom_boxplot(title= "Boxplot 24h")






# -------------------------------------------------------------------------

if (!requireNamespace('BiocManager', quietly = TRUE))
  install.packages('BiocManager')
BiocManager::install('EnhancedVolcano')


library(EnhancedVolcano)



row.names(alldt2h) <- alldt2h$Attributes2h
alldt2h <- alldt2h[,2:10]


row.names(alldt4h) <- alldt4h$Attributes4h
alldt4h <- alldt4h[,2:10]


row.names(alldt8h) <- alldt8h$Attributes8h
alldt8h <- alldt8h[,2:10]

row.names(alldt16h) <- alldt16h$Attributes16h
alldt16h <- alldt16h[,2:10]

row.names(alldt24h) <- alldt24h$Attributes24h
alldt24h <- alldt24h[,2:10]


# 2h ----------------------------------------------------------------------



EnhancedVolcano(alldt2h,
                lab = rownames(alldt2h),
                x = 'LogFC2h',
                y = 'Adj.P.Val2h',
                xlim=c(-6,6),
                xlab = bquote(~Log[2]~ 'fold change'),
                ylab = bquote(~-Log[10]~adjusted~italic(P)),
                pCutoff = 0.05,
                FCcutoff = 1.0,
                labSize = 4.0,
                colAlpha = 1,
                legend=c('NS','Log2 FC','Adjusted p-value',
                         'Adjusted p-value & Log2 FC'),
                legendPosition = 'bottom',
                legendLabSize = 10,
                legendIconSize = 3.0)







# 4h ----------------------------------------------------------------------




EnhancedVolcano(alldt4h,
                lab = rownames(alldt4h),
                x = 'LogFC4h',
                y = 'Adj.P.Val4h',
                xlim=c(-6,6),
                xlab = bquote(~Log[2]~ 'fold change'),
                ylab = bquote(~-Log[10]~adjusted~italic(P)),
                pCutoff = 0.05,
                FCcutoff = 1.0,
                labSize = 4.0,
                colAlpha = 1,
                legend=c('NS','Log2 FC','Adjusted p-value',
                         'Adjusted p-value & Log2 FC'),
                legendPosition = 'bottom',
                legendLabSize = 10,
                legendIconSize = 3.0)



# 8h ----------------------------------------------------------------------


EnhancedVolcano(alldt8h,
                lab = rownames(alldt8h),
                x = 'LogFC8h',
                y = 'Adj.P.Val8h',
                xlim=c(-6,6),
                xlab = bquote(~Log[2]~ 'fold change'),
                ylab = bquote(~-Log[10]~adjusted~italic(P)),
                pCutoff = 0.05,
                FCcutoff = 1.0,
                labSize = 4.0,
                colAlpha = 1,
                legend=c('NS','Log2 FC','Adjusted p-value',
                         'Adjusted p-value & Log2 FC'),
                legendPosition = 'bottom',
                legendLabSize = 10,
                legendIconSize = 3.0)

# 16h ---------------------------------------------------------------------


EnhancedVolcano(alldt16h,
                lab = rownames(alldt16h),
                x = 'LogFC16h',
                y = 'Adj.P.Val16h',
                xlim=c(-6,6),
                xlab = bquote(~Log[2]~ 'fold change'),
                ylab = bquote(~-Log[10]~adjusted~italic(P)),
                pCutoff = 0.05,
                FCcutoff = 1.0,
                labSize = 4.0,
                colAlpha = 1,
                legend=c('NS','Log2 FC','Adjusted p-value',
                         'Adjusted p-value & Log2 FC'),
                legendPosition = 'bottom',
                legendLabSize = 10,
                legendIconSize = 3.0)

# 24h ---------------------------------------------------------------------


EnhancedVolcano(alldt24h,
                lab = rownames(alldt24h),
                x = 'LogFC24h',
                y = 'Adj.P.Val24h',
                xlim=c(-6,6),
                xlab = bquote(~Log[2]~ 'fold change'),
                ylab = bquote(~-Log[10]~adjusted~italic(P)),
                pCutoff = 0.05,
                FCcutoff = 1.0,
                labSize = 4.0,
                colAlpha = 1,
                legend=c('NS','Log2 FC','Adjusted p-value',
                         'Adjusted p-value & Log2 FC'),
                legendPosition = 'bottom',
                legendLabSize = 10,
                legendIconSize = 3.0)




# -------------------------------------------------------------------------


# how to do Sets: union, intersect and setdiff ----------------------------

#if the column values are already TRUE-FALSE
alldt2h.up = alldt2h[alldt2h$Up2h, ]

#if the column values are as decimal/numeric
alldt2h.up = alldt2h[which(alldt2h$Up2h == TRUE), ]

#now APPLY SETDIFF for 2h
temp = setdiff(as.character(alldt2h.up$Attributes), as.character(alldt4h.up$Attributes4h))

alldt2h.diff.4h.up = alldt2h[which(alldt2h$Attributes2h %in% temp),]


# Temporal det diff-----------------------------------------------------------------------

  a <- read.csv(file = "2h_vs_0h.csv", header = TRUE)
  b <- read.csv(file = "4h_vs_2h.csv", header = TRUE)
  c <- read.csv(file = "8h_vs_4h.csv", header = TRUE)
  d <- read.csv(file = "16h_vs_8h.csv", header = TRUE)
  e <- read.csv(file = "24h_vs_16h.csv", header = TRUE)

  temp = setdiff(as.character(a$X), as.character(b$X))

 ab= a[which(a$X %in% temp),]
 