alldt2h <- read.csv(file = "Alldt2h.csv", header = TRUE)
alldt4h <- read.csv(file = "Alldt4h.csv", header = TRUE)
alldt8h <- read.csv(file = "Alldt8h.csv", header = TRUE)
alldt16h <- read.csv(file = "Alldt16h.csv", header = TRUE)
alldt24h <- read.csv(file = "Alldt24h.csv", header = TRUE)




# ## Plotting-------------------------------------------------------------------------



# -----------------------Boxplot ~ Gordana ------------------------------------------


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






# ---------------------ENHANCED VOLCANO PLOT----------------------------------------------------

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


# 2h ---------------------



EnhancedVolcano(alldt2h,
                lab = rownames(alldt2h),
                x = 'LogFC2h',
                y = 'Adj.P.Val2h',
                xlim=c(-6,6, 1),
                xlab = bquote(~Log[2]~ 'fold change'),
                ylab = bquote(~-Log[10]~adjusted~italic(P)),
                pCutoff = 0.05,
                shadeAlpha = .5,
                FCcutoff = 1.0,
                labSize = 4.0,
                colAlpha = 1,
                pointSize = 1.0,
                legend=c('NS','Log2 FC','Adjusted p-value',
                         'S'),
                legendPosition = 'bottom',
                legendLabSize = 10,
                legendIconSize = 3.0)







# 4h ---------------------




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



# 8h ----------------------


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

# 16h ------------------------


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

# 24h ---------------------------


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





# -------------------------HEATMAP------------------------------------------------



mean_expr <- cbind(alldt2h[,c(1,3)], alldt4h$AveExpr4h, alldt8h$AveExpr8h, alldt16h$AveExpr16h, alldt24h$AveExpr24h)
colnames(mean_expr)<- c("Attributes", "2h","4h","8h","16h", "24h")
class(mean_expr)
#[1] "data.frame"

mean_expr_matrix <-data.matrix(mean_expr)
class(mean_expr_matrix)
#[1] "matrix"

mean_heatmap <- heatmap(mean_expr_matrix, Rowv = NULL, Colv= NULL, col=cm.colors(256), scale= "column", margins=c(5,10))
#or, mean_heatmap <- heatmap(mean_expr_matrix, Rowv = NULL, Colv= NULL, col=heat.colors(256), scale= "column", margins=c(5,10))

heatmap3(mean_expr_matrix,showRowDendro=FALSE,colorCell=colorCell,
         highlightCell=highlightCell)