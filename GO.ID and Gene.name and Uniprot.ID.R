data <- read.csv("Bio_uni.csv")
genes <- unique(data$Gene.name)


tt <- sapply(1:length(genes), function(i) paste(gsub("[[:space:]]", "", data[which(data$Gene.name==genes[i]),]$Gene.Ontology.Terms),collapse=";"))
tt.1 <- sapply(1:3,
               function(i) paste(data[which(data$Gene.name==genes[i]),]$Gene.ontology.IDs))

tt.2 <- sapply(1:length(genes), function(i) as.character(data[which(data$Gene.name==genes[i]),]$Gene.Ontology.Terms))


tt.3 <- sapply(1:3, function(i) paste(gsub("[[:space:]]","", unique(c(unique(unlist(strsplit(tt.1[[i]],"; "))),tt.2[[i]]))),collapse=";"))



tt[1:3] <- ""
result = paste(tt,tt.3)

result<-as.data.frame(result)

df<-data.frame(Gene.name=unique(data$Gene.name),
               Uniprot.ID=unique(data$Uniprot.ID),
               GO.term=result)

write.csv(df,file ="DF_Bio_uni.csv", row.names = FALSE)
