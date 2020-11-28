GO_geneSetEnrichmentAnalysis <- function(moduleID, clusters = clusters, min.GeneSet=5, max.GeneSet=60,  pop.filepath = pop.filepath){
  # require(org.Hs.eg.db)
  require(stringi)
  require(AnnotationDbi)
  require(dplyr)
  require(data.table)
  require(clusterProfiler)
  
  # set of genes in the current cluster
  givenSet = clusters[[moduleID]]
  # current cluster name
  moduleName = names(clusters[moduleID])
  
  # process the GO term annotaion file [remove too specific or too generic terms]
  num.col = max(count.fields(pop.filepath, sep = "\t", blank.lines.skip = TRUE), na.rm = TRUE)
  GO_terms = read.table(pop.filepath, sep = "\t", header = F, stringsAsFactors = FALSE, fill = TRUE, col.names = 1:num.col, blank.lines.skip = TRUE, quote = "")
  newGO_terms = GO_terms[-which(rowSums(GO_terms != "", na.rm = T) <= min.GeneSet),] # throw away terms with minimum size
  newGO_terms = newGO_terms[-which(rowSums(newGO_terms != "", na.rm = T) >= max.GeneSet),] # throw away terms with maximum size
  
  all.GO.terms = newGO_terms
  row.names(all.GO.terms) = gsub(")","",lapply(lapply(stri_split_fixed(stri_reverse(all.GO.terms[,1]),"(",n = 2), FUN = stri_reverse), "[[",1))
  nTerms = nrow(all.GO.terms)
  pop.genes = unique(c(as.matrix(all.GO.terms[,-1])))
  
  N = length(pop.genes)
  
  # enrichment test using HyperGeometric test
  info.gene.overrep = data.frame(matrix(Inf, nrow = nrow(all.GO.terms), ncol = 7))
  
  K = length(givenSet)
  
  
  for (i in 1:nTerms) {
    # message(i)
    GO.term.genes.symbol = as.character(all.GO.terms[i,which(all.GO.terms[i,] != '')])[-1]
    M = length(GO.term.genes.symbol)
    # print(M)
    x.overlap.genes.symbol = intersect(GO.term.genes.symbol, givenSet)
    
    x = length(x.overlap.genes.symbol)
    # print(x)
    info.gene.overrep[i,1] = paste0("Cluster_", moduleName)
    info.gene.overrep[i,2] = all.GO.terms[i,1]
    if(x > 0){
      info.gene.overrep[i,3] = paste0(GO.term.genes.symbol, collapse = "/")
      info.gene.overrep[i,4] = paste0(x.overlap.genes.symbol, collapse = "/")
      info.gene.overrep[i,5] = x
      info.gene.overrep[i,6] = phyper(x, M, N-M, K, lower.tail = FALSE) #He wrote : overlap.genes-1
      #insert FDR val
      info.gene.overrep[i,7] = p.adjust(info.gene.overrep[i,6], method = "bonferroni", n = nTerms)
    }
  }
  colnames(info.gene.overrep) = c("ClusterID", "GO_Term", "GO_term.geneSymbol","Overlapping.geneSymbol", "Count", "pvalue", "p.adjust")
  # only return those overlapping rows that has at least one overlapping genes (i.e. x > 0)
  return(info.gene.overrep[which(info.gene.overrep$pvalue != Inf),])
}

# ----------------- 

# read GO Annotation data
library(data.table)
library(dplyr)
library(stringr)
dat <- as.data.frame(fread("Annotation.enrichment.csv"))
GO.dat <- dat[,c(2,4,6)] %>% group_by(Gene.Ontology.Terms) %>% 
            summarise_each(funs(paste0(., collapse = "\t"))) %>% 
            as.data.frame()
fwrite(GO.dat[,c(1,2)], file = "GO_uniprots_SL1344.csv", sep="\t", col.names = F,  quote = F)

# Read clusters and Preprocess
newClusters.dat <- as.data.frame(fread("Temporal cluster genelist 0.125.csv", header = T))
newClusters <- sapply(newClusters.dat[,2], function(x){
  return(str_split(x, ","))
})
names(newClusters) <- as.character(str_replace(newClusters.dat[,1], "Cluster",""))

#make list of lists for annotation
GO <- sapply(GO.dat[,2], function(x){
  return(str_split(x, "\t"))
})
names(GO) <- as.character(str_replace(GO.dat[,1], "Gene.Ontology.terms",""))


pop.filepath = "GO_uniprots_SL1344.csv"
enrichmentResult <- do.call("rbind", lapply(seq_len(length(newClusters)), function(i){
  return(GO_geneSetEnrichmentAnalysis(i, newClusters, pop.filepath = pop.filepath))
}))
fwrite(enrichmentResult, file="GO_enrichmentResults_TempCLusters_0.3.csv")
