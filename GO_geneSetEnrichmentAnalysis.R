GO_geneSetEnrichmentAnalysis <- function(givenSetFile = givenSetFile, min.GeneSet=15, max.GeneSet=100,  pop.filepath = pop.filepath, outFilePath = outFilePath){
  require(org.Hs.eg.db)
  require(stringi)
  require(AnnotationDbi)
  require(dplyr)
  require(data.table)
  require(clusterProfiler)
  
  
  givenSet = read.table(givenSetFile, sep = "\t", header = T)
  givenSet = as.character(unique(givenSet$Gene.Symbol))
  
  num.col = max(count.fields(pop.filepath, sep = "\t", blank.lines.skip = TRUE), na.rm = TRUE)
  GO_terms = read.table(pop.filepath, sep = "\t", header = FALSE, stringsAsFactors = FALSE, fill = TRUE, col.names = 1:num.col, blank.lines.skip = TRUE, quote = "")
  newGO_terms = GO_terms[-which(rowSums(GO_terms != "", na.rm = T)-1 <= min.GeneSet),]
  newGO_terms = newGO_terms[-which(rowSums(newGO_terms != "", na.rm = T)-1 >= max.GeneSet),]
  newGO_terms = dplyr::select(newGO_terms, -2) # throw away the the column with all NA (2nd-Column)
  
  all.GO.terms = newGO_terms
  row.names(all.GO.terms) = gsub(")","",lapply(lapply(stri_split_fixed(stri_reverse(all.GO.terms[,1]),"(",n = 2), FUN = stri_reverse), "[[",1))
  nTerms = nrow(all.GO.terms)
  pop.genes = unique(c(as.matrix(all.GO.terms[,-1])))
  N = length(pop.genes)
  
  # enrichment test using HyperGeometric test
  info.gene.overrep = data.frame(matrix(Inf, nrow = nrow(all.GO.terms), ncol = 8))
  
  K = length(givenSet)
  
  
  for (i in 1:nTerms) {
    
    GO.term.genes.symbol = as.character(all.GO.terms[i,which(all.GO.terms[i,] != '')])[-1]
    M = length(GO.term.genes.symbol)
    # print(M)
    x.overlap.genes.symbol = intersect(GO.term.genes.symbol, givenSet)
    # print(i)
    x = length(x.overlap.genes.symbol)
    # print(x)
    info.gene.overrep[i,1] = all.GO.terms[i,1]
    if(x > 0){
      info.gene.overrep[i,2] = paste0(GO.term.genes.symbol, collapse = "/")
      info.gene.overrep[i,3] = paste0(x.overlap.genes.symbol, collapse = "/")
      
      x.overlap.genes.df = bitr(x.overlap.genes.symbol, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
      x.overlap.genes.entrezID = as.character(x.overlap.genes.df$ENTREZID)
      
      info.gene.overrep[i,4] = paste0(x.overlap.genes.entrezID, collapse = "/")
      info.gene.overrep[i,5] = x
      info.gene.overrep[i,6] = paste0(x, "/", K)
      
      
      info.gene.overrep[i,7] = phyper(x, M, N-M, K, lower.tail = FALSE) #He wrote : overlap.genes-1
      #insert FDR val
      info.gene.overrep[i,8] = p.adjust(info.gene.overrep[i,7], method = "bonferroni", n = nTerms)
    }
  }
  colnames(info.gene.overrep) = c("Description", "Pathway.geneSymbol","Overlapping.geneSymbol", "Overlapping.geneID", "Count", "GeneRatio", "pvalue", "p.adjust")
  fwrite(info.gene.overrep[which(info.gene.overrep$pvalue != Inf),], file = outFilePath, sep = ",", col.names = T, quote = "auto")
}

# ----------------- GO_CC
pop.filepath = "Data/GO terms/GO_CC_2018.txt"
study = "NT"
givenSetFile = paste0("Data/Angela_HandPicked_DifferentialMarker/HandPicked_DifferentialMarkers_",study,".txt")
filename = tools::file_path_sans_ext(basename(pop.filepath))
outFilePath = paste0("Data/Angela_HandPicked_DifferentialMarker/GeneSet Enrichment results/",filename,"_ER_",study,".csv")
GO_geneSetEnrichmentAnalysis(givenSetFile, min.GeneSet=15, max.GeneSet=100,  pop.filepath, outFilePath)

study = "T"
givenSetFile = paste0("Data/Angela_HandPicked_DifferentialMarker/HandPicked_DifferentialMarkers_",study,".txt")
filename = tools::file_path_sans_ext(basename(pop.filepath))
outFilePath = paste0("Data/Angela_HandPicked_DifferentialMarker/GeneSet Enrichment results/",filename,"_ER_",study,".csv")
GO_geneSetEnrichmentAnalysis(givenSetFile, min.GeneSet=15, max.GeneSet=100,  pop.filepath, outFilePath)

# ----------------- GO_MF
pop.filepath = "Data/GO terms/GO_MF_2018.txt"
study = "NT"
givenSetFile = paste0("Data/Angela_HandPicked_DifferentialMarker/HandPicked_DifferentialMarkers_",study,".txt")
filename = tools::file_path_sans_ext(basename(pop.filepath))
outFilePath = paste0("Data/Angela_HandPicked_DifferentialMarker/GeneSet Enrichment results/",filename,"_ER_",study,".csv")
GO_geneSetEnrichmentAnalysis(givenSetFile, min.GeneSet=15, max.GeneSet=100,  pop.filepath, outFilePath)

study = "T"
givenSetFile = paste0("Data/Angela_HandPicked_DifferentialMarker/HandPicked_DifferentialMarkers_",study,".txt")
filename = tools::file_path_sans_ext(basename(pop.filepath))
outFilePath = paste0("Data/Angela_HandPicked_DifferentialMarker/GeneSet Enrichment results/",filename,"_ER_",study,".csv")
GO_geneSetEnrichmentAnalysis(givenSetFile, min.GeneSet=15, max.GeneSet=100,  pop.filepath, outFilePath)

# ----------------- GO_BP
pop.filepath = "Data/GO terms/GO_BP_2018.txt"
study = "NT"
givenSetFile = paste0("Data/Angela_HandPicked_DifferentialMarker/HandPicked_DifferentialMarkers_",study,".txt")
filename = tools::file_path_sans_ext(basename(pop.filepath))
outFilePath = paste0("Data/Angela_HandPicked_DifferentialMarker/GeneSet Enrichment results/",filename,"_ER_",study,".csv")
GO_geneSetEnrichmentAnalysis(givenSetFile, min.GeneSet=15, max.GeneSet=100,  pop.filepath, outFilePath)

study = "T"
givenSetFile = paste0("Data/Angela_HandPicked_DifferentialMarker/HandPicked_DifferentialMarkers_",study,".txt")
filename = tools::file_path_sans_ext(basename(pop.filepath))
outFilePath = paste0("Data/Angela_HandPicked_DifferentialMarker/GeneSet Enrichment results/",filename,"_ER_",study,".csv")
GO_geneSetEnrichmentAnalysis(givenSetFile, min.GeneSet=15, max.GeneSet=100,  pop.filepath, outFilePath)
