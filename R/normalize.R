normalize <- function(data, pheno){

  group <- pheno[colnames(data), "groups"]
  dge <- edgeR::DGEList(counts = data, group = group, genes = rownames(data))
  data.norm <- edgeR::cpm(edgeR::calcNormFactors(dge, method="TMM"), log = TRUE)

  return (data.norm)
}
