filter <- function(data, pheno){

  group <- pheno[colnames(data), "groups"]
  keep <- edgeR::filterByExpr(data, group = group, min.count = log2(10))
  data.fil <- data[keep, ]
  data.fil <- data.frame(data.fil)

  return (data.fil)
}
