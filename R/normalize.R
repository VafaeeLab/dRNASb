norm_methods <- list(
  list(
    norm_name = "log_TMM",
    norm_func = function(data, pheno) {
      group <- pheno[colnames(data), "groups"]
      dge <-
        edgeR::DGEList(counts = data,
                       group = group,
                       genes = rownames(data))
      data.norm <-
        edgeR::cpm(edgeR::calcNormFactors(dge, method = "TMM"), log = TRUE)
      return (data.norm)
    }
  ),
  list(
    norm_name = "TMM",
    norm_func = function(data, pheno) {
      group <- pheno[colnames(data), "groups"]
      dge <-
        edgeR::DGEList(counts = data,
                       group = group,
                       genes = rownames(data))
      data.norm <-
        edgeR::cpm(edgeR::calcNormFactors(dge, method = "TMM"), log = FALSE)
      return (data.norm)
    }
  ),
  list(
    norm_name = "CPM",
    norm_func = function(data, pheno) {
      data.norm <- edgeR::cpm(data, log = FALSE)
      return (data.norm)
    }
  ),
  list(
    norm_name = "log_CPM",
    norm_func = function(data, pheno) {
      data.norm <- edgeR::cpm(data, log = TRUE)
      return (data.norm)
    }
  ),
  list(
    norm_name = "znorm_log_CPM",
    norm_func = function(data, pheno) {
      data.norm <- edgeR::cpm(data, log = TRUE)
      data.norm <- t(as.matrix(data.norm))
      data.norm <- scale(data.norm)
      data.norm <- t(as.matrix(data.norm))
      return (data.norm)
    }
  ),
  list(
    norm_name = "none",
    norm_func = function(data, pheno) {
      return (as.matrix(data))
    }
  )
)


#' show_allowed_norm_methods
#' @export
show_allowed_norm_methods <- function(){
  for (nm in norm_methods) {
    print(nm[["norm_name"]])
  }
}
# show_allowed_norm_methods()

normalize <- function(data, pheno, norm_method = "log_TMM"){
  for(nm in norm_methods){
    if(nm[["norm_name"]] == norm_method){
      return (nm[["norm_func"]](data, pheno))
    }
  }
  print("unknown norm_method - no norm done !")
  return (data)
}

