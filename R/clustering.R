perform_clustering <- function(replicate_mean_hourly,
                               ann_fun,
                               timepoint = c(0, 2, 4, 8, 16, 24),
                               num_of_clust = 10,
                               cluer_pval_cutoff = 0.01,
                               output_dir_path = "Results/Mfuzz_Clustering/",
                               result_file_prefix = ""){
  # timepoint <- c(0, 2, 4, 8, 16, 24)

  # num_of_clust = 10

  print("performing clustering ...")

  y.dat <- as.matrix(replicate_mean_hourly)
  y.dat <- y.dat[which(apply(y.dat, 1, var) > 2 & apply(y.dat, 1, mean) > 2), 1:6]
  y.dat <- rbind(timepoint, y.dat)
  rownames(y.dat)[1]<- "time"
  tmp <- tempfile()
  write.table(y.dat, file = tmp, sep = '\t', quote = FALSE, col.names = NA)


  #problem : below line requires : library(Biobase)
  z.data <- Mfuzz::table2eset(tmp)


  data.z <- Mfuzz::standardise(z.data)
  # class(data.z)
  m1 <- Mfuzz::mestimate(data.z)

  #below line requires library(e1071)
  # Mfuzz::Dmin(data.z, m=m1, crange=seq(2,22,1), repeats = 3, visu = TRUE)

  c <- Mfuzz::mfuzz(data.z, c = num_of_clust, m = m1)

  # output_dir_path <- "Results/Mfuzz_Clustering/"
  check_and_create_directory(output_dir_path)
  output_file_name <- paste0(result_file_prefix, "mfuzz.plot.tiff")

  tiff(filename = paste0(output_dir_path, output_file_name), compression = "lzw")
  Mfuzz::mfuzz.plot(
    data.z,
    cl = c,
    mfrow = c(4, 4),
    min.mem = 0.5,
    time.labels = timepoint,
    new.window = FALSE
  )
  dev.off()

  output_file_name <- paste0(result_file_prefix, "mfuzz.plot.pdf")
  pdf(file = paste0(output_dir_path, output_file_name),
      width = 10,
      height = 10)
  Mfuzz::mfuzz.plot(
    data.z,
    cl = c,
    mfrow = c(4, 4),
    min.mem = 0.5,
    time.labels = timepoint,
    new.window = FALSE
  )
  dev.off()

  membership <- c$membership
  membership <- data.frame(membership)
  fd <- data.frame(cor(t(c[[1]])))
  acore <- Mfuzz::acore(data.z, c, min.acore = 0.5)
  acore_list <-
    do.call(rbind, lapply(seq_along(acore), function(i) {
      data.frame(CLUSTER = i, acore[[i]])
    }))
  colnames(acore_list)[2] <- "Gene.name"
  genelist <- Mfuzz::acore(data.z, cl = c, min.acore = 0.7)
  temp <- do.call("rbind", lapply(
    genelist,
    FUN = function(x) {
      return(paste0(as.character(x$NAME), collapse = ","))
    }
  ))
  Cluster_list <- as.data.frame(temp)
  colnames(Cluster_list) <- "Gene.name"
  Cluster_list <-
    stringr::str_split_fixed(Cluster_list$Gene.name, ",", n = Inf)
  Cluster_list <- t(Cluster_list)

  cluster_list_col_names <- paste0("Cluster", c(1:num_of_clust))
  colnames(Cluster_list) <- cluster_list_col_names

  output_file_name <- paste0(result_file_prefix, "cluster.acore_list.csv")
  write.csv(
    acore_list,
    file = paste0(output_dir_path, output_file_name),
    quote = F,
    row.names = FALSE
  )



  # Make list ---------------------------------------------------------------
  anno <- ann_fun
  GO <- unique(anno$Gene.ontology.ID)
  Uniprot.ID <- sapply(1:length(GO), function(i) paste(gsub("[[:space:]]", "", anno[which(anno$Gene.ontology.ID==GO[i]),]$Uniprot.ID),collapse=" "))
  Uniprot.ID <- as.data.frame(Uniprot.ID)
  GO_Pro_ID <- data.frame(GO.ID=unique(anno$Gene.ontology.ID),
                          Uniprot.ID=Uniprot.ID)


  # Make list of list -------------------------------------------------------
  Anno <- list()
  GO_IDs <- as.vector(GO_Pro_ID[,1])

  for (i in GO_IDs) {
    myindex <- which(GO_Pro_ID == i)
    Anno[i] <- strsplit(as.character(GO_Pro_ID[myindex, 2]), " ")
  }


  # Enrichment using "ClueR" ------------------------------------------------
  # cluer_pval_cutoff <- 0.01

  ce <- ClueR::clustEnrichment(c, annotation = Anno,
                               effectiveSize = c(2,100), pvalueCutoff = cluer_pval_cutoff)

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

  output_file_name <- paste0(result_file_prefix, "cluster.enrichment.csv")
  write.csv(
    out,
    file = paste0(output_dir_path, output_file_name),
    quote = F,
    row.names = FALSE
  )


  # Function using GO.term --------------------------------------------------
  fu <- ann_fun[, c(1,2)]
  u <- merge(out, fu, by="Gene.ontology.ID")

  output_file_name <- paste0(result_file_prefix, "cluster.enrichment.function.csv")
  write.csv(
    u,
    file = paste0(output_dir_path, output_file_name),
    quote = F,
    row.names = FALSE
  )

  #  Variable and  frequency of Function
  cl <- as.data.frame(out[, c(5)])
  colnames(cl) <- "Gene.name"
  gene.list <-
    cl %>% tidyr::separate_rows(Gene.name, sep = "\\|") %>%
    dplyr::group_by(Gene.name)
  freq.list.gene <-
    cl %>% tidyr::separate_rows(Gene.name, sep = "\\|") %>%
    dplyr::group_by(Gene.name) %>%
    dplyr::summarize(n = dplyr::n()) %>%
    dplyr::arrange(desc(n))
  output_file_name <- paste0(result_file_prefix, "Variable.enriched.genes.and.their.frequency.csv")
  write.csv(
    freq.list.gene,
    file = paste0(output_dir_path, output_file_name),
    quote = F,
    row.names = FALSE
  )
}
