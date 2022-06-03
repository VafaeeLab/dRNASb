#' dRNASb pipeline
#' @param data_file_path full file path of input data file containing read counts in (gene x samples) format
#' @param phenotype_file_path full file path of phenotype file - contains mapping of sample to groups with groups column
#' @param annotation_function_file_path full file path of annotation function file
#' @param ppi_file_path full file path of ppi file
#' @param result_file_prefix a prefix string to be added to all results file
#' @export
dRNASb <- function(data_file_path, phenotype_file_path,
                   annotation_function_file_path, ppi_file_path,
                   result_file_prefix = ""){

  # Read data---------------------------------------------------------------
  data <- read.csv(data_file_path, row.names = 1)
  pheno <- read.csv(phenotype_file_path, row.names = 1)
  ann_fun <- read.csv(annotation_function_file_path)
  ppi <- read.csv(ppi_file_path)
  # Read data end---------------------------------------------------------------


  data.norm <- normalize(data, pheno)
  data.norm <- filter(data.norm, pheno)


  # Differential gene expression analysis ---------------------

  DE <- de_analysis_hourwise(data = data.norm, pheno = pheno, method = "limma")

  #obtain hourwise results
  hour_mapping <- c("2h", "4h", "8h", "16h", "24h")
  de_results_dir_path <- "Results/Differential_gene_expression_analysis/"
  logFC_cutoff <- 1

  DE_selected <- list()
  DE_selected_upreg <- list()
  DE_selected_downreg <- list()

  check_and_create_directory(de_results_dir_path)
  for(i in c(1:length(DE))){
    DE_per_hour <- DE[[i]]
    DE_per_hour <- cbind("Gene.name" = row.names(DE_per_hour), DE_per_hour)

    logFC_col_name <- paste("logFC", hour_mapping[i], sep = ".")
    colnames(DE_per_hour)[2] <- logFC_col_name

    de_file_name <- paste0(result_file_prefix, "DE-", hour_mapping[i], ".csv")
    write.csv(DE_per_hour, file = paste0(de_results_dir_path, "/", de_file_name), row.names = FALSE)

    DE_selected[[i]] <- DE_per_hour[abs(DE_per_hour[[logFC_col_name]]) > logFC_cutoff, c(1, 2)]
    DE_selected_upreg[[i]] <- DE_per_hour[DE_per_hour[[logFC_col_name]] > logFC_cutoff, c(1, 2)]
    DE_selected_downreg[[i]] <- DE_per_hour[DE_per_hour[[logFC_col_name]] < -logFC_cutoff, c(1, 2)]
  }


# ################### 2h
# D<-DE[[1]]
# Gene.name<-as.data.frame(row.names(D))
# colnames(Gene.name)<-"Gene.name"
# D<-cbind(Gene.name,D)
# colnames(D)[2]<-"logFC.2h"
# Qp2<-subset(D,D$logFC.2h<(-1)|D$logFC.2h>1)
# Qp2<-Qp2[,c(1,2)]
# dpQ2<-D%>% dplyr::filter(logFC.2h<(-1))
# dpQ2<-dpQ2[,c(1,2)]
# upQ2<-D%>% dplyr::filter(logFC.2h>(1))
# upQ2<-upQ2[,c(1,2)]
# write.csv(D,file = paste0("./Results/","./Differential_gene_expression_analysis/","Pathogen.DE-2h.csv"), row.names = FALSE)
#
# ################### 4h
# D<-DE[[2]]
# Gene.name<-as.data.frame(row.names(D))
# colnames(Gene.name)<-"Gene.name"
# D<-cbind(Gene.name,D)
# colnames(D)[2]<-"logFC.4h"
# Qp4<-subset(D,D$logFC.4h<(-1)|D$logFC.4h>1)
# Qp4<-Qp4[,c(1,2)]
# dpQ4<-D%>% dplyr::filter(logFC.4h<(-1))
# dpQ4<-dpQ4[,c(1,2)]
# upQ4<-D%>% dplyr::filter(logFC.4h>(1))
# upQ4<-upQ4[,c(1,2)]
# write.csv(D,file = paste0("./Results/","./Differential_gene_expression_analysis/","Pathogen.DE-4h.csv"), row.names = FALSE)
#
# ################### 8h
# D<-DE[[3]]
# Gene.name<-as.data.frame(row.names(D))
# colnames(Gene.name)<-"Gene.name"
# D<-cbind(Gene.name,D)
# colnames(D)[2]<-"logFC.8h"
# Qp8<-subset(D,D$logFC.8h<(-1)|D$logFC.8h>1)
# Qp8<-Qp8[,c(1,2)]
# dpQ8<-D%>% dplyr::filter(logFC.8h<(-1))
# dpQ8<-dpQ8[,c(1,2)]
# upQ8<-D%>% dplyr::filter(logFC.8h>(1))
# upQ8<-upQ8[,c(1,2)]
# write.csv(D,file = paste0("./Results/","./Differential_gene_expression_analysis/","Pathogen.DE-8h.csv"), row.names = FALSE)
#
# ################### 16h
# D<-DE[[4]]
# Gene.name<-as.data.frame(row.names(D))
# colnames(Gene.name)<-"Gene.name"
# D<-cbind(Gene.name,D)
# colnames(D)[2]<-"logFC.16h"
# Qp16<-subset(D,D$logFC.16h<(-1)|D$logFC.16h>1)
# Qp16<-Qp16[,c(1,2)]
# dpQ16<-D%>% dplyr::filter(logFC.16h<(-1))
# dpQ16<-dpQ16[,c(1,2)]
# upQ16<-D%>% dplyr::filter(logFC.16h>(1))
# upQ16<-upQ16[,c(1,2)]
# write.csv(D,file = paste0("./Results/","./Differential_gene_expression_analysis/","Pathogen.DE-16h.csv"), row.names = FALSE)
#
# ################### 24h
# D<-DE[[5]]
# Gene.name<-as.data.frame(row.names(D))
# colnames(Gene.name)<-"Gene.name"
# D<-cbind(Gene.name,D)
# colnames(D)[2]<-"logFC.24h"
# Qp24<-subset(D,D$logFC.24h<(-1)|D$logFC.24h>1)
# Qp24<-Qp24[,c(1,2)]
# dpQ24<-D%>% dplyr::filter(logFC.24h<(-1))
# dpQ24<-dpQ24[,c(1,2)]
# upQ24<-D%>% dplyr::filter(logFC.24h>(1))
# upQ24<-upQ24[,c(1,2)]
# write.csv(D,file = paste0("./Results/","./Differential_gene_expression_analysis/","Pathogen.DE-24h.csv"), row.names = FALSE)


  # Average replicates across each time -------------------------------------
  replicates <- 3
  hour_mapping <- c("0h", "2h", "4h", "8h", "16h", "24h")

  replicate_mean_hourly <- data.frame()
  for(i in c(1: length(hour_mapping))){
    col_start <- (i-1)*replicates + 1
    col_end <- i*replicates
    if(i == 1){
      replicate_mean_hourly <- rowMeans(data[, c(col_start:col_end)], na.rm = TRUE)
    } else{
      replicate_mean_hourly <- cbind(replicate_mean_hourly,
                                        rowMeans(data[, c(col_start:col_end)], na.rm = TRUE))
    }
  }
  colnames(replicate_mean_hourly) <- paste("Mean", hour_mapping, sep = ".")

  output_dir_path <- "Results/Average_data/"
  output_file_name <- paste0("All.", result_file_prefix, "mean.data.csv")
  check_and_create_directory(output_dir_path)
  replicate_mean_hourly <- cbind(Gene.name = rownames(replicate_mean_hourly),
                                 replicate_mean_hourly)
  write.csv(replicate_mean_hourly,
            paste0(output_dir_path, output_file_name), row.names = FALSE)

  # -------------------------------------


  ### Select gene that shows DE atleast in one time point
  for(i in c(1: length(DE))){
    DE_selected[[i]][DE_selected[[i]][, 2] < 0, 2] <- 1
    if(i == 1){
      DE_selected_all <- DE_selected[[i]]
    } else{
      DE_selected_all <- plyr::rbind.fill(DE_selected_all, DE_selected[[i]])
    }
  }
  DE_selected_all[is.na(DE_selected_all)]<-0
  DE_selected_all <-
    DE_selected_all %>% tidyr::pivot_longer(dplyr::starts_with("logFC"),
                                            names_to = "timepoint",
                                            values_to = "sign") %>%
    dplyr::filter(sign == 1) %>%
    tidyr::pivot_wider(names_from = timepoint,
                       values_from = sign,
                       values_fill = 0)



  ### Obtain mean value of DE genes
  DE_mean_expr_value <- merge(DE_selected_all, replicate_mean_hourly, by = "Gene.name")
  DE_mean_expr_value <- DE_mean_expr_value %>%
    dplyr::select(-c(dplyr::starts_with("logFC"), "Mean.0h"))

  output_dir_path <- "Results/Average_data/"
  check_and_create_directory(output_dir_path)
  output_file_name <- paste0("All.", result_file_prefix, "DE.gene.csv")
  write.csv(DE_mean_expr_value, paste0(output_dir_path, output_file_name), row.names = FALSE)



  # Mfuzz Clustering   ------------------------------------------------
  y.dat<- as.matrix(replicate_mean_hourly)
  y.dat <- y.dat[which(apply(y.dat, 1, var)>2 & apply(y.dat,1,mean)>2), 1:6]
  timepoint <- c(0,2,4,8,16,24)
  y.dat <- rbind(timepoint, y.dat)
  rownames(y.dat)[1]<- "time"
  tmp<- tempfile()
  write.table(y.dat,file=tmp, sep='\t',quote=FALSE, col.names=NA)


  #problem : below line requires : library(Biobase)
  z.data <- Mfuzz::table2eset(tmp)


  data.z <-Mfuzz::standardise(z.data)
  class(data.z)
  m1 <-Mfuzz::mestimate(data.z)

  #below line requires library(e1071)
  # Mfuzz::Dmin(data.z, m=m1, crange=seq(2,22,1), repeats = 3, visu = TRUE)

  clust=10
  c<- Mfuzz::mfuzz(data.z, c=clust, m=m1)

  tiff(filename ="./Results/Mfuzz_Clustering/Pathogen.mfuzz.plot.tiff", compression = "lzw")
  Mfuzz::mfuzz.plot(data.z,cl=c,mfrow=c(4,4),min.mem=0.5,time.labels=c(0,2,4,8,16,24),new.window=FALSE)
  dev.off()

  pdf(file = "./Results/Mfuzz_Clustering/Pathogen.mfuzz.plot.pdf",width = 10, height = 10)
  Mfuzz::mfuzz.plot(data.z,cl=c,mfrow=c(4,4),min.mem=0.5,time.labels=c(0,2,4,8,16,24),new.window=FALSE)
  dev.off()

  membership<-c$membership
  membership<-data.frame(membership)
  fd<-data.frame(cor(t(c[[1]])))
  acore<-Mfuzz::acore(data.z,c,min.acore = 0.5)
  acore_list<-do.call(rbind,lapply(seq_along(acore), function(i){data.frame(CLUSTER=i, acore[[i]])}))
  colnames(acore_list)[2]<-"Gene.name"
  genelist<- Mfuzz::acore(data.z,cl=c,min.acore=0.7)
  temp <- do.call("rbind", lapply(genelist, FUN = function(x){
    return(paste0(as.character(x$NAME), collapse = ","))
  }))
  Cluster_list<-as.data.frame(temp)
  colnames(Cluster_list) <-"Gene.name"
  Cluster_list<-stringr::str_split_fixed(Cluster_list$Gene.name,",", n=Inf)
  Cluster_list<-t(Cluster_list)
  colnames(Cluster_list)<- c("Cluster1", "Cluster2","Cluster3","Cluster4","Cluster5","Cluster6","Cluster7","Cluster8","Cluster9","Cluster10")  ### ? how we can make it depanded to cluster number???

  write.csv(acore_list,file = paste0("./Results/","./Mfuzz_Clustering/","Pathogen.cluster.acore_list.csv"), quote = F, row.names = FALSE)



  # Make list ---------------------------------------------------------------
  anno<-ann_fun
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
  ce <- ClueR::clustEnrichment(c, annotation=Anno, effectiveSize=c(2,100), pvalueCutoff=0.01)

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
  fu<-ann_fun[,c(1,2)]
  u<-merge(out, fu, by="Gene.ontology.ID")
  write.csv(u, file= paste0 ("./Results/","./Mfuzz_Clustering/","Pathogen.cluster.enrichment.function.csv"),quote = F, row.names = FALSE)



  #  Variable and  frequency of Function
  cl<-as.data.frame(out[,c(5)])
  colnames(cl)<-"Gene.name"
  gene.list<-cl %>% tidyr::separate_rows(Gene.name, sep = "\\|") %>% group_by(Gene.name)
  freq.list.pathogen.gene<-cl %>% tidyr::separate_rows(Gene.name, sep = "\\|") %>% group_by(Gene.name)%>% dplyr::summarize(n = n()) %>% arrange(desc(n))
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



  grid::grid.newpage()
  venn.plot <- VennDiagram::draw.quintuple.venn(
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
  grid::grid.draw(venn.plot)
  dev.off()

  pdf(file ="./Results/Venn_diagram/Pathogen.downregulate.Venn.diagram.plot.pdf",width = 5, height = 5)
  grid::grid.draw(venn.plot)
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


  grid::grid.newpage()
  venn.plot <- VennDiagram::draw.quintuple.venn(
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
  grid::grid.draw(venn.plot)
  dev.off()

  pdf(file ="./Results//Venn_diagram/Pathogen.upregulate.Venn.diagram.plot.pdf",width = 5, height = 5)
  grid::grid.draw(venn.plot)
  dev.off()


  # Upset Plot --------------------------------------------------------------
  ### Select downregulated gene
  dpQ2$logFC.2h [dpQ2$logFC.2h<0]<-1
  dpQ4$logFC.4h [dpQ4$logFC.4h<0]<-1
  dpQ8$logFC.8h [dpQ8$logFC.8h<0]<-1
  dpQ16$logFC.16h [dpQ16$logFC.16h<0]<-1
  dpQ24$logFC.24h [dpQ24$logFC.24h<0]<-1
  dpQ<-plyr::rbind.fill(dpQ2,dpQ4,dpQ8,dpQ16,dpQ24)
  dpQ[is.na(dpQ[1:6])]<-0
  dpQ <- dpQ %>% tidyr::pivot_longer(logFC.2h:logFC.24h, names_to = "timepoint", values_to = "sign") %>%
    dplyr::filter(sign == 1) %>%
    tidyr::pivot_wider(names_from = timepoint, values_from = sign, values_fill = 0)%>% t()

  write.table(dpQ,file = paste0("./Inputs/","Pathogen.all.downregulated.ED.csv"),  sep=",",col.names= FALSE,row.names = TRUE)


  ### Select upregulated gene
  upQ2$logFC.2h [upQ2$logFC.2h>0]<-1
  upQ4$logFC.4h [upQ4$logFC.4h>0]<-1
  upQ8$logFC.8h [upQ8$logFC.8h>0]<-1
  upQ16$logFC.16h [upQ16$logFC.16h>0]<-1
  upQ24$logFC.24h [upQ24$logFC.24h>0]<-1
  upQ<-plyr::rbind.fill(upQ2,upQ4,upQ8,upQ16,upQ24)
  upQ[is.na(upQ[1:6])]<-0
  upQ <- upQ %>% tidyr::pivot_longer(logFC.2h:logFC.24h, names_to = "timepoint", values_to = "sign") %>%
    dplyr::filter(sign == 1) %>%
    tidyr::pivot_wider(names_from = timepoint, values_from = sign, values_fill = 0)%>% t()

  write.table(upQ,file = paste0("./Inputs/","Pathogen.all.upregulated.ED.csv"),  sep=",",col.names= FALSE,row.names = TRUE)



  ### Downregulated
  y <-t(read.csv("Inputs/Pathogen.all.downregulated.ED.csv")) %>% as.data.frame()
  colnames(y)<-y[1,]
  y <- y[-1,] %>% dplyr::mutate_if(is.character, as.numeric)
  UpSetR::upset(y)

  # Setting colors
  main_bar_col <- c("blue4")
  sets_bar_col <- c("coral1")
  matrix_col <- c("forestgreen")
  shade_col <- c("wheat4")


  # Setting Set Variables
  mb_ratio1 <- c(0.55,0.45)

  tiff(filename ="./Results//Upset_plot/Pathogen.downregulate.upset.plot.tiff", compression = "lzw")
  UpSetR::upset(y,
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
  UpSetR::upset(y,
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
  y <- y[-1,] %>% dplyr::mutate_if(is.character, as.numeric)
  UpSetR::upset(y)

  # Setting colors
  main_bar_col <- c("violetred4")
  sets_bar_col <- c("turquoise4")
  matrix_col <- c("slateblue4")
  shade_col <- c("wheat4")


  # Setting Set Variables
  mb_ratio1 <- c(0.55,0.45)

  tiff(filename ="./Results//Upset_plot/Pathogen.upregulate.upset.plot.tiff", compression = "lzw")
  UpSetR::upset(y,
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
  UpSetR::upset(y,
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
  net<-igraph::graph.data.frame(unique(ppi[,c(1,4)]),directed = FALSE)
  igraph::V(net)$label<-igraph::V(net)$name
  deg<-igraph::degree(net,v=igraph::V(net), mode = c("total"),
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
  b<-igraph::betweenness(net, v = igraph::V(net), directed = TRUE, weights = NULL,
                         nobigint = TRUE, normalized = FALSE)

  b.p<-as.data.frame(b)
  gene.name<-as.data.frame(row.names(b.p))
  b.p<-cbind(gene.name,b.p)
  colnames(b.p)[1]<-"Gene.name"
  colnames(b.p)[2]<-"Betweennes"
  Betweennes.pathogen<-b.p%>% dplyr::filter(Betweennes>100)
  write.csv(Betweennes.pathogen,file = paste0 ("./Results/","./Network_analysis/","Pathogen.betweennes.csv"),quote = F, row.names = FALSE)

  # Network short distance --------------------------------------------------
  d<-igraph::distances(
    net,
    v = igraph::V(net),
    to = igraph::V(net),
    mode = c("all", "out", "in"),
    weights = NULL,
    algorithm = c("automatic", "unweighted", "dijkstra", "bellman-ford", "johnson")
  )
  g<-as.data.frame(row.names(d))
  colnames(g)[1]<-"Gene.name"
  d<-cbind(g,d)
  write.csv(d,file = paste0 ("./Results/","./Network_analysis/","Pathogen.shortest.path.csv"),quote = F, row.names = FALSE)


  # Network Modules ---------------------------------------------------------
  cl<-igraph::cluster_louvain(net, weights = NULL)
  t<-as.data.frame(cl$membership)
  t1<-as.data.frame(cl$names)
  t2<-cbind(t1,t)
  colnames(t2)[1]<-"Gene.name"
  colnames(t2)[2]<-"membership"
  write.csv(t2,file = paste0 ("./Results/","./Network_analysis/","Pathogen.modules.in.ppi.network.csv"),quote = F, row.names = FALSE)


  # Plot Modules ------------------------------------------------------------
  g_grouped = net

  for(i in unique(igraph::V(net)$community)){
    groupV = which(igraph::V(net)$community == i)
    g_grouped = igraph::add_edges(g_grouped, combn(groupV, 2), attr=list(weight = 2))
  }

  l <- igraph::layout_nicely(g_grouped)

  # Writing to file
  pdf(file ="./Results/Network_analysis/Pathogen.modules.in.ppi.network.plot.pdf",width = 5, height = 5)
  plot(cl,net, layout = igraph::layout_with_fr,
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
       vertex.label=igraph::V(net)$v,
       vertex.label.color="black",
       vertex.label.cex=0.5,
       vertex.label.cex = 0.5 )
  dev.off()

  tiff(filename ="./Results/Network_analysis/Pathogen.modules.in.ppi.network.plot.tiff", compression = "lzw")
  plot(cl,net, layout = igraph::layout_with_fr,
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
       vertex.label=igraph::V(net)$v,
       vertex.label.color="black",
       vertex.label.cex=0.5,
       vertex.label.cex = 0.5 )

  dev.off()

}
