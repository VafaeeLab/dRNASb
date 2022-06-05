#' dRNASb pipeline
#' @param data_file_path full file path of input data file containing read counts in (gene x samples) format
#'   should have samples corresponding to different hours and replicates for each hour
#' @param phenotype_file_path full file path of phenotype file - contains mapping of sample to groups with groups column
#' @param annotation_function_file_path full file path of annotation function file
#' @param ppi_file_path full file path of ppi file
#' @param result_file_prefix a prefix string to be added to all results file
#' @param de_method differential expression method to be used
#' @param norm_method normalization method to be used - use show_allowed_norm_methods() for available norm methods
#' @param perform_filter Should filter preprocessing step be performed
#' @param hours_in_data Vector of different hours information present in the data
#' @param replicates_in_data Number of replicates for each hour in the data
#' @param logFC_cutoff log fold change cutoff used to select differentially expressed genes
#' @importFrom magrittr "%>%"
#' @importFrom Biobase AnnotatedDataFrame
#' @importFrom e1071 cmeans
#' @export
dRNASb_pipeline <- function(data_file_path,
                            phenotype_file_path,
                            annotation_function_file_path,
                            ppi_file_path,
                            result_file_prefix = "",
                            de_method = "limma",
                            norm_method = "log_TMM",
                            perform_filter = TRUE,
                            hours_in_data = c("0h", "2h", "4h", "8h", "16h", "24h"),
                            replicates_in_data = 3,
                            logFC_cutoff = 1) {


  # Read data---------------------------------------------------------------
  data <- read.csv(data_file_path, row.names = 1)
  pheno <- read.csv(phenotype_file_path, row.names = 1)
  ann_fun <- read.csv(annotation_function_file_path)
  ppi <- read.csv(ppi_file_path)
  # ---------------------------------------------------------------

  data.norm <- normalize(data, pheno, norm_method)

  if(perform_filter){
    data.norm <- filter(data.norm, pheno)
  }

  # Differential gene expression analysis ---------------------

  DE <- de_analysis_hourwise(data = data.norm, pheno = pheno, method = de_method)

  #obtain hourwise results - excluding that for 0th hour - assumed to be the first element
  hour_mapping <- hours_in_data[-1] #remove 0th hour
  de_results_dir_path <- "Results/Differential_gene_expression_analysis/"

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

  # Average replicates across each time -------------------------------------
  replicates <- replicates_in_data
  hour_mapping <- hours_in_data

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
  replicate_mean_hourly_with_genename <-
    cbind(data.frame("Gene.name" = row.names(replicate_mean_hourly)),
          replicate_mean_hourly)
  write.csv(replicate_mean_hourly_with_genename,
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
  DE_mean_expr_value <- merge(DE_selected_all, replicate_mean_hourly_with_genename, by = "Gene.name")
  DE_mean_expr_value <- DE_mean_expr_value %>%
    dplyr::select(-c(dplyr::starts_with("logFC"), "Mean.0h"))

  output_dir_path <- "Results/Average_data/"
  check_and_create_directory(output_dir_path)
  output_file_name <- paste0("All.", result_file_prefix, "DE.gene.csv")
  write.csv(DE_mean_expr_value, paste0(output_dir_path, output_file_name), row.names = FALSE)


  # Mfuzz Clustering   ------------------------------------------------
  perform_clustering(replicate_mean_hourly, ann_fun)



  # Venn diagram ------------------------------------------------------------------
  create_venn(DE_selected_upreg, DE_selected_downreg)


  # Upset Plot --------------------------------------------------------------
  create_upset_plot(DE_selected_upreg, DE_selected_downreg)


  # Network analysis using igraph --------------------------------
  perform_network_analysis(ppi)

}
