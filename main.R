library(dRNASb)


dRNASb_pipeline(data_file_path = "Data/Pathogen.data.csv",
                phenotype_file_path = "Data/Pheno.csv",
                annotation_function_file_path = "Data/Pathogen.annotation.function.csv",
                ppi_file_path = "Data/Pathogen.ppi.csv",
                result_file_prefix = "Pathogen.",
                de_method = "limma",
                norm_method = "log_TMM",
                perform_filter = TRUE,
                hours_in_data = c("0h", "2h", "4h", "8h", "16h", "24h"),
                replicates_in_data = 3,
                logFC_cutoff = 1)
