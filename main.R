library(dRNASb)

#pathogen
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
                logFC_cutoff = 1,
                num_of_clust = 10,
                hub_gene_cutoff = 10,
                betweenness_cutoff = 100)


#host (execution of below line takes approx 4 hours)
dRNASb_pipeline(data_file_path = "Data/Host.data.csv",
                phenotype_file_path = "Data/Pheno.csv",
                annotation_function_file_path = "Data/Host.annotation.function.csv",
                ppi_file_path = "Data/Host.ppi.csv",
                result_file_prefix = "Host.",
                de_method = "limma",
                norm_method = "log_TMM",
                perform_filter = TRUE,
                hours_in_data = c("0h", "2h", "4h", "8h", "16h", "24h"),
                replicates_in_data = 3,
                logFC_cutoff = 1,
                num_of_clust = 10,
                hub_gene_cutoff = 100,
                betweenness_cutoff = 1000)





#correlation analysis with all gene in study

R.corr <- correlation_analysis(
  geneset_file_path = "Data/Select.gene.set.for.correlation.study/All.mean.transpose.for.cor.csv",
  x_indices = c(1:293),
  y_indices = c(294:422),
  corr_method = "pearson",
  corr_adj_method = "holm",
  create_corrplot = FALSE,
  write_results = FALSE,
  output_dir_path = "Results/Correlation_analysis/"
)

negative <- R.corr %>% dplyr::filter(corr<(-0.7))
create_connected_genes_plot(corr_data = negative,
                            plot_file_name = "Negative.highly.corrplot.plot",
                            output_dir_path = "Results/Correlation_analysis/")

positive <- R.corr %>% dplyr::filter(corr>0.7, corr<1)
create_connected_genes_plot(corr_data = positive,
                            plot_file_name = "Positive.highly.corrplot.plot",
                            output_dir_path = "Results/Correlation_analysis/")


#correlation analysis with selected gene in study
#between host and pathogen

R.corr <- correlation_analysis(
  geneset_file_path = "Data/Select.gene.set.for.correlation.study/Transpose.both gene set.csv",
  x_indices = c(1:54),
  y_indices = c(1:54),
  corr_method = "pearson",
  corr_adj_method = "holm",
  create_corrplot = TRUE,
  write_results = TRUE,
  corrplot_file_name = "Corrplot.plot",
  output_file_name = "R.corr.host.pathogen.csv",
  output_dir_path = "Results/Correlation_analysis/"
)

negative <- R.corr %>% dplyr::filter(corr<(-0.7))
create_connected_genes_plot(corr_data = negative,
                            plot_file_name = "Negative.70.highly.corrplot.plot",
                            output_dir_path = "Results/Correlation_analysis/")

negative <- R.corr %>% dplyr::filter(corr<(-0.9))
create_connected_genes_plot(corr_data = negative,
                            plot_file_name = "Negative.90.highly.corrplot.plot",
                            output_dir_path = "Results/Correlation_analysis/")

positive <- R.corr %>% dplyr::filter(corr>0.7, corr<1)
create_connected_genes_plot(corr_data = positive,
                            plot_file_name = "Positive.70.highly.corrplot.plot",
                            output_dir_path = "Results/Correlation_analysis/")

positive <- R.corr%>% dplyr::filter(corr>0.9, corr<1)
create_connected_genes_plot(corr_data = positive,
                            plot_file_name = "Positive.90.highly.corrplot.plot",
                            output_dir_path = "Results/Correlation_analysis/")
