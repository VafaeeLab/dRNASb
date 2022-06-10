#' correlation analysis - returns correlation matrix as a dataframe
#' @param geneset_file_path full file path of geneset to perform correlation on
#' @param x_indices subset of columns from geneset to be used as x value to obtain correlation
#' @param y_indices subset of columns from geneset to be used as y value to obtain correlation
#' @param corr_method correlation method
#' @param corr_adj_method p value adjustment method
#' @param create_corrplot boolean indicating if corrplot should be created
#' @param write_results boolean indicating if results should be written to file
#' @param corrplot_file_name corrplot file name
#' @param output_file_name results file name
#' @param output_dir_path file path for all outputs - results and plots
#' @export
correlation_analysis <- function(geneset_file_path,
                                 x_indices,
                                 y_indices,
                                 corr_method = "pearson",
                                 corr_adj_method = "holm",
                                 create_corrplot = FALSE,
                                 write_results = FALSE,
                                 corrplot_file_name = "Corrplot.plot",
                                 output_file_name = "R.corr.csv",
                                 output_dir_path = "Results/Correlation_analysis/"
                                 ){

  print("performing correlation analysis ...")
  check_and_create_directory(output_dir_path)


  c <- read.csv(geneset_file_path,
              header = T, row.names = 1, check.names = F)
  All.cor <- psych::corr.test(c[x_indices], c[y_indices],
                            method = corr_method, adjust = corr_adj_method,
                            ci = FALSE)
  r.df <- as.data.frame(All.cor$r)
  R.corr <- r.df %>%
    dplyr::mutate(gene1 = row.names(r.df)) %>%
    tidyr::pivot_longer(-gene1,
                        names_to = "gene2", names_ptypes = list(gene2=character()),
                        values_to = "corr") %>%
    dplyr::mutate(gene1 = stringr::str_replace(gene1, "\\.\\.", " ("),
                  gene1 = stringr::str_replace(gene1, "\\.$", ")"),
                  gene2 = stringr::str_replace(gene2, "\\.\\.", " ("),
                  gene2 = stringr::str_replace(gene2, "\\.$", ")")) %>%
    dplyr::mutate(comb = paste(gene1, "-", gene2))

  if(create_corrplot){
    tiff(filename = paste0(output_dir_path, corrplot_file_name ,".tiff"), compression = "lzw")
    corrplot::corrplot(as.matrix(r.df), is.corr = FALSE, method = "circle", order = "hclust")
    dev.off()

    pdf(file = paste0(output_dir_path, corrplot_file_name ,".pdf"), width = 5, height = 5)
    corrplot::corrplot(as.matrix(r.df), is.corr = FALSE, method = "circle", order = "hclust")
    dev.off()
  }

  if(write_results){
    write.csv(R.corr, file = paste0(output_dir_path, output_file_name),
              quote = F, row.names = FALSE)
  }

  return (R.corr)

}



#' create connected genes plot - creates tiff and pdf plots
#' @param corr_data correlation data based on which to create network
#' @param plot_file_name plot file name
#' @param output_dir_path file path for all outputs - results and plots
#' @export
create_connected_genes_plot <- function(corr_data,
                                        plot_file_name,
                                        output_dir_path = "Results/Correlation_analysis/"){

  # output_dir_path = "Results/Correlation_analysis/"
  # plot_file_name = "Negative.highly.corrplot.plot"

  check_and_create_directory(output_dir_path)

  net <- igraph::graph.data.frame(unique(corr_data[,c(1,2)]),
                                directed = FALSE)

  cl <- igraph::cluster_louvain(net, weights = NULL)
  t <- as.data.frame(cl$membership)
  t1 <- as.data.frame(cl$names)
  t2 <- cbind(t1,t)
  colnames(t2)[1] <- "Gene.name"
  colnames(t2)[2] <- "membership"

  g_grouped = net

  for(i in unique(igraph::V(net)$community)){
    groupV = which(igraph::V(net)$community == i)
    g_grouped = igraph::add_edges(g_grouped, combn(groupV, 2), attr=list(weight = 2))
  }

  l <- igraph::layout_nicely(g_grouped)

  tiff(filename = paste0(output_dir_path, plot_file_name ,".tiff"),
       compression = "lzw")
  plot(cl, net, layout = igraph::layout_with_fr,
       vertex.size =10,
       edge.width = 1,
       vertex.label.dist=0.001,
       vertex.color ='gold',
       vertex.frame.color="#555555",
       edge.label=net$v,
       vertex.size=1,
       edge.color="gray",
       vertex.label.font=0.5,
       edge.label.font =0.5,
       edge.label.cex = 0.5,
       edge.arrow.size=0.2,
       edge.curved=0,
       vertex.label=igraph::V(net)$v,
       vertex.label.color="black",
       vertex.label.cex=0.5,
       vertex.label.cex = 0.5 )
  dev.off()

  pdf(file = paste0(output_dir_path, plot_file_name ,".pdf"),
      width = 5, height = 5)
  plot(cl,net, layout = igraph::layout_with_fr,
       vertex.size =10,
       edge.width = 1,
       vertex.label.dist=0.001,
       vertex.color ='gold',
       vertex.frame.color="#555555",
       edge.label=net$v,
       vertex.size=1,
       edge.color="gray",
       vertex.label.font=0.5,
       edge.label.font =0.5,
       edge.label.cex = 0.5,
       edge.arrow.size=0.2,
       edge.curved=0,
       vertex.label=igraph::V(net)$v,
       vertex.label.color="black",
       vertex.label.cex=0.5,
       vertex.label.cex = 0.5 )
  dev.off()
}


