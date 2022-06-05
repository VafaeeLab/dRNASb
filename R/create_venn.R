create_venn <- function(DE_selected_upreg, DE_selected_downreg,
                        hour_mapping = c("2h", "4h", "8h", "16h", "24h"),
                        output_dir_path = "Results/Venn_diagram/"){

  #########creating venn diagram with hourwise downregulated genes

  down_reg_list <- list()
  for(i in c(1:length(DE_selected_downreg))){
    down_reg_list[[hour_mapping[i]]] <- DE_selected_downreg[[i]]$Gene.name
  }

  grid::grid.newpage()
  venn.plot <- VennDiagram::venn.diagram(
    down_reg_list,
    fill = c("#84b3e7", "#317456", "#abcdef", "#ff99cc", "#bd0000"),
    cat.col = c("#84b3e7", "#317456", "#abcdef", "#ff99cc", "#bd0000"),
    cat.cex = 2,
    margin = 0.05,
    filename = NULL
  )

  check_and_create_directory(output_dir_path)

  output_file_name <- paste0(result_file_prefix, "downregulate.Venn.diagram.tiff")
  tiff(filename = paste0(output_dir_path, output_file_name), compression = "lzw")
  grid::grid.draw(venn.plot)
  dev.off()

  output_file_name <- paste0(result_file_prefix, "downregulate.Venn.diagram.plot.pdf")
  pdf(file = paste0(output_dir_path, output_file_name), width = 5, height = 5)
  grid::grid.draw(venn.plot)
  dev.off()


  #########creating venn diagram with hourwise upregulated genes

  up_reg_list <- list()
  for(i in c(1:length(DE_selected_upreg))){
    up_reg_list[[hour_mapping[i]]] <- DE_selected_upreg[[i]]$Gene.name
  }

  grid::grid.newpage()
  venn.plot <- VennDiagram::venn.diagram(
    up_reg_list,
    fill = c("#e1bebe", "darkgoldenrod1", "#c70000", "#ff99cc", "#9ecbff"),
    cat.col = c("#e1bebe", "darkgoldenrod1", "#c70000", "#ff99cc", "#9ecbff"),
    cat.cex = 2,
    margin = 0.05,
    filename = NULL
  )

  output_file_name <- paste0(result_file_prefix, "upregulate.Venn.diagram.tiff")
  tiff(filename = paste0(output_dir_path, output_file_name), compression = "lzw")
  grid::grid.draw(venn.plot)
  dev.off()

  output_file_name <- paste0(result_file_prefix, "upregulate.Venn.diagram.plot.pdf")
  pdf(file = paste0(output_dir_path, output_file_name), width = 5, height = 5)
  grid::grid.draw(venn.plot)
  dev.off()
}
