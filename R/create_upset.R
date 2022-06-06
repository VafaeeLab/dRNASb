create_upset_plot <- function(DE_selected_upreg, DE_selected_downreg,
                              hour_mapping = c("2h", "4h", "8h", "16h", "24h"),
                              output_dir_path = "Results/Upset_plot/",
                              result_file_prefix = ""){

  print("creating upset plots ...")
  check_and_create_directory(output_dir_path)

  ################ Downregulated ################
  down_reg_list <- list()
  for(i in c(1:length(DE_selected_downreg))){
    down_reg_list[[hour_mapping[i]]] <- DE_selected_downreg[[i]]$Gene.name
  }

  # Setting colors
  main_bar_col <- c("blue4")
  sets_bar_col <- c("coral1")
  matrix_col <- c("forestgreen")
  shade_col <- c("wheat4")

  # Set Variables
  mb_ratio1 <- c(0.55, 0.45)

  output_file_name <- paste0(result_file_prefix, "downregulate.upset.plot.tiff")
  tiff(filename = paste0(output_dir_path, output_file_name), compression = "lzw")
  print({
    UpSetR::upset(UpSetR::fromList(down_reg_list),
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
  })
  dev.off()

  output_file_name <- paste0(result_file_prefix, "downregulate.upset.plot.pdf")
  pdf(file = paste0(output_dir_path, output_file_name), width = 5, height = 5)
  print({
    UpSetR::upset(UpSetR::fromList(down_reg_list),
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
  })

  dev.off()



  ################ Upregulated ################
  up_reg_list <- list()
  for(i in c(1:length(DE_selected_upreg))){
    up_reg_list[[hour_mapping[i]]] <- DE_selected_upreg[[i]]$Gene.name
  }


  # Setting colors
  main_bar_col <- c("violetred4")
  sets_bar_col <- c("turquoise4")
  matrix_col <- c("slateblue4")
  shade_col <- c("wheat4")


  # Set Variables
  mb_ratio1 <- c(0.55,0.45)

  output_file_name <- paste0(result_file_prefix, "upregulate.upset.plot.tiff")
  tiff(filename = paste0(output_dir_path, output_file_name), compression = "lzw")
  print({
    UpSetR::upset(UpSetR::fromList(up_reg_list),
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
                  shade.color = shade_col)
  })
  dev.off()


  output_file_name <- paste0(result_file_prefix, "upregulate.upset.plot.pdf")
  pdf(file = paste0(output_dir_path, output_file_name), width = 5, height = 5)
  print({
    UpSetR::upset(UpSetR::fromList(up_reg_list),
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
                  shade.color = shade_col)
  })
  dev.off()
}
