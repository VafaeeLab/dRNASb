create_upset_plot <- function(DE_selected_upreg, DE_selected_downreg){

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
}
