create_venn <- function(DE_selected_upreg, DE_selected_downreg){

  grid::grid.newpage()
  venn.plot <- VennDiagram::venn.diagram(
    list(
      "2h" = dpQ2$Gene.name,
      "4h" = dpQ4$Gene.name,
      "8h" = dpQ8$Gene.name,
      "16h" = dpQ16$Gene.name,
      "24h" = dpQ24$Gene.name
    ),
    fill = c("#84b3e7", "#317456", "#abcdef", "#ff99cc", "#bd0000"),
    cat.col = c("#84b3e7", "#317456", "#abcdef", "#ff99cc", "#bd0000"),
    cat.cex = 2,
    margin = 0.05,
    filename = NULL
  )
  # Writing to file
  tiff(filename ="./Results/Venn_diagram/Pathogn.downregulate.Venn.diagram_new.tiff", compression = "lzw")
  grid::grid.draw(venn.plot)
  dev.off()

  pdf(file ="./Results/Venn_diagram/Pathogen.downregulate.Venn.diagram.plot_new.pdf",width = 5, height = 5)
  grid::grid.draw(venn.plot)
  dev.off()


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

}
