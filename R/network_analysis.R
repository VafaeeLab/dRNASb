perform_network_analysis <- function(ppi){
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
