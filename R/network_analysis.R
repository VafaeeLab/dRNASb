perform_network_analysis <- function(ppi, output_dir_path = "Results/Network_analysis/",
                                     hub_gene_cutoff = 10,
                                     betweenness_cutoff = 100,
                                     result_file_prefix = ""){

  print("performing network analysis ...")
  check_and_create_directory(output_dir_path)

  net <- igraph::graph.data.frame(unique(ppi[,c(1,4)]), directed = FALSE)
  igraph::V(net)$label <- igraph::V(net)$name
  deg <- igraph::degree(net, v = igraph::V(net), mode = c("total"),
                      loops = TRUE, normalized = FALSE)

  # Hub gene  ---------------------------------------------------
  hub <- as.data.frame(deg[which(deg >= hub_gene_cutoff)])
  gene.name <- as.data.frame(row.names(hub))
  hub <- cbind(gene.name, hub)
  colnames(hub)[1] <- "Gene.name"
  colnames(hub)[2] <- "edge.number"
  output_file_name <- paste0(result_file_prefix, "hub.more.", hub_gene_cutoff, ".csv")
  write.csv(hub, file = paste0(output_dir_path, output_file_name), quote = F, row.names = FALSE)


  # Network betweennes --------------------------------------------------
  b <- igraph::betweenness(net, v = igraph::V(net), directed = TRUE, weights = NULL,
                         normalized = FALSE)
  b <- as.data.frame(b)
  gene.name <- as.data.frame(row.names(b))
  b <- cbind(gene.name,b)
  colnames(b)[1] <- "Gene.name"
  colnames(b)[2] <- "Betweennes"
  betweenness <- b %>% dplyr::filter(Betweennes > betweenness_cutoff)
  output_file_name <- paste0(result_file_prefix, "betweennes.csv")
  write.csv(betweenness, file = paste0(output_dir_path, output_file_name),
            quote = F, row.names = FALSE)

  # Network short distance --------------------------------------------------
  d <- igraph::distances(
    net,
    v = igraph::V(net),
    to = igraph::V(net),
    mode = c("all", "out", "in"),
    weights = NULL,
    algorithm = c("automatic", "unweighted", "dijkstra", "bellman-ford", "johnson")
  )
  g <- as.data.frame(row.names(d))
  colnames(g)[1] <- "Gene.name"
  d <- cbind(g, d)
  output_file_name <- paste0(result_file_prefix, "shortest.path.csv")
  write.csv(d, file = paste0(output_dir_path, output_file_name), quote = F, row.names = FALSE)

  # Network Modules ---------------------------------------------------------
  cl <- igraph::cluster_louvain(net, weights = NULL)
  t <- as.data.frame(cl$membership)
  t1 <- as.data.frame(cl$names)
  t2 <- cbind(t1, t)
  colnames(t2)[1] <- "Gene.name"
  colnames(t2)[2] <- "membership"
  output_file_name <- paste0(result_file_prefix, "modules.in.ppi.network.csv")
  write.csv(t2, file = paste0(output_dir_path, output_file_name), quote = F, row.names = FALSE)


  # Plot Modules ------------------------------------------------------------
  g_grouped = net

  for(i in unique(igraph::V(net)$community)){
    groupV = which(igraph::V(net)$community == i)
    g_grouped = igraph::add_edges(g_grouped, combn(groupV, 2), attr = list(weight = 2))
  }

  l <- igraph::layout_nicely(g_grouped)

  # Writing to file
  output_file_name <- paste0(result_file_prefix, "modules.in.ppi.network.plot.pdf")
  pdf(file = paste0(output_dir_path, output_file_name), width = 5, height = 5)
  plot(
    cl,
    net,
    layout = igraph::layout_with_fr,
    vertex.size = 10,
    edge.width = 1,
    vertex.label.dist = 0.001,
    vertex.color = 'gold',
    vertex.frame.color = "#555555",
    edge.label = net$v,
    vertex.size = 1,
    edge.color = "gray",
    vertex.label.font = 1,
    edge.label.font = 1,
    edge.label.cex = 1,
    edge.arrow.size = 0.2,
    edge.curved = 0,
    vertex.label = igraph::V(net)$v,
    vertex.label.color = "black",
    vertex.label.cex = 0.5,
    vertex.label.cex = 0.5
  )
  dev.off()

  output_file_name <- paste0(result_file_prefix, "modules.in.ppi.network.plot.tiff")
  tiff(filename = paste0(output_dir_path, output_file_name), compression = "lzw")
  plot(
    cl,
    net,
    layout = igraph::layout_with_fr,
    vertex.size = 10,
    edge.width = 1,
    vertex.label.dist = 0.001,
    vertex.color = 'gold',
    vertex.frame.color = "#555555",
    edge.label = net$v,
    vertex.size = 1,
    edge.color = "gray",
    vertex.label.font = 1,
    edge.label.font = 1,
    edge.label.cex = 1,
    edge.arrow.size = 0.2,
    edge.curved = 0,
    vertex.label = igraph::V(net)$v,
    vertex.label.color = "black",
    vertex.label.cex = 0.5,
    vertex.label.cex = 0.5
  )
  dev.off()
}
