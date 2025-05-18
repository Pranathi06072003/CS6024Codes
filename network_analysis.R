library(SpiecEasi)
library(devtools)
library(igraph)
library(vegan)
library(Matrix)
library(reshape2)
library(plyr)
library(dplyr)
library(gridExtra)
library(grid)
library(igraph)
library(intergraph)
library(GGally)
library(network)
library(intergraph)
library(RColorBrewer)
library(dplyr)
library(microbiome)
library(NetCoMi)

library(SpiecEasi)
library(devtools)
library(igraph)
library(vegan)
library(Matrix)
library(reshape2)
library(plyr)
library(dplyr)
library(gridExtra)
library(grid)
library(future)
library(future.apply)
library(igraph)
library(intergraph)
library(GGally)
library(network)
library(intergraph)
library(RColorBrewer)
library(dplyr)
library(microbiome)


node_gen <- function(merge_data1, filepath_node, filepath_edge, metrics_path, degree_path){
  
  merge_data1.f <- microbiomeutilities::format_to_besthit(merge_data1)
  otu.table <- as(otu_table(merge_data1.f), "matrix")
  env.all <- as(sample_data(merge_data1.f), "data.frame")
  otu.table.all <- t(otu.table)
  row.names(env.all) == row.names(otu.table.all)
  
  if (nrow(otu.table) != nsamples(merge_data1.f)) {
    # If samples are not rows, transpose the OTU table
    otu.table <- t(otu.table)
    
  }
  num_cores <- parallel::detectCores() - 1  # Use one less than the total number of cores to avoid system overload
  run_spiec_easi <- function(otu_table, num_cores) {
    future.apply::future_lapply(1, function(x) {
      spiec.easi(
        otu_table,
        method = 'mb',
        lambda.min.ratio = 1e-2,
        nlambda = 20,
        icov.select.params = list(rep.num = 50, ncores = 1)  # ncores = 1 because future will handle parallelization
      )
    }, future.seed = TRUE)[[1]]  # Set future.seed = TRUE as an argument to future_lapply
  }
  
  # Set up the parallel backend
  spieceasi.net <-  run_spiec_easi(otu.table, num_cores)
  # The result is an object of class "spiec.easi"
  print(spieceasi.net)
  
  # Reset the plan to default
  plan(sequential)
  
  ##Generating the adjacency matrix
  adjacency_matrix <- symBeta(getOptBeta(spieceasi.net))
  colnames(adjacency_matrix) <- rownames(adjacency_matrix) <- colnames(otu.table)
  adjacency_matrix_return <- as.matrix(adjacency_matrix)
  
  return(adjacency_matrix_return)
}

pre_physeq <- subset_samples(physeq_entire, 
                                collection_point %in% c("pre"))
post_physeq <- subset_samples(physeq_entire, 
                                  collection_point %in% c("post"))
in_physeq <- subset_samples(physeq_entire, 
                            collection_point %in% c("in"))

complete_data1 <- core(pre_physeq, detection=50, prevalence = 0.01, include.lowest = FALSE)
non_zero_samples <- sample_sums(complete_data1) > 0
complete_data1 <- prune_samples(non_zero_samples, complete_data1)
complete_data1 <- prune_taxa(taxa_sums(complete_data1) > 0, complete_data1)
adj_complete <- node_gen(complete_data1, "C:/Desktop/PranathiR/IITM/ComputationalSystemsBiologyLab/QZA_Results/nodes_1.csv", "C:/Desktop/PranathiR/IITM/ComputationalSystemsBiologyLab/QZA_Results/edges_1.csv", "C:/Desktop/PranathiR/IITM/ComputationalSystemsBiologyLab/QZA_Results/metrics_1.csv", "C:/Desktop/PranathiR/IITM/ComputationalSystemsBiologyLab/QZA_Results/degree_5.png")
adj_pre <- adj_complete

in_physeq1 <- core(in_physeq, detection=50, prevalence = 0.01, include.lowest = FALSE)
non_zero_samples <- sample_sums(in_physeq1) > 0
marine_physeq1 <- prune_samples(non_zero_samples, in_physeq1)
marine_physeq1 <- prune_taxa(taxa_sums(in_physeq1) > 0, in_physeq1)
adj_in <- node_gen(in_physeq1, "C:/Desktop/PranathiR/IITM/ComputationalSystemsBiologyLab/QZA_Results/nodes_marine.csv", "C:/Desktop/PranathiR/IITM/ComputationalSystemsBiologyLab/QZA_Results/edges_marine.csv", "C:/Desktop/PranathiR/IITM/ComputationalSystemsBiologyLab/QZA_Results/metrics_marine.csv", "C:/Desktop/PranathiR/IITM/ComputationalSystemsBiologyLab/QZA_Results/degree_marine.png")

post_physeq1 <- core(post_physeq, detection=50, prevalence = 0.01, include.lowest = FALSE)
non_zero_samples <- sample_sums(post_physeq1) > 0
post_physeq1 <- prune_samples(non_zero_samples, post_physeq1)
post_physeq1 <- prune_taxa(taxa_sums(post_physeq1) > 0, post_physeq1)
adj_post <- node_gen(post_physeq1, "C:/Desktop/PranathiR/IITM/ComputationalSystemsBiologyLab/QZA_Results/nodes_sediment.csv", "C:/Desktop/PranathiR/IITM/ComputationalSystemsBiologyLab/QZA_Results/edges_sediment.csv", "C:/Desktop/PranathiR/IITM/ComputationalSystemsBiologyLab/QZA_Results/metrics_sediment.csv", "C:/Desktop/PranathiR/IITM/ComputationalSystemsBiologyLab/QZA_Results/degree_sediment.png")

net_pre <- netConstruct(data=adj_pre,
                         normMethod = "none", zeroMethod = "none",
                         sparsMethod = "none", dataType = "condDependence",
                         verbose = 3)

- 
net_in <- netConstruct(data=adj_in,
                             normMethod = "none", zeroMethod = "none",
                             sparsMethod = "none", dataType = "condDependence",
                             verbose = 3)

net_post <- netConstruct(data=adj_post,
                       normMethod = "none", zeroMethod = "none",
                       sparsMethod = "none", dataType = "condDependence",
                       verbose = 3)

in_analysis <- netAnalyze(net_in,centrLCC = FALSE,avDissIgnoreInf = TRUE,
                              hubPar = "degree", hubQuant = 0.95,
                              normDeg = TRUE, normBetw = TRUE, normEigen = TRUE, verbose=2)


pre_analysis <- netAnalyze(net_pre,centrLCC = FALSE,avDissIgnoreInf = TRUE,
                                hubPar = "degree", hubQuant = 0.95,
                                normDeg = TRUE, normBetw = TRUE, normEigen = TRUE, verbose=2)

post_analysis <- netAnalyze(net_post,centrLCC = FALSE,avDissIgnoreInf = TRUE,
                                    hubPar = "degree", hubQuant = 0.95,
                                    normDeg = TRUE, normBetw = TRUE, normEigen = TRUE, verbose=2)


plot_pre <- plot(pre_analysis,
                      layout = "spring",
                      sameLayout = FALSE,
                      # layout = "circle",
                      edgeFilter="threshold",
                      edgeFilterPar = 0.05,
                      layoutGroup = "union",
                      nodeColor = "cluster", 
                      rmSingles = "inboth",
                      nodeSize = "degree",
                      highlightHubs = TRUE,
                      title1 = "Pre-Flight", 
                      showTitle = TRUE,
                      cexNodes = 1,
                      cexHubs = 6,
                      cexHubLabels = 3,
                      cexTitle = 1.5,
                      mar =c(3, 4, 5, 2),
                      labelScale = TRUE,                            labelFont = 0,
                      labels = FALSE
)


plot_in <- plot(in_analysis,
                 layout = "spring",
                 sameLayout = FALSE,
                 # layout = "circle",
                 edgeFilter="threshold",
                 edgeFilterPar = 0.05,
                 layoutGroup = "union",
                 nodeColor = "cluster", 
                 rmSingles = "inboth",
                 nodeSize = "degree",
                 highlightHubs = TRUE,
                 title1 = "In-Flight", 
                 showTitle = TRUE,
                 cexNodes = 1,
                 cexHubs = 6,
                 cexHubLabels = 3,
                 cexTitle = 1.5,
                 mar =c(3, 4, 5, 2),
                 labelScale = TRUE,                            labelFont = 0,
                 labels = FALSE
)


plot_post <- plot(post_analysis,
                 layout = "spring",
                 sameLayout = FALSE,
                 # layout = "circle",
                 edgeFilter="threshold",
                 edgeFilterPar = 0.05,
                 layoutGroup = "union",
                 nodeColor = "cluster", 
                 rmSingles = "inboth",
                 nodeSize = "degree",
                 highlightHubs = TRUE,
                 title1 = "Post-Flight", 
                 showTitle = TRUE,
                 cexNodes = 1,
                 cexHubs = 6,
                 cexHubLabels = 3,
                 cexTitle = 1.5,
                 mar =c(3, 4, 5, 2),
                 labelScale = TRUE,                            labelFont = 0,
                 labels = FALSE
)

adj_in_mod <- adj_post  # Make a copy
adj_in_mod[adj_in_mod < 0] <- 0  # Replace negatives with 0

g_in <- graph_from_adjacency_matrix(adj_in_mod, mode = "undirected", weighted = TRUE, diag = FALSE)
g <- g_in
library(igraph)
# Number of nodes and edges
vcount(g)
ecount(g)

# Density of the graph
graph.density(g)

# Degree (number of connections per node)
mean(degree(g))

# Clustering coefficient
transitivity(g, type = "global")  # Global clustering coefficient
transitivity(g, type = "local")   # Local clustering per node

# Betweenness centrality
mean(betweenness(g))

# Closeness centrality
closeness(g)

# Eigenvector centrality
eigen_centrality(g)$vector

# Average path length
average.path.length(g)

# Diameter (longest shortest path)
diameter(g)

# Assortativity (degree correlation)
assortativity_degree(g)

# Modularity (community structure)
cluster <- cluster_fast_greedy(g)
modularity(cluster)

deg <- degree(g)

# Basic histogram of degrees
hist(deg,
     breaks = seq(0, max(deg) + 1, by = 1),
     main = "Degree Distribution",
     xlab = "Degree",
     ylab = "Frequency",
     col = "skyblue",
     border = "white")

