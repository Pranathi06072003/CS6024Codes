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

marine_physeq <- subset_samples(physeq_entire, 
                                collection_point %in% c("in"))
sediment_physeq <- subset_samples(physeq_entire, 
                                  collection_point %in% c("post"))
coral_sponge_physeq <- subset_samples(physeq_entire, 
                                      collection_point %in% c("pre"))

complete_data1 <- core(physeq_combined_Algo, detection=50, prevalence = 0.01, include.lowest = FALSE)
non_zero_samples <- sample_sums(complete_data1) > 0
complete_data1 <- prune_samples(non_zero_samples, complete_data1)
complete_data1 <- prune_taxa(taxa_sums(complete_data1) > 0, complete_data1)
adj_complete <- node_gen(complete_data1, "C:/Desktop/PranathiR/IITM/ComputationalSystemsBiologyLab/QZA_Results/nodes_1.csv", "C:/Desktop/PranathiR/IITM/ComputationalSystemsBiologyLab/QZA_Results/edges_1.csv", "C:/Desktop/PranathiR/IITM/ComputationalSystemsBiologyLab/QZA_Results/metrics_1.csv", "C:/Desktop/PranathiR/IITM/ComputationalSystemsBiologyLab/QZA_Results/degree_5.png")

marine_physeq1 <- core(marine_physeq, detection=50, prevalence = 0.01, include.lowest = FALSE)
non_zero_samples <- sample_sums(marine_physeq1) > 0
marine_physeq1 <- prune_samples(non_zero_samples, marine_physeq1)
marine_physeq1 <- prune_taxa(taxa_sums(marine_physeq1) > 0, marine_physeq1)
adj_marine <- node_gen(marine_physeq1, "C:/Desktop/PranathiR/IITM/project2/nodes_in.csv", "C:/Desktop/PranathiR/IITM/project2/edges_in.csv", "C:/Desktop/PranathiR/IITM/project2/metrics_in.csv", "C:/Desktop/PranathiR/IITM/project2/degree_in.png")

sediment_physeq1 <- core(sediment_physeq, detection=50, prevalence = 0.01, include.lowest = FALSE)
non_zero_samples <- sample_sums(sediment_physeq1) > 0
sediment_physeq1 <- prune_samples(non_zero_samples, sediment_physeq1)
sediment_physeq1 <- prune_taxa(taxa_sums(sediment_physeq1) > 0, sediment_physeq1)
adj_sediment <-  node_gen(marine_physeq1, "C:/Desktop/PranathiR/IITM/project2/nodes_out.csv", "C:/Desktop/PranathiR/IITM/project2/edges_out.csv", "C:/Desktop/PranathiR/IITM/project2/metrics_out.csv", "C:/Desktop/PranathiR/IITM/project2/degree_out.png")

coral_sponge_physeq1 <- core(coral_sponge_physeq, detection=50, prevalence = 0.01, include.lowest = FALSE)
non_zero_samples <- sample_sums(coral_sponge_physeq1) > 0
coral_sponge_physeq1 <- prune_samples(non_zero_samples, coral_sponge_physeq1)
coral_sponge_physeq1 <- prune_taxa(taxa_sums(coral_sponge_physeq1) > 0, coral_sponge_physeq1)
adj_coral <- node_gen(coral_sponge_physeq1, "C:/Desktop/PranathiR/IITM/ComputationalSystemsBiologyLab/QZA_Results/nodes_coral.csv", "C:/Desktop/PranathiR/IITM/ComputationalSystemsBiologyLab/QZA_Results/edges_coral.csv", "C:/Desktop/PranathiR/IITM/ComputationalSystemsBiologyLab/QZA_Results/metrics_coral.csv", "C:/Desktop/PranathiR/IITM/ComputationalSystemsBiologyLab/QZA_Results/degree_coral.png")


combined <- netConstruct(data=adj_complete,
                         normMethod = "none", zeroMethod = "none",
                         sparsMethod = "none", dataType = "condDependence",
                         verbose = 3)

- 
  marine_net <- netConstruct(data=adj_marine,
                             normMethod = "none", zeroMethod = "none",
                             sparsMethod = "none", dataType = "condDependence",
                             verbose = 3)



marine_analysis <- netAnalyze(marine_net,centrLCC = FALSE,avDissIgnoreInf = TRUE,
                              hubPar = "degree", hubQuant = 0.95,
                              normDeg = TRUE, normBetw = TRUE, normEigen = TRUE, verbose=2)

sediment_net <- netConstruct(data=adj_sediment,
                             normMethod = "none", zeroMethod = "none",
                             sparsMethod = "none", dataType = "condDependence",
                             verbose = 3)

sediment_analysis <- netAnalyze(sediment_net,centrLCC = FALSE,avDissIgnoreInf = TRUE,
                                hubPar = "degree", hubQuant = 0.95,
                                normDeg = TRUE, normBetw = TRUE, normEigen = TRUE, verbose=2)

coral_sponge_net <- netConstruct(data=adj_coral,
                                 normMethod = "none", zeroMethod = "none",
                                 sparsMethod = "none", dataType = "condDependence",
                                 verbose = 3)

coral_sponge_analysis <- netAnalyze(coral_sponge_net,centrLCC = FALSE,avDissIgnoreInf = TRUE,
                                    hubPar = "degree", hubQuant = 0.95,
                                    normDeg = TRUE, normBetw = TRUE, normEigen = TRUE, verbose=2)


##Performing PERMANOVA
library(phyloseq)
library(vegan)
library(limma)

remove_prefixes <- function(physeq_1) {
  tax_table1 = tax_table(physeq_1)
  taxa_names <- gsub("^d_", "", tax_table1) # Remove 'k_' prefix
  taxa_names <- gsub("^p_", "", taxa_names)
  taxa_names <- gsub("^c_", "", taxa_names)
  taxa_names <- gsub("^o_", "", taxa_names)
  taxa_names <- gsub("^f_", "", taxa_names)
  taxa_names <- gsub("^s_", "", taxa_names) # Remove 's_' prefix
  taxa_names <- gsub("^_", "", taxa_names) # Remove 's_' prefix
  tax_names2 <- tax_table(taxa_names)
  otu_table1 <- otu_table(physeq_1)
  metadata <- sam_data(physeq_1)
  physeq2 <- phyloseq(otu_table1, tax_names2, metadata)
  
  
  return(physeq2)
}

combined_data1 <- remove_prefixes(physeq_entire)
#combined_data1 <- physeq_combined_Alg
combined_genus_agglomeration1 <- combined_data1
complete_data1 <- prune_taxa(taxa_sums(combined_genus_agglomeration1) > 0, combined_genus_agglomeration1)

complete_data1 <- core(complete_data1, detection=50, prevalence = 0.01, include.lowest = FALSE)
non_zero_samples <- sample_sums(complete_data1) > 0
complete_data1 <- prune_samples(non_zero_samples, complete_data1)

combined_N <- transform_sample_counts(complete_data1, function(x) x / sum(x))

dist_matrix <- phyloseq::distance(combined_N, method = "bray")
metadata <- data.frame(sample_data(combined_N))

result_1 <- adonis2(dist_matrix ~ collection_point, data = metadata, permutations = 999)
