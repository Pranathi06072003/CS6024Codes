library(phyloseq)
library(ggplot2)
library(RColorBrewer)
library(pheatmap)
library(ggpubr)
library(dplyr)

combined_data <- physeq_entire
combined_data <- subset_samples(combined_data, !(collection_point %in% c("post control 1", "post control 2")))

# Convert sample data to a dataframe
sample_data_df <- as(sample_data(combined_data), "data.frame")
unique_sample_origins <- unique(sample_data_df$collection_point)

# Shannon Diversity Analysis
alpha_div <- estimate_richness(combined_data, measures = c("Shannon"))
alpha_div_df <- as.data.frame(alpha_div)
metadata_df <- as(sample_data(combined_data), "data.frame")
alpha_div_df <- cbind(metadata_df, alpha_div_df)

library(dplyr)
color_palette <- scale_fill_brewer(palette = "Set3") # Or customize colors as you prefer

alpha_div <- ggplot(alpha_div_df, aes(x = collection_point, y = Shannon, fill = collection_point)) +
  geom_boxplot() +
  theme_minimal() +
  labs(title = "Shannon Diversity Index by Sample Origin",
       x = "Sample Origin",
       y = "Shannon Diversity Index",
       fill = "Sample Origin") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + # Rotate x-axis labels for readability
  color_palette

# Load necessary libraries
library(phyloseq)
library(ggplot2)
library(vegan)  # For performing PCA
library(dplyr)
library(NetCoMi)
library(microbiome)

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

combined_data1 <- remove_prefixes(combined_data)
# Load necessary libraries
library(phyloseq)
library(ggplot2)
library(vegan)  # For performing PCA
library(dplyr)
library(NetCoMi)
library(microbiome)

#combined_data1 <- remove_prefixes(physeq_combined_Algo)
#combined_data1 <- physeq_combined_Alg
combined_genus_agglomeration1 <- combined_data1
complete_data1 <- prune_taxa(taxa_sums(combined_genus_agglomeration1) > 0, combined_genus_agglomeration1)

complete_data1 <- core(complete_data1, detection=50, prevalence = 0.01, include.lowest = FALSE)
non_zero_samples <- sample_sums(complete_data1) > 0
complete_data1 <- prune_samples(non_zero_samples, complete_data1)

combined_N <- transform_sample_counts(complete_data1, function(x) x / sum(x))

#merged_NMDS <- ordinate(combined_data1, method = "NMDS", distance = "bray")
#merged_PCoA <- ordinate(complete_data1, method = "PCoA", distance = "bray")
merged_PCoA_N <- ordinate(complete_data1, method = "PCoA", distance = "bray")
#combined_data1 <- physeq_combined_Alg

pc1pc2 <- plot_ordination(complete_data1, merged_PCoA_N, axes= c(1,2), type= "samples", color="collection_point", title = "B")


dist_matrix <- phyloseq::distance(combined_N, method = "bray")
metadata <- data.frame(sample_data(combined_N))

result_1 <- adonis2(dist_matrix ~ collection_point, data = metadata, permutations = 999)
