
# Load necessary libraries
library(qiime2R)
library(phyloseq)
library(dplyr)
library(tidyr)
library(qiime2R)
library(Biostrings)
library(phyloseq)


otu_table_path = "C:/Desktop/PranathiR/IITM/project2/table.qza"
taxonomy_path = "C:/Desktop/PranathiR/IITM/project2/taxonomy.qza"
table_qza <- read_qza(otu_table_path)
taxonomy_qza <- read_qza(taxonomy_path)
feature_table <- as.data.frame(as.matrix(table_qza$data))  # Convert to a matrix or data frame

taxonomy <- as.data.frame(taxonomy_qza$data)

# Step 4: Split the 'Taxon' column into separate columns for Kingdom, Phylum, Class, Order, Family, Genus, Species
taxonomy_separated <- taxonomy %>%
  separate(Taxon, into = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"), 
           sep = ";", fill = "right", extra = "drop")

# Step 5: Remove the 'Confidence' column (if you don't want it)
taxonomy_cleaned <- taxonomy_separated[, -ncol(taxonomy_separated)]

# Step 6: Trim leading and trailing white spaces from all taxonomic columns
taxonomy_cleaned[] <- lapply(taxonomy_cleaned, trimws)

# Step 7: Set 'Feature.ID' as row names and remove it from the dataframe
rownames(taxonomy_cleaned) <- taxonomy$Feature.ID
taxonomy_cleaned <- taxonomy_cleaned[, -which(names(taxonomy_cleaned) == "Feature.ID")]

# Step 8: Convert the cleaned taxonomy data frame to a matrix for phyloseq
taxonomy_matrix <- as.matrix(taxonomy_cleaned)

# Step 9: Create the phyloseq taxonomy object
taxonomy_ps <- tax_table(taxonomy_matrix)

# Step 10: Create the OTU table for phyloseq
otu_table_ps <- otu_table(feature_table, taxa_are_rows = TRUE)

physeq_entire <- phyloseq(otu_table_ps, taxonomy_ps)

sample_names_vector <- sample_names(physeq_entire)
updated_sample_names <- gsub("_1.fastq.gz$", "", sample_names_vector)
sample_names(physeq_entire) <- updated_sample_names
metadata <- read.csv("C:/Desktop/PranathiR/IITM/project2/metadata.csv", row.names = 1, stringsAsFactors = FALSE)
sample_data(physeq_entire) <- sample_data(metadata)

