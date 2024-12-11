library(tidyverse)

meta <- read_delim("gastric_cancer_metadata.tsv", delim = "\t")
#### Format sample metadata ####
# Save everything except sampleid as new data frame
metadata <- as.data.frame(meta[,-1])
# Make sampleids the rownames
rownames(metadata)<- meta$'sample-id'

otu_table <- read.delim(file="table.txt", skip=1, row.names = 1)

taxonomy <- read.delim("taxonomy.tsv", header = TRUE, row.names = 1, sep="\t")

# Step 1: Identify the OTU corresponding to Fusobacterium in your taxonomy
fuso_otu <- rownames(taxonomy)[grepl("g__Fusobacterium", taxonomy$Taxon)]

# Step 2: Initialize all samples as "No" in the metadata
metadata$fuso_detected <- "Low"

# Step 3: If Fusobacterium OTU(s) are found, proceed with percentage-based detection
if (length(fuso_otu) > 0) {
  # Get the absolute abundance of Fusobacterium for each sample
  fuso_abundance <- otu_table[fuso_otu, , drop = FALSE]
  
  # Get the total abundance per sample
  total_abundance_per_sample <- colSums(otu_table)
  
  # Convert Fusobacterium abundance to percentage
  fuso_percentage <- colSums(fuso_abundance) / total_abundance_per_sample * 100
  
  # Step 4: Set detection threshold (e.g., > 1% abundance)
  detected_samples <- fuso_percentage > 1
  
  
  # Set "Yes" for samples with detected Fusobacterium (i.e., percentage > 1%)
  metadata$fuso_detected[detected_samples] <- "High"
}
# Rename the column after all processing
colnames(metadata)[colnames(metadata) == "fuso_detected"] <- "Fusobacterium_abundance"

metadata <- metadata %>%
  mutate(Stage_FusoAbundance = paste(Group, Fusobacterium_abundance, sep = "_"))

metadata <- metadata %>%
  mutate(Stage_Pylori = paste(Group, `H. pylori status                                 13 C-urea breath test`, sep = "_"))

# Save the updated metadata to a new TSV file
write.table(metadata, file = "LH_gastric_cancer_metadata.tsv", 
            sep = "\t", quote = FALSE, col.names = TRUE, row.names = TRUE)