# Load necessary libraries
library(tidyverse)
library(phyloseq)
library(DESeq2)

# Load the rarefied phyloseq object and metadata
load("ps_rare.RData")
metadata <- read_tsv("LH_gastric_cancer_metadata.tsv")

# Transform the sample counts to avoid zero counts
ps_plus1 <- transform_sample_counts(ps_rare, function(x) x + 1)

# Set "Healthy control (HC)" as the reference level for the Group variable
sample_data(ps_plus1)$Group <- factor(sample_data(ps_plus1)$Group, levels = c("Healthy control (HC)", "Gastric cancer (GC)"))

# Subset to retain only "Healthy control (HC)" and "Gastric cancer (GC)" samples
ps_plus1 <- subset_samples(ps_plus1, Group %in% c("Healthy control (HC)", "Gastric cancer (GC)"))

# Conduct DESeq analysis with 'Group' as the independent variable
ps_deseq <- phyloseq_to_deseq2(ps_plus1, ~ Group)
ps_fecal <- DESeq(ps_deseq)

# Extract results, setting Healthy control (HC) as the reference group
res <- results(ps_fecal, tidy = TRUE)

# Filter significant ASVs based on p-adjusted value and log2FoldChange thresholds
sigASVs <- res %>%
  filter(padj < 0.01 & abs(log2FoldChange) > 2) %>%
  dplyr::rename(ASV = row)

# Get a vector of ASV names for significant ASVs
sigASVs_vec <- sigASVs %>%
  pull(ASV)

# Prune the phyloseq object to retain only the significant ASVs
ps_filt <- prune_taxa(sigASVs_vec, ps_rare)

# Add taxonomy to the DESeq results table
merged_results <- tax_table(ps_rare) %>%
  as.data.frame() %>%
  rownames_to_column(var = "ASV") %>%
  right_join(sigASVs, by = "ASV") %>%
  arrange(desc(abs(log2FoldChange))) %>%
  mutate(Genus = make.unique(Genus)) %>%
  mutate(Genus = factor(Genus, levels = unique(Genus)))

# Select the top 50 most significant ASVs by absolute log2FoldChange
top_50_ASVs <- merged_results %>%
  slice_max(order_by = abs(log2FoldChange), n = 50) %>%
  mutate(Genus = factor(Genus, levels = unique(Genus)))

# Create bar plot of top 50 significant ASVs with error bars
plot <- ggplot(top_50_ASVs) +
  geom_bar(aes(x=Genus, y=log2FoldChange), stat="identity")+
  geom_errorbar(aes(x = Genus, ymin = log2FoldChange - lfcSE, ymax = log2FoldChange + lfcSE), width = 0.2) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 10),
        axis.title = element_text(size = 12)) +
  labs(x = "Genus", y = "Log2 Fold Change", title = "Top 50 Significant ASVs: Gastric Cancer vs. Healthy Control")


# Save the bar plot as a PNG file in the project folder
ggsave("ps_deseq_plot_top50_GC.png", plot = plot, height = 6, width = 10)

# Print plot
print(plot)

