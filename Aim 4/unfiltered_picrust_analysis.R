#### Load packages ####
# Load all necessary libraries
library(readr)
library(ggpicrust2)
library(tibble)
library(tidyverse)
library(ggprism)
library(patchwork)
library(DESeq2)
library(ggh4x)
library(dplyr)


#### Import files and preparing tables ####
#Importing the pathway PICrsut2
abundance_file <- "pathway_abundance.tsv"
abundance_data <- read_delim(abundance_file, delim = "\t", skip = 1, col_names = TRUE, trim_ws = TRUE)
abundance_data  =as.data.frame(abundance_data)

#Import your metadata file, no need to filter yet
metadata <- as.data.frame(read_delim("LH_gastric_cancer_metadata.tsv"))

#Remove NAs for your column of interest in this case subject
metadata = metadata[!is.na(metadata$Stage_FusoAbundance),]

#Filter your metadata to include only 2 at a time
healthyFusoLowvsHigh_metadata <- filter(metadata, Stage_FusoAbundance %in% c("Healthy control (HC)_Low", "Healthy control (HC)_High"))
cgFusoLowvsHigh_metadata <- filter(metadata, Stage_FusoAbundance %in% c("Chronic gastritis (CG)_Low", "Chronic gastritis (CG)_High"))
imFusoLowvsHigh_metadata <- filter(metadata, Stage_FusoAbundance %in% c("Intestinal metaplasia (IM）_Low", "Intestinal metaplasia (IM）_High"))
ineoFusoLowvsHigh_metadata <- filter(metadata, Stage_FusoAbundance %in% c("Intraepithelial neoplasia (IN)_Low", "Intraepithelial neoplasia (IN)_High"))
cancerFusoLowvsHigh_metadata <-filter(metadata, Stage_FusoAbundance %in% c("Gastric cancer (GC)_Low", "Gastric cancer (GC)_High"))


#######Comparing Healthy Fuso Low vs High#######

#Filtering the abundance table to only include samples that are in the filtered metadata
healthyFusoLowvsHigh_sample_names = healthyFusoLowvsHigh_metadata$'sample-id'
healthyFusoLowvsHigh_sample_names = append(healthyFusoLowvsHigh_sample_names, "#OTU ID")
healthyFusoLowvsHigh_abundance_data_filtered = abundance_data[, colnames(abundance_data) %in% healthyFusoLowvsHigh_sample_names] #This step is the actual filtering

#Removing individuals with no data that caused a problem for pathways_daa()
healthyFusoLowvsHigh_abundance_data_filtered =  healthyFusoLowvsHigh_abundance_data_filtered[, colSums(healthyFusoLowvsHigh_abundance_data_filtered != 0) > 0]

#Ensuring the rownames for the abundance_data_filtered is empty. This is required for their functions to run.
rownames(healthyFusoLowvsHigh_abundance_data_filtered) = NULL

#verify samples in metadata match samples in abundance_data
healthyFusoLowvsHigh_abun_samples = rownames(t(healthyFusoLowvsHigh_abundance_data_filtered[,-1])) #Getting a list of the sample names in the newly filtered abundance data
healthyFusoLowvsHigh_metadata = healthyFusoLowvsHigh_metadata[healthyFusoLowvsHigh_metadata$`sample-id` %in% healthyFusoLowvsHigh_abun_samples,] #making sure the filtered metadata only includes these samples

#### DESEq Healthy Fuso LOW vs HIGH ####
#Perform pathway DAA using DESEQ2 method
healthyFusoLowvsHigh_abundance_daa_results_df <- pathway_daa(abundance = healthyFusoLowvsHigh_abundance_data_filtered %>% column_to_rownames("#OTU ID"), 
                                        metadata = healthyFusoLowvsHigh_metadata, group = "Stage_FusoAbundance", daa_method = "DESeq2")

# Annotate MetaCyc pathway so they are more descriptive
healthyFusoLowvsHigh_metacyc_daa_annotated_results_df <- pathway_annotation(pathway = "MetaCyc", 
                                                                            daa_results_df = healthyFusoLowvsHigh_abundance_daa_results_df, ko_to_kegg = FALSE)

# Filter p-values to only significant ones
healthyFusoLowvsHigh_feature_with_p_0.05 <- healthyFusoLowvsHigh_abundance_daa_results_df %>% filter(p_values < 0.05)

#Changing the pathway column to description for the results 
healthyFusoLowvsHigh_feature_desc = inner_join(healthyFusoLowvsHigh_feature_with_p_0.05,healthyFusoLowvsHigh_metacyc_daa_annotated_results_df, by = "feature")
healthyFusoLowvsHigh_feature_desc$feature = healthyFusoLowvsHigh_feature_desc$description
healthyFusoLowvsHigh_feature_desc = healthyFusoLowvsHigh_feature_desc[,c(1:7)]
colnames(healthyFusoLowvsHigh_feature_desc) = colnames(healthyFusoLowvsHigh_feature_with_p_0.05)

#Changing the pathway column to description for the abundance table
healthyFusoLowvsHigh_abundance = healthyFusoLowvsHigh_abundance_data_filtered %>% filter(`#OTU ID` %in% healthyFusoLowvsHigh_feature_with_p_0.05$feature)
colnames(healthyFusoLowvsHigh_abundance)[1] = "feature"
healthyFusoLowvsHigh_abundance_desc = inner_join(healthyFusoLowvsHigh_abundance,healthyFusoLowvsHigh_metacyc_daa_annotated_results_df, by = "feature")
healthyFusoLowvsHigh_abundance_desc$feature = healthyFusoLowvsHigh_abundance_desc$description
#this line will change for each dataset. 34 represents the number of samples in the filtered abundance table
healthyFusoLowvsHigh_abundance_desc = healthyFusoLowvsHigh_abundance_desc[,-c(14:ncol(healthyFusoLowvsHigh_abundance_desc))] 

# Generate a heatmap
healthyFusoLowvsHigh_pathway_heatmap <- pathway_heatmap(abundance = healthyFusoLowvsHigh_abundance_desc %>% column_to_rownames("feature"), metadata = healthyFusoLowvsHigh_metadata, group = "Stage_FusoAbundance")

# Generate pathway PCA plot
healthyFusoLowvsHigh_pathway_pca <- pathway_pca(abundance = healthyFusoLowvsHigh_abundance_data_filtered %>% column_to_rownames("#OTU ID"), metadata = healthyFusoLowvsHigh_metadata, group = "Stage_FusoAbundance")

# Generating a bar plot representing log2FC from the custom deseq2 function

# Go to the Deseq2 function script and update the metadata category of interest

# Lead the function in
source("DESeq2_function.R")

# Run the function on your own data
healthyFusoLowvsHigh_res =  DEseq2_function(healthyFusoLowvsHigh_abundance_data_filtered, metadata, "Stage_FusoAbundance")
healthyFusoLowvsHigh_res$feature =rownames(healthyFusoLowvsHigh_res)
healthyFusoLowvsHigh_res_desc = inner_join(healthyFusoLowvsHigh_res,healthyFusoLowvsHigh_metacyc_daa_annotated_results_df, by = "feature")
healthyFusoLowvsHigh_res_desc = healthyFusoLowvsHigh_res_desc[, -c(8:13)]
View(healthyFusoLowvsHigh_res_desc)

# Filter to only include significant pathways
sig_res = res_desc %>%
  filter(pvalue < 0.05)
# You can also filter by Log2fold change (like include only +-2 log fold change)

sig_res <- sig_res[order(sig_res$log2FoldChange),]
ggplot(data = sig_res, aes(y = reorder(description, sort(as.numeric(log2FoldChange))), x= log2FoldChange, fill = pvalue))+
  geom_bar(stat = "identity")+ 
  theme_bw()+
  labs(x = "Log2FoldChange", y="Pathways")
