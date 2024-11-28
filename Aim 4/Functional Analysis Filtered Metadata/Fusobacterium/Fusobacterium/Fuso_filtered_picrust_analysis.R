##### Install packages #####
# Start by installing all necessary packages when asked if you want to install
# from source, please just type Yes in the terminal below

# If you don't have BiocManager, here is the code to install it
# A lot of you probably already have this so you can skip
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

# Create a list of all the packages you need to install
pkgs <- c("ALDEx2", "SummarizedExperiment", "Biobase", "devtools", 
          "ComplexHeatmap", "BiocGenerics", "metagenomeSeq", 
          "Maaslin2", "edgeR", "lefser", "limma", "KEGGREST", "DESeq2")

# Use the above list to install all the packages using a for loop
for (pkg in pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE))
    BiocManager::install(pkg)
}
# when asked if you want to update all, some or none please type "n" for none

# After installing all of its above dependencies, install ggpicrust2
install.packages("ggpicrust2")

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
source("Fuso_DESeq2_function.R")

# Run the function on your own data
healthyFusoLowvsHigh_res =  DEseq2_function(healthyFusoLowvsHigh_abundance_data_filtered,healthyFusoLowvsHigh_metadata, "Stage_FusoAbundance")
healthyFusoLowvsHigh_res$feature =rownames(healthyFusoLowvsHigh_res)
healthyFusoLowvsHigh_res_desc = inner_join(healthyFusoLowvsHigh_res,healthyFusoLowvsHigh_metacyc_daa_annotated_results_df, by = "feature")
healthyFusoLowvsHigh_res_desc = healthyFusoLowvsHigh_res_desc[, -c(8:13)]
View(healthyFusoLowvsHigh_res_desc)

# Filter to only include significant pathways
healthyFusoLowvsHigh_sig_res = healthyFusoLowvsHigh_res_desc %>%
  filter(pvalue < 0.05)
# You can also filter by Log2fold change (like include only +-2 log fold change)
#healthyFusoLowvsHigh_sig_res <- sig_res[order(healthyFusoLowvsHigh_sig_res$log2FoldChange),]

Healthy_FilteredFusoLowVSHigh_DESeq2 <- ggplot(data = healthyFusoLowvsHigh_sig_res, aes(y = reorder(description, sort(as.numeric(log2FoldChange))), x= log2FoldChange, fill = pvalue))+
  geom_bar(stat = "identity")+ 
  theme_bw()+
  labs(x = "Log2FoldChange", y="Pathways")

ggsave("Healthy_FilteredFusoLowVSHigh_DESeq2.png", Healthy_FilteredFusoLowVSHigh_DESeq2)



#######Comparing Chronic Gastritis Fuso Low vs High#######

#Filtering the abundance table to only include samples that are in the filtered metadata
cgFusoLowvsHigh_sample_names = cgFusoLowvsHigh_metadata$'sample-id'
cgFusoLowvsHigh_sample_names = append(cgFusoLowvsHigh_sample_names, "#OTU ID")
cgFusoLowvsHigh_abundance_data_filtered = abundance_data[, colnames(abundance_data) %in% cgFusoLowvsHigh_sample_names] #This step is the actual filtering

#Removing individuals with no data that caused a problem for pathways_daa()
cgFusoLowvsHigh_abundance_data_filtered =  cgFusoLowvsHigh_abundance_data_filtered[, colSums(cgFusoLowvsHigh_abundance_data_filtered != 0) > 0]

#Ensuring the rownames for the abundance_data_filtered is empty. This is required for their functions to run.
rownames(cgFusoLowvsHigh_abundance_data_filtered) = NULL

#verify samples in metadata match samples in abundance_data
cgFusoLowvsHigh_abun_samples = rownames(t(cgFusoLowvsHigh_abundance_data_filtered[,-1])) #Getting a list of the sample names in the newly filtered abundance data
cgFusoLowvsHigh_metadata = cgFusoLowvsHigh_metadata[cgFusoLowvsHigh_metadata$`sample-id` %in% cgFusoLowvsHigh_abun_samples,] #making sure the filtered metadata only includes these samples

#### DESEq Chronic Gastritis Fuso LOW vs HIGH ####
#Perform pathway DAA using DESEQ2 method
cgFusoLowvsHigh_abundance_daa_results_df <- pathway_daa(abundance = cgFusoLowvsHigh_abundance_data_filtered %>% column_to_rownames("#OTU ID"), 
                                                        metadata = cgFusoLowvsHigh_metadata, group = "Stage_FusoAbundance", daa_method = "DESeq2")

# Annotate MetaCyc pathway so they are more descriptive
cgFusoLowvsHigh_metacyc_daa_annotated_results_df <- pathway_annotation(pathway = "MetaCyc", 
                                                                       daa_results_df = cgFusoLowvsHigh_abundance_daa_results_df, ko_to_kegg = FALSE)

# Filter p-values to only significant ones
cgFusoLowvsHigh_feature_with_p_0.05 <- cgFusoLowvsHigh_abundance_daa_results_df %>% filter(p_values < 0.05)

#Changing the pathway column to description for the results 
cgFusoLowvsHigh_feature_desc = inner_join(cgFusoLowvsHigh_feature_with_p_0.05,cgFusoLowvsHigh_metacyc_daa_annotated_results_df, by = "feature")
cgFusoLowvsHigh_feature_desc$feature = cgFusoLowvsHigh_feature_desc$description
cgFusoLowvsHigh_feature_desc = cgFusoLowvsHigh_feature_desc[,c(1:7)]
colnames(cgFusoLowvsHigh_feature_desc) = colnames(cgFusoLowvsHigh_feature_with_p_0.05)

#Changing the pathway column to description for the abundance table
cgFusoLowvsHigh_abundance = cgFusoLowvsHigh_abundance_data_filtered %>% filter(`#OTU ID` %in% cgFusoLowvsHigh_feature_with_p_0.05$feature)
colnames(cgFusoLowvsHigh_abundance)[1] = "feature"
cgFusoLowvsHigh_abundance_desc = inner_join(cgFusoLowvsHigh_abundance,cgFusoLowvsHigh_metacyc_daa_annotated_results_df, by = "feature")
cgFusoLowvsHigh_abundance_desc$feature = cgFusoLowvsHigh_abundance_desc$description

#this line will change for each dataset. 34 represents the number of samples in the filtered abundance table
cgFusoLowvsHigh_abundance_desc = cgFusoLowvsHigh_abundance_desc[,-c(38:ncol(cgFusoLowvsHigh_abundance_desc))] 

# Generate a heatmap
cgFusoLowvsHigh_pathway_heatmap <- pathway_heatmap(abundance = cgFusoLowvsHigh_abundance_desc %>% column_to_rownames("feature"), metadata = cgFusoLowvsHigh_metadata, group = "Stage_FusoAbundance")

# Generate pathway PCA plot
cgFusoLowvsHigh_pathway_pca <- pathway_pca(abundance = cgFusoLowvsHigh_abundance_data_filtered %>% column_to_rownames("#OTU ID"), metadata = cgFusoLowvsHigh_metadata, group = "Stage_FusoAbundance")

# Generating a bar plot representing log2FC from the custom deseq2 function

# Go to the Deseq2 function script and update the metadata category of interest

# Lead the function in (already done so commented out)
# source("Fuso_DESeq2_function.R")

# Run the function on your own data
cgFusoLowvsHigh_res =  DEseq2_function(cgFusoLowvsHigh_abundance_data_filtered,cgFusoLowvsHigh_metadata, "Stage_FusoAbundance")
cgFusoLowvsHigh_res$feature =rownames(cgFusoLowvsHigh_res)
cgFusoLowvsHigh_res_desc = inner_join(cgFusoLowvsHigh_res,cgFusoLowvsHigh_metacyc_daa_annotated_results_df, by = "feature")
cgFusoLowvsHigh_res_desc = cgFusoLowvsHigh_res_desc[, -c(8:13)]
View(cgFusoLowvsHigh_res_desc)

# Filter to only include significant pathways
cgFusoLowvsHigh_sig_res = cgFusoLowvsHigh_res_desc %>%
  filter(pvalue < 0.05)
# You can also filter by Log2fold change (like include only +-2 log fold change)
#cgFusoLowvsHigh_sig_res <- sig_res[order(cgFusoLowvsHigh_sig_res$log2FoldChange),]

CG_FilteredFusoLowVSHigh_DESeq2 <- ggplot(data = cgFusoLowvsHigh_sig_res, aes(y = reorder(description, sort(as.numeric(log2FoldChange))), x= log2FoldChange, fill = pvalue))+
  geom_bar(stat = "identity")+ 
  theme_bw()+
  labs(x = "Log2FoldChange", y="Pathways")

ggsave("CG_FilteredFusoLowVSHigh_DESeq2.png", CG_FilteredFusoLowVSHigh_DESeq2)


#######Comparing Intestinal Metaplasia Fuso Low vs High#######

#Filtering the abundance table to only include samples that are in the filtered metadata
imFusoLowvsHigh_sample_names = imFusoLowvsHigh_metadata$'sample-id'
imFusoLowvsHigh_sample_names = append(imFusoLowvsHigh_sample_names, "#OTU ID")
imFusoLowvsHigh_abundance_data_filtered = abundance_data[, colnames(abundance_data) %in% imFusoLowvsHigh_sample_names] #This step is the actual filtering

#Removing individuals with no data that caused a problem for pathways_daa()
imFusoLowvsHigh_abundance_data_filtered =  imFusoLowvsHigh_abundance_data_filtered[, colSums(imFusoLowvsHigh_abundance_data_filtered != 0) > 0]

#Ensuring the rownames for the abundance_data_filtered is empty. This is required for their functions to run.
rownames(imFusoLowvsHigh_abundance_data_filtered) = NULL

#verify samples in metadata match samples in abundance_data
imFusoLowvsHigh_abun_samples = rownames(t(imFusoLowvsHigh_abundance_data_filtered[,-1])) #Getting a list of the sample names in the newly filtered abundance data
imFusoLowvsHigh_metadata = imFusoLowvsHigh_metadata[imFusoLowvsHigh_metadata$`sample-id` %in% imFusoLowvsHigh_abun_samples,] #making sure the filtered metadata only includes these samples

#### DESEq Intestinal Metaplasia Fuso LOW vs HIGH ####
#Perform pathway DAA using DESEQ2 method
imFusoLowvsHigh_abundance_daa_results_df <- pathway_daa(abundance = imFusoLowvsHigh_abundance_data_filtered %>% column_to_rownames("#OTU ID"), 
                                                        metadata = imFusoLowvsHigh_metadata, group = "Stage_FusoAbundance", daa_method = "DESeq2")

# Annotate MetaCyc pathway so they are more descriptive
imFusoLowvsHigh_metacyc_daa_annotated_results_df <- pathway_annotation(pathway = "MetaCyc", 
                                                                       daa_results_df = imFusoLowvsHigh_abundance_daa_results_df, ko_to_kegg = FALSE)

# Filter p-values to only significant ones
imFusoLowvsHigh_feature_with_p_0.05 <- imFusoLowvsHigh_abundance_daa_results_df %>% filter(p_values < 0.05)

#Changing the pathway column to description for the results 
imFusoLowvsHigh_feature_desc = inner_join(imFusoLowvsHigh_feature_with_p_0.05,imFusoLowvsHigh_metacyc_daa_annotated_results_df, by = "feature")
imFusoLowvsHigh_feature_desc$feature = imFusoLowvsHigh_feature_desc$description
imFusoLowvsHigh_feature_desc = imFusoLowvsHigh_feature_desc[,c(1:7)]
colnames(imFusoLowvsHigh_feature_desc) = colnames(imFusoLowvsHigh_feature_with_p_0.05)

#Changing the pathway column to description for the abundance table
imFusoLowvsHigh_abundance = imFusoLowvsHigh_abundance_data_filtered %>% filter(`#OTU ID` %in% imFusoLowvsHigh_feature_with_p_0.05$feature)
colnames(imFusoLowvsHigh_abundance)[1] = "feature"
imFusoLowvsHigh_abundance_desc = inner_join(imFusoLowvsHigh_abundance,imFusoLowvsHigh_metacyc_daa_annotated_results_df, by = "feature")
imFusoLowvsHigh_abundance_desc$feature = imFusoLowvsHigh_abundance_desc$description

#this line will change for each dataset. 34 represents the number of samples in the filtered abundance table
imFusoLowvsHigh_abundance_desc = imFusoLowvsHigh_abundance_desc[,-c(54:ncol(imFusoLowvsHigh_abundance_desc))] 

# Generate a heatmap
imFusoLowvsHigh_pathway_heatmap <- pathway_heatmap(abundance = imFusoLowvsHigh_abundance_desc %>% column_to_rownames("feature"), metadata = imFusoLowvsHigh_metadata, group = "Stage_FusoAbundance")

# Generate pathway PCA plot
imFusoLowvsHigh_pathway_pca <- pathway_pca(abundance = imFusoLowvsHigh_abundance_data_filtered %>% column_to_rownames("#OTU ID"), metadata = imFusoLowvsHigh_metadata, group = "Stage_FusoAbundance")

# Generating a bar plot representing log2FC from the custom deseq2 function

# Go to the Deseq2 function script and update the metadata category of interest

# Lead the function in (already done so commented out)
#source("Fuso_DESeq2_function.R")

# Run the function on your own data
imFusoLowvsHigh_res =  DEseq2_function(imFusoLowvsHigh_abundance_data_filtered,imFusoLowvsHigh_metadata, "Stage_FusoAbundance")
imFusoLowvsHigh_res$feature =rownames(imFusoLowvsHigh_res)
imFusoLowvsHigh_res_desc = inner_join(imFusoLowvsHigh_res,imFusoLowvsHigh_metacyc_daa_annotated_results_df, by = "feature")
imFusoLowvsHigh_res_desc = imFusoLowvsHigh_res_desc[, -c(8:13)]
View(imFusoLowvsHigh_res_desc)

# Filter to only include significant pathways
imFusoLowvsHigh_sig_res = imFusoLowvsHigh_res_desc %>%
  filter(pvalue < 0.05)
# You can also filter by Log2fold change (like include only +-2 log fold change)
#imFusoLowvsHigh_sig_res <- sig_res[order(imFusoLowvsHigh_sig_res$log2FoldChange),]

IM_FilteredFusoLowVSHigh_DESeq2 <- ggplot(data = imFusoLowvsHigh_sig_res, aes(y = reorder(description, sort(as.numeric(log2FoldChange))), x= log2FoldChange, fill = pvalue))+
  geom_bar(stat = "identity")+ 
  theme_bw()+
  labs(x = "Log2FoldChange", y="Pathways")


ggsave("IM_FilteredFusoLowVSHigh_DESeq2.png", IM_FilteredFusoLowVSHigh_DESeq2)


#######Comparing Intestinal Neoplasia Fuso Low vs High#######

#Filtering the abundance table to only include samples that are in the filtered metadata
ineoFusoLowvsHigh_sample_names = ineoFusoLowvsHigh_metadata$'sample-id'
ineoFusoLowvsHigh_sample_names = append(ineoFusoLowvsHigh_sample_names, "#OTU ID")
ineoFusoLowvsHigh_abundance_data_filtered = abundance_data[, colnames(abundance_data) %in% ineoFusoLowvsHigh_sample_names] #This step is the actual filtering

#Removing individuals with no data that caused a problem for pathways_daa()
ineoFusoLowvsHigh_abundance_data_filtered =  ineoFusoLowvsHigh_abundance_data_filtered[, colSums(ineoFusoLowvsHigh_abundance_data_filtered != 0) > 0]

#Ensuring the rownames for the abundance_data_filtered is empty. This is required for their functions to run.
rownames(ineoFusoLowvsHigh_abundance_data_filtered) = NULL

#verify samples in metadata match samples in abundance_data
ineoFusoLowvsHigh_abun_samples = rownames(t(ineoFusoLowvsHigh_abundance_data_filtered[,-1])) #Getting a list of the sample names in the newly filtered abundance data
ineoFusoLowvsHigh_metadata = ineoFusoLowvsHigh_metadata[ineoFusoLowvsHigh_metadata$`sample-id` %in% ineoFusoLowvsHigh_abun_samples,] #making sure the filtered metadata only includes these samples

#### DESEq Intestinal Metaplasia Fuso LOW vs HIGH ####
#Perform pathway DAA using DESEQ2 method
ineoFusoLowvsHigh_abundance_daa_results_df <- pathway_daa(abundance = ineoFusoLowvsHigh_abundance_data_filtered %>% column_to_rownames("#OTU ID"), 
                                                          metadata = ineoFusoLowvsHigh_metadata, group = "Stage_FusoAbundance", daa_method = "DESeq2")

# Annotate MetaCyc pathway so they are more descriptive
ineoFusoLowvsHigh_metacyc_daa_annotated_results_df <- pathway_annotation(pathway = "MetaCyc", 
                                                                         daa_results_df = ineoFusoLowvsHigh_abundance_daa_results_df, ko_to_kegg = FALSE)

# Filter p-values to only significant ones
ineoFusoLowvsHigh_feature_with_p_0.05 <- ineoFusoLowvsHigh_abundance_daa_results_df %>% filter(p_values < 0.05)

#Changing the pathway column to description for the results 
ineoFusoLowvsHigh_feature_desc = inner_join(ineoFusoLowvsHigh_feature_with_p_0.05,ineoFusoLowvsHigh_metacyc_daa_annotated_results_df, by = "feature")
ineoFusoLowvsHigh_feature_desc$feature = ineoFusoLowvsHigh_feature_desc$description
ineoFusoLowvsHigh_feature_desc = ineoFusoLowvsHigh_feature_desc[,c(1:7)]
colnames(ineoFusoLowvsHigh_feature_desc) = colnames(ineoFusoLowvsHigh_feature_with_p_0.05)

#Changing the pathway column to description for the abundance table
ineoFusoLowvsHigh_abundance = ineoFusoLowvsHigh_abundance_data_filtered %>% filter(`#OTU ID` %in% ineoFusoLowvsHigh_feature_with_p_0.05$feature)
colnames(ineoFusoLowvsHigh_abundance)[1] = "feature"
ineoFusoLowvsHigh_abundance_desc = inner_join(ineoFusoLowvsHigh_abundance,ineoFusoLowvsHigh_metacyc_daa_annotated_results_df, by = "feature")
ineoFusoLowvsHigh_abundance_desc$feature = ineoFusoLowvsHigh_abundance_desc$description

#this line will change for each dataset. 34 represents the number of samples in the filtered abundance table
ineoFusoLowvsHigh_abundance_desc = ineoFusoLowvsHigh_abundance_desc[,-c(76:ncol(imFusoLowvsHigh_abundance_desc))] 

# Generate a heatmap
ineoFusoLowvsHigh_pathway_heatmap <- pathway_heatmap(abundance = ineoFusoLowvsHigh_abundance_desc %>% column_to_rownames("feature"), metadata = ineoFusoLowvsHigh_metadata, group = "Stage_FusoAbundance")

# Generate pathway PCA plot
ineoFusoLowvsHigh_pathway_pca <- pathway_pca(abundance = ineoFusoLowvsHigh_abundance_data_filtered %>% column_to_rownames("#OTU ID"), metadata = ineoFusoLowvsHigh_metadata, group = "Stage_FusoAbundance")

# Generating a bar plot representing log2FC from the custom deseq2 function

# Go to the Deseq2 function script and update the metadata category of interest

# Lead the function in (already done so commented out)
#source("DESeq2_function.R")

# Run the function on your own data
ineoFusoLowvsHigh_res =  DEseq2_function(ineoFusoLowvsHigh_abundance_data_filtered,ineoFusoLowvsHigh_metadata, "Stage_FusoAbundance")
ineoFusoLowvsHigh_res$feature =rownames(ineoFusoLowvsHigh_res)
ineoFusoLowvsHigh_res_desc = inner_join(ineoFusoLowvsHigh_res,ineoFusoLowvsHigh_metacyc_daa_annotated_results_df, by = "feature")
ineoFusoLowvsHigh_res_desc = ineoFusoLowvsHigh_res_desc[, -c(8:13)]
View(ineoFusoLowvsHigh_res_desc)

# Filter to only include significant pathways
ineoFusoLowvsHigh_sig_res = ineoFusoLowvsHigh_res_desc %>%
  filter(pvalue < 0.05)
# You can also filter by Log2fold change (like include only +-2 log fold change)
# ineoFusoLowvsHigh_sig_res <- sig_res[order(ineoFusoLowvsHigh_sig_res$log2FoldChange),]

INeo_FilteredFusoLowVSHigh_DESeq2 <- ggplot(data = ineoFusoLowvsHigh_sig_res, aes(y = reorder(description, sort(as.numeric(log2FoldChange))), x= log2FoldChange, fill = pvalue))+
  geom_bar(stat = "identity")+ 
  theme_bw()+
  labs(x = "Log2FoldChange", y="Pathways")

ggsave("INeo_FilteredFusoLowVSHigh_DESeq2.png", INeo_FilteredFusoLowVSHigh_DESeq2)


#######Comparing Gastric Cancer Fuso Low vs High#######

#Filtering the abundance table to only include samples that are in the filtered metadata
cancerFusoLowvsHigh_sample_names = cancerFusoLowvsHigh_metadata$'sample-id'
cancerFusoLowvsHigh_sample_names = append(cancerFusoLowvsHigh_sample_names, "#OTU ID")
cancerFusoLowvsHigh_abundance_data_filtered = abundance_data[, colnames(abundance_data) %in% cancerFusoLowvsHigh_sample_names] #This step is the actual filtering

#Removing individuals with no data that caused a problem for pathways_daa()
cancerFusoLowvsHigh_abundance_data_filtered =  cancerFusoLowvsHigh_abundance_data_filtered[, colSums(cancerFusoLowvsHigh_abundance_data_filtered != 0) > 0]

#Ensuring the rownames for the abundance_data_filtered is empty. This is required for their functions to run.
rownames(cancerFusoLowvsHigh_abundance_data_filtered) = NULL

#verify samples in metadata match samples in abundance_data
cancerFusoLowvsHigh_abun_samples = rownames(t(cancerFusoLowvsHigh_abundance_data_filtered[,-1])) #Getting a list of the sample names in the newly filtered abundance data
cancerFusoLowvsHigh_metadata = cancerFusoLowvsHigh_metadata[cancerFusoLowvsHigh_metadata$`sample-id` %in% cancerFusoLowvsHigh_abun_samples,] #making sure the filtered metadata only includes these samples

#### DESEq Gastric Cancer Fuso LOW vs HIGH ####
#Perform pathway DAA using DESEQ2 method
cancerFusoLowvsHigh_abundance_daa_results_df <- pathway_daa(abundance = cancerFusoLowvsHigh_abundance_data_filtered %>% column_to_rownames("#OTU ID"), 
                                                            metadata = cancerFusoLowvsHigh_metadata, group = "Stage_FusoAbundance", daa_method = "DESeq2")

# Annotate MetaCyc pathway so they are more descriptive
cancerFusoLowvsHigh_metacyc_daa_annotated_results_df <- pathway_annotation(pathway = "MetaCyc", 
                                                                           daa_results_df = cancerFusoLowvsHigh_abundance_daa_results_df, ko_to_kegg = FALSE)

# Filter p-values to only significant ones
cancerFusoLowvsHigh_feature_with_p_0.05 <- cancerFusoLowvsHigh_abundance_daa_results_df %>% filter(p_values < 0.05)

#Changing the pathway column to description for the results 
cancerFusoLowvsHigh_feature_desc = inner_join(cancerFusoLowvsHigh_feature_with_p_0.05,cancerFusoLowvsHigh_metacyc_daa_annotated_results_df, by = "feature")
cancerFusoLowvsHigh_feature_desc$feature = cancerFusoLowvsHigh_feature_desc$description
cancerFusoLowvsHigh_feature_desc = cancerFusoLowvsHigh_feature_desc[,c(1:7)]
colnames(cancerFusoLowvsHigh_feature_desc) = colnames(cancerFusoLowvsHigh_feature_with_p_0.05)

#Changing the pathway column to description for the abundance table
cancerFusoLowvsHigh_abundance = cancerFusoLowvsHigh_abundance_data_filtered %>% filter(`#OTU ID` %in% cancerFusoLowvsHigh_feature_with_p_0.05$feature)
colnames(cancerFusoLowvsHigh_abundance)[1] = "feature"
cancerFusoLowvsHigh_abundance_desc = inner_join(cancerFusoLowvsHigh_abundance,cancerFusoLowvsHigh_metacyc_daa_annotated_results_df, by = "feature")
cancerFusoLowvsHigh_abundance_desc$feature = cancerFusoLowvsHigh_abundance_desc$description

#this line will change for each dataset. 34 represents the number of samples in the filtered abundance table
cancerFusoLowvsHigh_abundance_desc = cancerFusoLowvsHigh_abundance_desc[,-c(85:ncol(cancerFusoLowvsHigh_abundance_desc))] 

# Generate a heatmap
cancerFusoLowvsHigh_pathway_heatmap <- pathway_heatmap(abundance = cancerFusoLowvsHigh_abundance_desc %>% column_to_rownames("feature"), metadata = cancerFusoLowvsHigh_metadata, group = "Stage_FusoAbundance")

# Generate pathway PCA plot
cancerFusoLowvsHigh_pathway_pca <- pathway_pca(abundance = cancerFusoLowvsHigh_abundance_data_filtered %>% column_to_rownames("#OTU ID"), metadata = cancerFusoLowvsHigh_metadata, group = "Stage_FusoAbundance")

# Generating a bar plot representing log2FC from the custom deseq2 function

# Go to the Deseq2 function script and update the metadata category of interest

# Lead the function in (already done so commented out)
# source("Fuso_DESeq2_function.R")

# Run the function on your own data
cancerFusoLowvsHigh_res =  DEseq2_function(cancerFusoLowvsHigh_abundance_data_filtered,cancerFusoLowvsHigh_metadata, "Stage_FusoAbundance")
cancerFusoLowvsHigh_res$feature =rownames(cancerFusoLowvsHigh_res)
cancerFusoLowvsHigh_res_desc = inner_join(cancerFusoLowvsHigh_res,cancerFusoLowvsHigh_metacyc_daa_annotated_results_df, by = "feature")
cancerFusoLowvsHigh_res_desc = cancerFusoLowvsHigh_res_desc[, -c(8:13)]
View(cancerFusoLowvsHigh_res_desc)

# Filter to only include significant pathways
cancerFusoLowvsHigh_sig_res = cancerFusoLowvsHigh_res_desc %>%
  filter(pvalue < 0.05)
# You can also filter by Log2fold change (like include only +-2 log fold change)
#cancerFusoLowvsHigh_sig_res <- sig_res[order(cancerFusoLowvsHigh_sig_res$log2FoldChange),]

cancer_FilteredFusoLowVSHigh_DESeq2 <- ggplot(data = cancerFusoLowvsHigh_sig_res, aes(y = reorder(description, sort(as.numeric(log2FoldChange))), x= log2FoldChange, fill = pvalue))+
  geom_bar(stat = "identity")+ 
  theme_bw()+
  labs(x = "Log2FoldChange", y="Pathways")

ggsave("cancer_FilteredFusoLowVSHigh_DESeq2.png", cancer_FilteredFusoLowVSHigh_DESeq2)

