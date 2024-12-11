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

install.packages("ggokabeito")

install.packages("viridis")

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
library(ggplot2)
library(tidyr)
library(RColorBrewer)
library(ggokabeito)
library(viridis)


#### Import files and preparing tables ####
#Importing the pathway PICrsut2
abundance_file <- "pathway_abundance.tsv"
abundance_data <- read_delim(abundance_file, delim = "\t", skip = 1, col_names = TRUE, trim_ws = TRUE)
abundance_data  =as.data.frame(abundance_data)

#Import your metadata file, no need to filter yet
metadata <- as.data.frame(read_delim("LH_gastric_cancer_metadata.tsv"))

#Remove NAs for your column of interest in this case subject
metadata = metadata[!is.na(metadata$Stage_Pylori),]

#Filter your metadata to include only 2 at a time
cgPyloriPlusvsMinus_metadata <- filter(metadata, Stage_Pylori %in% c("Chronic gastritis (CG)_Negative", "Chronic gastritis (CG)_Positive"))
imPyloriPlusvsMinus_metadata <- filter(metadata, Stage_Pylori %in% c("Intestinal metaplasia (IM）_Negative", "Intestinal metaplasia (IM）_Positive"))
ineoPyloriPlusvsMinus_metadata <- filter(metadata, Stage_Pylori %in% c("Intraepithelial neoplasia (IN)_Negative", "Intraepithelial neoplasia (IN)_Positive"))
cancerPyloriPlusvsMinus_metadata <-filter(metadata, Stage_Pylori %in% c("Gastric cancer (GC)_Negative", "Gastric cancer (GC)_Positive"))


#######Comparing Chronic Gastritis H. Pylori Yes vs No #######

#Filtering the abundance table to only include samples that are in the filtered metadata
cgPyloriPlusvsMinus_sample_names = cgPyloriPlusvsMinus_metadata$'sample-id'
cgPyloriPlusvsMinus_sample_names = append(cgPyloriPlusvsMinus_sample_names, "#OTU ID")
cgPyloriPlusvsMinus_abundance_data_filtered = abundance_data[, colnames(abundance_data) %in% cgPyloriPlusvsMinus_sample_names] #This step is the actual filtering

#Removing individuals with no data that caused a problem for pathways_daa()
cgPyloriPlusvsMinus_abundance_data_filtered =  cgPyloriPlusvsMinus_abundance_data_filtered[, colSums(cgPyloriPlusvsMinus_abundance_data_filtered != 0) > 0]

#Ensuring the rownames for the abundance_data_filtered is empty. This is required for their functions to run.
rownames(cgPyloriPlusvsMinus_abundance_data_filtered) = NULL

#verify samples in metadata match samples in abundance_data
cgPyloriPlusvsMinus_abun_samples = rownames(t(cgPyloriPlusvsMinus_abundance_data_filtered[,-1])) #Getting a list of the sample names in the newly filtered abundance data
cgPyloriPlusvsMinus_metadata = cgPyloriPlusvsMinus_metadata[cgPyloriPlusvsMinus_metadata$`sample-id` %in% cgPyloriPlusvsMinus_abun_samples,] #making sure the filtered metadata only includes these samples

#### DESEq Chronic Gastritis H. Pylori Yes vs No ####
#Perform pathway DAA using DESEQ2 method
cgPyloriPlusvsMinus_abundance_daa_results_df <- pathway_daa(abundance = cgPyloriPlusvsMinus_abundance_data_filtered %>% column_to_rownames("#OTU ID"), 
                                                        metadata = cgPyloriPlusvsMinus_metadata, group = "Stage_Pylori", daa_method = "DESeq2")

# Annotate MetaCyc pathway so they are more descriptive
cgPyloriPlusvsMinus_metacyc_daa_annotated_results_df <- pathway_annotation(pathway = "MetaCyc", 
                                                                       daa_results_df = cgPyloriPlusvsMinus_abundance_daa_results_df, ko_to_kegg = FALSE)

# Filter p-values to only significant ones
cgPyloriPlusvsMinus_feature_with_p_0.05 <- cgPyloriPlusvsMinus_abundance_daa_results_df %>% filter(p_values < 0.05)

#Changing the pathway column to description for the results 
cgPyloriPlusvsMinus_feature_desc = inner_join(cgPyloriPlusvsMinus_feature_with_p_0.05,cgPyloriPlusvsMinus_metacyc_daa_annotated_results_df, by = "feature")
cgPyloriPlusvsMinus_feature_desc$feature = cgPyloriPlusvsMinus_feature_desc$description
cgPyloriPlusvsMinus_feature_desc = cgPyloriPlusvsMinus_feature_desc[,c(1:7)]
colnames(cgPyloriPlusvsMinus_feature_desc) = colnames(cgPyloriPlusvsMinus_feature_with_p_0.05)

#Changing the pathway column to description for the abundance table
cgPyloriPlusvsMinus_abundance = cgPyloriPlusvsMinus_abundance_data_filtered %>% filter(`#OTU ID` %in% cgPyloriPlusvsMinus_feature_with_p_0.05$feature)
colnames(cgPyloriPlusvsMinus_abundance)[1] = "feature"
cgPyloriPlusvsMinus_abundance_desc = inner_join(cgPyloriPlusvsMinus_abundance,cgPyloriPlusvsMinus_metacyc_daa_annotated_results_df, by = "feature")
cgPyloriPlusvsMinus_abundance_desc$feature = cgPyloriPlusvsMinus_abundance_desc$description

#this line will change for each dataset. 34 represents the number of samples in the filtered abundance table
cgPyloriPlusvsMinus_abundance_desc = cgPyloriPlusvsMinus_abundance_desc[,-c(39:ncol(cgPyloriPlusvsMinus_abundance_desc))] 

# Generate a heatmap
cgPyloriPlusvsMinus_pathway_heatmap <- pathway_heatmap(abundance = cgPyloriPlusvsMinus_abundance_desc %>% column_to_rownames("feature"), metadata = cgPyloriPlusvsMinus_metadata, group = "Stage_Pylori")

# Generate pathway PCA plot
cgPyloriPlusvsMinus_pathway_pca <- pathway_pca(abundance = cgPyloriPlusvsMinus_abundance_data_filtered %>% column_to_rownames("#OTU ID"), metadata = cgPyloriPlusvsMinus_metadata, group = "Stage_Pylori")

# Generating a bar plot representing log2FC from the custom deseq2 function

# Go to the Deseq2 function script and update the metadata category of interest

# Lead the function in
source("H_Pylori_DESeq2_function.R")

# Run the function on your own data
cgPyloriPlusvsMinus_res =  DEseq2_function(cgPyloriPlusvsMinus_abundance_data_filtered,cgPyloriPlusvsMinus_metadata, "Stage_Pylori")
cgPyloriPlusvsMinus_res$feature =rownames(cgPyloriPlusvsMinus_res)
cgPyloriPlusvsMinus_res_desc = inner_join(cgPyloriPlusvsMinus_res,cgPyloriPlusvsMinus_metacyc_daa_annotated_results_df, by = "feature")
cgPyloriPlusvsMinus_res_desc = cgPyloriPlusvsMinus_res_desc[, -c(8:13)]
# View(cgPyloriPlusvsMinus_res_desc)

# Filter to only include significant pathways
cgPyloriPlusvsMinus_sig_res = cgPyloriPlusvsMinus_res_desc %>%
  filter(pvalue < 0.05)

# You can also filter by Log2fold change (like include only +-2 log fold change)
  #cgPyloriPlusvsMinus_sig_res <- sig_res[order(cgPyloriPlusvsMinus_sig_res$log2FoldChange),]

cg_Filtered_HPylori_YesVSNo_DESeq2 <- ggplot(data = cgPyloriPlusvsMinus_sig_res, aes(y = reorder(description, sort(as.numeric(log2FoldChange))), x= log2FoldChange, fill = pvalue))+
  geom_bar(stat = "identity")+ 
  theme_bw()+
  labs(x = "Log2FoldChange", y="Pathways")

ggsave("cg_Filtered_HPylori_YesVSNo_DESeq2.png", cg_Filtered_HPylori_YesVSNo_DESeq2)


cg_Ordered_Filtered_HPylori_YesVSNo_DESeq2 <- ggplot(data = cgPyloriPlusvsMinus_sig_res, aes(y = reorder(description, -log2FoldChange), x= log2FoldChange, fill = pvalue))+
  geom_bar(stat = "identity") +
  theme_test()+
  theme(plot.margin = margin(t = 10, r = 10, b = 30, l = 10), 
        axis.text.y = element_text(size = 12), 
        axis.text.x = element_text(size = 12), 
        axis.title.x = element_text(size = 16, vjust = -5), 
        axis.title.y = element_text(size = 16), 
        legend.title = element_text(size = 16)) +
  scale_fill_viridis() +
  labs(x = "Log2FoldChange", y="Pathways", fill = "p value")

ggsave("cg_Ordered_Filtered_HPylori_YesVSNo_DESeq2.png", cg_Ordered_Filtered_HPylori_YesVSNo_DESeq2, width = 13, height = 8.5)




cg_summary_table <- cgPyloriPlusvsMinus_sig_res %>%
    summarize(
    Positive_Log2FC = sum(log2FoldChange > 0),  # Count positive values
    Negative_Log2FC = sum(log2FoldChange < 0)   # Count negative values
  ) %>%
  mutate(Group = "Chronic Gastritis")


#######Comparing Intestinal Metaplasia H. Pylori Yes vs No #######

#Filtering the abundance table to only include samples that are in the filtered metadata
imPyloriPlusvsMinus_sample_names = imPyloriPlusvsMinus_metadata$'sample-id'
imPyloriPlusvsMinus_sample_names = append(imPyloriPlusvsMinus_sample_names, "#OTU ID")
imPyloriPlusvsMinus_abundance_data_filtered = abundance_data[, colnames(abundance_data) %in% imPyloriPlusvsMinus_sample_names] #This step is the actual filtering

#Removing individuals with no data that caused a problem for pathways_daa()
imPyloriPlusvsMinus_abundance_data_filtered =  imPyloriPlusvsMinus_abundance_data_filtered[, colSums(imPyloriPlusvsMinus_abundance_data_filtered != 0) > 0]

#Ensuring the rownames for the abundance_data_filtered is empty. This is required for their functions to run.
rownames(imPyloriPlusvsMinus_abundance_data_filtered) = NULL

#verify samples in metadata match samples in abundance_data
imPyloriPlusvsMinus_abun_samples = rownames(t(imPyloriPlusvsMinus_abundance_data_filtered[,-1])) #Getting a list of the sample names in the newly filtered abundance data
imPyloriPlusvsMinus_metadata = imPyloriPlusvsMinus_metadata[imPyloriPlusvsMinus_metadata$`sample-id` %in% imPyloriPlusvsMinus_abun_samples,] #making sure the filtered metadata only includes these samples

#### DESEq Intestinal Metaplasia H. Pylori Yes vs No ####
#Perform pathway DAA using DESEQ2 method
imPyloriPlusvsMinus_abundance_daa_results_df <- pathway_daa(abundance = imPyloriPlusvsMinus_abundance_data_filtered %>% column_to_rownames("#OTU ID"), 
                                                            metadata = imPyloriPlusvsMinus_metadata, group = "Stage_Pylori", daa_method = "DESeq2")

# Annotate MetaCyc pathway so they are more descriptive
imPyloriPlusvsMinus_metacyc_daa_annotated_results_df <- pathway_annotation(pathway = "MetaCyc", 
                                                                           daa_results_df = imPyloriPlusvsMinus_abundance_daa_results_df, ko_to_kegg = FALSE)

# Filter p-values to only significant ones
imPyloriPlusvsMinus_feature_with_p_0.05 <- imPyloriPlusvsMinus_abundance_daa_results_df %>% filter(p_values < 0.05)

#Changing the pathway column to description for the results 
imPyloriPlusvsMinus_feature_desc = inner_join(imPyloriPlusvsMinus_feature_with_p_0.05,imPyloriPlusvsMinus_metacyc_daa_annotated_results_df, by = "feature")
imPyloriPlusvsMinus_feature_desc$feature = imPyloriPlusvsMinus_feature_desc$description
imPyloriPlusvsMinus_feature_desc = imPyloriPlusvsMinus_feature_desc[,c(1:7)]
colnames(imPyloriPlusvsMinus_feature_desc) = colnames(imPyloriPlusvsMinus_feature_with_p_0.05)

#Changing the pathway column to description for the abundance table
imPyloriPlusvsMinus_abundance = imPyloriPlusvsMinus_abundance_data_filtered %>% filter(`#OTU ID` %in% imPyloriPlusvsMinus_feature_with_p_0.05$feature)
colnames(imPyloriPlusvsMinus_abundance)[1] = "feature"
imPyloriPlusvsMinus_abundance_desc = inner_join(imPyloriPlusvsMinus_abundance,imPyloriPlusvsMinus_metacyc_daa_annotated_results_df, by = "feature")
imPyloriPlusvsMinus_abundance_desc$feature = imPyloriPlusvsMinus_abundance_desc$description

#this line will change for each dataset. 34 represents the number of samples in the filtered abundance table
imPyloriPlusvsMinus_abundance_desc = imPyloriPlusvsMinus_abundance_desc[,-c(55:ncol(imPyloriPlusvsMinus_abundance_desc))] 

# Generate a heatmap
imPyloriPlusvsMinus_pathway_heatmap <- pathway_heatmap(abundance = imPyloriPlusvsMinus_abundance_desc %>% column_to_rownames("feature"), metadata = imPyloriPlusvsMinus_metadata, group = "Stage_Pylori")

# Generate pathway PCA plot
imPyloriPlusvsMinus_pathway_pca <- pathway_pca(abundance = imPyloriPlusvsMinus_abundance_data_filtered %>% column_to_rownames("#OTU ID"), metadata = imPyloriPlusvsMinus_metadata, group = "Stage_Pylori")

# Generating a bar plot representing log2FC from the custom deseq2 function

# Go to the Deseq2 function script and update the metadata category of interest

# Lead the function in
#source("H_Pylori_DESeq2_function.R")

# Run the function on your own data
imPyloriPlusvsMinus_res =  DEseq2_function(imPyloriPlusvsMinus_abundance_data_filtered,imPyloriPlusvsMinus_metadata, "Stage_Pylori")
imPyloriPlusvsMinus_res$feature =rownames(imPyloriPlusvsMinus_res)
imPyloriPlusvsMinus_res_desc = inner_join(imPyloriPlusvsMinus_res,imPyloriPlusvsMinus_metacyc_daa_annotated_results_df, by = "feature")
imPyloriPlusvsMinus_res_desc = imPyloriPlusvsMinus_res_desc[, -c(8:13)]
# View(imPyloriPlusvsMinus_res_desc)

# Filter to only include significant pathways
imPyloriPlusvsMinus_sig_res = imPyloriPlusvsMinus_res_desc %>%
  filter(pvalue < 0.05)
# You can also filter by Log2fold change (like include only +-2 log fold change)
  #imPyloriPlusvsMinus_sig_res <- sig_res[order(imPyloriPlusvsMinus_sig_res$log2FoldChange),]

im_Filtered_HPylori_YesVSNo_DESeq2 <- ggplot(data = imPyloriPlusvsMinus_sig_res, aes(y = reorder(description, sort(as.numeric(log2FoldChange))), x= log2FoldChange, fill = pvalue))+
  geom_bar(stat = "identity")+ 
  theme_bw()+
  labs(x = "Log2FoldChange", y="Pathways")

ggsave("im_Filtered_HPylori_YesVSNo_DESeq2.png", im_Filtered_HPylori_YesVSNo_DESeq2)


im_Ordered_Filtered_HPylori_YesVSNo_DESeq2 <- ggplot(data = imPyloriPlusvsMinus_sig_res, aes(y = reorder(description, -log2FoldChange), x= log2FoldChange, fill = pvalue))+
  geom_bar(stat = "identity") +
  theme_test()+
  theme(plot.margin = margin(t = 10, r = 10, b = 30, l = 10), 
        axis.text.y = element_text(size = 12), 
        axis.text.x = element_text(size = 12), 
        axis.title.x = element_text(size = 16, vjust = -5), 
        axis.title.y = element_text(size = 16), 
        legend.title = element_text(size = 16)) +
  scale_fill_viridis() +
  scale_y_discrete(labels = c(
    "superpathway of N-acetylglucosamine, N-acetylmannosamine and N-acetylneuraminate degradation" = "superpathway of N-acetylglucosamine, N-acetylmannosamine \nand N-acetylneuraminate degradation"
  )) +
  labs(x = "Log2FoldChange", y="Pathways", fill = "p value")

ggsave("im_Ordered_Filtered_HPylori_YesVSNo_DESeq2.png", im_Ordered_Filtered_HPylori_YesVSNo_DESeq2, width = 8.5, height = 8.5)


im_summary_table <- imPyloriPlusvsMinus_sig_res %>%
  summarize(
    Positive_Log2FC = sum(log2FoldChange > 0),  # Count positive values
    Negative_Log2FC = sum(log2FoldChange < 0)   # Count negative values
  )%>%
  mutate(Group = "Intestinal Metaplasia")



#######Comparing Intraepithelial Neoplasia H. Pylori Yes vs No #######

#Filtering the abundance table to only include samples that are in the filtered metadata
ineoPyloriPlusvsMinus_sample_names = ineoPyloriPlusvsMinus_metadata$'sample-id'
ineoPyloriPlusvsMinus_sample_names = append(ineoPyloriPlusvsMinus_sample_names, "#OTU ID")
ineoPyloriPlusvsMinus_abundance_data_filtered = abundance_data[, colnames(abundance_data) %in% ineoPyloriPlusvsMinus_sample_names] #This step is the actual filtering

#Removing individuals with no data that caused a problem for pathways_daa()
ineoPyloriPlusvsMinus_abundance_data_filtered =  ineoPyloriPlusvsMinus_abundance_data_filtered[, colSums(ineoPyloriPlusvsMinus_abundance_data_filtered != 0) > 0]

#Ensuring the rownames for the abundance_data_filtered is empty. This is required for their functions to run.
rownames(ineoPyloriPlusvsMinus_abundance_data_filtered) = NULL

#verify samples in metadata match samples in abundance_data
ineoPyloriPlusvsMinus_abun_samples = rownames(t(ineoPyloriPlusvsMinus_abundance_data_filtered[,-1])) #Getting a list of the sample names in the newly filtered abundance data
ineoPyloriPlusvsMinus_metadata = ineoPyloriPlusvsMinus_metadata[ineoPyloriPlusvsMinus_metadata$`sample-id` %in% ineoPyloriPlusvsMinus_abun_samples,] #making sure the filtered metadata only includes these samples

#### DESEq Intraepithelial Neoplasia H. Pylori Yes vs No ####
#Perform pathway DAA using DESEQ2 method
ineoPyloriPlusvsMinus_abundance_daa_results_df <- pathway_daa(abundance = ineoPyloriPlusvsMinus_abundance_data_filtered %>% column_to_rownames("#OTU ID"), 
                                                            metadata = ineoPyloriPlusvsMinus_metadata, group = "Stage_Pylori", daa_method = "DESeq2")

# Annotate MetaCyc pathway so they are more descriptive
ineoPyloriPlusvsMinus_metacyc_daa_annotated_results_df <- pathway_annotation(pathway = "MetaCyc", 
                                                                           daa_results_df = ineoPyloriPlusvsMinus_abundance_daa_results_df, ko_to_kegg = FALSE)

# Filter p-values to only significant ones
ineoPyloriPlusvsMinus_feature_with_p_0.05 <- ineoPyloriPlusvsMinus_abundance_daa_results_df %>% filter(p_values < 0.05)

#Changing the pathway column to description for the results 
ineoPyloriPlusvsMinus_feature_desc = inner_join(ineoPyloriPlusvsMinus_feature_with_p_0.05,ineoPyloriPlusvsMinus_metacyc_daa_annotated_results_df, by = "feature")
ineoPyloriPlusvsMinus_feature_desc$feature = ineoPyloriPlusvsMinus_feature_desc$description
ineoPyloriPlusvsMinus_feature_desc = ineoPyloriPlusvsMinus_feature_desc[,c(1:7)]
colnames(ineoPyloriPlusvsMinus_feature_desc) = colnames(ineoPyloriPlusvsMinus_feature_with_p_0.05)

#Changing the pathway column to description for the abundance table
ineoPyloriPlusvsMinus_abundance = ineoPyloriPlusvsMinus_abundance_data_filtered %>% filter(`#OTU ID` %in% ineoPyloriPlusvsMinus_feature_with_p_0.05$feature)
colnames(ineoPyloriPlusvsMinus_abundance)[1] = "feature"
ineoPyloriPlusvsMinus_abundance_desc = inner_join(ineoPyloriPlusvsMinus_abundance,ineoPyloriPlusvsMinus_metacyc_daa_annotated_results_df, by = "feature")
ineoPyloriPlusvsMinus_abundance_desc$feature = ineoPyloriPlusvsMinus_abundance_desc$description

#this line will change for each dataset. 34 represents the number of samples in the filtered abundance table
ineoPyloriPlusvsMinus_abundance_desc = ineoPyloriPlusvsMinus_abundance_desc[,-c(76:ncol(ineoPyloriPlusvsMinus_abundance_desc))]

# Generate a heatmap
ineoPyloriPlusvsMinus_pathway_heatmap <- pathway_heatmap(abundance = ineoPyloriPlusvsMinus_abundance_desc %>% column_to_rownames("feature"), metadata = ineoPyloriPlusvsMinus_metadata, group = "Stage_Pylori")

# Generate pathway PCA plot
ineoPyloriPlusvsMinus_pathway_pca <- pathway_pca(abundance = ineoPyloriPlusvsMinus_abundance_data_filtered %>% column_to_rownames("#OTU ID"), metadata = ineoPyloriPlusvsMinus_metadata, group = "Stage_Pylori")

# Generating a bar plot representing log2FC from the custom deseq2 function

# Go to the Deseq2 function script and update the metadata category of interest

# Lead the function in
# source("H_Pylori_DESeq2_function.R")

# Run the function on your own data
ineoPyloriPlusvsMinus_res =  DEseq2_function(ineoPyloriPlusvsMinus_abundance_data_filtered,ineoPyloriPlusvsMinus_metadata, "Stage_Pylori")
ineoPyloriPlusvsMinus_res$feature =rownames(ineoPyloriPlusvsMinus_res)
ineoPyloriPlusvsMinus_res_desc = inner_join(ineoPyloriPlusvsMinus_res,ineoPyloriPlusvsMinus_metacyc_daa_annotated_results_df, by = "feature")
ineoPyloriPlusvsMinus_res_desc = ineoPyloriPlusvsMinus_res_desc[, -c(8:13)]
# View(ineoPyloriPlusvsMinus_res_desc)

# Filter to only include significant pathways
ineoPyloriPlusvsMinus_sig_res = ineoPyloriPlusvsMinus_res_desc %>%
  filter(pvalue < 0.05)
# You can also filter by Log2fold change (like include only +-2 log fold change)
  #ineoPyloriPlusvsMinus_sig_res <- sig_res[order(ineoPyloriPlusvsMinus_sig_res$log2FoldChange),]

INeo_Filtered_HPylori_YesVSNo_DESeq2 <- ggplot(data = ineoPyloriPlusvsMinus_sig_res, aes(y = reorder(description, sort(as.numeric(log2FoldChange))), x= log2FoldChange, fill = pvalue))+
  geom_bar(stat = "identity")+ 
  theme_bw()+
  labs(x = "Log2FoldChange", y="Pathways")

ggsave("INeo_Filtered_HPylori_YesVSNo_DESeq2.png", INeo_Filtered_HPylori_YesVSNo_DESeq2)


INeo_Ordered_Filtered_HPylori_YesVSNo_DESeq2 <- ggplot(data = ineoPyloriPlusvsMinus_sig_res, aes(y = reorder(description, -log2FoldChange), x= log2FoldChange, fill = pvalue))+
  geom_bar(stat = "identity") +
  theme_test()+
  theme(plot.margin = margin(t = 10, r = 10, b = 30, l = 10), 
        axis.text.y = element_text(size = 12), 
        axis.text.x = element_text(size = 12), 
        axis.title.x = element_text(size = 16, vjust = -5), 
        axis.title.y = element_text(size = 16), 
        legend.title = element_text(size = 16)) +
  scale_fill_viridis() +
  labs(x = "Log2FoldChange", y="Pathways", fill = "p value")

ggsave("INeo_Ordered_Filtered_HPylori_YesVSNo_DESeq2.png", INeo_Ordered_Filtered_HPylori_YesVSNo_DESeq2, width = 13, height = 8.5)


ineo_summary_table <- ineoPyloriPlusvsMinus_sig_res %>%
  summarize(
    Positive_Log2FC = sum(log2FoldChange > 0),  # Count positive values
    Negative_Log2FC = sum(log2FoldChange < 0)   # Count negative values
  )%>%
  mutate(Group = "Intraepithelial Neoplasia")


#######Comparing gastric cancer H. Pylori Yes vs No #######

#Filtering the abundance table to only include samples that are in the filtered metadata
cancerPyloriPlusvsMinus_sample_names = cancerPyloriPlusvsMinus_metadata$'sample-id'
cancerPyloriPlusvsMinus_sample_names = append(cancerPyloriPlusvsMinus_sample_names, "#OTU ID")
cancerPyloriPlusvsMinus_abundance_data_filtered = abundance_data[, colnames(abundance_data) %in% cancerPyloriPlusvsMinus_sample_names] #This step is the actual filtering

#Removing individuals with no data that caused a problem for pathways_daa()
cancerPyloriPlusvsMinus_abundance_data_filtered =  cancerPyloriPlusvsMinus_abundance_data_filtered[, colSums(cancerPyloriPlusvsMinus_abundance_data_filtered != 0) > 0]

#Ensuring the rownames for the abundance_data_filtered is empty. This is required for their functions to run.
rownames(cancerPyloriPlusvsMinus_abundance_data_filtered) = NULL

#verify samples in metadata match samples in abundance_data
cancerPyloriPlusvsMinus_abun_samples = rownames(t(cancerPyloriPlusvsMinus_abundance_data_filtered[,-1])) #Getting a list of the sample names in the newly filtered abundance data
cancerPyloriPlusvsMinus_metadata = cancerPyloriPlusvsMinus_metadata[cancerPyloriPlusvsMinus_metadata$`sample-id` %in% cancerPyloriPlusvsMinus_abun_samples,] #making sure the filtered metadata only includes these samples

#### DESEq Intestinal Neoplasia H. Pylori Yes vs No ####
#Perform pathway DAA using DESEQ2 method
cancerPyloriPlusvsMinus_abundance_daa_results_df <- pathway_daa(abundance = cancerPyloriPlusvsMinus_abundance_data_filtered %>% column_to_rownames("#OTU ID"), 
                                                              metadata = cancerPyloriPlusvsMinus_metadata, group = "Stage_Pylori", daa_method = "DESeq2")

# Annotate MetaCyc pathway so they are more descriptive
cancerPyloriPlusvsMinus_metacyc_daa_annotated_results_df <- pathway_annotation(pathway = "MetaCyc", 
                                                                             daa_results_df = cancerPyloriPlusvsMinus_abundance_daa_results_df, ko_to_kegg = FALSE)

# Filter p-values to only significant ones
cancerPyloriPlusvsMinus_feature_with_p_0.05 <- cancerPyloriPlusvsMinus_abundance_daa_results_df %>% filter(p_values < 0.05)

#Changing the pathway column to description for the results 
cancerPyloriPlusvsMinus_feature_desc = inner_join(cancerPyloriPlusvsMinus_feature_with_p_0.05,cancerPyloriPlusvsMinus_metacyc_daa_annotated_results_df, by = "feature")
cancerPyloriPlusvsMinus_feature_desc$feature = cancerPyloriPlusvsMinus_feature_desc$description
cancerPyloriPlusvsMinus_feature_desc = cancerPyloriPlusvsMinus_feature_desc[,c(1:7)]
colnames(cancerPyloriPlusvsMinus_feature_desc) = colnames(cancerPyloriPlusvsMinus_feature_with_p_0.05)

#Changing the pathway column to description for the abundance table
cancerPyloriPlusvsMinus_abundance = cancerPyloriPlusvsMinus_abundance_data_filtered %>% filter(`#OTU ID` %in% cancerPyloriPlusvsMinus_feature_with_p_0.05$feature)
colnames(cancerPyloriPlusvsMinus_abundance)[1] = "feature"
cancerPyloriPlusvsMinus_abundance_desc = inner_join(cancerPyloriPlusvsMinus_abundance,cancerPyloriPlusvsMinus_metacyc_daa_annotated_results_df, by = "feature")
cancerPyloriPlusvsMinus_abundance_desc$feature = cancerPyloriPlusvsMinus_abundance_desc$description

#this line will change for each dataset. 34 represents the number of samples in the filtered abundance table
cancerPyloriPlusvsMinus_abundance_desc = cancerPyloriPlusvsMinus_abundance_desc[,-c(85:ncol(cancerPyloriPlusvsMinus_abundance_desc))]

# Generate a heatmap
cancerPyloriPlusvsMinus_pathway_heatmap <- pathway_heatmap(abundance = cancerPyloriPlusvsMinus_abundance_desc %>% column_to_rownames("feature"), metadata = cancerPyloriPlusvsMinus_metadata, group = "Stage_Pylori")

# Generate pathway PCA plot
cancerPyloriPlusvsMinus_pathway_pca <- pathway_pca(abundance = cancerPyloriPlusvsMinus_abundance_data_filtered %>% column_to_rownames("#OTU ID"), metadata = cancerPyloriPlusvsMinus_metadata, group = "Stage_Pylori")

# Generating a bar plot representing log2FC from the custom deseq2 function

# Go to the Deseq2 function script and update the metadata category of interest

# Lead the function in
# source("H_Pylori_DESeq2_function.R")

# Run the function on your own data
cancerPyloriPlusvsMinus_res =  DEseq2_function(cancerPyloriPlusvsMinus_abundance_data_filtered,cancerPyloriPlusvsMinus_metadata, "Stage_Pylori")
cancerPyloriPlusvsMinus_res$feature =rownames(cancerPyloriPlusvsMinus_res)
cancerPyloriPlusvsMinus_res_desc = inner_join(cancerPyloriPlusvsMinus_res,cancerPyloriPlusvsMinus_metacyc_daa_annotated_results_df, by = "feature")
cancerPyloriPlusvsMinus_res_desc = cancerPyloriPlusvsMinus_res_desc[, -c(8:13)]
# View(cancerPyloriPlusvsMinus_res_desc)

# Filter to only include significant pathways
cancerPyloriPlusvsMinus_sig_res = cancerPyloriPlusvsMinus_res_desc %>%
  filter(pvalue < 0.05)
# You can also filter by Log2fold change (like include only +-2 log fold change)
#ineoPyloriPlusvsMinus_sig_res <- sig_res[order(ineoPyloriPlusvsMinus_sig_res$log2FoldChange),]

cancer_Filtered_HPylori_YesVSNo_DESeq2 <- ggplot(data = cancerPyloriPlusvsMinus_sig_res, aes(y = reorder(description, sort(as.numeric(log2FoldChange))), x= log2FoldChange, fill = pvalue))+
  geom_bar(stat = "identity")+ 
  theme_bw()+
  labs(x = "Log2FoldChange", y="Pathways")

ggsave("cancer_Filtered_HPylori_YesVSNo_DESeq2.png", cancer_Filtered_HPylori_YesVSNo_DESeq2)


cancer_Ordered_Filtered_HPylori_YesVSNo_DESeq2 <- ggplot(data = cancerPyloriPlusvsMinus_sig_res, aes(y = reorder(description, -log2FoldChange), x= log2FoldChange, fill = pvalue))+
  geom_bar(stat = "identity") +
  theme_test()+
  theme(plot.margin = margin(t = 10, r = 10, b = 30, l = 10), 
        axis.text.y = element_text(size = 12), 
        axis.text.x = element_text(size = 12), 
        axis.title.x = element_text(size = 16, vjust = -5), 
        axis.title.y = element_text(size = 16), 
        legend.title = element_text(size = 16)) +
  scale_fill_viridis() +
  labs(x = "Log2FoldChange", y="Pathways", fill = "p value")

ggsave("cancer_Ordered_Filtered_HPylori_YesVSNo_DESeq2.png", cancer_Ordered_Filtered_HPylori_YesVSNo_DESeq2, width = 13.4, height = 8.5)



cancer_summary_table <- cancerPyloriPlusvsMinus_sig_res %>%
  summarize(
    Positive_Log2FC = sum(log2FoldChange > 0),  # Count positive values
    Negative_Log2FC = sum(log2FoldChange < 0)   # Count negative values
  )%>%
  mutate(Group = "Gastric Cancer")

combined_summary <- bind_rows(cg_summary_table, im_summary_table, ineo_summary_table, cancer_summary_table)


#### Stacked Bar Plot ####

combined_summary$Group <- factor(combined_summary$Group, levels = c("Chronic Gastritis", "Intestinal Metaplasia", "Intraepithelial Neoplasia", "Gastric Cancer"))
  
summary_long <- combined_summary %>% 
   pivot_longer(cols = c(Positive_Log2FC, Negative_Log2FC), 
               names_to = "Log2FC_Type", 
               values_to = "Count")

Stacked_HPylori_YesVSNo_Counts_Summary <- ggplot(summary_long, aes(x = Group, y = Count, fill = Log2FC_Type)) +
  geom_bar(stat = "identity") +
  labs(title = "Log2FoldChange Counts by Group", x = "Group", y = "# of Pathways", fill = " ") +
  scale_fill_manual(values = c("Positive_Log2FC" = "lightblue", "Negative_Log2FC" = "salmon"), 
                    labels = c("Uprepresented in H. Pylori Positive", "Uprepresented in H. Pylori Negative")) +  # Custom colors
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), 
        plot.title = element_text(hjust = 0.5))  # Rotate x-axis labels if needed

OneColour_Stacked_HPylori_YesVSNo_Counts_Summary <- ggplot(summary_long, aes(x = Group, y = Count, fill = Group)) +
  geom_bar(stat = "identity") +
  labs(title = expression(atop("Comparison of Significantly Changed Pathways between ", italic(H.~pylori)~" Status Across Disease Stages")),
       x = " ", y = "Number of Pathways") +
  theme_minimal() +
  scale_fill_okabe_ito() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), 
        plot.title = element_text(hjust = 0.5, margin = margin(b = 10)), 
        legend.position = "none")  # Rotate x-axis labels if needed



ggsave("OnrColour_Stacked_Filtered_HPylori_YesVSNo_SummaryBarGraph.png", OneColour_Stacked_HPylori_YesVSNo_Counts_Summary, height = 7, width = 7)

#### Side by side Bar Plot ####

sideBySide_HPylori_YesVSNo_Counts_Summary <- ggplot (summary_long, aes (x = Group, y = Count, fill = Log2FC_Type)) + 
  geom_bar(stat = "identity", position = "dodge") + 
  labs(title = "Log2FoldChange Counts by Group", x = "Group", y = "# of Pathways", fill = " ") +
  scale_fill_manual(values = c("Positive_Log2FC" = "lightblue", "Negative_Log2FC" = "salmon"), 
                    labels = c("Uprepresented in H. Pylori Positive", "Uprepresented in H. Pylori Negative")) +  # Custom colors
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), 
        plot.title = element_text(hjust = 0.5))

ggsave("sideBySide_Filtered_HPylori_YesVSNo_SummaryBarGraph.png", sideBySide_HPylori_YesVSNo_Counts_Summary)


