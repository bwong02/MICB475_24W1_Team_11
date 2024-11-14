#### Load Packages ####
library(tidyverse)
library(phyloseq)
library(indicspecies)

#### Load Data ####
load("3AB_ps_filt.RData")

#### Prepare Data ####
ps_genus <- tax_glom(ps_filt, "Genus", NArm = FALSE) ##group data to the Genus level
ps_genus_RA <- transform_sample_counts(ps_genus, fun=function(x) x/sum(x)) ##Convert counts in OTU table to relative abundance
sample_data_df <- as.data.frame(sample_data(ps_filt)) #create dataframe 

#### Fusobacterium Indicator Species Analysis ####
fuso_isa <- multipatt(t(otu_table(ps_genus_RA)), cluster = sample_data_df$Stage_FusoAbundance) #conduct indicator species analysis using fusobacterium abundance at each stage
fuso_taxtable <- tax_table(ps_filt) %>% as.data.frame() %>% rownames_to_column(var="ASV") #extract taxonomy table
#Merge taxonomy table with phyloseq object and filter by significant p-value
fuso_isa_sum <- fuso_isa$sign %>%
  rownames_to_column(var = "ASV") %>%
  left_join(fuso_taxtable, by = "ASV") %>%
  filter(p.value < 0.05) 

#### H. Pylori Indicator Species Analysis ####
#conduct indicator species analysis
hp_isa <- multipatt(t(otu_table(ps_genus_RA)), cluster = sample_data_df$Stage_Pylori) #conduct indicator species analysis using H. Pylori presence at each stage
hp_taxtable <- tax_table(ps_filt) %>% as.data.frame() %>% rownames_to_column(var="ASV") #extract taxonomy table
#Merge taxonomy table with phyloseq object and filter by significant p-value
hp_isa_sum <- hp_isa$sign %>%
  rownames_to_column(var = "ASV") %>%
  left_join(hp_taxtable, by = "ASV") %>%
  filter(p.value < 0.05) 

#### Cancer Stage Indicator Species Analysis ####
stage_isa <- multipatt(t(otu_table(ps_genus_RA)), cluster = sample_data_df$Group) #conduct indicator species analysis using gastric cancer stage
stage_taxtable <- tax_table(ps_filt) %>% as.data.frame() %>% rownames_to_column(var="ASV") #extract taxonomy table
#Merge taxonomy table with phyloseq object and filter by significant p-value
stage_isa_sum <- stage_isa$sign %>%
  rownames_to_column(var = "ASV") %>%
  left_join(stage_taxtable, by = "ASV") %>%
  filter(p.value < 0.05)

#### Save Data #### 
save(fuso_isa_sum, file = "fuso_isa.RData") #save fusobacterium ISA
save(hp_isa_sum, file = "hp_isa.RData") #save H.pylori
save(stage_isa_sum, file = "stage_isa.RData") #GC stage ISA
