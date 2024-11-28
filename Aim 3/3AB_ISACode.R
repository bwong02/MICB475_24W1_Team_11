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

#### Save Raw Data #### 
save(fuso_isa_sum, file = "fuso_isa.RData") #save fusobacterium ISA
save(hp_isa_sum, file = "hp_isa.RData") #save H.pylori
save(stage_isa_sum, file = "stage_isa.RData") #GC stage ISA

#### Filter and save Fusobacterium indicator species analysis data ####
columns_to_check_fuso <- c("s.Chronic gastritis (CG)_High", "s.Chronic gastritis (CG)_Low",
                      "s.Gastric cancer (GC)_High", "s.Gastric cancer (GC)_Low",
                      "s.Healthy control (HC)_High", "s.Healthy control (HC)_Low",
                      "s.Intestinal metaplasia (IM）_High", "s.Intestinal metaplasia (IM）_Low",
                      "s.Intraepithelial neoplasia (IN)_High", "s.Intraepithelial neoplasia (IN)_Low") #specify columns to check for 1s
final_fuso_isa <- fuso_isa_sum[rowSums(fuso_isa_sum[columns_to_check_fuso] == 1, na.rm = TRUE) <= 1, ] #remove rows with 1s in multiple columns
save(final_fuso_isa, file = "final_fuso_isa.RData") #filtered fuso abundance ISA

#### Filter and save H. pylori indicator species analysis data ####
columns_to_check_hp <- c("s.Chronic gastritis (CG)_Negative", "s.Chronic gastritis (CG)_Positive",
                           "s.Gastric cancer (GC)_Negative", "s.Gastric cancer (GC)_Positive",
                           "s.Healthy control (HC)_Negative",
                           "s.Intestinal metaplasia (IM）_Negative", "s.Intestinal metaplasia (IM）_Positive",
                           "s.Intraepithelial neoplasia (IN)_Negative", "s.Intraepithelial neoplasia (IN)_Positive") #specify columns to check for 1s
final_hp_isa <- hp_isa_sum[rowSums(hp_isa_sum[columns_to_check_hp] == 1, na.rm = TRUE) <= 1, ] #remove rows with 1s in multiple columns
save(final_hp_isa, file = "final_hp_isa.RData") #filtered H. pylori abundance ISA

#### Filter and save cancer stage indicator species analysis data ####
columns_to_check_group <- c("s.Chronic gastritis (CG)", "s.Gastric cancer (GC)", 
                      "s.Healthy control (HC)", "s.Intestinal metaplasia (IM）",
                      "s.Intraepithelial neoplasia (IN)") #specify columns to check for 1s
final_stages_isa <- stage_isa_sum[rowSums(stage_isa_sum[columns_to_check_group] == 1, na.rm = TRUE) <= 1, ] #remove rows with 1s in multiple columns
save(final_stages_isa, file = "final_stage_isa.RData") #filtered GC stage ISA
