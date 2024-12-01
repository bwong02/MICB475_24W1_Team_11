library(phyloseq)
library(ape)
library(tidyverse)
library(picante)

load("ps_filt.RData")
load("ps_rare.RData")

#### Beta diversity analysis #####

#Calculate the unweighted unifrac distance between samples in the rarefied phyloseq object.
unweighted_unifrac <- distance(ps_rare, method="unifrac", weighted = FALSE)

#Perform principal coordinates analysis (PCoA) on the unifrac distance matrix.
pcoa_unifrac <- ordinate(ps_rare, method="PCoA", distance=unweighted_unifrac)

###Fusobacterium Abundance and Disease Stage PCoA Plot###

#Reorder factor levels in Stage_FusoAbundance to make color gradient and legend follow biological disease progression
ps_rare@sam_data$Stage_FusoAbundance <- factor(ps_rare@sam_data$Stage_FusoAbundance, levels = c("Healthy control (HC)_Low", "Healthy control (HC)_High", "Chronic gastritis (CG)_Low", "Chronic gastritis (CG)_High", "Intestinal metaplasia (IM）_Low", "Intestinal metaplasia (IM）_High", "Intraepithelial neoplasia (IN)_Low", "Intraepithelial neoplasia (IN)_High", "Gastric cancer (GC)_Low", "Gastric cancer (GC)_High"))

#Check if variables were ordered correctly.
table(ps_rare@sam_data$Stage_FusoAbundance)

#Generate visualization of PCoA results based on disease stage and Fusobacterium abundance.
gg_pcoa_fuso <- plot_ordination(ps_rare, pcoa_unifrac, color = "Stage_FusoAbundance") +
  labs(col = expression(italic("Fusobacterium") ~ "Abundance and Disease Stage"),
       title = expression(italic("Fusobacterium") ~ "Abundance and Disease Stage"))
gg_pcoa_fuso

#Save visualization file.
ggsave("gg_pcoa_fuso_final.png"
       , gg_pcoa_fuso
       , height=4, width=8)

###H. pylori Abundance and Disease Stage PCoA Plot###
table(ps_rare@sam_data$Stage_Pylori)

#Reorder factor levels in Stage_Pylori to make color gradient and legend follow biological disease progression
ps_rare@sam_data$Stage_Pylori <- factor(ps_rare@sam_data$Stage_Pylori, levels = c("Healthy control (HC)_Negative", "Chronic gastritis (CG)_Negative", "Chronic gastritis (CG)_Positive", "Intestinal metaplasia (IM）_Negative", "Intestinal metaplasia (IM）_Positive", "Intraepithelial neoplasia (IN)_Negative", "Intraepithelial neoplasia (IN)_Positive", "Gastric cancer (GC)_Negative", "Gastric cancer (GC)_Positive"))

#Check if variables were ordered correctly.
table(ps_rare@sam_data$Stage_Pylori)

#Generate visualization of PCoA results based on disease stage and H. pylori status.
gg_pcoa_hpyl <- plot_ordination(ps_rare, pcoa_unifrac, color = "Stage_Pylori") +
  labs(col = expression(italic("H. pylori") ~ "Status and Disease Stage"),
       title = expression(italic("Helicobacter pylori") ~ "Status and Disease Stage"))
gg_pcoa_hpyl

#Save visualization file.
ggsave("gg_pcoa_hpyl_final.png"
       , gg_pcoa_hpyl
       , height=4, width=8)
