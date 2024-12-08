##PERMANOVA

library(phyloseq)
library(ape)
library(tidyverse)
library(picante)
library(vegan)
library(RColorBrewer) 
library(ggplot2)

load("ps_rare.RData")
##ASK HANS (why estimate richness, isn't that alpha diversity?)
samp_data_wdiv <- data.frame(sample_data(ps_rare), estimate_richness(ps_rare))

#### Beta diversity analysis #####

#calculate weighted unifrac distance matrix
unifrac <- UniFrac(ps_rare, weighted = FALSE)

#Perform principal coordinates analysis (PCoA) on the unifrac distance matrix.
pcoa_unifrac <- ordinate(ps_rare, method="PCoA", distance = "unifrac", weighted = FALSE)

## Fusobacterium **

#run PERMANOVA 
adonis2(unifrac ~ Group, data = samp_data_wdiv, permutations = 999)

#Reorder factor levels in Stage_FusoAbundance to make color gradient and legend follow biological disease progression
ps_rare@sam_data$Group <- factor(ps_rare@sam_data$Group, levels = c("Healthy control (HC)", "Chronic gastritis (CG)", "Intestinal metaplasia (IM)",  "Intraepithelial neoplasia (IN)", "Gastric cancer (GC)"))

#plot PCoA with ellipses + significant differences 
gg_pcoa_fuso_ellipse <- plot_ordination(ps_rare, pcoa_unifrac, color = "Group", shape = "Fusobacterium_abundance") +
  stat_ellipse(type = "norm") +
  labs(col = expression(italic("Fusobacterium") ~ "Abundance and Disease Stage"))+
  scale_fill_manual(labels = c(
    "Healthy control (HC)" = "Healthy Control",
    "Chronic gastritis (CG)" = "Chronic Gastritis",
    "Intestinal metaplasia (IM）" = "Intestinal Metaplasia",
    "Intraepithelial neoplasia (IN)" = "Intraepithelial Neoplasia",
    "Gastric cancer (GC)" = "Gastric Cancer"
  ))
gg_pcoa_fuso_ellipse

ggsave("gg_pcoa_fuso_ellipse.png"
       , gg_pcoa_fuso
       , height=4, width=8)


## H. pylori **

#run PERMANOVA 
adonis2(unifrac ~ Stage_Pylori, data = samp_data_wdiv, permutations = 999)

#Reorder factor levels in Stage_Pylori to make color gradient and legend follow biological disease progression
ps_rare@sam_data$Stage_Pylori <- factor(ps_rare@sam_data$Stage_Pylori, levels = c("Healthy control (HC)_Negative", "Chronic gastritis (CG)_Negative", "Chronic gastritis (CG)_Positive", "Intestinal metaplasia (IM）_Negative", "Intestinal metaplasia (IM）_Positive", "Intraepithelial neoplasia (IN)_Negative", "Intraepithelial neoplasia (IN)_Positive", "Gastric cancer (GC)_Negative", "Gastric cancer (GC)_Positive"))

#Generate visualization of PCoA results based on disease stage and H. pylori status.
gg_pcoa_hpyl <- plot_ordination(ps_rare, pcoa_unifrac, color = "Group", shape = "H. pylori status                                 13 C-urea breath test") +
  stat_ellipse(type = "norm") +
  labs(col = expression(italic("H. pylori") ~ "Status and Disease Stage")) +
  theme_test() +
  scale_x_discrete(labels = c("Healthy control (HC)" = 
    "Chronic gastritis (CG)" = "Chronic gastritis",
    "Intestinal metaplasia (IM）" = "Intestinal metaplasia",
    "Intraepithelial neoplasia (IN)" = "Intraepithelial neoplasia",
    "Gastric cancer (GC)" = "Gastric cancer"
  )) 
gg_pcoa_hpyl

#Save visualization file.
ggsave("gg_pcoa_hpyl_ellipse.png"
       , gg_pcoa_hpyl
       , height=4, width=8)

