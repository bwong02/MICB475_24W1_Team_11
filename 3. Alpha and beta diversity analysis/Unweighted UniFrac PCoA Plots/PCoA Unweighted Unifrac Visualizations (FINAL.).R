library(tidyverse)
library(phyloseq)
library(vegan)
library(ggokabeito)

load("ps_rare_renamed.RData")
samp_data_wdiv <- data.frame(sample_data(ps_rare), estimate_richness(ps_rare))

#### Beta diversity analysis #####

#calculate weighted unifrac distance matrix
unifrac <- UniFrac(ps_rare, weighted = FALSE)

#Perform principal coordinates analysis (PCoA) on the unifrac distance matrix.
pcoa_unifrac <- ordinate(ps_rare, method="PCoA", distance = "unifrac", weighted = FALSE)

#Reorder factor levels in Group to make legend follow biological disease progression
ps_rare@sam_data$Group <- factor(ps_rare@sam_data$Group, levels = c("Healthy control (HC)", "Chronic gastritis (CG)", "Intestinal metaplasia (IM）",  "Intraepithelial neoplasia (IN)", "Gastric cancer (GC)"))

#### Fusobacterium ####

#plot PCoA with ellipses + significant differences

gg_pcoa_fuso_simple <- gg_pcoa_fuso_simple <- plot_ordination(ps_rare, pcoa_unifrac, color = "Group", shape = "Fusobacterium_abundance") +
                      stat_ellipse(type = "norm", level = 0.95, aes(group = Group)) +
                      scale_color_okabe_ito(
                        labels = c(
                          "Healthy control (HC)" = "Healthy Control",
                          "Chronic gastritis (CG)" = "Chronic Gastritis",
                          "Intestinal metaplasia (IM）" = "Intestinal Metaplasia",
                          "Intraepithelial neoplasia (IN)" = "Intraepithelial Neoplasia",
                          "Gastric cancer (GC)" = "Gastric Cancer")) + 
                      guides(color = guide_legend(title = "Disease stage"), shape = guide_legend(title = expression(italic("Fusobacterium") ~ "abundance"))) +
                      theme_test()

gg_pcoa_fuso_simple

#run PERMANOVA 
adonis2(unifrac ~ Stage_FusoAbundance, data = samp_data_wdiv, permutations = 999)


#Save visualization file.
ggsave("gg_pcoa_fuso_simple.png"
       , gg_pcoa_fuso_simple
       , height=4, width=8)

#### H. pylori ####

#plot PCoA with ellipses + significant differences

gg_pcoa_hpyl_simple <- gg_pcoa_fuso_simple <- plot_ordination(ps_rare, pcoa_unifrac, color = "Group", shape = "H.pylori_Status") +
                        stat_ellipse(type = "norm", level = 0.95, aes(group = Group)) +
                        scale_color_okabe_ito(
                          labels = c(
                            "Healthy control (HC)" = "Healthy Control",
                            "Chronic gastritis (CG)" = "Chronic Gastritis",
                            "Intestinal metaplasia (IM）" = "Intestinal Metaplasia",
                            "Intraepithelial neoplasia (IN)" = "Intraepithelial Neoplasia",
                            "Gastric cancer (GC)" = "Gastric Cancer")) + 
                        guides(color = guide_legend(title = "Disease stage"), shape = guide_legend(title = expression(italic("H. pylori") ~ "status"))) +
                        theme_test()

gg_pcoa_hpyl_simple

#run PERMANOVA 
adonis2(unifrac ~ Stage_Pylori, data = samp_data_wdiv, permutations = 999)

#Save visualization file.
ggsave("gg_pcoa_hpyl_simple.png"
       , gg_pcoa_hpyl_simple
       , height=4, width=7.5)

#run PERMANOVA based on disease stage 
adonis2(unifrac ~ Group, data = samp_data_wdiv, permutations = 999)
