#### Load Packages ####
library(tidyverse)
library(phyloseq)
library(indicspecies)

#### Load data ####
load("3AB_ps_filt.RData")

#### Indicator Species/Taxa Analysis ####
# glom to Genus
ps_genus <- tax_glom(ps_filt, "Genus", NArm = FALSE)
ps_genus_RA <- transform_sample_counts(ps_genus, fun=function(x) x/sum(x))

#ISA
sample_data_df <- as.data.frame(sample_data(ps_filt))
isa_ps <- multipatt(t(otu_table(ps_genus_RA)), cluster = sample_data_df$Stage_FusoAbundance)
summary(isa_ps)
taxtable <- tax_table(ps_filt) %>% as.data.frame() %>% rownames_to_column(var="ASV")

### Combine and Save ISA and Taxatable ####
isa_ps_sum <- isa_ps$sign %>%
  rownames_to_column(var = "ASV") %>%
  left_join(taxtable, by = "ASV") %>%
  filter(p.value < 0.05)

save(isa_ps_sum, file = "isa_ps_sum.RData")
