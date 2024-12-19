#Load Packages
library(tidyverse)
library(phyloseq)
library(microbiome)
library(ggVennDiagram)

#Load unrarefied, filtered pyhloseq object
load("3AB_ps_filt.RData")

#convert to relative abundance
gc_RA <- transform_sample_counts(ps_filt, fun=function(x) x/sum(x))

####Core Microbiome Analysis of Groups ####

#Filter by group
in_group <- subset_samples(gc_RA, `Group`=="Intraepithelial neoplasia (IN)") #samples in intraepithelial neoplasia group
hc_group <- subset_samples(gc_RA, `Group`=="Healthy control (HC)") #samples in healthy control group
im_group <- subset_samples(gc_RA, `Group`=="Intestinal metaplasia (IM）") #samples in intestinal metaplasia group
gc_group <- subset_samples(gc_RA, `Group`=="Gastric cancer (GC)") #samples within gastric cancer group
cg_group <- subset_samples(gc_RA, `Group`=="Chronic gastritis (CG)") #samples within chronic gastritis

#Determine which ASVs are present in more than 90% of samples in each group
in_group_ASVs <- core_members(in_group, detection=0, prevalence = 0.9) #ASVs in 90% of samples in intraepithelial neoplasia group
hc_group_ASVs <- core_members(hc_group, detection=0, prevalence = 0.9) #ASVs in 90% of samples in healthy control group
im_group_ASVs <- core_members(im_group, detection=0, prevalence = 0.9) #ASVs in 90% of samples in intestinal metaplasia group
gc_group_ASVs <- core_members(gc_group, detection=0, prevalence = 0.9) #ASVs in 90% of samples within gastric cancer group
cg_group_ASVs <- core_members(cg_group, detection=0, prevalence = 0.9) #ASVs in 90% of samples within the chronic gastritis group

#Prune taxa tables and assign ASVs to genus and species level to create table with only assigned ASVs from each group
tax_table(prune_taxa(in_group_ASVs,ps_filt))
tax_table(prune_taxa(hc_group_ASVs,ps_filt))
tax_table(prune_taxa(im_group_ASVs,ps_filt))
tax_table(prune_taxa(gc_group_ASVs,ps_filt))
tax_table(prune_taxa(cg_group_ASVs,ps_filt))

#Create venn diagram comparing the ASVs in each group
group_list_full <- list(intneo = in_group_ASVs, hc = hc_group_ASVs,im = im_group_ASVs, gc = gc_group_ASVs, cg = cg_group_ASVs)
group_venn <- ggVennDiagram(x = group_list_full)

#Save the plot
ggsave("group_venn.png", group_venn)

####Core Microbiome Analysis of Fuso Abundance Low vs High by Groups ####

#Filter by high and low fuso abundance within each group
in_fuso_high <- subset_samples(gc_RA, `Stage_FusoAbundance`=="Intraepithelial neoplasia (IN)_High") #samples with high fuso abundance within intraepithalial neoplasia group
in_fuso_low <- subset_samples(gc_RA, `Stage_FusoAbundance`=="Intraepithelial neoplasia (IN)_Low") #samples with low fuso abundance within intraepithalial neoplasia group

hc_fuso_high <- subset_samples(gc_RA, `Stage_FusoAbundance`=="Healthy control (HC)_High") #samples with high fuso abundance within healthy control group
hc_fuso_low <- subset_samples(gc_RA, `Stage_FusoAbundance`=="Healthy control (HC)_Low") #samples with low fuso abundance within healthy control group

im_fuso_high <- subset_samples(gc_RA, `Stage_FusoAbundance`=="Intestinal metaplasia (IM）_High") #samples with high fuso abundance within intestinal metaplasia group
im_fuso_low <- subset_samples(gc_RA, `Stage_FusoAbundance`=="Intestinal metaplasia (IM）_Low") #samples with low fuso abundance within intestinal metaplasia group

gc_fuso_high <- subset_samples(gc_RA, `Stage_FusoAbundance`=="Gastric cancer (GC)_High") #samples with high fuso abundance within gastric cancer group
gc_fuso_low <- subset_samples(gc_RA, `Stage_FusoAbundance`=="Gastric cancer (GC)_Low") #samples with low fuso abundance within gastric cancer group

cg_fuso_high <- subset_samples(gc_RA, `Stage_FusoAbundance`=="Chronic gastritis (CG)_High") #samples with high fuso abundance within chronic gastritis group
cg_fuso_low <- subset_samples(gc_RA, `Stage_FusoAbundance`=="Chronic gastritis (CG)_Low") #samples with low fuso abundance within chronic gastritis group


#Determine which ASVs are present in more than 90% of samples in each group
in_fuso_high_ASVs <- core_members(in_fuso_high, detection=0, prevalence = 0.9) #ASVs in 90% of with high fuso abundance within intraepithalial neoplasia group
in_fuso_low_ASVs <- core_members(in_fuso_low, detection=0, prevalence = 0.9) #ASVs in 90% of with low fuso abundance within intraepithalial neoplasia group

hc_fuso_high_ASVs <- core_members(hc_fuso_high, detection=0, prevalence = 0.9) #ASVs in 90% of samples with high fuso abundance within healthy control group
hc_fuso_low_ASVs <- core_members(hc_fuso_low, detection=0, prevalence = 0.9) #ASVs in 90% of samples with low fuso abundance within healthy control group

im_fuso_high_ASVs <- core_members(im_fuso_high, detection=0, prevalence = 0.9) #ASVs in 90% of samples with high fuso abundance within intestinal metaplasia group
im_fuso_low_ASVs <- core_members(im_fuso_low, detection=0, prevalence = 0.9) #ASVs in 90% of samples with low fuso abundance within intestinal metaplasia group

gc_fuso_high_ASVs <- core_members(gc_fuso_high, detection=0, prevalence = 0.9) #ASVs in 90% of samples with high fuso abundance within gastric cancer group
gc_fuso_low_ASVs <- core_members(gc_fuso_low, detection=0, prevalence = 0.9) #ASVs in 90% of samples with low fuso abundance within gastric cancer group

cg_fuso_high_ASVs <- core_members(cg_fuso_high, detection=0, prevalence = 0.9) #ASVs in 90% of samples with high fuso abundance within chronic gastritis group
cg_fuso_low_ASVs <- core_members(cg_fuso_low, detection=0, prevalence = 0.9) #ASVs in 90% of samples low high fuso abundance within chronic gastritis group

#Prune taxa tables and assign ASVs to genus and species level to create table with only assigned ASVs from each group
tax_table(prune_taxa(in_fuso_high_ASVs,ps_filt))
tax_table(prune_taxa(in_fuso_low_ASVs,ps_filt))

tax_table(prune_taxa(hc_fuso_high_ASVs,ps_filt))
tax_table(prune_taxa(hc_fuso_low_ASVs,ps_filt))

tax_table(prune_taxa(im_fuso_high_ASVs,ps_filt))
tax_table(prune_taxa(im_fuso_low_ASVs,ps_filt))

tax_table(prune_taxa(gc_fuso_high_ASVs,ps_filt))
tax_table(prune_taxa(gc_fuso_low_ASVs,ps_filt))

tax_table(prune_taxa(cg_fuso_high_ASVs,ps_filt))
tax_table(prune_taxa(cg_fuso_low_ASVs,ps_filt))

#Create venn diagram comparing the ASVs in the intraepithalial neoplasia fuso high vs low groups
group_list_full <- list(intneo_fuso_high = in_fuso_high_ASVs, intneo_fuso_low = in_fuso_low_ASVs)
in_fuso_venn <- ggVennDiagram(x = group_list_full)
ggsave("in_fuso_venn.png", in_fuso_venn)

#Create venn diagram comparing the ASVs in the healthy control fuso high vs low groups
group_list_full <- list(hc_fuso_high = hc_fuso_high_ASVs, hc_fuso_low = hc_fuso_low_ASVs)
hc_fuso_venn <- ggVennDiagram(x = group_list_full)
ggsave("hc_fuso_venn.png", hc_fuso_venn)

#Create venn diagram comparing the ASVs in the intestinal metaplasia fuso high vs low groups
group_list_full <- list(im_fuso_high = im_fuso_high_ASVs, im_fuso_low = im_fuso_low_ASVs)
im_fuso_venn <- ggVennDiagram(x = group_list_full)
ggsave("im_fuso_venn.png", im_fuso_venn)

#Create venn diagram comparing the ASVs in the gastric cancer fuso high vs low groups
group_list_full <- list(gc_fuso_high = gc_fuso_high_ASVs, gc_fuso_low = gc_fuso_low_ASVs)
gc_fuso_venn <- ggVennDiagram(x = group_list_full)
ggsave("gc_fuso_venn.png", gc_fuso_venn)

#Create venn diagram comparing the ASVs in the chronic gastritis fuso high vs low groups
group_list_full <- list(cg_fuso_high = cg_fuso_high_ASVs, cg_fuso_low = cg_fuso_low_ASVs)
cg_fuso_venn <- ggVennDiagram(x = group_list_full)
ggsave("cg_fuso_venn.png", cg_fuso_venn)

####Core Microbiome Analysis of H. Pylori Presence Low vs High by Groups ####

#Filter by high and low H. pylori presence within each group 
#No healthy control because there are no H. Pylori positive healthy control samples
in_hp_high <- subset_samples(gc_RA, `Stage_Pylori`=="Intraepithelial neoplasia (IN)_Positive") #samples with high H. pylori presence within intraepithalial neoplasia group
in_hp_low <- subset_samples(gc_RA, `Stage_Pylori`=="Intraepithelial neoplasia (IN)_Negative") #samples with low H. pylori presence within intraepithalial neoplasia group

im_hp_high <- subset_samples(gc_RA, `Stage_Pylori`=="Intestinal metaplasia (IM）_Positive") #samples with high H. pylori presence within intestinal metaplasia group
im_hp_low <- subset_samples(gc_RA, `Stage_Pylori`=="Intestinal metaplasia (IM）_Negative") #samples with low H. pylori presence within intestinal metaplasia group

gc_hp_high <- subset_samples(gc_RA, `Stage_Pylori`=="Gastric cancer (GC)_Positive") #samples with high H. pylori presence within gastric cancer group
gc_hp_low <- subset_samples(gc_RA, `Stage_Pylori`=="Gastric cancer (GC)_Negative") #samples with low H. pylori presence within gastric cancer group

cg_hp_high <- subset_samples(gc_RA, `Stage_Pylori`=="Chronic gastritis (CG)_Positive") #samples with high H. pylori presence within chronic gastritis group
cg_hp_low <- subset_samples(gc_RA, `Stage_Pylori`=="Chronic gastritis (CG)_Negative") #samples with low H. pylori presence within chronic gastritis group


#Determine which ASVs are present in more than 90% of samples in each group
in_hp_high_ASVs <- core_members(in_hp_high, detection=0, prevalence = 0.9) #ASVs in 90% of with high H. pylori presence within intraepithalial neoplasia group
in_hp_low_ASVs <- core_members(in_hp_low, detection=0, prevalence = 0.9) #ASVs in 90% of with low H. pyloripresence within intraepithalial neoplasia group

im_hp_high_ASVs <- core_members(im_hp_high, detection=0, prevalence = 0.9) #ASVs in 90% of samples with high H. pylori presence within intestinal metaplasia group
im_hp_low_ASVs <- core_members(im_hp_low, detection=0, prevalence = 0.9) #ASVs in 90% of samples with low H. pylori presence within intestinal metaplasia group

gc_hp_high_ASVs <- core_members(gc_hp_high, detection=0, prevalence = 0.9) #ASVs in 90% of samples with high H. pylori presence within gastric cancer group
gc_hp_low_ASVs <- core_members(gc_hp_low, detection=0, prevalence = 0.9) #ASVs in 90% of samples with low H. pylori presence within gastric cancer group

cg_hp_high_ASVs <- core_members(cg_hp_high, detection=0, prevalence = 0.9) #ASVs in 90% of samples with high H. pylori presence within chronic gastritis group
cg_hp_low_ASVs <- core_members(cg_hp_low, detection=0, prevalence = 0.9) #ASVs in 90% of samples low high H. pylori presence within chronic gastritis group

#Prune taxa tables and assign ASVs to genus and species level to create table with only assigned ASVs from each group
tax_table(prune_taxa(in_hp_high_ASVs,ps_filt))
tax_table(prune_taxa(in_hp_low_ASVs,ps_filt))

tax_table(prune_taxa(im_hp_high_ASVs,ps_filt))
tax_table(prune_taxa(im_hp_low_ASVs,ps_filt))

tax_table(prune_taxa(gc_hp_high_ASVs,ps_filt))
tax_table(prune_taxa(gc_hp_low_ASVs,ps_filt))

tax_table(prune_taxa(cg_hp_high_ASVs,ps_filt))
tax_table(prune_taxa(cg_hp_low_ASVs,ps_filt))

#Create venn diagram comparing the ASVs in the intraepithalial neoplasia H. pylori pos vs neg groups
group_list_full <- list(intneo_hp_high = in_hp_high_ASVs, intneo_hp_low = in_hp_low_ASVs)
in_hp_venn <- ggVennDiagram(x = group_list_full)
ggsave("in_hp_venn.png", in_hp_venn)

#Create venn diagram comparing the ASVs in the intestinal metaplasia H. pylori pos vs neg groups
group_list_full <- list(im_hp_high = im_hp_high_ASVs, im_hp_low = im_hp_low_ASVs)
im_hp_venn <- ggVennDiagram(x = group_list_full)
ggsave("im_hp_venn.png", im_hp_venn)

#Create venn diagram comparing the ASVs in the gastric cancer H. pylori pos vs neg groups
group_list_full <- list(gc_hp_high = gc_hp_high_ASVs, gc_hp_low = gc_hp_low_ASVs)
gc_hp_venn <- ggVennDiagram(x = group_list_full)
ggsave("gc_hp_venn.png", gc_hp_venn)

#Create venn diagram comparing the ASVs in the chronic gastritis H. pylori pos vs neg groups
group_list_full <- list(cg_hp_high = cg_hp_high_ASVs, cg_hp_low = cg_hp_low_ASVs)
cg_hp_venn <- ggVennDiagram(x = group_list_full)
ggsave("cg_hp_venn.png", cg_hp_venn)
