#### Install Packages ####
library(phyloseq)
library(ape) # importing trees
library(tidyverse)
library(vegan)

#### Load Data ####
#load metadata
metafp <- "LH_gastric_cancer_metadata.tsv"
meta <- read_delim(metafp, delim="\t")
#load feature table
otufp <- "feature-table.txt" #LOAD THIS
otu <- read_delim(file = asvfp, delim="\t", skip=1) 
#load taxonomy data
taxfp <- "taxonomy.tsv" #LOAD THIS
tax <- read_delim(taxfp, delim="\t")
#load phylogenetic tree
phylotreefp <- "tree.nwk" #LOAD THIS
phylotree <- read.tree(phylotreefp)

#### Format Feature Table ####
#save everything except first column (OTU ID) into a matrix
otu_mat <- as.matrix(otu[,-1])
# Make first column (#OTU ID) the rownames of the new matrix
rownames(otu_mat) <- otu$`#OTU ID`
# Use the "otu_table" function to make an OTU table
OTU <- otu_table(otu_mat, taxa_are_rows = TRUE) 
class(OTU)

#### Format Sample Metadata ####
#Save everything except sampleid as new data frame
samp_df <- as.data.frame(meta[,-1])
# Make sample ids the rownames
rownames(samp_df)<- meta$'sample-id'
# Make phyloseq sample data with sample_data() function
SAMP <- sample_data(samp_df)
class(SAMP)

#### Format Taxonomy Data ####
# Convert taxon strings to a table with separate taxa rank columns
tax_mat <- tax %>% select(-Confidence)%>%
  separate(col=Taxon, sep="; "
           , into = c("Domain","Phylum","Class","Order","Family","Genus","Species")) %>%
  as.matrix() # Saving as a matrix
# Save everything except feature IDs 
tax_mat <- tax_mat[,-1]
# Make sample ids the rownames
rownames(tax_mat) <- tax$`Feature ID`
# Make taxa table
TAX <- tax_table(tax_mat)
class(TAX)

#### Create Phyloseq Object ####
# Merge all into a phyloseq object
ps <- phyloseq(OTU, SAMP, TAX, phylotree)

#### Filter Data ####
# Remove non-bacterial sequences, if any
ps_filt <- subset_taxa(mpt,  Domain == "d__Bacteria" & Class!="c__Chloroplast" & Family !="f__Mitochondria")

#### Rarefy samples ####
rarecurve(t(as.data.frame(otu_table(ps_filt))), cex=0.1)
ps_rare <- rarefy_even_depth(ps_filt, rngseed = 1, sample.size = 1000)

#### Save Data ####
save(ps_filt, file="ps_filt.RData")
save(ps_rare, file="ps_rare.RData")