#Load in libraries
library(tidyverse)
library(ggplot2)
library(stringr)
library(viridis)


newMetadataFP <- "LH_gastric_cancer_metadata.tsv"
newMetadata <- read_delim(newMetadataFP)

#Create group names vector
group_names <- c("Healthy control (HC)" = "Healthy Control", 
            "Chronic gastritis (CG)" = "Chronic Gastritis", 
            "Intestinal metaplasia (IM）" = "Intestinal Metaplasia", 
            "Intraepithelial neoplasia (IN)" = "Intraepithelial Neoplasia", 
            "Gastric cancer (GC)" = "Gastric Cancer")

#Create new columns in metadata with: updated group name, H.pylori status, updated site name (capitalized and extra spaces removed)
newMetadata <- newMetadata %>%
  mutate(`Group Name` = recode(`Group`, !!!group_names)) %>%
  mutate(`H. Pylori Status` = `H. pylori status                                 13 C-urea breath test`) %>%
  mutate(`Site Name` = str_squish(str_to_title(newMetadata$Site)))


#summarise number of + and - H.pylori cases by group and site
Hpylori_groupSite_summarised <- newMetadata %>%
  group_by(`Group Name`, `Site Name`, `H. Pylori Status`) %>%
  summarise(count=n()) %>%
  mutate(total_count = sum(count)) %>%
  mutate(`relative_percentage` = (count/total_count) *100)


# Create a heat map using ggplot2
hPylori_gg_heatmap <- ggplot(Hpylori_groupSite_summarised, aes(x = `Site Name`, y = `H. Pylori Status`, fill = relative_percentage)) +
  geom_tile(color = "white") +                      # Add borders to the tiles
  scale_fill_viridis(limits = c(10, 100),name = "Relative Percentage") +  # Color gradient
  labs(
    x = "Site", 
    y = "Status") +
  theme_test() + # Cleaner theme
  theme(axis.text.x = element_text(angle=90),
        plot.title = element_text(hjust = 0.5)) +
  facet_wrap(~ `Group Name`)
  
ggsave("hPyloris_initial_heat_map.png", hPylori_gg_heatmap)


#### Fusobacterium site map ####

#Summarise number of low and high fusobacterium cases by group 
fuso_groupSite_summarised <- newMetadata %>%
  group_by(`Group Name`, `Site Name`, `Fusobacterium_abundance`) %>%
  summarise(count=n()) %>%
  mutate(total_count = sum(count)) %>%
  mutate(`relative_percentage` = (count/total_count) *100)

# Create a heat map using ggplot2
fuso_gg_heatmap <- ggplot(fuso_groupSite_summarised, aes(x = `Site Name`, y = Fusobacterium_abundance, fill = relative_percentage)) +
  geom_tile(color = "white") +                      # Add borders to the tiles
  scale_fill_viridis(limits = c(10, 100), name = "Relative Percentage") +  # Color gradient
  labs(
    x = "Site", 
    y = "Status") +
  theme_test() + # Cleaner theme
  theme(axis.text.x = element_text(angle=90), 
        plot.title = element_text(hjust = 0.5)) +
  facet_wrap(~ `Group Name`)

ggsave("fuso_initial_heat_map.png", fuso_gg_heatmap)
