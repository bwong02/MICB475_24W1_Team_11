#Load in libraries
library(tidyverse)
library(ggplot2)
library(stringr)


#Load in metadata
metadatafp <- "gastric_cancer_metadata.tsv"
metadata <- read_delim(metadatafp)

#Create group names vector
group_names <- c("Healthy control (HC)" = "Healthy Control", 
            "Chronic gastritis (CG)" = "Chronic Gastritis", 
            "Intestinal metaplasia (IM）" = "Intestinal Metaplasia", 
            "Intraepithelial neoplasia (IN)" = "Intraepithelial neoplasia", 
            "Gastric cancer (GC)" = "Gastric Cancer")

#Create new columns in metadata with: updated group name, H.pylori status, updated site name (capitalized and extra spaces removed)
metadata <- metadata %>%
  mutate(`Group Name` = recode(`Group`, !!!group_names)) %>%
  mutate(`Status` = `H. pylori status                                 13 C-urea breath test`) %>%
  mutate(`Site Name` = str_squish(str_to_title(metadata$Site)))

#summarise number of + and - H.pylori cases by group and site
groupSite_summarised <- metadata %>%
  group_by(`Group Name`, `Site Name`, `Status`) %>%
  summarise(count=n()) %>%
  mutate(total_count = sum(count)) %>%
  mutate(`relative_percentage` = (count/total_count) *100)


# Create a heat map using ggplot2
gg_heatmap <- ggplot(groupSite_summarised, aes(x = `Site Name`, y = Status, fill = relative_percentage)) +
  geom_tile(color = "white") +                      # Add borders to the tiles
  scale_fill_gradient(low = "white", high = "red", name = "Relative Percentage") +  # Color gradient
  labs(
    title = expression("Heatmap of Group and Site by " * italic("H. pylori") * " Status Count"), 
    x = "Group and Site", 
    y = "Status") +
  theme_minimal() + # Cleaner theme
  theme(axis.text.x = element_text(angle=90),
        plot.title = element_text(hjust = 0.5)) +
  facet_wrap(~ `Group Name`)
  
ggsave("hPyloris_initial_heat_map.png", gg_heatmap)

length(unique(metadata$'sample-id'))
min(metadata$'Age')
max(metadata$'Age')
mean(metadata$'Age')

female <- tibble(metadata$'Gender ' == "female")
sum(female$`metadata$"Gender " == "female"`==TRUE)
sum(female$`metadata$"Gender " == "female"`==FALSE)



#### Fusobacterium site map ####
newMetadataFP <- "LH_gastric_cancer_metadata.tsv"
newMetadata <- read_delim(newMetadataFP)

#Create new columns in metadata with updated site name (capitalized and extra spaces removed)
newMetadata <- newMetadata %>%
  mutate(`Group Name` = recode(`Group`, !!!group_names)) %>%
  mutate(`Site Name` = str_squish(str_to_title(newMetadata$Site)))

#Summarise number of low and high fusobacterium cases by group 
fuso_groupSite_summarised <- newMetadata %>%
  group_by(`Group Name`, `Site Name`, `Fusobacterium_abundance`) %>%
  summarise(count=n()) %>%
  mutate(total_count = sum(count)) %>%
  mutate(`relative_percentage` = (count/total_count) *100)

# Create a heat map using ggplot2
fuso_gg_heatmap <- ggplot(fuso_groupSite_summarised, aes(x = `Site Name`, y = Fusobacterium_abundance, fill = relative_percentage)) +
  geom_tile(color = "white") +                      # Add borders to the tiles
  scale_fill_gradient(low = "white", high = "red", name = "Relative Percentage") +  # Color gradient
  labs(
    title = expression("Heatmap of Group and Site by " * italic("Fusobacterium") * " Status Count"), 
    x = "Group and Site", 
    y = "Status") +
  theme_minimal() + # Cleaner theme
  theme(axis.text.x = element_text(angle=90), 
        plot.title = element_text(hjust = 0.5)) +
  facet_wrap(~ `Group Name`)

ggsave("fuso_initial_heat_map.png", fuso_gg_heatmap)
