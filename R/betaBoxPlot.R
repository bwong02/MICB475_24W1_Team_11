library(phyloseq)
library(tidyverse)
library(ape)
library(vegan)
library(picante)
library(ggplot2)
library(ggsignif)

### H. PYLORI

# Load the data
data <- read.delim("raw_data_uunifrac_stagePylori.tsv", header = TRUE, sep = "\t")

# Extract the stage from Group1 and Group2, and assign variant type (Positive/Negative)
data$Stage1 <- sub("_.*", "", data$Group1)  # Extract stage from Group1
data$Stage2 <- sub("_.*", "", data$Group2)  # Extract stage from Group2

# Identify the variant type (Positive/Negative) in both Group1 and Group2
data$Variant1 <- ifelse(grepl("Positive", data$Group1), "Positive", "Negative")
data$Variant2 <- ifelse(grepl("Positive", data$Group2), "Positive", "Negative")

# Filter for rows where:
# - Stage1 == Stage2 (both samples are from the same stage)
# - Group1 is Negative and Group2 is Positive
valid_data <- data[data$Stage1 == data$Stage2 & 
                   data$Variant1 == "Negative" & 
                   data$Variant2 == "Positive", ]

# Prepare the data for plotting
valid_data_filtered$Stage <- valid_data_filtered$Stage1  # Use Stage1 as the final stage column

# Calculate the median distance for Negative/Negative group for each stage
negative_negative_data <- data[data$Variant1 == "Negative" & data$Variant2 == "Negative", ]

# Compute the median for each stage
negative_medians <- negative_negative_data %>%
  group_by(Stage1) %>%
  summarise(median_distance = median(Distance, na.rm = TRUE))

# Filter for Negative/Positive group (ensure it's the right comparison)
valid_data_filtered <- data[data$Variant1 == "Negative" & data$Variant2 == "Positive", ]
# Remove Healthy control as no positive samples exist
valid_data_filtered <- valid_data_filtered[valid_data_filtered$Stage1 != "Healthy control (HC)", ]
valid_data_filtered <- valid_data_filtered[valid_data_filtered$Stage2 != "Healthy control (HC)", ]

# Merge the Negative/Negative medians with the Negative/Positive data
valid_data_filtered <- merge(valid_data_filtered, negative_medians, by.x = "Stage1", by.y = "Stage1")

# Calculate the change in distance (subtract median Negative/Negative distance)
valid_data_filtered$Distance_change <- valid_data_filtered$Distance - valid_data_filtered$median_distance

# Flip the sign of the Distance_change to make the values positive
valid_data_filtered$Distance_change <- -valid_data_filtered$Distance_change

# Recalculate the median changes with the new flipped values
median_changes <- valid_data_filtered %>%
  group_by(Stage1) %>%
  summarise(median_distance_change = median(Distance_change, na.rm = TRUE))

# Reorder Stage1 based on the new median change in distance
valid_data_filtered$Stage1 <- factor(valid_data_filtered$Stage1, 
                                     levels = median_changes$Stage1[order(median_changes$median_distance_change)])

# Replot with significance markers for IN and GC stages, without visual outliers
ggplot(valid_data_filtered, aes(x = Stage1, y = Distance_change)) +
  geom_boxplot(aes(fill = Stage1), outlier.shape = NA) +  # Remove visual outliers
  labs(x = "Disease Stage", y = "Unweighted Unifrac Distance", 
       title = "Distance of H. pylori positive samples relative to the median of negative samples across stages") +
  scale_fill_discrete(name = "Disease Stage") +  # Rename legend title
  scale_x_discrete(labels = c(
    "Chronic gastritis (CG)" = "CG",
    "Intestinal metaplasia (IM）" = "IM",
    "Intraepithelial neoplasia (IN)" = "IN",
    "Gastric cancer (GC)" = "GC"
  )) +  # Change x-axis names to acronyms
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +  # Rotate x-axis labels for readability
  annotate("text", x = "Intraepithelial neoplasia (IN)", y = max(valid_data_filtered$Distance_change), label = "*", size = 6, color = "red") + 
  annotate("text", x = "Gastric cancer (GC)", y = max(valid_data_filtered$Distance_change), label = "*", size = 6, color = "red")


### FUSOBACTERIUM

# Load the data
data <- read.delim("raw_data_uunifrac_stageFuso.tsv", header = TRUE, sep = "\t")

# Extract the stage from Group1 and Group2, and assign variant type
data$Stage1 <- sub("_.*", "", data$Group1)  # Extract stage from Group1
data$Stage2 <- sub("_.*", "", data$Group2)  # Extract stage from Group2

data$Variant1 <- ifelse(grepl("High", data$Group1), "High", "Low")
data$Variant2 <- ifelse(grepl("High", data$Group2), "High", "Low")

valid_data <- data[data$Stage1 == data$Stage2 & 
                   data$Variant1 == "Low" & 
                   data$Variant2 == "High", ]

# Prepare the data for plotting
valid_data_filtered$Stage <- valid_data_filtered$Stage1  # Use Stage1 as the final stage column

# Calculate the median distance
negative_negative_data <- data[data$Variant1 == "Low" & data$Variant2 == "Low", ]

# Compute the median for each stage
negative_medians <- negative_negative_data %>%
  group_by(Stage1) %>%
  summarise(median_distance = median(Distance, na.rm = TRUE))

valid_data_filtered <- data[data$Variant1 == "Low" & data$Variant2 == "High", ]

valid_data_filtered <- merge(valid_data_filtered, negative_medians, by.x = "Stage1", by.y = "Stage1")

# Calculate the change in distance 
valid_data_filtered$Distance_change <- valid_data_filtered$Distance - valid_data_filtered$median_distance

# Flip the sign of the Distance_change to make the values positive
valid_data_filtered$Distance_change <- -valid_data_filtered$Distance_change

# Recalculate the median changes with the new flipped values
median_changes <- valid_data_filtered %>%
  group_by(Stage1) %>%
  summarise(median_distance_change = median(Distance_change, na.rm = TRUE))

# Reorder Stage1 based on the new median change in distance
valid_data_filtered$Stage1 <- factor(valid_data_filtered$Stage1, 
                                     levels = median_changes$Stage1[order(median_changes$median_distance_change)])

# Replot with significance markers for IN and GC stages, without visual outliers
ggplot(valid_data_filtered, aes(x = Stage1, y = Distance_change)) +
  geom_boxplot(aes(fill = Stage1), outlier.shape = NA) +  # Remove visual outliers
  labs(x = "Disease Stage", y = "Unweighted Unifrac Distance", 
       title = "Distance of Fusobacterium high samples relative to the median of low samples across stages") +
  scale_fill_discrete(name = "Disease Stage") +  # Rename legend title
  scale_x_discrete(labels = c(
    "Chronic gastritis (CG)" = "CG",
    "Intestinal metaplasia (IM）" = "IM",
    "Intraepithelial neoplasia (IN)" = "IN",
    "Gastric cancer (GC)" = "GC",
    "Healthy control (HC)" = "HC"
  )) +  # Change x-axis names to acronyms
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) # Rotate x-axis labels for readability
