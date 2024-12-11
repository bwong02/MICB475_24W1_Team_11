library(tidyverse)
library(ggplot2)
library(ggokabeito)

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
  labs(
    x = expression("" * italic(H.~pylori) * ""), 
    y = "Unweighted Unifrac Distance", 
    title = "" #expression(atop("Distance of " * italic(H.~pylori) * " Positive Samples", "Relative to Negative Samples Across Stages"))
  ) +
  scale_fill_okabe_ito(name = "Disease Stage",
  labels = c(
    "Chronic gastritis (CG)" = "Chronic gastritis",
    "Intestinal metaplasia (IM）" = "Intestinal metaplasia",
    "Intraepithelial neoplasia (IN)" = "Intraepithelial neoplasia",
    "Gastric cancer (GC)" = "Gastric cancer",
    "Healthy control (HC)" = "Healthy control"
  )) +  # Rename legend title
  scale_x_discrete(labels = c(
    "Chronic gastritis (CG)" = "",
    "Intestinal metaplasia (IM）" = "",
    "Intraepithelial neoplasia (IN)" = "",
    "Gastric cancer (GC)" = "",
    "Healthy control (HC)" = ""
  )) +  # Change x-axis names to acronyms
  theme_test() +
  theme(
    axis.text.x = element_text(angle = 0, size = 24),
    axis.text.y = element_text(size = 16),
    aspect.ratio = 1,  # Adjust the ratio to make it wider
    plot.title = element_text(hjust = 0.5, size = 28),
    #legend.position = "none",  # Remove legend
    legend.text = element_text(size = 18),
    legend.title = element_text(size = 20),
    axis.title = element_text(size = 20)
  ) +
  coord_cartesian(ylim = c(min(valid_data_filtered$Distance_change), max(valid_data_filtered$Distance_change) - 0.05)) +  # Adjust y-axis to remove empty space at the top
  annotate("text", x = "Intraepithelial neoplasia (IN)", y = max(valid_data_filtered$Distance_change) - 0.05, label = "**", size = 6, color = "black") + 
  annotate("text", x = "Gastric cancer (GC)", y = max(valid_data_filtered$Distance_change) - 0.05, label = "**", size = 6, color = "black")


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

# Remove "Healthy control (HC)" from the data
valid_data_filtered <- valid_data_filtered[valid_data_filtered$Stage1 != "Healthy control (HC)", ]

# Reorder Stage1 based on the desired order
valid_data_filtered$Stage1 <- factor(valid_data_filtered$Stage1, 
                                     levels = c("Chronic gastritis (CG)", "Intestinal metaplasia (IM）", 
                                                "Intraepithelial neoplasia (IN)", "Gastric cancer (GC)"))

# Replot with significance markers for IN and GC stages, without visual outliers
ggplot(valid_data_filtered, aes(x = Stage1, y = Distance_change)) +
  geom_boxplot(aes(fill = Stage1), outlier.shape = NA) +  # Remove visual outliers
  labs(x = expression("" * italic(Fusobacterium) * ""), y = "Unweighted Unifrac Distance", 
       title = "") + #"Distance of Fusobacterium High Abundance Samples \nRelative to Low Abundance Samples Across Stages") +
  scale_fill_okabe_ito(name = "Disease Stage",
  labels = c(
    "Chronic gastritis (CG)" = "Chronic \ngastritis",
    "Intestinal metaplasia (IM）" = "Intestinal \nmetaplasia",
    "Intraepithelial neoplasia (IN)" = "Intraepithelial \nneoplasia",
    "Gastric cancer (GC)" = "Gastric \ncancer",
    "Healthy control (HC)" = "Healthy \ncontrol"
  )) +  # Rename legend title
  scale_x_discrete(labels = c(
    "Chronic gastritis (CG)" = "",
    "Intestinal metaplasia (IM）" = "",
    "Intraepithelial neoplasia (IN)" = "",
    "Gastric cancer (GC)" = "",
    "Healthy control (HC)" = ""
  )) +  # Change x-axis names to acronyms
  theme_test() +
  theme(axis.text.x = element_text(angle = 0, size = 24),
  axis.text.y = element_text(size = 16),
  aspect.ratio = 1,  # Adjust the ratio to make it wider
  plot.title = element_text(hjust = 0.5, size = 28),
  legend.position = "none",
  axis.title = element_text(size=20)) + # Rotate x-axis labels for readability
  coord_cartesian(ylim = c(min(valid_data_filtered$Distance_change), max(valid_data_filtered$Distance_change) - 0.05))


###weighted


### H. PYLORI

# Load the data
data <- read.delim("raw_data_wUnifrac_stagePylori.tsv", header = TRUE, sep = "\t")

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

# Reorder Stage1 based on the desired order
valid_data_filtered$Stage1 <- factor(valid_data_filtered$Stage1, 
                                     levels = c("Chronic gastritis (CG)", "Intestinal metaplasia (IM）", 
                                                "Intraepithelial neoplasia (IN)", "Gastric cancer (GC)"))


# Replot with significance markers for IN and GC stages, without visual outliers
ggplot(valid_data_filtered, aes(x = Stage1, y = Distance_change)) +
  geom_boxplot(aes(fill = Stage1), outlier.shape = NA) +  # Remove visual outliers
  labs(
    x = expression("Disease Stage of " * italic(H.~pylori) * " Patients"), 
    y = "Weighted Unifrac Distance", 
    title = ""
  ) +
  scale_fill_okabe_ito(name = "Disease Stage") +  # Rename legend title
  scale_x_discrete(labels = c(
    "Chronic gastritis (CG)" = "Chronic \ngastritis",
    "Intestinal metaplasia (IM）" = "Intestinal \nmetaplasia",
    "Intraepithelial neoplasia (IN)" = "Intraepithelial \nneoplasia",
    "Gastric cancer (GC)" = "Gastric \ncancer"
  )) +  # Change x-axis names to acronyms
  theme_test() +
  theme(
    axis.text.x = element_text(angle = 0, size = 24),
    axis.text.y = element_text(size = 16),
    aspect.ratio = 1,  # Adjust the ratio to make it wider
    plot.title = element_text(hjust = 0.5, size = 28),
    axis.title = element_text(size = 20),
    legend.position = "none"  # Remove legend
  ) +
  coord_cartesian(ylim = c(min(valid_data_filtered$Distance_change), max(valid_data_filtered$Distance_change) - 0.05)) +  # Adjust y-axis to remove empty space at the top
    annotate("text", x = "Chronic gastritis (CG)", y = max(valid_data_filtered$Distance_change) - 0.05, label = "**", size = 6, color = "black") + 
  annotate("text", x = "Intraepithelial neoplasia (IN)", y = max(valid_data_filtered$Distance_change) - 0.05, label = "**", size = 6, color = "black") + 
  annotate("text", x = "Gastric cancer (GC)", y = max(valid_data_filtered$Distance_change) - 0.05, label = "**", size = 6, color = "black")


### FUSOBACTERIUM

# Load the data
data <- read.delim("raw_data_wUnifrac_stageFuso.tsv", header = TRUE, sep = "\t")

# Extract the stage from Group1 and Group2, and assign variant type
data$Stage1 <- sub("_.*", "", data$Group1)  # Extract stage from Group1
data$Stage2 <- sub("_.*", "", data$Group2)  # Extract stage from Group2

data$Variant1 <- ifelse(grepl("High", data$Group1), "High", "Low")
data$Variant2 <- ifelse(grepl("High", data$Group2), "High", "Low")

valid_data <- data[data$Stage1 == data$Stage2 & 
                   data$Variant1 == "Low" & 
                   data$Variant2 == "High", ]

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

# Remove "Healthy control (HC)" from the data
valid_data_filtered <- valid_data_filtered[valid_data_filtered$Stage1 != "Healthy control (HC)", ]

# Reorder Stage1 based on the desired order
valid_data_filtered$Stage1 <- factor(valid_data_filtered$Stage1, 
                                     levels = c("Chronic gastritis (CG)", "Intestinal metaplasia (IM）", 
                                                "Intraepithelial neoplasia (IN)", "Gastric cancer (GC)"))


# Replot with significance markers for IN and GC stages, without visual outliers
ggplot(valid_data_filtered, aes(x = Stage1, y = Distance_change)) +
  geom_boxplot(aes(fill = Stage1), outlier.shape = NA) +  # Remove visual outliers
  labs(x = expression("Disease Stage of " * italic(Fusobacterium) * " Patients"), y = "Weighted Unifrac Distance", 
       title = "") +
  scale_fill_okabe_ito(name = "Disease Stage") +  # Rename legend title
  scale_x_discrete(labels = c(
    "Chronic gastritis (CG)" = "Chronic \ngastritis",
    "Intestinal metaplasia (IM）" = "Intestinal \nmetaplasia",
    "Intraepithelial neoplasia (IN)" = "Intraepithelial \nneoplasia",
    "Gastric cancer (GC)" = "Gastric \ncancer"
  )) +  # Change x-axis names to acronyms
  theme_test() +
  theme(axis.text.x = element_text(angle = 0, size = 24),
  axis.text.y = element_text(size = 16),
  aspect.ratio = 1,  # Adjust the ratio to make it wider
  plot.title = element_text(hjust = 0.5, size = 28),
  axis.title = element_text(size=20),
  legend.position = "none") + # Rotate x-axis labels for readability
  coord_cartesian(ylim = c(min(valid_data_filtered$Distance_change), max(valid_data_filtered$Distance_change) - 0.05)) +
  annotate("text", x = "Chronic gastritis (CG)", y = max(valid_data_filtered$Distance_change) - 0.05, label = "*", size = 6, color = "black")
  
### Jaccard

### H. PYLORI

# Load the data
data <- read.delim("raw_data_jaccard_pylori.tsv", header = TRUE, sep = "\t")

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
  labs(
    x = expression("Disease Stage of " * italic(H.~pylori) * " Patients"), 
    y = "Jaccard Distance", 
    title = "" #expression(atop("Distance of " * italic(H.~pylori) * " Positive Samples", "Relative to Negative Samples Across Stages"))
  ) +
  scale_fill_okabe_ito(name = "Disease Stage",
  labels = c(
    "Chronic gastritis (CG)" = "Chronic gastritis",
    "Intestinal metaplasia (IM）" = "Intestinal metaplasia",
    "Intraepithelial neoplasia (IN)" = "Intraepithelial neoplasia",
    "Gastric cancer (GC)" = "Gastric cancer",
    "Healthy control (HC)" = "Healthy control"
  )) +  # Rename legend title
  scale_x_discrete(labels = c(
    "Chronic gastritis (CG)" = "Chronic \ngastritis",
    "Intestinal metaplasia (IM）" = "Intestinal \nmetaplasia",
    "Intraepithelial neoplasia (IN)" = "Intraepithelial \nneoplasia",
    "Gastric cancer (GC)" = "Gastric \ncancer"
  )) +  # Change x-axis names to acronyms
  theme_test() +
  theme(
    axis.text.x = element_text(angle = 0, size = 24),
    axis.text.y = element_text(size = 16),
    aspect.ratio = 1,  # Adjust the ratio to make it wider
    plot.title = element_text(hjust = 0.5, size = 28),
    legend.position = "none",  # Remove legend
    legend.text = element_text(size = 18),
    legend.title = element_text(size = 20),
    axis.title = element_text(size = 20)
  ) +
  coord_cartesian(ylim = c(min(valid_data_filtered$Distance_change), max(valid_data_filtered$Distance_change) - 0.07)) +  # Adjust y-axis to remove empty space at the top
  annotate("text", x = "Intraepithelial neoplasia (IN)", y = max(valid_data_filtered$Distance_change) - 0.07, label = "**", size = 6, color = "black") + 
  annotate("text", x = "Gastric cancer (GC)", y = max(valid_data_filtered$Distance_change) - 0.07, label = "**", size = 6, color = "black")


### Fuso

# Load the data
data <- read.delim("raw_data_jaccard_fuso.tsv", header = TRUE, sep = "\t")

# Extract the stage from Group1 and Group2, and assign variant type
data$Stage1 <- sub("_.*", "", data$Group1)  # Extract stage from Group1
data$Stage2 <- sub("_.*", "", data$Group2)  # Extract stage from Group2

data$Variant1 <- ifelse(grepl("High", data$Group1), "High", "Low")
data$Variant2 <- ifelse(grepl("High", data$Group2), "High", "Low")

valid_data <- data[data$Stage1 == data$Stage2 & 
                   data$Variant1 == "Low" & 
                   data$Variant2 == "High", ]

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

# Remove "Healthy control (HC)" from the data
valid_data_filtered <- valid_data_filtered[valid_data_filtered$Stage1 != "Healthy control (HC)", ]

# Reorder Stage1 based on the desired order
valid_data_filtered$Stage1 <- factor(valid_data_filtered$Stage1, 
                                     levels = c("Chronic gastritis (CG)", "Intestinal metaplasia (IM）", 
                                                "Intraepithelial neoplasia (IN)", "Gastric cancer (GC)"))


# Replot with significance markers for IN and GC stages, without visual outliers
ggplot(valid_data_filtered, aes(x = Stage1, y = Distance_change)) +
  geom_boxplot(aes(fill = Stage1), outlier.shape = NA) +  # Remove visual outliers
  labs(x = expression("Disease Stage of " * italic(Fusobacterium) * " Patients"), y = "Jaccard Distance", 
       title = "") +
  scale_fill_okabe_ito(name = "Disease Stage") +  # Rename legend title
  scale_x_discrete(labels = c(
    "Chronic gastritis (CG)" = "Chronic \ngastritis",
    "Intestinal metaplasia (IM）" = "Intestinal \nmetaplasia",
    "Intraepithelial neoplasia (IN)" = "Intraepithelial \nneoplasia",
    "Gastric cancer (GC)" = "Gastric \ncancer"
  )) +  # Change x-axis names to acronyms
  theme_test() +
  theme(axis.text.x = element_text(angle = 0, size = 24),
  axis.text.y = element_text(size = 16),
  aspect.ratio = 1,  # Adjust the ratio to make it wider
  plot.title = element_text(hjust = 0.5, size = 28),
  axis.title = element_text(size=20),
  legend.position = "none") +  # Rotate x-axis labels for readability
coord_cartesian(ylim = c(-0.25, max(valid_data_filtered$Distance_change)))  # Adjust y-axis to remove empty space at the top


### Bray curtis

### H. PYLORI

# Load the data
data <- read.delim("raw_data_brayC_pylori.tsv", header = TRUE, sep = "\t")

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

# Reorder Stage1 based on the desired order
valid_data_filtered$Stage1 <- factor(valid_data_filtered$Stage1, 
                                     levels = c("Chronic gastritis (CG)", "Intestinal metaplasia (IM）", 
                                                "Intraepithelial neoplasia (IN)", "Gastric cancer (GC)"))


# Replot with significance markers for IN and GC stages, without visual outliers
ggplot(valid_data_filtered, aes(x = Stage1, y = Distance_change)) +
  geom_boxplot(aes(fill = Stage1), outlier.shape = NA) +  # Remove visual outliers
  labs(
    x = expression("Disease Stage of " * italic(H.~pylori) * " Patients"), 
    y = "Bray Curtis Distance", 
    title = ""
  ) +
  scale_fill_okabe_ito(name = "Disease Stage") +  # Rename legend title
  scale_x_discrete(labels = c(
    "Chronic gastritis (CG)" = "Chronic \ngastritis",
    "Intestinal metaplasia (IM）" = "Intestinal \nmetaplasia",
    "Intraepithelial neoplasia (IN)" = "Intraepithelial \nneoplasia",
    "Gastric cancer (GC)" = "Gastric \ncancer"
  )) +  # Change x-axis names to acronyms
  theme_test() +
  theme(
    axis.text.x = element_text(angle = 0, size = 24),
    axis.text.y = element_text(size = 16),
    aspect.ratio = 1,  # Adjust the ratio to make it wider
    plot.title = element_text(hjust = 0.5, size = 28),
    axis.title = element_text(size = 20),
    legend.position = "none"  # Remove legend
  ) +
  coord_cartesian(ylim = c(min(valid_data_filtered$Distance_change) + 0.07, max(valid_data_filtered$Distance_change))) +  # Adjust y-axis to remove empty space at the top
    annotate("text", x = "Chronic gastritis (CG)", y = max(valid_data_filtered$Distance_change) + 0.01, label = "**", size = 6, color = "black") + 
  annotate("text", x = "Intraepithelial neoplasia (IN)", y = max(valid_data_filtered$Distance_change) + 0.01, label = "**", size = 6, color = "black") + 
  annotate("text", x = "Gastric cancer (GC)", y = max(valid_data_filtered$Distance_change) + 0.01, label = "**", size = 6, color = "black")


### FUSOBACTERIUM

# Load the data
data <- read.delim("raw_data_brayC_fuso.tsv", header = TRUE, sep = "\t")

# Extract the stage from Group1 and Group2, and assign variant type
data$Stage1 <- sub("_.*", "", data$Group1)  # Extract stage from Group1
data$Stage2 <- sub("_.*", "", data$Group2)  # Extract stage from Group2

data$Variant1 <- ifelse(grepl("High", data$Group1), "High", "Low")
data$Variant2 <- ifelse(grepl("High", data$Group2), "High", "Low")

valid_data <- data[data$Stage1 == data$Stage2 & 
                   data$Variant1 == "Low" & 
                   data$Variant2 == "High", ]

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

# Remove "Healthy control (HC)" from the data
valid_data_filtered <- valid_data_filtered[valid_data_filtered$Stage1 != "Healthy control (HC)", ]

# Reorder Stage1 based on the desired order
valid_data_filtered$Stage1 <- factor(valid_data_filtered$Stage1, 
                                     levels = c("Chronic gastritis (CG)", "Intestinal metaplasia (IM）", 
                                                "Intraepithelial neoplasia (IN)", "Gastric cancer (GC)"))


# Replot with significance markers for IN and GC stages, without visual outliers
ggplot(valid_data_filtered, aes(x = Stage1, y = Distance_change)) +
  geom_boxplot(aes(fill = Stage1), outlier.shape = NA) +  # Remove visual outliers
  labs(x = expression("Disease Stage of " * italic(Fusobacterium) * " Patients"), y = "Bray Curtis Distance", 
       title = "") +
  scale_fill_okabe_ito(name = "Disease Stage") +  # Rename legend title
  scale_x_discrete(labels = c(
    "Chronic gastritis (CG)" = "Chronic \ngastritis",
    "Intestinal metaplasia (IM）" = "Intestinal \nmetaplasia",
    "Intraepithelial neoplasia (IN)" = "Intraepithelial \nneoplasia",
    "Gastric cancer (GC)" = "Gastric \ncancer"
  )) +  # Change x-axis names to acronyms
  theme_test() +
  theme(axis.text.x = element_text(angle = 0, size = 24),
  axis.text.y = element_text(size = 16),
  aspect.ratio = 1,  # Adjust the ratio to make it wider
  plot.title = element_text(hjust = 0.5, size = 28),
  axis.title = element_text(size=20),
  legend.position = "none") + # Rotate x-axis labels for readability
  coord_cartesian(ylim = c(min(valid_data_filtered$Distance_change) + 0.2, max(valid_data_filtered$Distance_change) + 0.05)) +
  annotate("text", x = "Intraepithelial neoplasia (IN)", y = max(valid_data_filtered$Distance_change) + 0.05, label = "*", size = 6, color = "black")
  