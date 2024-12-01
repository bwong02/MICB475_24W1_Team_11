#Load Packages
library(tidyverse)
library(phyloseq)
library(microbiome)
library(ggVennDiagram)
library(ggsignif)
library(ggplot2)

#Load unrarefied, filtered pyhloseq object
load("3AB_ps_filt.RData")


# Load your phyloseq object
physeq <-  ps_filt

# Estimate richness, including observed features
richness <- estimate_richness(physeq, measures = "Observed")

# Extract sample data
sample_data(physeq)$ObservedFeatures <- richness$Observed


# Create a data frame combining sample metadata and observed features
plot_data <- data.frame(sample_data(physeq), ObservedFeatures = richness$Observed)

custom_labels <- c("Chronic gastritis (CG)_Negative" = "Chronic gastritis",
                   "Chronic gastritis (CG)_Positive" = "",
                   "Gastric cancer (GC)_Negative" = "Gastric cancer",
                   "Gastric cancer (GC)_Positive" = "",
                   "Healthy control (HC)_Negative" = "HC-",
                   "Intestinal metaplasia (IM）_Negative" = "Intestinal metaplasia",
                   "Intestinal metaplasia (IM）_Positive" = "",
                   "Intraepithelial neoplasia (IN)_Negative" = "Intraepithelial neoplasia",
                   "Intraepithelial neoplasia (IN)_Positive" = "")

# Filter out the "Healthy control (HC)_Negative" group since it doesn't have a Positive counterpart
plot_data_filtered <- plot_data[plot_data$Stage_Pylori != "Healthy control (HC)_Negative", ]

# Create a new variable to distinguish positive/negative within each group
plot_data_filtered$PosNeg <- factor(
  ifelse(grepl("Positive", plot_data_filtered$Stage_Pylori), "Positive", "Negative"),
  levels = c("Negative", "Positive")
)

# Reorder the Stage_Pylori factor levels to your desired order
plot_data_filtered$Stage_Pylori <- factor(plot_data_filtered$Stage_Pylori,
                                          levels = c("Chronic gastritis (CG)_Negative", "Chronic gastritis (CG)_Positive",
                                                     "Intestinal metaplasia (IM）_Negative", "Intestinal metaplasia (IM）_Positive",
                                                     "Intraepithelial neoplasia (IN)_Negative", "Intraepithelial neoplasia (IN)_Positive",
                                                     "Gastric cancer (GC)_Negative", "Gastric cancer (GC)_Positive"))

# Create the plot with reordered groups and significance testing
a_pyloriPlot = ggplot(plot_data_filtered, aes(x = Stage_Pylori, y = ObservedFeatures, fill = PosNeg)) + 
  geom_boxplot() +
  labs(
    title = "Alpha Diversity Comparison of Observed Features by H. pylori Status Across Disease Stages",
    x = "Disease Stage",
    y = "Number of Observed Features"
  ) +
  scale_x_discrete(labels = custom_labels) +  # Apply custom labels
  scale_fill_manual(
    values = c("Negative" = "lightblue", "Positive" = "salmon"),  # Colors for positive and negative
    name = "H. pylori status"  # Change the legend title here
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0), # Rotate x-axis labels
    plot.title = element_text(hjust = 0.5, size = 16)
  ) +
  geom_signif(
    comparisons = list(
      c("Chronic gastritis (CG)_Negative", "Chronic gastritis (CG)_Positive"),
      c("Gastric cancer (GC)_Negative", "Gastric cancer (GC)_Positive"),
      c("Intestinal metaplasia (IM）_Negative", "Intestinal metaplasia (IM）_Positive"),
      c("Intraepithelial neoplasia (IN)_Negative", "Intraepithelial neoplasia (IN)_Positive")
    ),
    map_signif_level = TRUE,
    step_increase = 0.1,
    test = "wilcox.test",  # Specify test type if needed
    test.args = list(exact = FALSE),  # Allow approximate p-value to avoid tie issues
    y_position = c(900, 810, 730, 640)
  )













custom_labels <- c("Chronic gastritis (CG)_Low" = "Chronic gastritis",
                   "Chronic gastritis (CG)_High" = "",
                   "Gastric cancer (GC)_Low" = "Gastric cancer",
                   "Gastric cancer (GC)_High" = "",
                   "Healthy control (HC)_Low" = "HC-",
                   "Intestinal metaplasia (IM）_Low" = "Intestinal metaplasia",
                   "Intestinal metaplasia (IM）_High" = "",
                   "Intraepithelial neoplasia (IN)_Low" = "Intraepithelial neoplasia",
                   "Intraepithelial neoplasia (IN)_High" = "")

# Create a new variable to distinguish positive/negative within each group
plot_data_filtered$HighLow <- factor(
  ifelse(grepl("High", plot_data_filtered$Stage_FusoAbundance), "High", "Low"),
  levels = c("Low", "High")
)

# Reorder the Stage_FusoAbundance factor levels to your desired order
plot_data_filtered$Stage_FusoAbundance <- factor(plot_data_filtered$Stage_FusoAbundance,
                                          levels = c("Chronic gastritis (CG)_Low", "Chronic gastritis (CG)_High",
                                                     "Intestinal metaplasia (IM）_Low", "Intestinal metaplasia (IM）_High",
                                                     "Intraepithelial neoplasia (IN)_Low", "Intraepithelial neoplasia (IN)_High",
                                                     "Gastric cancer (GC)_Low", "Gastric cancer (GC)_High"))

# Create the plot with reordered groups and significance testing
a_fusoPlot = ggplot(plot_data_filtered, aes(x = Stage_FusoAbundance, y = ObservedFeatures, fill = HighLow)) + 
  geom_boxplot() +
  labs(
    title = "Alpha Diversity Comparison of Observed Features by Fusobacterium Abundance Across Disease Stages",
    x = "Disease Stage",
    y = "Number of Observed Features"
  ) +
  scale_x_discrete(labels = custom_labels) +  # Apply custom labels
  scale_fill_manual(
    values = c("Low" = "lightblue", "High" = "salmon"),  # Colors for positive and negative
    name = "Fusobacterium\n abundance"  # Change the legend title here
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0), # Rotate x-axis labels
    plot.title = element_text(hjust = 0.1, size = 16)
  ) +
  geom_signif(
    comparisons = list(
      c("Chronic gastritis (CG)_Low", "Chronic gastritis (CG)_High"),
      c("Gastric cancer (GC)_Low", "Gastric cancer (GC)_High"),
      c("Intestinal metaplasia (IM）_Low", "Intestinal metaplasia (IM）_High"),
      c("Intraepithelial neoplasia (IN)_Low", "Intraepithelial neoplasia (IN)_High")
    ),
    map_signif_level = TRUE,
    step_increase = 0.1,
    test = "wilcox.test",  # Specify test type if needed
    test.args = list(exact = FALSE),  # Allow approximate p-value to avoid tie issues
    y_position = c(900, 810, 730, 640)
  )


ggsave("alpha_pylori_features.png", a_pyloriPlot)
ggsave("alpha_fuso_features.png", a_fusoPlot)
