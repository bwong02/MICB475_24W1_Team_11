#Load Packages
library(tidyverse)
library(phyloseq)
library(ggsignif)
library(ggplot2)

#Load unrarefied, filtered pyhloseq object
load("3AB_ps_filt.RData")


# Load phyloseq object
physeq <-  ps_filt

# Estimate richness, including observed features
richness <- estimate_richness(physeq, measures = "Observed")

# Extract sample data
sample_data(physeq)$ObservedFeatures <- richness$Observed


# Create a data frame combining sample metadata and observed features
plot_data <- data.frame(sample_data(physeq), ObservedFeatures = richness$Observed)

custom_labels <- c("Chronic gastritis (CG)_Negative" = " Chronic \n gastritis",
                   "Chronic gastritis (CG)_Positive" = "",
                   "Gastric cancer (GC)_Negative" = "   Gastric \n   cancer",
                   "Gastric cancer (GC)_Positive" = "",
                   "Healthy control (HC)_Negative" = "HC-",
                   "Intestinal metaplasia (IM）_Negative" = "  Intestinal \n metaplasia",
                   "Intestinal metaplasia (IM）_Positive" = "",
                   "Intraepithelial neoplasia (IN)_Negative" = "Intraepithelial \n  neoplasia",
                   "Intraepithelial neoplasia (IN)_Positive" = "")

# Filter out the "Healthy control (HC)_Negative" group since it doesn't have a Positive counterpart
plot_data_filtered <- plot_data[plot_data$Stage_Pylori != "Healthy control (HC)_Negative", ]

# Create a new variable to distinguish positive/negative within each group
plot_data_filtered$PosNeg <- factor(
  ifelse(grepl("Positive", plot_data_filtered$Stage_Pylori), "Positive", "Negative"),
  levels = c("Negative", "Positive")
)

# Reorder the Stage_Pylori factor levels to order
plot_data_filtered$Stage_Pylori <- factor(plot_data_filtered$Stage_Pylori,
                                          levels = c("Chronic gastritis (CG)_Negative", "Chronic gastritis (CG)_Positive",
                                                     "Intestinal metaplasia (IM）_Negative", "Intestinal metaplasia (IM）_Positive",
                                                     "Intraepithelial neoplasia (IN)_Negative", "Intraepithelial neoplasia (IN)_Positive",
                                                     "Gastric cancer (GC)_Negative", "Gastric cancer (GC)_Positive"))

# Create the plot with reordered groups and significance testing
#a_pyloriPlot = 
ggplot(plot_data_filtered, aes(x = Stage_Pylori, y = ObservedFeatures, fill = PosNeg)) + 
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(
    aes(color = "black"),  
    width = 0.1,          # Adjust the horizontal spread of jitter
    size = 3,             # Adjust point size
    alpha = 0.3           # Make points slightly transparent
  ) +
  labs(
    title = "", #"Alpha Diversity Comparison of Observed Features by \nH. pylori Status Across Disease Stages",
    x = "Disease Stage",
    y = "Number of Observed Features"
  ) +
  scale_x_discrete(labels = custom_labels) +  # Apply custom labels
  scale_fill_manual(
    values = c("Negative" = "lightblue", "Positive" = "salmon"),  # Colors for positive and negative
    name = expression("" * italic(H.~pylori) * " Status")  # Change the legend title here
  ) +
  scale_color_manual(
    values = c("Negative" = "lightblue", "Positive" = "salmon"),  # Match jitter colors with boxplot fill
    name = expression("" * italic(H.~pylori) * " Status")
  ) +
  theme_test() +
  theme(
    axis.text.x = element_text(angle = 0, size = 22, hjust = 0.2), # Rotate x-axis labels
    plot.title = element_text(hjust = 0.5, size = 28),
    axis.text.y = element_text(size = 16),
    aspect.ratio = 1,  # Adjust the ratio to make it wider
    axis.title = element_text(size=20),
    #legend.position = "none", #removes legend
    legend.text = element_text(size = 16),  # Increase legend text size
    legend.title = element_text(size = 16)  # Increase legend title size
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


### Fusobacterium

# Initialize a data frame to store results
comparison_results <- data.frame(
  Comparison = character(),
  P_Value = numeric(),
  stringsAsFactors = FALSE
)

# Define the comparisons
comparisons <- list(
  c("Chronic gastritis (CG)_Negative", "Chronic gastritis (CG)_Positive"),
  c("Gastric cancer (GC)_Negative", "Gastric cancer (GC)_Positive"),
  c("Intestinal metaplasia (IM）_Negative", "Intestinal metaplasia (IM）_Positive"),
  c("Intraepithelial neoplasia (IN)_Negative", "Intraepithelial neoplasia (IN)_Positive")
)

# Perform the Wilcoxon test for each comparison
for (comparison in comparisons) {
  group1 <- comparison[1]
  group2 <- comparison[2]
  
  # Subset data for the two groups
  group1_data <- plot_data_filtered$ObservedFeatures[plot_data_filtered$Stage_Pylori == group1]
  group2_data <- plot_data_filtered$ObservedFeatures[plot_data_filtered$Stage_Pylori == group2]
  
  # Perform the test
  test_result <- wilcox.test(group1_data, group2_data, exact = FALSE)
  
  # Add to results data frame
  comparison_results <- rbind(
    comparison_results,
    data.frame(Comparison = paste(group1, "vs", group2),
               P_Value = test_result$p.value)
  )
}


custom_labels <- c("Chronic gastritis (CG)_Low" = " Chronic \n gastritis",
                   "Chronic gastritis (CG)_High" = "",
                   "Gastric cancer (GC)_Low" = "   Gastric \n   cancer",
                   "Gastric cancer (GC)_High" = "",
                   "Healthy control (HC)_Low" = "HC-",
                   "Intestinal metaplasia (IM）_Low" = "  Intestinal \n metaplasia",
                   "Intestinal metaplasia (IM）_High" = "",
                   "Intraepithelial neoplasia (IN)_Low" = "Intraepithelial \n  neoplasia",
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
#a_fusoPlot = 
ggplot(plot_data_filtered, aes(x = Stage_FusoAbundance, y = ObservedFeatures, fill = HighLow)) + 
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(
    aes(color = "black"),  
    width = 0.1,          # Adjust the horizontal spread of jitter
    size = 3,             # Adjust point size
    alpha = 0.3           # Make points slightly transparent
  ) +
  labs(
    title = "",
    x = "Disease Stage",
    y = "Number of Observed Features"
  ) +
  scale_x_discrete(labels = custom_labels) +  # Apply custom labels
  scale_fill_manual(
    values = c("Low" = "lightblue", "High" = "salmon"),  # Colors for positive and negative
    name = expression(atop("" * italic(Fusobacterium) * "", "Abundance"))  # Change the legend title
  ) +
  scale_color_manual(
    values = c("Low" = "lightblue", "High" = "salmon"),  # Match jitter colors with boxplot fill
    name = expression(atop("" * italic(Fusobacterium) * "", "Abundance")) 
  ) +
  theme_test() +
  theme(
    axis.text.x = element_text(angle = 0, size = 22, hjust = 0.2), # Rotate x-axis labels
    plot.title = element_text(hjust = 0.5, size = 28),
    axis.text.y = element_text(size = 16),
    aspect.ratio = 1,  # Adjust the ratio to make it wider
    axis.title = element_text(size=20),
    #legend.position = "none", #removes legend
    legend.text = element_text(size = 16),  # Increase legend text size
    legend.title = element_text(size = 16)  # Increase legend title size
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
    test = "wilcox.test",  # Specify test type
    test.args = list(exact = FALSE),  # Allow approximate p-value to avoid tie issues
    y_position = c(900, 810, 730, 640)
  )


ggsave("alpha_pylori_features.png", a_pyloriPlot)
ggsave("alpha_fuso_features.png", a_fusoPlot)
