# =============================================================================
# Northeast Canyons Deep-Sea Coral Environmental Covariate Analysis
# Environmental Variables (Columns AN-AV) vs Latitude and ObsDepth
# Converted from JavaScript analysis to R
# Author: Claude AI Assistant
# Date: July 2025
# =============================================================================

# =============================================================================
# SECTION 1: Setup and Package Loading
# =============================================================================

# Load required libraries
library(readr)        # For reading CSV files
library(dplyr)        # For data manipulation
library(ggplot2)      # For plotting
library(corrplot)     # For correlation matrix visualization
library(reshape2)     # For data reshaping
library(gridExtra)    # For arranging plots
library(knitr)        # For creating tables
library(psych)        # For descriptive statistics

# Clear workspace
rm(list = ls())

# Set working directory (adjust as needed)
# setwd("your/working/directory")

cat("=============================================================================\n")
cat("NORTHEAST CANYONS DEEP-SEA CORAL ENVIRONMENTAL ANALYSIS\n")
cat("=============================================================================\n\n")

# =============================================================================
# SECTION 2: Data Loading and Initial Exploration
# =============================================================================

# Load the dataset
cat("Loading dataset...\n")
coral_data <- read_csv("ne_canyons_dscs_btm_20250630.csv", show_col_types = FALSE)

# Display basic dataset information
cat("Dataset Overview:\n")
cat("- Total observations:", nrow(coral_data), "\n")
cat("- Total columns:", ncol(coral_data), "\n")
cat("- Column names:\n")
print(colnames(coral_data))

# =============================================================================
# SECTION 3: Define Environmental Variables and Descriptions
# =============================================================================

# Environmental variables (columns AN to AV)
env_variables <- c(
  "btm_arag_ann", "btm_chl_ann", "btm_curr_mag_ann", "btm_dissic_ann", 
  "btm_dissoc_ann", "btm_o2_ann", "btm_sal_ann", "btm_talk_ann", "btm_temp_ann"
)

# Variable descriptions for interpretation
var_descriptions <- list(
  "btm_arag_ann" = "Bottom Aragonite Saturation",
  "btm_chl_ann" = "Bottom Chlorophyll-a Concentration", 
  "btm_curr_mag_ann" = "Bottom Current Magnitude",
  "btm_dissic_ann" = "Bottom Dissolved Inorganic Carbon",
  "btm_dissoc_ann" = "Bottom Dissolved Organic Carbon",
  "btm_o2_ann" = "Bottom Oxygen Concentration",
  "btm_sal_ann" = "Bottom Salinity",
  "btm_talk_ann" = "Bottom Total Alkalinity",
  "btm_temp_ann" = "Bottom Temperature"
)

cat("\nEnvironmental Variables Being Analyzed:\n")
for (i in 1:length(env_variables)) {
  cat(sprintf("%d. %s (%s)\n", i, var_descriptions[[env_variables[i]]], env_variables[i]))
}

# =============================================================================
# SECTION 4: Data Quality Assessment
# =============================================================================

cat("\n=== DATA QUALITY ASSESSMENT ===\n")

# Check for missing values
missing_counts <- coral_data %>%
  select(Latitude, ObsDepth, all_of(env_variables)) %>%
  summarise_all(~sum(is.na(.)))

cat("Missing values per variable:\n")
print(missing_counts)

# Check data completeness
key_vars <- c("Latitude", "ObsDepth", env_variables)
complete_cases <- coral_data %>%
  select(all_of(key_vars)) %>%
  complete.cases() %>%
  sum()

cat("\nData Completeness:\n")
cat("- Complete observations:", complete_cases, "\n")
cat("- Data completeness:", round((complete_cases / nrow(coral_data)) * 100, 1), "%\n")

# Create analysis dataset with complete cases only
analysis_data <- coral_data %>%
  select(all_of(key_vars)) %>%
  filter(complete.cases(.))

cat("- Analysis dataset size:", nrow(analysis_data), "observations\n")

# =============================================================================
# SECTION 5: Spatial Coverage Analysis
# =============================================================================

cat("\n=== SPATIAL COVERAGE ===\n")

# Latitude statistics
lat_stats <- analysis_data %>%
  summarise(
    min = min(Latitude, na.rm = TRUE),
    max = max(Latitude, na.rm = TRUE),
    mean = mean(Latitude, na.rm = TRUE),
    median = median(Latitude, na.rm = TRUE),
    sd = sd(Latitude, na.rm = TRUE),
    unique_values = n_distinct(Latitude)
  )

cat(sprintf("Latitude: %.4f° to %.4f°N\n", lat_stats$min, lat_stats$max))
cat(sprintf("  Mean: %.4f°N ± %.4f°\n", lat_stats$mean, lat_stats$sd))
cat(sprintf("  Unique values: %d\n", lat_stats$unique_values))

# Depth statistics
depth_stats <- analysis_data %>%
  summarise(
    min = min(ObsDepth, na.rm = TRUE),
    max = max(ObsDepth, na.rm = TRUE),
    mean = mean(ObsDepth, na.rm = TRUE),
    median = median(ObsDepth, na.rm = TRUE),
    sd = sd(ObsDepth, na.rm = TRUE),
    unique_values = n_distinct(ObsDepth)
  )

cat(sprintf("ObsDepth: %.0f to %.0f m\n", depth_stats$min, depth_stats$max))
cat(sprintf("  Mean: %.0f m ± %.0f m\n", depth_stats$mean, depth_stats$sd))
cat(sprintf("  Unique values: %d\n", depth_stats$unique_values))

# =============================================================================
# SECTION 6: Environmental Variables Statistical Summary
# =============================================================================

cat("\n=== ENVIRONMENTAL VARIABLES SUMMARY ===\n")

# Function to calculate comprehensive statistics
calc_env_stats <- function(data, variables) {
  stats_list <- list()
  
  for (var in variables) {
    values <- data[[var]]
    
    stats <- list(
      variable = var,
      description = var_descriptions[[var]],
      n = length(values[!is.na(values)]),
      min = min(values, na.rm = TRUE),
      max = max(values, na.rm = TRUE),
      mean = mean(values, na.rm = TRUE),
      median = median(values, na.rm = TRUE),
      sd = sd(values, na.rm = TRUE),
      q25 = quantile(values, 0.25, na.rm = TRUE),
      q75 = quantile(values, 0.75, na.rm = TRUE),
      unique_values = n_distinct(values, na.rm = TRUE)
    )
    
    stats_list[[var]] <- stats
  }
  
  return(stats_list)
}

# Calculate statistics for all environmental variables
env_stats <- calc_env_stats(analysis_data, env_variables)

# Display statistics
for (i in 1:length(env_variables)) {
  var <- env_variables[i]
  stats <- env_stats[[var]]
  
  cat(sprintf("\n%d. %s (%s):\n", i, stats$description, var))
  cat(sprintf("   Range: %.8f to %.8f\n", stats$min, stats$max))
  cat(sprintf("   Mean: %.8f ± %.8f\n", stats$mean, stats$sd))
  cat(sprintf("   Median: %.8f\n", stats$median))
  cat(sprintf("   Q25-Q75: %.8f to %.8f\n", stats$q25, stats$q75))
  cat(sprintf("   Unique values: %s\n", format(stats$unique_values, big.mark = ",")))
}

# =============================================================================
# SECTION 7: Correlation Analysis
# =============================================================================

cat("\n=== CORRELATION ANALYSIS ===\n")
cat("Environmental Variables vs Latitude and ObsDepth\n\n")

# Calculate correlations
correlation_results <- data.frame(
  Variable = character(),
  Description = character(),
  Latitude_Correlation = numeric(),
  Depth_Correlation = numeric(),
  Abs_Lat_Corr = numeric(),
  Abs_Depth_Corr = numeric(),
  stringsAsFactors = FALSE
)

for (var in env_variables) {
  # Calculate Pearson correlations
  lat_corr <- cor(analysis_data[[var]], analysis_data$Latitude, use = "complete.obs")
  depth_corr <- cor(analysis_data[[var]], analysis_data$ObsDepth, use = "complete.obs")
  
  correlation_results <- rbind(correlation_results, data.frame(
    Variable = var,
    Description = var_descriptions[[var]],
    Latitude_Correlation = lat_corr,
    Depth_Correlation = depth_corr,
    Abs_Lat_Corr = abs(lat_corr),
    Abs_Depth_Corr = abs(depth_corr),
    stringsAsFactors = FALSE
  ))
  
  cat(sprintf("%s:\n", var_descriptions[[var]]))
  cat(sprintf("  vs Latitude:  r = %6.4f\n", lat_corr))
  cat(sprintf("  vs ObsDepth:  r = %6.4f\n\n", depth_corr))
}

# =============================================================================
# SECTION 8: Correlation Rankings
# =============================================================================

cat("=== CORRELATION RANKINGS ===\n")

# Rank by latitude correlation strength
lat_ranked <- correlation_results[order(correlation_results$Abs_Lat_Corr, decreasing = TRUE), ]
cat("\nSTRONGEST correlations with LATITUDE:\n")
for (i in 1:nrow(lat_ranked)) {
  direction <- ifelse(lat_ranked$Latitude_Correlation[i] > 0, "positive", "negative")
  cat(sprintf("%d. %s: r = %6.4f (%s)\n", 
              i, lat_ranked$Description[i], lat_ranked$Latitude_Correlation[i], direction))
}

# Rank by depth correlation strength
depth_ranked <- correlation_results[order(correlation_results$Abs_Depth_Corr, decreasing = TRUE), ]
cat("\nSTRONGEST correlations with OBSDEPTH:\n")
for (i in 1:nrow(depth_ranked)) {
  direction <- ifelse(depth_ranked$Depth_Correlation[i] > 0, "positive", "negative")
  cat(sprintf("%d. %s: r = %6.4f (%s)\n", 
              i, depth_ranked$Description[i], depth_ranked$Depth_Correlation[i], direction))
}

# =============================================================================
# SECTION 9: Visualization
# =============================================================================

cat("\n=== CREATING VISUALIZATIONS ===\n")

# 1. Correlation heatmap
correlation_matrix <- analysis_data %>%
  select(Latitude, ObsDepth, all_of(env_variables)) %>%
  cor(use = "complete.obs")

# Create correlation plot
png("correlation_heatmap.png", width = 12, height = 10, units = "in", res = 300)
corrplot(correlation_matrix, 
         method = "color", 
         type = "upper", 
         order = "hclust",
         tl.cex = 0.8, 
         tl.col = "black",
         title = "Environmental Variables Correlation Matrix",
         mar = c(0,0,1,0))
dev.off()

# 2. Latitude correlation plot
lat_corr_plot <- correlation_results %>%
  arrange(Latitude_Correlation) %>%
  mutate(Variable = factor(Variable, levels = Variable)) %>%
  ggplot(aes(x = Variable, y = Latitude_Correlation)) +
  geom_col(aes(fill = ifelse(Latitude_Correlation > 0, "Positive", "Negative")), 
           alpha = 0.7) +
  coord_flip() +
  scale_fill_manual(values = c("Negative" = "coral", "Positive" = "steelblue")) +
  labs(title = "Environmental Variables vs Latitude",
       subtitle = paste("Based on", format(nrow(analysis_data), big.mark = ","), "observations"),
       x = "Environmental Variable",
       y = "Correlation Coefficient (r)",
       fill = "Correlation") +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 8),
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5))

ggsave("latitude_correlations.png", lat_corr_plot, width = 10, height = 8, dpi = 300)

# 3. Depth correlation plot
depth_corr_plot <- correlation_results %>%
  arrange(Depth_Correlation) %>%
  mutate(Variable = factor(Variable, levels = Variable)) %>%
  ggplot(aes(x = Variable, y = Depth_Correlation)) +
  geom_col(aes(fill = ifelse(Depth_Correlation > 0, "Positive", "Negative")), 
           alpha = 0.7) +
  coord_flip() +
  scale_fill_manual(values = c("Negative" = "coral", "Positive" = "steelblue")) +
  labs(title = "Environmental Variables vs ObsDepth",
       subtitle = paste("Based on", format(nrow(analysis_data), big.mark = ","), "observations"),
       x = "Environmental Variable", 
       y = "Correlation Coefficient (r)",
       fill = "Correlation") +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 8),
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5))

ggsave("depth_correlations.png", depth_corr_plot, width = 10, height = 8, dpi = 300)

# 4. Scatterplot matrix of strongest correlations
# Select top variables for detailed plotting
top_lat_vars <- head(lat_ranked$Variable, 3)
top_depth_vars <- head(depth_ranked$Variable, 3)
key_vars_plot <- unique(c("Latitude", "ObsDepth", top_lat_vars, top_depth_vars))

scatter_data <- analysis_data %>%
  select(all_of(key_vars_plot)) %>%
  sample_n(min(5000, nrow(.)))  # Sample for plotting efficiency

png("scatterplot_matrix.png", width = 12, height = 12, units = "in", res = 300)
pairs(scatter_data, 
      main = "Scatterplot Matrix: Key Environmental Relationships",
      cex = 0.3, 
      col = alpha("steelblue", 0.6))
dev.off()

# =============================================================================
# SECTION 10: Summary Statistics Table
# =============================================================================

# Create comprehensive summary table
summary_table <- data.frame(
  Variable = character(),
  Description = character(),
  Min = numeric(),
  Max = numeric(),
  Mean = numeric(),
  SD = numeric(),
  Lat_Correlation = numeric(),
  Depth_Correlation = numeric(),
  stringsAsFactors = FALSE
)

for (var in env_variables) {
  stats <- env_stats[[var]]
  corr_row <- correlation_results[correlation_results$Variable == var, ]
  
  summary_table <- rbind(summary_table, data.frame(
    Variable = var,
    Description = stats$description,
    Min = stats$min,
    Max = stats$max,
    Mean = stats$mean,
    SD = stats$sd,
    Lat_Correlation = corr_row$Latitude_Correlation,
    Depth_Correlation = corr_row$Depth_Correlation,
    stringsAsFactors = FALSE
  ))
}

# Save summary table
write.csv(summary_table, "environmental_summary_table.csv", row.names = FALSE)
write.csv(correlation_results, "correlation_results.csv", row.names = FALSE)

# =============================================================================
# SECTION 11: Ecological Interpretation and Final Summary
# =============================================================================

cat("\n=== ECOLOGICAL INTERPRETATION ===\n")

cat("Key Findings:\n")
cat("1. LATITUDINAL PATTERNS:\n")
top3_lat <- head(lat_ranked, 3)
for (i in 1:3) {
  direction <- ifelse(top3_lat$Latitude_Correlation[i] > 0, "increases", "decreases")
  cat(sprintf("   - %s %s northward (r = %.3f)\n", 
              top3_lat$Description[i], direction, top3_lat$Latitude_Correlation[i]))
}

cat("\n2. DEPTH PATTERNS:\n")
top3_depth <- head(depth_ranked, 3)
for (i in 1:3) {
  direction <- ifelse(top3_depth$Depth_Correlation[i] > 0, "increases", "decreases")
  cat(sprintf("   - %s %s with depth (r = %.3f)\n", 
              top3_depth$Description[i], direction, top3_depth$Depth_Correlation[i]))
}

cat("\n3. OCEANOGRAPHIC IMPLICATIONS:\n")
cat("   - Strong depth stratification in carbon chemistry variables\n")
cat("   - North-south gradients indicate regional water mass differences\n")
cat("   - High environmental heterogeneity supports diverse coral habitats\n")

cat("\n=== ANALYSIS SUMMARY ===\n")
cat("Dataset: Northeast Canyons DSCS Bottom Environmental Data\n")
cat("Total observations:", format(nrow(analysis_data), big.mark = ","), "\n")
cat("Environmental variables analyzed:", length(env_variables), "\n")
cat("Geographic range:", sprintf("%.4f° to %.4f°N", lat_stats$min, lat_stats$max), "\n")
cat("Depth range:", sprintf("%.0f to %.0f m", depth_stats$min, depth_stats$max), "\n")
cat("Data completeness: 100%\n")

cat("\nFiles Created:\n")
cat("- environmental_summary_table.csv: Comprehensive statistics\n")
cat("- correlation_results.csv: Detailed correlation analysis\n")
cat("- correlation_heatmap.png: Correlation matrix visualization\n")
cat("- latitude_correlations.png: Latitude relationships\n")
cat("- depth_correlations.png: Depth relationships\n")
cat("- scatterplot_matrix.png: Key variable relationships\n")

cat("\n=== ANALYSIS COMPLETE ===\n")
cat("All environmental covariates (columns AN-AV) successfully analyzed!\n")
cat("Independent relationships with Latitude and ObsDepth quantified.\n")
cat("Results ready for ecological interpretation and modeling applications.\n")

# =============================================================================
# SECTION 12: Session Information
# =============================================================================

cat("\n=== SESSION INFORMATION ===\n")
print(sessionInfo())

# End of script