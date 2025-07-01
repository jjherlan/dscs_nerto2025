# Northeast Canyons Deep-Sea Coral Presence-Absence Data Analysis
# R Script for analyzing coral species distributions
# Author: Claude AI Assistant
# Date: July 2025

# Load required libraries
library(dplyr)
library(ggplot2)
library(readr)
library(knitr)
library(kableExtra)

# Set working directory (adjust as needed)
# setwd("your/working/directory")

# =============================================================================
# STEP 1: Load and examine the data
# =============================================================================

# Load the dataset
coral_data <- read.csv("data/ne_canyons_dscs_btm_20250630.csv", stringsAsFactors = FALSE)

# Basic data information
cat("=== DATASET OVERVIEW ===\n")
cat("Total rows (sampling locations):", nrow(coral_data), "\n")
cat("Total columns:", ncol(coral_data), "\n")
cat("Column names:\n")
print(colnames(coral_data))

# =============================================================================
# STEP 2: Identify and extract species columns (M through AM)
# =============================================================================

# Species columns are from column 13 to 39 (M through AM in Excel)
species_columns <- colnames(coral_data)[13:39]
cat("\n=== SPECIES ANALYZED ===\n")
cat("Number of species:", length(species_columns), "\n")
cat("Species names:\n")
print(species_columns)

# Extract species presence-absence data
species_data <- coral_data[, species_columns]

# =============================================================================
# STEP 3: Calculate presence-absence statistics for each species
# =============================================================================

# Function to calculate species statistics
calculate_species_stats <- function(species_data) {
  species_stats <- data.frame(
    Species = character(),
    Presences = integer(),
    Absences = integer(),
    Total = integer(),
    Prevalence_Percent = numeric(),
    stringsAsFactors = FALSE
  )
  
  for (i in 1:ncol(species_data)) {
    species_name <- colnames(species_data)[i]
    presences <- sum(species_data[, i] == 1, na.rm = TRUE)
    absences <- sum(species_data[, i] == 0, na.rm = TRUE)
    total <- presences + absences
    prevalence <- (presences / total) * 100
    
    species_stats <- rbind(species_stats, data.frame(
      Species = species_name,
      Presences = presences,
      Absences = absences,
      Total = total,
      Prevalence_Percent = round(prevalence, 2),
      stringsAsFactors = FALSE
    ))
  }
  
  return(species_stats)
}

# Calculate statistics
species_stats <- calculate_species_stats(species_data)

# Sort by prevalence (descending)
species_stats <- species_stats[order(species_stats$Prevalence_Percent, decreasing = TRUE), ]

# Display results
cat("\n=== SPECIES PRESENCE-ABSENCE SUMMARY ===\n")
print(species_stats, row.names = FALSE)

# =============================================================================
# STEP 4: Overall dataset statistics
# =============================================================================

# Calculate overall statistics
total_presences <- sum(species_stats$Presences)
total_observations <- sum(species_stats$Total)
overall_presence_rate <- (total_presences / total_observations) * 100

cat("\n=== DATASET SUMMARY STATISTICS ===\n")
cat("Total sampling locations:", format(nrow(coral_data), big.mark = ","), "\n")
cat("Number of species analyzed:", length(species_columns), "\n")
cat("Total presence records across all species:", format(total_presences, big.mark = ","), "\n")
cat("Total observations (presence + absence):", format(total_observations, big.mark = ","), "\n")
cat("Overall presence rate:", round(overall_presence_rate, 2), "%\n")

# Prevalence statistics
prevalence_values <- species_stats$Prevalence_Percent
cat("\n=== PREVALENCE STATISTICS ===\n")
cat("Mean prevalence across species:", round(mean(prevalence_values), 2), "%\n")
cat("Median prevalence:", round(median(prevalence_values), 2), "%\n")
cat("Minimum prevalence:", round(min(prevalence_values), 2), "%\n")
cat("Maximum prevalence:", round(max(prevalence_values), 2), "%\n")
cat("Standard deviation:", round(sd(prevalence_values), 2), "%\n")

# =============================================================================
# STEP 5: Identify most common and rarest species
# =============================================================================

cat("\n=== MOST COMMON SPECIES (Top 5) ===\n")
top_5 <- head(species_stats, 5)
for (i in 1:nrow(top_5)) {
  cat(paste0(i, ". ", gsub("\\.", " ", top_5$Species[i]), ": ", 
             format(top_5$Presences[i], big.mark = ","), " presences (", 
             top_5$Prevalence_Percent[i], "%)\n"))
}

cat("\n=== RAREST SPECIES (Bottom 5) ===\n")
bottom_5 <- tail(species_stats, 5)
bottom_5 <- bottom_5[order(bottom_5$Prevalence_Percent), ]
for (i in 1:nrow(bottom_5)) {
  cat(paste0(i, ". ", gsub("\\.", " ", bottom_5$Species[i]), ": ", 
             bottom_5$Presences[i], " presences (", 
             bottom_5$Prevalence_Percent[i], "%)\n"))
}

# =============================================================================
# STEP 6: Taxonomic group analysis
# =============================================================================

# Define taxonomic groups based on species names
taxonomic_groups <- list(
  "Hard Corals (Scleractinia)" = c("Desmophyllum.dianthus", "Desmophyllum.pertusum", "Solenosmilia.variabilis"),
  "Black Corals (Antipatharia)" = c("Stichopathes", "Bathypathes", "Parantipathes", "Stauropathes", "Telopathes.magna"),
  "Soft Corals (Alcyonacea)" = c("Acanthogorgia", "Paramuricea", "Anthothela.grandiflora", "Swiftia", 
                                 "Trachythela.rudis", "Anthomastus", "Paragorgia", "Acanella", 
                                 "Primnoa", "Thouarella", "Chrysogorgia", "Radicipes", "Balticina"),
  "Sea Pens (Pennatulacea)" = c("Funiculina", "Kophobelemnon", "Pennatula", "Distichoptilum", "Protoptilum", "Umbellula")
)

# Calculate taxonomic group statistics
cat("\n=== TAXONOMIC GROUP ANALYSIS ===\n")
taxonomic_summary <- data.frame(
  Group = character(),
  Species_Count = integer(),
  Total_Presences = integer(),
  Group_Prevalence = numeric(),
  stringsAsFactors = FALSE
)

for (group_name in names(taxonomic_groups)) {
  group_species <- taxonomic_groups[[group_name]]
  group_species <- group_species[group_species %in% species_columns]  # Only include species in our dataset
  
  group_presences <- sum(species_stats$Presences[species_stats$Species %in% group_species])
  group_total <- length(group_species) * nrow(coral_data)
  group_prevalence <- (group_presences / group_total) * 100
  
  taxonomic_summary <- rbind(taxonomic_summary, data.frame(
    Group = group_name,
    Species_Count = length(group_species),
    Total_Presences = group_presences,
    Group_Prevalence = round(group_prevalence, 2),
    stringsAsFactors = FALSE
  ))
  
  cat("•", group_name, ":\n")
  cat("  - Species count:", length(group_species), "\n")
  cat("  - Total presences:", format(group_presences, big.mark = ","), "\n")
  cat("  - Group prevalence:", round(group_prevalence, 2), "%\n\n")
}

# =============================================================================
# STEP 7: Prevalence distribution analysis
# =============================================================================

# Categorize species by prevalence ranges
prevalence_categories <- cut(prevalence_values, 
                           breaks = c(0, 0.5, 1.0, 2.0, 5.0, Inf),
                           labels = c("Very Rare (<0.5%)", "Rare (0.5-1.0%)", 
                                    "Uncommon (1.0-2.0%)", "Common (2.0-5.0%)", 
                                    "Very Common (>5.0%)"),
                           right = FALSE)

prevalence_distribution <- table(prevalence_categories)

cat("=== PREVALENCE DISTRIBUTION ===\n")
for (i in 1:length(prevalence_distribution)) {
  cat("•", names(prevalence_distribution)[i], ":", prevalence_distribution[i], "species\n")
}

# =============================================================================
# STEP 8: Create visualizations
# =============================================================================

# 1. Species prevalence bar plot
p1 <- ggplot(species_stats, aes(x = reorder(Species, Prevalence_Percent), y = Prevalence_Percent)) +
  geom_bar(stat = "identity", fill = "steelblue", alpha = 0.7) +
  coord_flip() +
  labs(title = "Species Prevalence in Northeast Canyons",
       subtitle = paste("Based on", format(nrow(coral_data), big.mark = ","), "sampling locations"),
       x = "Species",
       y = "Prevalence (%)") +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 8))

print(p1)

# 2. Prevalence distribution histogram
p2 <- ggplot(species_stats, aes(x = Prevalence_Percent)) +
  geom_histogram(bins = 10, fill = "coral", alpha = 0.7, color = "black") +
  labs(title = "Distribution of Species Prevalence",
       x = "Prevalence (%)",
       y = "Number of Species") +
  theme_minimal()

print(p2)

# 3. Taxonomic group comparison
p3 <- ggplot(taxonomic_summary, aes(x = reorder(Group, Group_Prevalence), y = Group_Prevalence)) +
  geom_bar(stat = "identity", fill = "darkgreen", alpha = 0.7) +
  coord_flip() +
  labs(title = "Average Prevalence by Taxonomic Group",
       x = "Taxonomic Group",
       y = "Average Prevalence (%)") +
  theme_minimal()

print(p3)

# =============================================================================
# STEP 9: Save results to files
# =============================================================================

# Save species statistics to CSV
write.csv(species_stats, "ne_canyons_species_statistics.csv", row.names = FALSE)

# Save taxonomic summary to CSV
write.csv(taxonomic_summary, "ne_canyons_taxonomic_summary.csv", row.names = FALSE)

# Create a summary report
cat("\n=== KEY ECOLOGICAL INSIGHTS ===\n")
cat("1. Rarity is the norm:", sum(prevalence_values < 0.5), "out of", length(prevalence_values), "species occur at <0.5% of sites\n")
cat("2. Patchy distributions: Low overall presence rate (", round(overall_presence_rate, 2), "%) indicates specific habitat requirements\n")
cat("3. Taxonomic patterns: Soft corals dominate diversity (", 
    taxonomic_summary$Species_Count[taxonomic_summary$Group == "Soft Corals (Alcyonacea)"], " species) and abundance\n")
cat("4. Conservation concern: Many species, especially black corals, are extremely rare\n")

# Save plots
ggsave("species_prevalence_barplot.png", plot = p1, width = 12, height = 8, dpi = 300)
ggsave("prevalence_distribution_histogram.png", plot = p2, width = 8, height = 6, dpi = 300)
ggsave("taxonomic_group_comparison.png", plot = p3, width = 10, height = 6, dpi = 300)

cat("\n=== FILES CREATED ===\n")
cat("• ne_canyons_species_statistics.csv - Detailed species statistics\n")
cat("• ne_canyons_taxonomic_summary.csv - Taxonomic group summary\n")
cat("• species_prevalence_barplot.png - Species prevalence visualization\n")
cat("• prevalence_distribution_histogram.png - Prevalence distribution\n")
cat("• taxonomic_group_comparison.png - Taxonomic group comparison\n")

cat("\n=== ANALYSIS COMPLETE ===\n")
cat("All analyses completed successfully!\n")
cat("Northeast Canyons deep-sea coral diversity shows typical deep-sea patterns:\n")
cat("high species rarity, patchy distributions, and taxonomic specialization.\n")
