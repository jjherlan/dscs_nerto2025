# Modern GAM Analysis for Species Distribution Modeling
# Comprehensive version integrating stepwise selection and advanced plotting
# Updated to address deprecated functions and provide coherent analysis

# Load required libraries
library(biomod2)
library(MASS)
library(ggplot2)
library(gridExtra)  # For arranging multiple plots
library(dplyr)      # For data manipulation
library(viridis)    # For better color palettes
library(lattice)    # For additional plotting options

# Load the dataset
#mammals_data <- read.csv("data/tabular/species/mammals_and_bioclim_table.csv", row.names = 1)
ddia_btm <- read.csv("data/ne_canyons_dscs_btm_20250630.csv", row.names = 1)

# Data exploration and summary
# cat("=== Dataset Summary ===\n")
# cat("Total observations:", nrow(mammals_data), "\n")
# cat("Red fox presence records:", sum(mammals_data$VulpesVulpes == 1), "\n")
# cat("Red fox absence records:", sum(mammals_data$VulpesVulpes == 0), "\n")
# cat("Prevalence:", round(mean(mammals_data$VulpesVulpes), 3), "\n\n")

# Data exploration and summary
# cat("=== Dataset Summary ===\n")
# cat("Total observations:", nrow(ddia_adjsd_bpi), "\n")
# cat("Desmophyllum dianthus presence records:", sum(ddia_adjsd_bpi$Desmophyllum.dianthus == 1), "\n")
# cat("Desmophylum dianthus absence records:", sum(ddia_adjsd_bpi$Desmophyllum.dianthus == 0), "\n")
# cat("Prevalence:", round(mean(ddia_adjsd_bpi$Desmophyllum.dianthus), 3), "\n\n")

# Data exploration and summary
cat("=== Dataset Summary ===\n")
cat("Total observations:", nrow(ddia_btm), "\n")
cat("Desmophyllum dianthus presence records:", sum(ddia_btm$Desmophyllum.dianthus == 1), "\n")
cat("Desmophylum dianthus absence records:", sum(ddia_btm$Desmophyllum.dianthus == 0), "\n")
cat("Prevalence:", round(mean(ddia_btm$Desmophyllum.dianthus), 3), "\n\n")

# Check for missing values
# missing_summary <- colSums(is.na(mammals_data))
# if(any(missing_summary > 0)) {
#   cat("Missing values found:\n")
#   print(missing_summary[missing_summary > 0])
# } else {
#   cat("No missing values detected.\n")
# }

# Check for missing values
# missing_summary <- colSums(is.na(ddia_adjsd_bpi))
# if(any(missing_summary > 0)) {
#   cat("Missing values found:\n")
#   print(missing_summary[missing_summary > 0])
# } else {
#   cat("No missing values detected.\n")
# }

missing_summary <- colSums(is.na(ddia_btm))
if(any(missing_summary > 0)) {
  cat("Missing values found:\n")
  print(missing_summary[missing_summary > 0])
} else {
  cat("No missing values detected.\n")
}

### 10.3	Generalized Additive Models


# r code_10.3_Generalized_Additive_Models_b1

if(is.element("package:mgcv", search())) detach("package:mgcv") ## make sure the mgcv package is not loaded to avoid conflicts

library(gam)

#gam1 = gam(VulpesVulpes ~ s(bio3,2) + s(bio7,2) + s(bio11,2) + s(bio12,2), data=mammals_data, family="binomial")

#gam2 = gam(VulpesVulpes ~ s(bio3,4) + s(bio7,4) + s(bio11,4) + s(bio12,4), data=mammals_data, family="binomial")

gam1 = gam(Desmophyllum.dianthus ~ s(btm_arag_ann, 2) + s(btm_chl_ann, 2) + s(btm_curr_mag_ann, 2) + s(btm_dissic_ann, 2), 
           s(btm_dissoc_ann, 2) + s(btm_o2_ann, 2) + s(btm_sal_ann, 2) + s(btm_talk_ann, 2), s(btm_temp_ann, 2),
           data = ddia_btm, family = "binomial")

gam2 = gam(Desmophyllum.dianthus ~ s(btm_arag_ann, 4) + s(btm_chl_ann, 4) + s(btm_curr_mag_ann, 4) + s(btm_dissic_ann, 4), 
           s(btm_dissoc_ann, 4) + s(btm_o2_ann, 4) + s(btm_sal_ann, 4) + s(btm_talk_ann, 4), s(btm_temp_ann, 4),
           data = ddia_btm, family = "binomial")