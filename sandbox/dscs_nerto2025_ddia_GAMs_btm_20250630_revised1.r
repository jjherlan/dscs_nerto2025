# Modern GAM Analysis for Species Distribution Modeling
# Comprehensive version integrating stepwise selection and advanced plotting
# Updated to address deprecated functions and provide coherent analysis
# REVISED VERSION - Fixed GAM syntax errors and improved workflow

# Load required libraries
library(biomod2)
library(MASS)
library(ggplot2)
library(gridExtra)  # For arranging multiple plots
library(dplyr)      # For data manipulation
library(viridis)    # For better color palettes
library(lattice)    # For additional plotting options

# Load the dataset
ddia_btm <- read.csv("ne_canyons_dscs_btm_20250630.csv", row.names = 1)

# Data exploration and summary
cat("=== Dataset Summary ===\n")
cat("Total observations:", nrow(ddia_btm), "\n")
cat("Desmophyllum dianthus presence records:", sum(ddia_btm$Desmophyllum.dianthus == 1), "\n")
cat("Desmophyllum dianthus absence records:", sum(ddia_btm$Desmophyllum.dianthus == 0), "\n")
cat("Prevalence:", round(mean(ddia_btm$Desmophyllum.dianthus), 3), "\n\n")

# Check for missing values
missing_summary <- colSums(is.na(ddia_btm))
if(any(missing_summary > 0)) {
  cat("Missing values found:\n")
  print(missing_summary[missing_summary > 0])
} else {
  cat("No missing values detected.\n")
}

# Check unique values in environmental variables
env_vars <- c("btm_arag_ann", "btm_chl_ann", "btm_curr_mag_ann", "btm_dissic_ann", 
              "btm_dissoc_ann", "btm_o2_ann", "btm_sal_ann", "btm_talk_ann", "btm_temp_ann")

cat("\n=== Environmental Variables Summary ===\n")
for(var in env_vars) {
  unique_vals <- length(unique(ddia_btm[[var]], na.rm = TRUE))
  cat(sprintf("%s: %d unique values\n", var, unique_vals))
  
  if(unique_vals <= 3) {
    cat(sprintf("  WARNING: %s has only %d unique values - insufficient for GAM smoothing\n", var, unique_vals))
  }
}

### 10.3 Generalized Additive Models

# Ensure mgcv package is not loaded to avoid conflicts
if("package:mgcv" %in% search()) detach("package:mgcv")

library(gam)

# CORRECTED GAM SYNTAX - The original code had syntax errors
# Each s() function should be separate, connected with + not commas within s()

# GAM with degree 2 smoothing
cat("\n=== Fitting GAM with degree 2 smoothing ===\n")
gam1 <- gam(Desmophyllum.dianthus ~ s(btm_arag_ann, 2) + s(btm_chl_ann, 2) + 
            s(btm_curr_mag_ann, 2) + s(btm_dissic_ann, 2) + s(btm_dissoc_ann, 2) + 
            s(btm_o2_ann, 2) + s(btm_sal_ann, 2) + s(btm_talk_ann, 2) + s(btm_temp_ann, 2),
            data = ddia_btm, family = "binomial")

# GAM with degree 4 smoothing  
cat("=== Fitting GAM with degree 4 smoothing ===\n")
gam2 <- gam(Desmophyllum.dianthus ~ s(btm_arag_ann, 4) + s(btm_chl_ann, 4) + 
            s(btm_curr_mag_ann, 4) + s(btm_dissic_ann, 4) + s(btm_dissoc_ann, 4) + 
            s(btm_o2_ann, 4) + s(btm_sal_ann, 4) + s(btm_talk_ann, 4) + s(btm_temp_ann, 4),
            data = ddia_btm, family = "binomial")

cat("GAM models fitted successfully!\n")

# Model summaries
cat("\n=== GAM1 Summary (degree 2) ===\n")
print(summary(gam1))

cat("\n=== GAM2 Summary (degree 4) ===\n")
print(summary(gam2))

# =============================================================================
# MODERN PLOTTING FUNCTIONS (replacing deprecated level.plot and response.plot2)
# =============================================================================

# Function to create spatial distribution plots
create_spatial_plots <- function(data, response_vars, titles, coords_x = "Longitude", coords_y = "Latitude") {
  
  plot_list <- list()
  
  for(i in 1:length(response_vars)) {
    
    # Prepare data
    plot_data <- data.frame(
      x = data[[coords_x]],
      y = data[[coords_y]],
      response = response_vars[[i]]
    )
    
    # Remove any infinite or NA values
    plot_data <- plot_data[is.finite(plot_data$response), ]
    
    # Create plot
    p <- ggplot(plot_data, aes(x = x, y = y, color = response)) +
      geom_point(size = 0.5, alpha = 0.6) +
      scale_color_viridis_c(name = "Probability") +
      labs(title = titles[i],
           x = "Longitude", 
           y = "Latitude") +
      theme_minimal() +
      theme(
        plot.title = element_text(hjust = 0.5, size = 12),
        axis.text = element_text(size = 8),
        legend.position = "bottom"
      ) +
      coord_fixed()
    
    plot_list[[i]] <- p
  }
  
  return(plot_list)
}

# Function to create response curves (modern replacement for response.plot2)
create_response_curves <- function(model, data, env_vars) {
  
  response_data <- list()
  
  for(var in env_vars) {
    
    # Create sequence of values for this variable
    var_range <- range(data[[var]], na.rm = TRUE)
    var_seq <- seq(var_range[1], var_range[2], length.out = 100)
    
    # Create prediction data with other variables at their median
    pred_data <- data.frame(
      matrix(rep(sapply(data[env_vars], median, na.rm = TRUE), each = 100), 
             nrow = 100, ncol = length(env_vars))
    )
    colnames(pred_data) <- env_vars
    pred_data[[var]] <- var_seq
    
    # Get predictions
    predictions <- predict(model, newdata = pred_data, type = "response")
    
    # Store results
    response_data[[var]] <- data.frame(
      variable = var,
      value = var_seq,
      response = predictions
    )
  }
  
  # Combine all data
  response_df <- do.call(rbind, response_data)
  
  return(response_df)
}

# Create spatial distribution plots
cat("\n=== Creating spatial distribution plots ===\n")

spatial_plots <- create_spatial_plots(
  data = ddia_btm,
  response_vars = list(
    ddia_btm$Desmophyllum.dianthus,
    fitted(gam1),
    fitted(gam2)
  ),
  titles = c(
    "Original Data",
    "GAM1 (degree 2)",
    "GAM2 (degree 4)"
  )
)

# Display spatial plots
grid.arrange(grobs = spatial_plots, ncol = 2, nrow = 2)

# Create response curves
cat("=== Creating response curves ===\n")

response_data_gam1 <- create_response_curves(gam1, ddia_btm, env_vars)
response_data_gam2 <- create_response_curves(gam2, ddia_btm, env_vars)

# Add model identifier
response_data_gam1$model <- "GAM1 (degree 2)"
response_data_gam2$model <- "GAM2 (degree 4)"

# Combine data
combined_response_data <- rbind(response_data_gam1, response_data_gam2)

# Create response curve plot
response_plot <- ggplot(combined_response_data, aes(x = value, y = response, color = model)) +
  geom_line(size = 1) +
  facet_wrap(~ variable, scales = "free_x", nrow = 3, ncol = 3) +
  labs(title = "GAM Response Curves",
       x = "Environmental Variable Value",
       y = "Probability of Occurrence",
       color = "Model") +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14),
    strip.text = element_text(size = 10),
    axis.text = element_text(size = 8),
    legend.position = "bottom"
  ) +
  scale_color_manual(values = c("GAM1 (degree 2)" = "blue", "GAM2 (degree 4)" = "red"))

print(response_plot)

# =============================================================================
# MODEL EVALUATION AND COMPARISON
# =============================================================================

cat("\n=== Model Evaluation ===\n")

# AIC comparison
cat("AIC Comparison:\n")
cat("GAM1 (degree 2):", AIC(gam1), "\n")
cat("GAM2 (degree 4):", AIC(gam2), "\n")

# Deviance explained
cat("\nDeviance Explained:\n")
cat("GAM1:", round(((gam1$null.deviance - gam1$deviance) / gam1$null.deviance) * 100, 2), "%\n")
cat("GAM2:", round(((gam2$null.deviance - gam2$deviance) / gam2$null.deviance) * 100, 2), "%\n")

# Predictions summary
pred_gam1 <- fitted(gam1)
pred_gam2 <- fitted(gam2)

cat("\nPrediction Summaries:\n")
cat("GAM1 - Mean predicted probability:", round(mean(pred_gam1), 4), "\n")
cat("GAM1 - Range:", round(min(pred_gam1), 4), "to", round(max(pred_gam1), 4), "\n")
cat("GAM2 - Mean predicted probability:", round(mean(pred_gam2), 4), "\n")
cat("GAM2 - Range:", round(min(pred_gam2), 4), "to", round(max(pred_gam2), 4), "\n")

# =============================================================================
# ALTERNATIVE: Using mgcv package (more modern GAM implementation)
# =============================================================================

cat("\n=== Alternative: Using mgcv package ===\n")

# Detach gam package and load mgcv
detach("package:gam", unload = TRUE)
library(mgcv)

# Fit GAMs using mgcv (more flexible and modern)
cat("Fitting GAM using mgcv package...\n")

gam_mgcv1 <- gam(Desmophyllum.dianthus ~ s(btm_arag_ann, k = 4) + s(btm_chl_ann, k = 4) + 
                 s(btm_curr_mag_ann, k = 4) + s(btm_dissic_ann, k = 4) + s(btm_dissoc_ann, k = 4) + 
                 s(btm_o2_ann, k = 4) + s(btm_sal_ann, k = 4) + s(btm_talk_ann, k = 4) + s(btm_temp_ann, k = 4),
                 data = ddia_btm, family = binomial, method = "REML")

cat("mgcv GAM fitted successfully!\n")

# Summary
print(summary(gam_mgcv1))

# Plot smooths
cat("Creating mgcv diagnostic plots...\n")
par(mfrow = c(3, 3))
plot(gam_mgcv1, pages = 1, residuals = TRUE, pch = 19, cex = 0.3)
par(mfrow = c(1, 1))

cat("\n=== Analysis Complete ===\n")
cat("All GAM models fitted and visualized successfully!\n")
cat("Key improvements made:\n")
cat("1. Fixed GAM syntax errors\n")
cat("2. Replaced deprecated plotting functions\n")
cat("3. Added comprehensive model evaluation\n")
cat("4. Included modern mgcv alternative\n")
cat("5. Enhanced visualization with ggplot2\n")

# r GAMb2 10.6, opts.label = 'fig_half_page', fig.cap = 'Figure 10. 6: Response curves of model gam1 expressed in logit scale (function plot.gam() from the gam package).'}

par(mfrow = c(2,2))

plot(gam1, se = T)

# r GAMb3 10.7, message=FALSE,warning=FALSE, opts.label = 'fig_quarter_page', fig.cap = 'Figure 10. 7: Response curves of the gam1 (degree of smoothing = 2) and gam2 (degree of smoothing = 4) models.'}

# rp <- response.plot2(models = c('gam1', 'gam2'),
#                      Data = mammals_data[,c("bio3", "bio7", "bio11", "bio12")],
#                      show.variables = c("bio3",  "bio7", "bio11", "bio12"),
#                      fixed.var.metric = 'mean', plot = FALSE, use.formal.names = TRUE)

rp <- response.plot2(models = c('gam1', 'gam2'),
                     Data = mammals_data[,c("btm_arag_ann", "btm_chl_ann", "btm_curr_mag_ann", "btm_dissic_ann", 
                                            "btm_dissoc_ann", "btm_o2_ann", "btm_sal_ann", "btm_talk_ann", "btm_temp_ann")],
                     show.variables = c("bio3",  "bio7", "bio11", "bio12"),
                     fixed.var.metric = 'mean', plot = FALSE, use.formal.names = TRUE)

gam1 <- gam(Desmophyllum.dianthus ~ s(btm_arag_ann, 2) + s(btm_chl_ann, 2) + 
              s(btm_curr_mag_ann, 2) + s(btm_dissic_ann, 2) + s(btm_dissoc_ann, 2) + 
              s(btm_o2_ann, 2) + s(btm_sal_ann, 2) + s(btm_talk_ann, 2) + s(btm_temp_ann, 2),
            data = ddia_btm, family = "binomial")

gg.rp <- ggplot(rp, aes(x = expl.val, y = pred.val, lty = pred.name)) +
  geom_line() + ylab("prob of occ") + xlab("") + 
  rp.gg.theme + 
  facet_grid(~ expl.name, scales = 'free_x')

print(gg.rp)

# r code_10.3_Generalized_Additive_Models_b4

# gamStart <- gam(VulpesVulpes~1, data=mammals_data, family=binomial)

gamStart <- gam(Desmophyllum.dianthus ~ 1, data = ddia_btm , family = binomial)

# gamModAIC <- step.gam(gamStart, biomod2:::.scope(mammals_data[1:3,c("bio3",  "bio7", "bio11", "bio12")], "s", 4), trace=F, direction = "both")        

gamModAIC <- step.gam(gamStart, biomod2:::.scope(ddia_btm[1:3,c("btm_arag_ann", "btm_chl_ann", "btm_curr_mag_ann", "btm_dissic_ann", 
                                                            "btm_dissoc_ann", "btm_o2_ann", "btm_sal_ann", "btm_talk_ann", "btm_temp_ann")], "s", 4), trace=F, direction = "both")

# r code_10.3_Generalized_Additive_Models_b4b, echo=FALSE, eval = FALSE
## test of the multiple smoother case

#biomod2:::.scope(mammals_data[1:3,c("bio3",  "bio7", "bio11", "bio12")], "s", 2:4)

biomod2:::.scope(ddia_btm[1:3,c("btm_arag_ann", "btm_chl_ann", "btm_curr_mag_ann", "btm_dissic_ann", 
                                                            "btm_dissoc_ann", "btm_o2_ann", "btm_sal_ann", "btm_talk_ann", "btm_temp_ann")], "s", 2:4)

# gamModAIC.test <- step.gam(gamStart, biomod2:::.scope(mammals_data[1:3,c("bio3",  "bio7", "bio11", "bio12")], "s", 2:4), trace=F, direction = "both")

gamModAIC.test <- step.gam(gamStart, biomod2:::.scope(ddia_btm[1:3,c("btm_arag_ann", "btm_chl_ann", "btm_curr_mag_ann", "btm_dissic_ann", 
                                                                         "btm_dissoc_ann", "btm_o2_ann", "btm_sal_ann", "btm_talk_ann", "btm_temp_ann")], "s", 2:4), trace = F, direction = "both")

# r GAMb5 10.8, opts.label = 'fig_quarter_page', fig.cap = 'Figure 10. 8. Observed (black=presence, light gray= absence) and potential distribution of Vulpes vulpes extracted from gamModAIC. The gray scale of prediction illustratesshows habitat suitability values between 0 (light, unsuitable) and 1 (dark, highly suitable)'}

par(mfrow=c(1,2))

#level.plot(mammals_data$VulpesVulpes, XY=mammals_data[,c("X_WGS84", "Y_WGS84")], color.gradient = "grey", cex=0.3,level.range=c(0,1), show.scale=F, title="Original data")

level.plot(ddia_btm$Desmophyllum.dianthus, XY=ddia_btm[,c("Latitude", "Longitude")], color.gradient = "grey", cex=0.3,level.range=c(0,1), show.scale=F, title="Original data")

#level.plot(fitted(gamModAIC), XY=mammals_data[,c("X_WGS84", "Y_WGS84")], color.gradient = "grey", cex=0.3, level.range=c(0,1), show.scale=F, title="Stepwise GAM with AIC")

level.plot(fitted(gamModAIC), XY=ddia_btm[,c("Longitude", "Latitude")], color.gradient = "grey", cex=0.3, level.range=c(0,1), show.scale=F, title="Stepwise GAM with AIC")

# r code_10.3_Generalized_Additive_Models_b6, fig.keep = FALSE

if(is.element("package:gam", search())) detach("package:gam") ## make sure the gam package is not loaded to avoid conflicts

library(mgcv)

#gam_mgcv <- gam(VulpesVulpes~s(bio3)+s(bio7)+s(bio11)+s(bio12),data = mammals_data, family = "binomial")

gam_mgcv <- gam(Desmophyllum.dianthus ~ s(btm_arag_ann) + s(btm_chl_ann) + 
                  s(btm_curr_mag_ann) + s(btm_dissic_ann) + s(btm_dissoc_ann) + 
                  s(btm_o2_ann) + s(btm_sal_ann) + s(btm_talk_ann) + s(btm_temp_ann), data = ddia_btm, family = "binomial")


## see a range of summary statistics

summary(gam_mgcv)

gam.check(gam_mgcv)


# r GAMb7 10.9, opts.label = 'fig_half_page', fig.cap = 'Figure 10.9. Response curves of model gam_mgcv plotted using the internal function of mgcv().'}

plot(gam_mgcv,pages=1, seWithMean=TRUE)  

# r GAMb8 10.10, opts.label = 'fig_quarter_page', fig.cap = 'Figure 10.10. Response curves from the model calibrated with the mgcv package (gam_mgcv).'}

# rp <- response.plot2(models = c('gam_mgcv'),
#                      Data = mammals_data[,c("bio3", "bio7", "bio11", "bio12")],
#                      show.variables = c("bio3",  "bio7", "bio11", "bio12"),
#                      fixed.var.metric = 'mean', plot = FALSE, use.formal.names = TRUE)

rp <- response.plot2(models = c('gam_mgcv'),
                     Data = mammals_data[,c("btm_arag_ann", "btm_chl_ann", "btm_curr_mag_ann", "btm_dissic_ann", 
                                            "btm_dissoc_ann", "btm_o2_ann", "btm_sal_ann", "btm_talk_ann", "btm_temp_ann")],
                     show.variables = c("btm_arag_ann", "btm_chl_ann", "btm_curr_mag_ann", "btm_dissic_ann", 
                                        "btm_dissoc_ann", "btm_o2_ann", "btm_sal_ann", "btm_talk_ann", "btm_temp_ann"),
                     fixed.var.metric = 'mean', plot = FALSE, use.formal.names = TRUE)

gg.rp <- ggplot(rp, aes(x = expl.val, y = pred.val, lty = pred.name)) +
  geom_line() + ylab("prob of occ") + xlab("") + 
  rp.gg.theme + 
  facet_grid(~ expl.name, scales = 'free_x')

print(gg.rp)

# r GAMb9 10.11, opts.label = 'fig_quarter_page', fig.cap = 'Figure 10.11. Observed (black=presence, light gray= absence) and potential distribution of Vulpes vulpes extracted from the gam_mgcv object. The gray scale of predictions illustrates shows habitat suitability values between 0 (light, unsuitable) and 1 (dark, highly suitable).'}

par(mfrow=c(1,2))

#level.plot(mammals_data$VulpesVulpes, XY=mammals_data[,c("X_WGS84", "Y_WGS84")], color.gradient = "grey", cex=0.3,level.range=c(0,1), show.scale=F, title="Original data")

level.plot(ddia_btm$Desmophyllum.dianthus, XY=ddia_btm[,c("Latitude", "Longitude")], color.gradient = "grey", cex=0.3,level.range=c(0,1), show.scale=F, title="Original data")

#level.plot(fitted(gam_mgcv), XY=mammals_data[,c("X_WGS84", "Y_WGS84")], color.gradient = "grey", cex=0.3, level.range=c(0,1), show.scale=F, title="GAM with mgcv")

level.plot(fitted(gam_mgcv), XY=ddia_btm[,c("Latitude", "Longitude")], color.gradient = "grey", cex=0.3, level.range=c(0,1), show.scale=F, title="GAM with mgcv")

par(mfrow = c(1, 1))
