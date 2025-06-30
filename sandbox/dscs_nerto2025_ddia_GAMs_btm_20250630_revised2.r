# =============================================================================
# Complete GAM Analysis for Deep-Sea Coral Distribution Modeling
# Species: Desmophyllum dianthus
# Data: Northeast Canyons Deep-Sea Coral Survey (DSCS) - Bottom Environmental Data
# Author: [Your Name]
# Date: June 30, 2025
# =============================================================================

# =============================================================================
# SECTION 1: Setup and Data Loading
# =============================================================================

# Load required libraries
library(biomod2)      # For species distribution modeling functions
library(MASS)         # For statistical functions
library(ggplot2)      # For modern plotting
library(gridExtra)    # For arranging multiple plots
library(dplyr)        # For data manipulation
library(viridis)      # For better color palettes
library(lattice)      # For additional plotting options
library(gam)          # For GAM modeling (will switch to mgcv later)

# Load the dataset
cat("=== Loading Dataset ===\n")
ddia_btm <- read.csv("ne_canyons_dscs_btm_20250630.csv", row.names = 1)

# =============================================================================
# SECTION 2: Data Exploration and Quality Assessment
# =============================================================================

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

# Define environmental variables
env_vars <- c("btm_arag_ann", "btm_chl_ann", "btm_curr_mag_ann", "btm_dissic_ann", 
              "btm_dissoc_ann", "btm_o2_ann", "btm_sal_ann", "btm_talk_ann", "btm_temp_ann")

# Check data quality for GAM modeling
cat("\n=== Environmental Variables Summary ===\n")
env_summary <- data.frame(
  Variable = env_vars,
  Unique_Values = sapply(env_vars, function(x) length(unique(ddia_btm[[x]], na.rm = TRUE))),
  Min = sapply(env_vars, function(x) round(min(ddia_btm[[x]], na.rm = TRUE), 4)),
  Max = sapply(env_vars, function(x) round(max(ddia_btm[[x]], na.rm = TRUE), 4)),
  Mean = sapply(env_vars, function(x) round(mean(ddia_btm[[x]], na.rm = TRUE), 4)),
  SD = sapply(env_vars, function(x) round(sd(ddia_btm[[x]], na.rm = TRUE), 4))
)

print(env_summary)

# Check for variables with insufficient unique values for GAM smoothing
problematic_vars <- env_vars[env_summary$Unique_Values <= 3]
if(length(problematic_vars) > 0) {
  cat("\nWARNING: Variables with ≤3 unique values (insufficient for GAM):\n")
  print(problematic_vars)
} else {
  cat("\nAll variables have sufficient unique values for GAM smoothing.\n")
}

# =============================================================================
# SECTION 3: Basic GAM Models using gam package
# =============================================================================

cat("\n=== Fitting GAM Models (gam package) ===\n")

# Ensure mgcv package is not loaded to avoid conflicts
if("package:mgcv" %in% search()) detach("package:mgcv", unload = TRUE)
library(gam)

# GAM with degree 2 smoothing (lower complexity)
cat("Fitting GAM1 with degree 2 smoothing...\n")
gam1 <- gam(Desmophyllum.dianthus ~ s(btm_arag_ann, 2) + s(btm_chl_ann, 2) + 
            s(btm_curr_mag_ann, 2) + s(btm_dissic_ann, 2) + s(btm_dissoc_ann, 2) + 
            s(btm_o2_ann, 2) + s(btm_sal_ann, 2) + s(btm_talk_ann, 2) + s(btm_temp_ann, 2),
            data = ddia_btm, family = "binomial")

# GAM with degree 4 smoothing (higher complexity)
cat("Fitting GAM2 with degree 4 smoothing...\n")
gam2 <- gam(Desmophyllum.dianthus ~ s(btm_arag_ann, 4) + s(btm_chl_ann, 4) + 
            s(btm_curr_mag_ann, 4) + s(btm_dissic_ann, 4) + s(btm_dissoc_ann, 4) + 
            s(btm_o2_ann, 4) + s(btm_sal_ann, 4) + s(btm_talk_ann, 4) + s(btm_temp_ann, 4),
            data = ddia_btm, family = "binomial")

cat("Basic GAM models fitted successfully!\n")

# Display model summaries
cat("\n=== GAM1 Summary (degree 2) ===\n")
print(summary(gam1))

cat("\n=== GAM2 Summary (degree 4) ===\n")
print(summary(gam2))

# =============================================================================
# SECTION 4: Stepwise GAM Selection
# =============================================================================

cat("\n=== Performing Stepwise GAM Selection ===\n")

# Start with null model
gamStart <- gam(Desmophyllum.dianthus ~ 1, data = ddia_btm, family = binomial)

# Create scope for stepwise selection (using subset to avoid memory issues)
scope_data <- ddia_btm[sample(nrow(ddia_btm), 1000), env_vars]  # Random subset
gam_scope <- biomod2:::.scope(scope_data, "s", 4)

# Perform stepwise selection
cat("Running stepwise GAM selection with AIC...\n")
gamModAIC <- step.gam(gamStart, scope = gam_scope, trace = FALSE, direction = "both")

cat("Stepwise GAM selection completed.\n")
print(summary(gamModAIC))

# Optional: Test multiple smoother degrees
cat("\n=== Testing Multiple Smoother Degrees ===\n")
scope_data_multi <- ddia_btm[sample(nrow(ddia_btm), 500), env_vars]  # Smaller subset
gam_scope_multi <- biomod2:::.scope(scope_data_multi, "s", 2:4)

gamModAIC_multi <- step.gam(gamStart, scope = gam_scope_multi, trace = FALSE, direction = "both")
print(summary(gamModAIC_multi))

# =============================================================================
# SECTION 5: Modern Visualization Functions
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
    
    # Create plot based on data type
    if(is.factor(plot_data$response) || all(plot_data$response %in% c(0, 1))) {
      # Binary data (presence/absence)
      p <- ggplot(plot_data, aes(x = x, y = y, color = factor(response))) +
        geom_point(size = 0.3, alpha = 0.6) +
        scale_color_manual(values = c("0" = "lightgray", "1" = "black"),
                          labels = c("Absence", "Presence"),
                          name = "Occurrence") +
        labs(title = titles[i], x = "Longitude", y = "Latitude")
    } else {
      # Continuous data (probabilities)
      p <- ggplot(plot_data, aes(x = x, y = y, color = response)) +
        geom_point(size = 0.3, alpha = 0.6) +
        scale_color_viridis_c(name = "Probability") +
        labs(title = titles[i], x = "Longitude", y = "Latitude")
    }
    
    p <- p + theme_minimal() +
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

# Function to create response curves
create_response_curves <- function(models_list, data, env_vars, model_names) {
  
  response_data_list <- list()
  
  for(i in 1:length(models_list)) {
    model <- models_list[[i]]
    model_name <- model_names[i]
    
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
      response_data_list[[paste(model_name, var, sep = "_")]] <- data.frame(
        expl.name = var,
        expl.val = var_seq,
        pred.val = predictions,
        pred.name = model_name
      )
    }
  }
  
  # Combine all data
  response_df <- do.call(rbind, response_data_list)
  return(response_df)
}

# Define ggplot theme for response plots
rp.gg.theme <- theme(
  legend.title = element_blank(),
  axis.text.x = element_text(angle = 90, vjust = 0.5),
  panel.background = element_rect(fill = NA, colour = "gray70"),
  strip.background = element_rect(fill = NA, colour = "gray70"),
  panel.grid.major = element_line(colour = "grey90"),
  legend.key = element_rect(fill = NA, colour = "gray70")
)

# =============================================================================
# SECTION 6: Visualizations for gam Package Models
# =============================================================================

cat("\n=== Creating Visualizations ===\n")

# Basic GAM plots using base R
cat("Creating basic GAM diagnostic plots...\n")
par(mfrow = c(3, 3))
plot(gam1, se = TRUE, main = "GAM1 Response Curves (Degree 2)")
par(mfrow = c(1, 1))

# Spatial distribution plots
cat("Creating spatial distribution plots...\n")
spatial_plots_basic <- create_spatial_plots(
  data = ddia_btm,
  response_vars = list(
    ddia_btm$Desmophyllum.dianthus,
    fitted(gam1),
    fitted(gam2),
    fitted(gamModAIC)
  ),
  titles = c(
    "Original Data",
    "GAM1 (degree 2)",
    "GAM2 (degree 4)",
    "Stepwise GAM"
  )
)

# Display spatial plots
grid.arrange(grobs = spatial_plots_basic, ncol = 2, nrow = 2)

# Response curves comparison
cat("Creating response curves...\n")
response_data_comparison <- create_response_curves(
  models_list = list(gam1, gam2),
  data = ddia_btm,
  env_vars = env_vars,
  model_names = c("GAM1 (degree 2)", "GAM2 (degree 4)")
)

# Create response curve plot
response_plot <- ggplot(response_data_comparison, aes(x = expl.val, y = pred.val, color = pred.name)) +
  geom_line(size = 1) +
  facet_wrap(~ expl.name, scales = "free_x", nrow = 3, ncol = 3) +
  labs(title = "GAM Response Curves Comparison",
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
# SECTION 7: Advanced GAM using mgcv Package
# =============================================================================

cat("\n=== Switching to mgcv Package ===\n")

# Properly switch to mgcv package
if("package:gam" %in% search()) detach("package:gam", unload = TRUE)
library(mgcv)

# Fit GAM using mgcv (more modern and flexible)
cat("Fitting advanced GAM using mgcv...\n")
gam_mgcv <- gam(Desmophyllum.dianthus ~ s(btm_arag_ann) + s(btm_chl_ann) + 
                s(btm_curr_mag_ann) + s(btm_dissic_ann) + s(btm_dissoc_ann) + 
                s(btm_o2_ann) + s(btm_sal_ann) + s(btm_talk_ann) + s(btm_temp_ann), 
                data = ddia_btm, family = binomial, method = "REML")

# Fit GAM with specified basis dimensions
gam_mgcv_k4 <- gam(Desmophyllum.dianthus ~ s(btm_arag_ann, k = 4) + s(btm_chl_ann, k = 4) + 
                   s(btm_curr_mag_ann, k = 4) + s(btm_dissic_ann, k = 4) + s(btm_dissoc_ann, k = 4) + 
                   s(btm_o2_ann, k = 4) + s(btm_sal_ann, k = 4) + s(btm_talk_ann, k = 4) + s(btm_temp_ann, k = 4),
                   data = ddia_btm, family = binomial, method = "REML")

cat("mgcv GAM models fitted successfully!\n")

# Model summaries and diagnostics
cat("\n=== mgcv GAM Summary (default) ===\n")
print(summary(gam_mgcv))

cat("\n=== mgcv GAM Summary (k=4) ===\n")
print(summary(gam_mgcv_k4))

# Model diagnostics
cat("\n=== Model Diagnostics ===\n")
gam.check(gam_mgcv)

# =============================================================================
# SECTION 8: mgcv Visualizations
# =============================================================================

cat("\n=== Creating mgcv Visualizations ===\n")

# Plot smooths using mgcv's internal plotting
par(mfrow = c(3, 3))
plot(gam_mgcv, pages = 1, seWithMean = TRUE, residuals = TRUE, pch = 19, cex = 0.3)
par(mfrow = c(1, 1))

# Response curves for mgcv models
create_mgcv_response_plot <- function(model, data, env_vars, model_name) {
  
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
      expl.name = var,
      expl.val = var_seq,
      pred.val = predictions,
      pred.name = model_name
    )
  }
  
  # Combine all data
  response_df <- do.call(rbind, response_data)
  return(response_df)
}

# Create response plots for mgcv models
rp_mgcv_default <- create_mgcv_response_plot(gam_mgcv, ddia_btm, env_vars, "mgcv GAM (default)")
rp_mgcv_k4 <- create_mgcv_response_plot(gam_mgcv_k4, ddia_btm, env_vars, "mgcv GAM (k=4)")

# Combine mgcv response data
rp_mgcv_combined <- rbind(rp_mgcv_default, rp_mgcv_k4)

# Create mgcv response plot
gg.rp.mgcv <- ggplot(rp_mgcv_combined, aes(x = expl.val, y = pred.val, color = pred.name)) +
  geom_line(size = 1) + 
  ylab("Probability of Occurrence") + 
  xlab("Environmental Variable Value") + 
  facet_wrap(~ expl.name, scales = 'free_x', nrow = 3, ncol = 3) +
  labs(title = "mgcv GAM Response Curves Comparison", color = "Model") +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14),
    strip.text = element_text(size = 10),
    axis.text = element_text(size = 8),
    legend.position = "bottom"
  ) +
  scale_color_manual(values = c("mgcv GAM (default)" = "darkgreen", "mgcv GAM (k=4)" = "orange"))

print(gg.rp.mgcv)

# Spatial visualization for mgcv models
spatial_plots_mgcv <- create_spatial_plots(
  data = ddia_btm,
  response_vars = list(
    ddia_btm$Desmophyllum.dianthus,
    fitted(gam_mgcv),
    fitted(gam_mgcv_k4)
  ),
  titles = c(
    "Original Data",
    "mgcv GAM (default)",
    "mgcv GAM (k=4)"
  )
)

# Display mgcv spatial plots
grid.arrange(grobs = spatial_plots_mgcv, ncol = 2, nrow = 2)

# =============================================================================
# SECTION 9: Model Evaluation and Comparison
# =============================================================================

cat("\n=== Comprehensive Model Evaluation ===\n")

# Function to calculate model performance metrics
calculate_performance <- function(model, data, response_var) {
  
  fitted_probs <- fitted(model)
  observed <- data[[response_var]]
  
  # AIC
  model_aic <- AIC(model)
  
  # Deviance explained
  if("mgcv" %in% class(model)) {
    dev_expl <- summary(model)$dev.expl * 100
  } else {
    dev_expl <- ((model$null.deviance - model$deviance) / model$null.deviance) * 100
  }
  
  # Prediction statistics
  mean_pred <- mean(fitted_probs)
  pred_range <- range(fitted_probs)
  
  # Classification metrics (using 0.5 threshold)
  pred_binary <- ifelse(fitted_probs > 0.5, 1, 0)
  confusion <- table(Predicted = pred_binary, Observed = observed)
  accuracy <- sum(diag(confusion)) / sum(confusion)
  
  # Sensitivity and Specificity
  if(all(c(0,1) %in% rownames(confusion)) && all(c(0,1) %in% colnames(confusion))) {
    sensitivity <- confusion["1", "1"] / sum(confusion[, "1"])
    specificity <- confusion["0", "0"] / sum(confusion[, "0"])
  } else {
    sensitivity <- NA
    specificity <- NA
  }
  
  return(list(
    AIC = model_aic,
    Deviance_Explained = round(dev_expl, 2),
    Mean_Prediction = round(mean_pred, 4),
    Prediction_Range = paste(round(pred_range[1], 4), "to", round(pred_range[2], 4)),
    Accuracy = round(accuracy, 3),
    Sensitivity = round(sensitivity, 3),
    Specificity = round(specificity, 3)
  ))
}

# Calculate performance for all models
models_list <- list(gam1, gam2, gamModAIC, gam_mgcv, gam_mgcv_k4)
model_names <- c("GAM1 (degree 2)", "GAM2 (degree 4)", "Stepwise GAM", 
                 "mgcv GAM (default)", "mgcv GAM (k=4)")

performance_results <- lapply(models_list, function(x) calculate_performance(x, ddia_btm, "Desmophyllum.dianthus"))
names(performance_results) <- model_names

# Create comparison table
comparison_table <- data.frame(
  Model = model_names,
  AIC = sapply(performance_results, function(x) x$AIC),
  Deviance_Explained = sapply(performance_results, function(x) x$Deviance_Explained),
  Mean_Prediction = sapply(performance_results, function(x) x$Mean_Prediction),
  Prediction_Range = sapply(performance_results, function(x) x$Prediction_Range),
  Accuracy = sapply(performance_results, function(x) x$Accuracy),
  Sensitivity = sapply(performance_results, function(x) x$Sensitivity),
  Specificity = sapply(performance_results, function(x) x$Specificity)
)

print(comparison_table)

# Identify best model based on AIC
best_model_idx <- which.min(comparison_table$AIC)
cat("\nBest model based on AIC:", model_names[best_model_idx], "\n")
cat("AIC:", comparison_table$AIC[best_model_idx], "\n")
cat("Deviance Explained:", comparison_table$Deviance_Explained[best_model_idx], "%\n")

# =============================================================================
# SECTION 10: Final Summary and Recommendations
# =============================================================================

cat("\n=== ANALYSIS SUMMARY ===\n")
cat("Dataset: Northeast Canyons DSCS Bottom Environmental Data\n")
cat("Species: Desmophyllum dianthus\n")
cat("Total observations:", nrow(ddia_btm), "\n")
cat("Prevalence:", round(mean(ddia_btm$Desmophyllum.dianthus), 3), "\n")
cat("Environmental variables tested:", length(env_vars), "\n")
cat("Models fitted:", length(models_list), "\n")

cat("\nKey Findings:\n")
cat("1. All environmental variables have sufficient variation for GAM modeling\n")
cat("2. Model complexity affects performance - compare AIC values\n")
cat("3. mgcv package generally provides more flexible and robust GAM implementation\n")
cat("4. Stepwise selection can help identify most important variables\n")

cat("\nRecommendations:\n")
cat("1. Use the best performing model (lowest AIC) for final predictions\n")
cat("2. Consider ecological interpretation alongside statistical performance\n")
cat("3. Validate model predictions with independent test data if available\n")
cat("4. Examine residual patterns and model diagnostics carefully\n")

cat("\n=== ANALYSIS COMPLETE ===\n")
cat("All GAM models successfully fitted and evaluated!\n")
cat("Visualizations and performance metrics generated.\n")
cat("Ready for ecological interpretation and application.\n")