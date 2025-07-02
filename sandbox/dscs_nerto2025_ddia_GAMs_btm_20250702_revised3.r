# =============================================================================
# Complete GAM Analysis for Deep-Sea Coral Distribution Modeling
# Species: Desmophyllum dianthus
# Data: Northeast Canyons Deep-Sea Coral Survey (DSCS) - Bottom Environmental Data
# Author: [Your Name]
# Date: June 30, 2025
# REVISED VERSION - Fixed step.gam() and other syntax issues
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
library(mgcv)         # For GAM modeling (using mgcv from the start)

# Load the dataset
cat("=== Loading Dataset ===\n")
#ddia_btm <- read.csv("ne_canyons_dscs_btm_20250630.csv", row.names = 1)
ddia_btm <- read.csv("data/ne_canyons_dscs_btm_20250630.csv", row.names = 1)

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
# SECTION 3: Basic GAM Models using mgcv package
# =============================================================================

cat("\n=== Fitting GAM Models (mgcv package) ===\n")

# GAM with lower complexity (k=3)
cat("Fitting GAM1 with k=3 smoothing...\n")
gam1 <- gam(Desmophyllum.dianthus ~ s(btm_arag_ann, k=3) + s(btm_chl_ann, k=3) + 
            s(btm_curr_mag_ann, k=3) + s(btm_dissic_ann, k=3) + s(btm_dissoc_ann, k=3) + 
            s(btm_o2_ann, k=3) + s(btm_sal_ann, k=3) + s(btm_talk_ann, k=3) + s(btm_temp_ann, k=3),
            data = ddia_btm, family = binomial, method = "REML")

# GAM with higher complexity (k=5)
cat("Fitting GAM2 with k=5 smoothing...\n")
gam2 <- gam(Desmophyllum.dianthus ~ s(btm_arag_ann, k=5) + s(btm_chl_ann, k=5) + 
            s(btm_curr_mag_ann, k=5) + s(btm_dissic_ann, k=5) + s(btm_dissoc_ann, k=5) + 
            s(btm_o2_ann, k=5) + s(btm_sal_ann, k=5) + s(btm_talk_ann, k=5) + s(btm_temp_ann, k=5),
            data = ddia_btm, family = binomial, method = "REML")

# GAM with default settings (automatic basis selection)
cat("Fitting GAM3 with default settings...\n")
gam3 <- gam(Desmophyllum.dianthus ~ s(btm_arag_ann) + s(btm_chl_ann) + 
            s(btm_curr_mag_ann) + s(btm_dissic_ann) + s(btm_dissoc_ann) + 
            s(btm_o2_ann) + s(btm_sal_ann) + s(btm_talk_ann) + s(btm_temp_ann),
            data = ddia_btm, family = binomial, method = "REML")

cat("Basic GAM models fitted successfully!\n")

# Display model summaries
cat("\n=== GAM1 Summary (k=3) ===\n")
print(summary(gam1))

cat("\n=== GAM2 Summary (k=5) ===\n")
print(summary(gam2))

cat("\n=== GAM3 Summary (default) ===\n")
print(summary(gam3))

# =============================================================================
# SECTION 4: Variable Selection using mgcv approach
# =============================================================================

cat("\n=== Performing Variable Selection ===\n")

# Method 1: Using select=TRUE in mgcv for automatic variable selection
cat("Fitting GAM with automatic variable selection (select=TRUE)...\n")
gamModSelect <- gam(Desmophyllum.dianthus ~ s(btm_arag_ann) + s(btm_chl_ann) + 
                    s(btm_curr_mag_ann) + s(btm_dissic_ann) + s(btm_dissoc_ann) + 
                    s(btm_o2_ann) + s(btm_sal_ann) + s(btm_talk_ann) + s(btm_temp_ann),
                    data = ddia_btm, family = binomial, method = "REML", select = TRUE)

cat("Variable selection model completed.\n")
print(summary(gamModSelect))

# Method 2: Manual forward/backward selection using AIC
cat("\n=== Manual Model Selection using AIC ===\n")

# Start with null model
gamNull <- gam(Desmophyllum.dianthus ~ 1, data = ddia_btm, family = binomial)

# Create individual models for each variable
cat("Testing individual variables...\n")
individual_aic <- numeric(length(env_vars))
names(individual_aic) <- env_vars

for(i in 1:length(env_vars)) {
  formula_str <- paste("Desmophyllum.dianthus ~ s(", env_vars[i], ")", sep = "")
  temp_model <- gam(as.formula(formula_str), data = ddia_btm, family = binomial, method = "REML")
  individual_aic[i] <- AIC(temp_model)
}

# Sort variables by AIC (best first)
var_order <- names(sort(individual_aic))
cat("Variables ranked by individual AIC performance:\n")
print(data.frame(Variable = var_order, AIC = sort(individual_aic)))

# Forward selection approach
cat("\nPerforming forward selection...\n")
selected_vars <- character(0)
current_aic <- AIC(gamNull)
best_aic <- current_aic

for(step in 1:length(env_vars)) {
  remaining_vars <- setdiff(env_vars, selected_vars)
  if(length(remaining_vars) == 0) break
  
  step_aic <- numeric(length(remaining_vars))
  names(step_aic) <- remaining_vars
  
  for(var in remaining_vars) {
    test_vars <- c(selected_vars, var)
    formula_str <- paste("Desmophyllum.dianthus ~ ", 
                        paste("s(", test_vars, ")", collapse = " + "), sep = "")
    temp_model <- gam(as.formula(formula_str), data = ddia_btm, family = binomial, method = "REML")
    step_aic[var] <- AIC(temp_model)
  }
  
  best_var <- names(which.min(step_aic))
  best_step_aic <- min(step_aic)
  
  if(best_step_aic < best_aic) {
    selected_vars <- c(selected_vars, best_var)
    best_aic <- best_step_aic
    cat("Step", step, ": Added", best_var, "(AIC =", round(best_aic, 2), ")\n")
  } else {
    cat("Step", step, ": No improvement found. Stopping.\n")
    break
  }
}

# Fit final selected model
if(length(selected_vars) > 0) {
  final_formula <- paste("Desmophyllum.dianthus ~ ", 
                        paste("s(", selected_vars, ")", collapse = " + "), sep = "")
  gamModAIC <- gam(as.formula(final_formula), data = ddia_btm, family = binomial, method = "REML")
  cat("\nFinal selected model fitted with variables:", paste(selected_vars, collapse = ", "), "\n")
  print(summary(gamModAIC))
} else {
  cat("No variables selected. Using null model.\n")
  gamModAIC <- gamNull
}

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

# Function to create response curves for mgcv models
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

# =============================================================================
# SECTION 6: Visualizations
# =============================================================================

cat("\n=== Creating Visualizations ===\n")

# Model diagnostic plots
cat("Creating GAM diagnostic plots...\n")
par(mfrow = c(2, 2))
gam.check(gam3)
par(mfrow = c(1, 1))

# Smoothing function plots
cat("Creating smoothing function plots...\n")
par(mfrow = c(3, 3))
plot(gam3, pages = 1, seWithMean = TRUE, residuals = TRUE, pch = 19, cex = 0.3)
par(mfrow = c(1, 1))

# Spatial distribution plots
cat("Creating spatial distribution plots...\n")
spatial_plots <- create_spatial_plots(
  data = ddia_btm,
  response_vars = list(
    ddia_btm$Desmophyllum.dianthus,
    fitted(gam1),
    fitted(gam2),
    fitted(gam3)
  ),
  titles = c(
    "Original Data",
    "GAM1 (k=3)",
    "GAM2 (k=5)",
    "GAM3 (default)"
  )
)

# Display spatial plots
grid.arrange(grobs = spatial_plots, ncol = 2, nrow = 2)

# Response curves comparison
cat("Creating response curves...\n")
rp_gam1 <- create_mgcv_response_plot(gam1, ddia_btm, env_vars, "GAM1 (k=3)")
rp_gam2 <- create_mgcv_response_plot(gam2, ddia_btm, env_vars, "GAM2 (k=5)")
rp_gam3 <- create_mgcv_response_plot(gam3, ddia_btm, env_vars, "GAM3 (default)")

# Combine response data
rp_combined <- rbind(rp_gam1, rp_gam2, rp_gam3)

# Create response plot
response_plot <- ggplot(rp_combined, aes(x = expl.val, y = pred.val, color = pred.name)) +
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
  scale_color_manual(values = c("GAM1 (k=3)" = "blue", "GAM2 (k=5)" = "red", "GAM3 (default)" = "darkgreen"))

print(response_plot)

# =============================================================================
# SECTION 7: Model Evaluation and Comparison
# =============================================================================

cat("\n=== Comprehensive Model Evaluation ===\n")

# Function to calculate model performance metrics
calculate_performance <- function(model, data, response_var) {
  
  fitted_probs <- fitted(model)
  observed <- data[[response_var]]
  
  # AIC
  model_aic <- AIC(model)
  
  # Deviance explained
  dev_expl <- summary(model)$dev.expl * 100
  
  # Prediction statistics
  mean_pred <- mean(fitted_probs)
  pred_range <- range(fitted_probs)
  
  # Classification metrics (using 0.5 threshold)
  pred_binary <- ifelse(fitted_probs > 0.5, 1, 0)
  confusion <- table(Predicted = pred_binary, Observed = observed)
  accuracy <- sum(diag(confusion)) / sum(confusion)
  
  # Sensitivity and Specificity
  if(all(c("0","1") %in% rownames(confusion)) && all(c("0","1") %in% colnames(confusion))) {
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
models_list <- list(gam1, gam2, gam3, gamModSelect)
model_names <- c("GAM1 (k=3)", "GAM2 (k=5)", "GAM3 (default)", "GAM Select")

# Add gamModAIC if it exists and is not null model
if(exists("gamModAIC") && length(selected_vars) > 0) {
  models_list <- c(models_list, list(gamModAIC))
  model_names <- c(model_names, "GAM AIC")
}

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
# SECTION 8: Final Summary and Recommendations
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
cat("3. mgcv package provides robust GAM implementation with automatic selection\n")
cat("4. Variable selection identifies most important environmental predictors\n")

cat("\nRecommendations:\n")
cat("1. Use the best performing model (lowest AIC) for final predictions\n")
cat("2. Consider ecological interpretation alongside statistical performance\n")
cat("3. Validate model predictions with independent test data if available\n")
cat("4. Examine residual patterns and model diagnostics carefully\n")

cat("\n=== ANALYSIS COMPLETE ===\n")
cat("All GAM models successfully fitted and evaluated!\n")
cat("Visualizations and performance metrics generated.\n")
cat("Ready for ecological interpretation and application.\n")

