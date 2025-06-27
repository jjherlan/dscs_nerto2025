# Modern GLM Analysis for Species Distribution Modeling
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
mammals_data <- read.csv("data/tabular/species/mammals_and_bioclim_table.csv", row.names = 1)

# Data exploration and summary
cat("=== Dataset Summary ===\n")
cat("Total observations:", nrow(mammals_data), "\n")
cat("Red fox presence records:", sum(mammals_data$VulpesVulpes == 1), "\n")
cat("Red fox absence records:", sum(mammals_data$VulpesVulpes == 0), "\n")
cat("Prevalence:", round(mean(mammals_data$VulpesVulpes), 3), "\n\n")

# Check for missing values
missing_summary <- colSums(is.na(mammals_data))
if(any(missing_summary > 0)) {
  cat("Missing values found:\n")
  print(missing_summary[missing_summary > 0])
} else {
  cat("No missing values detected.\n")
}

# =============================================================================
# CHAPTER 10.2: Generalized Linear Models (GLM) - Base Models
# =============================================================================

# Model 1: GLM with linear terms
glm1 <- glm(VulpesVulpes ~ bio3 + bio7 + bio11 + bio12, 
            data = mammals_data, 
            family = "binomial")

# Model 2: GLM with quadratic terms
glm2 <- glm(VulpesVulpes ~ poly(bio3,2) + poly(bio7,2) + poly(bio11,2) + poly(bio12,2), 
            data = mammals_data, 
            family = "binomial")

# Model 3: GLM with interactions
glm3 <- glm(VulpesVulpes ~ poly(bio3,2) + poly(bio7,2) + poly(bio11,2) + poly(bio12,2) + 
            bio3:bio7 + bio3:bio11 + bio3:bio12 + bio7:bio11 + bio7:bio12 + bio11:bio12,
            data = mammals_data, 
            family = "binomial")

# =============================================================================
# STEPWISE MODEL SELECTION (Integrating lines 365-422)
# =============================================================================

# Starting model for stepwise selection
glmStart <- glm(VulpesVulpes ~ 1, data = mammals_data, family = binomial)

# Define the full formula for stepwise selection
glm.formula <- formula(VulpesVulpes ~ poly(bio3,2) + poly(bio7,2) + poly(bio11,2) + poly(bio12,2) + 
                      bio3:bio7 + bio3:bio11 + bio3:bio12 + bio7:bio11 + bio7:bio12 + bio11:bio12)

cat("=== Stepwise Model Selection ===\n")
cat("Full formula for selection:\n")
print(glm.formula)

# Stepwise selection with AIC
cat("\nPerforming stepwise selection with AIC...\n")
glmModAIC <- stepAIC(glmStart, 
                     glm.formula,
                     data = mammals_data,
                     direction = "both", 
                     trace = FALSE, 
                     k = 2, 
                     control = glm.control(maxit = 100))

# Stepwise selection with BIC
cat("Performing stepwise selection with BIC...\n")
glmModBIC <- stepAIC(glmStart, 
                     glm.formula, 
                     direction = "both", 
                     trace = FALSE,
                     k = log(nrow(mammals_data)),
                     control = glm.control(maxit = 100))

# Model summaries including stepwise models
cat("\n=== Extended Model Summaries ===\n")
cat("GLM1 (linear) AIC:", AIC(glm1), "\n")
cat("GLM2 (quadratic) AIC:", AIC(glm2), "\n")
cat("GLM3 (quadratic + interactions) AIC:", AIC(glm3), "\n")
cat("GLM Stepwise AIC AIC:", AIC(glmModAIC), "\n")
cat("GLM Stepwise BIC AIC:", AIC(glmModBIC), "\n\n")

# Display selected variables
cat("Variables selected by AIC stepwise:\n")
print(names(coef(glmModAIC)))
cat("\nVariables selected by BIC stepwise:\n")
print(names(coef(glmModBIC)))

# =============================================================================
# MODERN PLOTTING FUNCTIONS (replacing deprecated level.plot)
# =============================================================================

# Function to create distribution plots (replaces level.plot)
modern_distribution_plot <- function(values, coords, title, 
                                   plot_type = "continuous", 
                                   point_size = 0.3) {
  
  # Prepare data
  plot_data <- data.frame(
    X = coords[,"X_WGS84"],
    Y = coords[,"Y_WGS84"],
    Values = values
  )
  
  if(plot_type == "binary") {
    # For binary data (presence/absence)
    p <- ggplot(plot_data, aes(x = X, y = Y)) +
      geom_point(aes(color = factor(Values)), size = point_size, alpha = 0.7) +
      scale_color_manual(values = c("0" = "lightgray", "1" = "black"),
                        labels = c("Absence", "Presence"),
                        name = "") +
      theme_minimal()
  } else {
    # For continuous data (probabilities)
    p <- ggplot(plot_data, aes(x = X, y = Y)) +
      geom_point(aes(color = Values), size = point_size, alpha = 0.7) +
      scale_color_gradient(low = "lightgray", high = "black",
                          name = "Probability") +
      theme_minimal()
  }
  
  p <- p +
    labs(title = title,
         x = "Longitude (WGS84)",
         y = "Latitude (WGS84)") +
    theme(
      plot.title = element_text(hjust = 0.5, size = 11),
      axis.text = element_text(size = 8),
      legend.position = "right",
      legend.title = element_text(size = 9),
      legend.text = element_text(size = 8)
    ) +
    coord_fixed(ratio = 1)
  
  return(p)
}

# =============================================================================
# CREATE DISTRIBUTION PLOTS - Extended to include stepwise models
# =============================================================================

# Generate plots for all models including stepwise
plot1 <- modern_distribution_plot(
  values = mammals_data$VulpesVulpes,
  coords = mammals_data[,c("X_WGS84", "Y_WGS84")],
  title = "Original Data",
  plot_type = "binary"
)

plot2 <- modern_distribution_plot(
  values = fitted(glm1),
  coords = mammals_data[,c("X_WGS84", "Y_WGS84")],
  title = "GLM Linear",
  plot_type = "continuous"
)

plot3 <- modern_distribution_plot(
  values = fitted(glm2),
  coords = mammals_data[,c("X_WGS84", "Y_WGS84")],
  title = "GLM Quadratic",
  plot_type = "continuous"
)

plot4 <- modern_distribution_plot(
  values = fitted(glmModAIC),
  coords = mammals_data[,c("X_WGS84", "Y_WGS84")],
  title = "Stepwise GLM (AIC)",
  plot_type = "continuous"
)

# Display first set of plots
cat("\n=== Displaying Distribution Plots (Set 1) ===\n")
grid.arrange(plot1, plot2, plot3, plot4, ncol = 2, nrow = 2)

# Additional plot comparing stepwise models
plot5 <- modern_distribution_plot(
  values = fitted(glmModBIC),
  coords = mammals_data[,c("X_WGS84", "Y_WGS84")],
  title = "Stepwise GLM (BIC)",
  plot_type = "continuous"
)

plot6 <- modern_distribution_plot(
  values = fitted(glm3),
  coords = mammals_data[,c("X_WGS84", "Y_WGS84")],
  title = "GLM Interactions",
  plot_type = "continuous"
)

# Display stepwise comparison
cat("\n=== Displaying Stepwise Model Comparison ===\n")
grid.arrange(plot1, plot4, plot5, plot6, ncol = 2, nrow = 2)

# =============================================================================
# RESPONSE CURVES (Modern approach replacing response.plot2)
# =============================================================================

# Create response curves manually since response.plot2 may be deprecated
create_response_curves <- function(model, data, var_names, model_name) {
  
  response_data <- data.frame()
  
  for(var in var_names) {
    # Create sequence for the variable
    var_seq <- seq(min(data[[var]], na.rm = TRUE), 
                   max(data[[var]], na.rm = TRUE), 
                   length.out = 100)
    
    # Create prediction data with other variables at mean
    pred_data <- data
    for(other_var in var_names) {
      if(other_var != var) {
        pred_data[[other_var]] <- mean(data[[other_var]], na.rm = TRUE)
      }
    }
    
    # Replace the focal variable with the sequence
    pred_data_seq <- pred_data[rep(1, length(var_seq)), ]
    pred_data_seq[[var]] <- var_seq
    
    # Make predictions with error handling
    tryCatch({
      predictions <- predict(model, pred_data_seq, type = "response")
      
      # Store results
      temp_data <- data.frame(
        Variable = var,
        Value = var_seq,
        Prediction = predictions,
        Model = model_name
      )
      
      response_data <- rbind(response_data, temp_data)
    }, error = function(e) {
      cat("Warning: Could not generate response curve for", var, "in model", model_name, "\n")
    })
  }
  
  return(response_data)
}

# Generate response curves for all models
var_names <- c("bio3", "bio7", "bio11", "bio12")

cat("\n=== Generating Response Curves ===\n")
response_glm1 <- create_response_curves(glm1, mammals_data, var_names, "GLM Linear")
response_glm2 <- create_response_curves(glm2, mammals_data, var_names, "GLM Quadratic")
response_glmAIC <- create_response_curves(glmModAIC, mammals_data, var_names, "Stepwise AIC")
response_glmBIC <- create_response_curves(glmModBIC, mammals_data, var_names, "Stepwise BIC")

# Combine all response data
all_responses <- rbind(response_glm1, response_glm2, response_glmAIC, response_glmBIC)

# Define custom ggplot theme
rp.gg.theme <- theme(
  legend.title = element_blank(),
  axis.text.x = element_text(angle = 90, vjust = .5),
  panel.background = element_rect(fill = NA, colour = "gray70"),
  strip.background = element_rect(fill = NA, colour = "gray70"),
  panel.grid.major = element_line(colour = "grey90"),
  legend.key = element_rect(fill = NA, colour = "gray70")
)

# Create response curve plot
response_plot <- ggplot(all_responses, aes(x = Value, y = Prediction, color = Model)) +
  geom_line(size = 1) +
  facet_wrap(~ Variable, scales = "free_x", ncol = 2) +
  labs(title = "Response Curves for All GLM Models",
       x = "Environmental Variable Value",
       y = "Probability of Occurrence") +
  rp.gg.theme +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14),
    strip.text = element_text(size = 10),
    legend.position = "bottom"
  ) +
  scale_color_viridis_d(name = "Model Type")

print(response_plot)

# =============================================================================
# BIVARIATE RESPONSE CURVES (Modern approach)
# =============================================================================

# Function to create bivariate response surfaces
create_bivariate_response <- function(model, data, var_names, model_name) {
  
  bivariate_data <- data.frame()
  
  # Create all combinations of variables
  var_combinations <- combn(var_names, 2, simplify = FALSE)
  
  for(combo in var_combinations) {
    var1 <- combo[1]
    var2 <- combo[2]
    
    # Create sequences for both variables
    var1_seq <- seq(min(data[[var1]], na.rm = TRUE), 
                    max(data[[var1]], na.rm = TRUE), 
                    length.out = 20)  # Reduced for computational efficiency
    var2_seq <- seq(min(data[[var2]], na.rm = TRUE), 
                    max(data[[var2]], na.rm = TRUE), 
                    length.out = 20)
    
    # Create grid
    grid_data <- expand.grid(var1_val = var1_seq, var2_val = var2_seq)
    
    # Create prediction data with other variables at median
    pred_data <- data[rep(1, nrow(grid_data)), ]
    pred_data[[var1]] <- grid_data$var1_val
    pred_data[[var2]] <- grid_data$var2_val
    
    # Set other variables to median
    for(other_var in var_names) {
      if(!other_var %in% c(var1, var2)) {
        pred_data[[other_var]] <- median(data[[other_var]], na.rm = TRUE)
      }
    }
    
    # Make predictions
    tryCatch({
      predictions <- predict(model, pred_data, type = "response")
      
      # Store results
      temp_data <- data.frame(
        expl1.name = var1,
        expl2.name = var2,
        expl1.val = grid_data$var1_val,
        expl2.val = grid_data$var2_val,
        pred.val = predictions,
        Model = model_name
      )
      
      bivariate_data <- rbind(bivariate_data, temp_data)
    }, error = function(e) {
      cat("Warning: Could not generate bivariate response for", var1, "vs", var2, "in model", model_name, "\n")
    })
  }
  
  return(bivariate_data)
}

# Generate bivariate response for the best stepwise model (AIC)
cat("\n=== Generating Bivariate Response Surfaces ===\n")
bivariate_response <- create_bivariate_response(glmModAIC, mammals_data, var_names, "Stepwise AIC")

# Create bivariate response plot
if(nrow(bivariate_response) > 0) {
  bivariate_plot <- ggplot(bivariate_response, aes(x = expl1.val, y = expl2.val, fill = pred.val)) +
    geom_raster() + 
    rp.gg.theme + 
    ylab("") + xlab("") + 
    theme(legend.title = element_text()) +
    scale_fill_gradient(name = "Prob of Occurrence", low = "#f0f0f0", high = "#000000") + 
    facet_grid(expl2.name ~ expl1.name, scales = 'free') +
    labs(title = "Bivariate Response Surfaces (Stepwise AIC Model)")
  
  print(bivariate_plot)
} else {
  cat("No bivariate response data generated.\n")
}

# =============================================================================
# MODEL EVALUATION AND COMPARISON
# =============================================================================

# Check if pROC package is available
if (!requireNamespace("pROC", quietly = TRUE)) {
  cat("pROC package not found. Installing it now...\n")
  install.packages("pROC")
}
library(pROC)

# Function to calculate model metrics with error handling
evaluate_model <- function(model, data) {
  tryCatch({
    predictions <- predict(model, type = "response")
    
    # AUC calculation
    roc_obj <- roc(data$VulpesVulpes, predictions, quiet = TRUE)
    auc_value <- as.numeric(auc(roc_obj))
    
    # Deviance explained (pseudo R-squared)
    null_deviance <- model$null.deviance
    residual_deviance <- model$deviance
    dev_explained <- (null_deviance - residual_deviance) / null_deviance
    
    return(list(
      AIC = AIC(model),
      AUC = auc_value,
      DevExplained = dev_explained,
      df = model$df.residual,
      Success = TRUE
    ))
  }, error = function(e) {
    cat("Error evaluating model:", e$message, "\n")
    return(list(
      AIC = AIC(model),
      AUC = NA,
      DevExplained = NA,
      df = model$df.residual,
      Success = FALSE
    ))
  })
}

# Evaluate all models
eval_glm1 <- evaluate_model(glm1, mammals_data)
eval_glm2 <- evaluate_model(glm2, mammals_data)
eval_glm3 <- evaluate_model(glm3, mammals_data)
eval_glmAIC <- evaluate_model(glmModAIC, mammals_data)
eval_glmBIC <- evaluate_model(glmModBIC, mammals_data)

# Create comprehensive comparison table
comparison_table <- data.frame(
  Model = c("GLM Linear", "GLM Quadratic", "GLM Interactions", "Stepwise AIC", "Stepwise BIC"),
  AIC = c(eval_glm1$AIC, eval_glm2$AIC, eval_glm3$AIC, eval_glmAIC$AIC, eval_glmBIC$AIC),
  AUC = c(eval_glm1$AUC, eval_glm2$AUC, eval_glm3$AUC, eval_glmAIC$AUC, eval_glmBIC$AUC),
  DevExplained = c(eval_glm1$DevExplained, eval_glm2$DevExplained, eval_glm3$DevExplained, 
                   eval_glmAIC$DevExplained, eval_glmBIC$DevExplained),
  df = c(eval_glm1$df, eval_glm2$df, eval_glm3$df, eval_glmAIC$df, eval_glmBIC$df),
  stringsAsFactors = FALSE
)

# Display results
cat("\n=== Comprehensive Model Comparison ===\n")

# Create properly formatted table for display
comparison_display <- comparison_table
comparison_display$AIC <- round(comparison_display$AIC, 2)
comparison_display$AUC <- round(comparison_display$AUC, 4)
comparison_display$DevExplained <- round(comparison_display$DevExplained, 4)

print(comparison_display, row.names = FALSE)

# Determine best models
best_aic_idx <- which.min(comparison_table$AIC)
best_auc_idx <- which.max(comparison_table$AUC)

cat("\nBest model based on AIC:", comparison_table$Model[best_aic_idx], "\n")
cat("Best model based on AUC:", comparison_table$Model[best_auc_idx], "\n")

# =============================================================================
# ANOVA AND MODEL DIAGNOSTICS
# =============================================================================

cat("\n=== ANOVA for Stepwise AIC Model ===\n")
print(anova(glmModAIC, test = "Chisq"))

cat("\n=== Model Formulas ===\n")
cat("Stepwise AIC model formula:\n")
print(formula(glmModAIC))
cat("\nStepwise BIC model formula:\n")
print(formula(glmModBIC))

# Additional model diagnostics
cat("\n=== Additional Model Diagnostics ===\n")
cat("GLM Linear - Residual Deviance:", round(glm1$deviance, 2), 
    "on", glm1$df.residual, "degrees of freedom\n")
cat("GLM Quadratic - Residual Deviance:", round(glm2$deviance, 2), 
    "on", glm2$df.residual, "degrees of freedom\n")
cat("GLM Interactions - Residual Deviance:", round(glm3$deviance, 2), 
    "on", glm3$df.residual, "degrees of freedom\n")
cat("Stepwise AIC - Residual Deviance:", round(glmModAIC$deviance, 2), 
    "on", glmModAIC$df.residual, "degrees of freedom\n")
cat("Stepwise BIC - Residual Deviance:", round(glmModBIC$deviance, 2), 
    "on", glmModBIC$df.residual, "degrees of freedom\n")

# =============================================================================
# ALTERNATIVE: Base R plotting (for comparison)
# =============================================================================

# Function that mimics level.plot behavior
base_distribution_plot <- function(values, coords, title, is_binary = FALSE) {
  
  if(is_binary) {
    colors <- c("lightgray", "black")[values + 1]
  } else {
    # Create grayscale based on values
    gray_levels <- gray(1 - values)  # Invert so high values are dark
    colors <- gray_levels
  }
  
  plot(coords[,"X_WGS84"], coords[,"Y_WGS84"],
       col = colors,
       pch = 19, cex = 0.3,
       main = title,
       xlab = "Longitude (WGS84)",
       ylab = "Latitude (WGS84)",
       asp = 1)
  
  if(is_binary) {
    legend("topright", 
           legend = c("Absence", "Presence"),
           col = c("lightgray", "black"),
           pch = 19,
           cex = 0.7)
  }
}

# Create base R plots for stepwise models
cat("\n=== Creating Base R Plots for Final Comparison ===\n")
par(mfrow = c(2, 2))

base_distribution_plot(mammals_data$VulpesVulpes,
                      mammals_data[,c("X_WGS84", "Y_WGS84")],
                      "Original data", is_binary = TRUE)

base_distribution_plot(fitted(glmModAIC),
                      mammals_data[,c("X_WGS84", "Y_WGS84")],
                      "Stepwise GLM with AIC")

base_distribution_plot(fitted(glmModBIC),
                      mammals_data[,c("X_WGS84", "Y_WGS84")],
                      "Stepwise GLM with BIC")

base_distribution_plot(fitted(glm2),
                      mammals_data[,c("X_WGS84", "Y_WGS84")],
                      "GLM with quadratic terms")

par(mfrow = c(1, 1))

cat("\n=== Analysis Complete ===\n")
cat("Comprehensive GLM analysis completed successfully!\n")
cat("All modern plotting functions implemented to replace deprecated level.plot and response.plot2\n")
cat("Stepwise model selection integrated with base model comparison\n")
cat("Best performing model:", comparison_table$Model[best_aic_idx], "\n")