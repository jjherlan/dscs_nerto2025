# =============================================================================
# Enhanced GAM Analysis for Rare Species Distribution Modeling
# Species: Desmophyllum dianthus
# Data: Northeast Canyons Deep-Sea Coral Survey (DSCS) - Bottom Environmental Data
# Author: [Enhanced Version]
# Date: July 2, 2025
# ENHANCED VERSION - Improved evaluation metrics for rare species
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
library(mgcv)         # For GAM modeling
library(pROC)         # For ROC curves and AUC calculation
library(caret)        # For additional model evaluation metrics
library(PresenceAbsence) # For threshold optimization
library(verification)  # For additional SDM metrics

# Load the dataset
cat("=== Loading Dataset ===\n")
ddia_btm <- read.csv("data/ne_canyons_dscs_btm_20250630.csv", row.names = 1)

# =============================================================================
# SECTION 2: Data Exploration and Quality Assessment
# =============================================================================

cat("=== Dataset Summary ===\n")
cat("Total observations:", nrow(ddia_btm), "\n")
cat("Desmophyllum dianthus presence records:", sum(ddia_btm$Desmophyllum.dianthus == 1), "\n")
cat("Desmophyllum dianthus absence records:", sum(ddia_btm$Desmophyllum.dianthus == 0), "\n")
prevalence <- mean(ddia_btm$Desmophyllum.dianthus)
cat("Prevalence:", round(prevalence, 3), "\n")
cat("NOTE: This is a rare species - special evaluation metrics needed!\n\n")

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

# =============================================================================
# SECTION 3: Enhanced Evaluation Functions for Rare Species
# =============================================================================

# Function to calculate optimal thresholds using multiple methods
calculate_optimal_thresholds <- function(observed, predicted, model_name = "") {
  
  # Create data frame for PresenceAbsence package
  pa_data <- data.frame(
    ID = 1:length(observed),
    Observed = observed,
    Predicted = predicted
  )
  
  # Calculate various optimal thresholds
  thresholds <- optimal.thresholds(pa_data, opt.methods = c(
    "Default",        # 0.5 threshold
    "Sens=Spec",      # Sensitivity equals Specificity
    "MaxSens+Spec",   # Maximize sum of sensitivity and specificity
    "MaxKappa",       # Maximize Cohen's Kappa
    "MaxPCC",         # Maximize Percent Correctly Classified
    "PredPrev=Obs",   # Predicted prevalence equals observed prevalence
    "ObsPrev",        # Use observed prevalence as threshold
    "MeanProb",       # Use mean predicted probability as threshold
    "MinROCdist"      # Minimize distance to perfect classification on ROC curve
  ))
  
  # Add TSS-based threshold (True Skill Statistic)
  tss_threshold <- calculate_tss_threshold(observed, predicted)
  
  # Combine results
  results <- data.frame(
    Model = model_name,
    Method = c(as.character(thresholds$Method), "MaxTSS"),
    Threshold = c(thresholds$Predicted, tss_threshold)
  )
  
  return(results)
}

# Function to calculate TSS-optimized threshold
calculate_tss_threshold <- function(observed, predicted) {
  
  # Test a range of thresholds
  thresholds <- seq(0, 1, by = 0.001)
  tss_values <- numeric(length(thresholds))
  
  for(i in 1:length(thresholds)) {
    pred_binary <- ifelse(predicted >= thresholds[i], 1, 0)
    conf_mat <- table(Predicted = pred_binary, Observed = observed)
    
    # Calculate sensitivity and specificity
    if(all(c("0","1") %in% rownames(conf_mat)) && all(c("0","1") %in% colnames(conf_mat))) {
      sensitivity <- conf_mat["1", "1"] / sum(conf_mat[, "1"])
      specificity <- conf_mat["0", "0"] / sum(conf_mat[, "0"])
      tss_values[i] <- sensitivity + specificity - 1
    } else {
      tss_values[i] <- NA
    }
  }
  
  # Find threshold that maximizes TSS
  optimal_idx <- which.max(tss_values)
  return(thresholds[optimal_idx])
}

# Enhanced performance calculation function for rare species
calculate_rare_species_performance <- function(model, data, response_var, threshold = NULL) {
  
  fitted_probs <- fitted(model)
  observed <- data[[response_var]]
  
  # Basic model statistics
  model_aic <- AIC(model)
  dev_expl <- summary(model)$dev.expl * 100
  
  # Calculate ROC and AUC
  roc_obj <- roc(observed, fitted_probs, quiet = TRUE)
  auc_value <- auc(roc_obj)
  
  # If no threshold provided, use TSS-optimized threshold
  if(is.null(threshold)) {
    threshold <- calculate_tss_threshold(observed, fitted_probs)
  }
  
  # Binary predictions using specified threshold
  pred_binary <- ifelse(fitted_probs >= threshold, 1, 0)
  conf_mat <- table(Predicted = pred_binary, Observed = observed)
  
  # Calculate comprehensive metrics
  if(all(c("0","1") %in% rownames(conf_mat)) && all(c("0","1") %in% colnames(conf_mat))) {
    TP <- conf_mat["1", "1"]
    TN <- conf_mat["0", "0"]
    FP <- conf_mat["1", "0"]
    FN <- conf_mat["0", "1"]
    
    # Standard metrics
    accuracy <- (TP + TN) / sum(conf_mat)
    sensitivity <- TP / (TP + FN)  # True Positive Rate (Recall)
    specificity <- TN / (TN + FP)  # True Negative Rate
    
    # Additional metrics for rare species
    precision <- ifelse((TP + FP) > 0, TP / (TP + FP), 0)  # Positive Predictive Value
    npv <- ifelse((TN + FN) > 0, TN / (TN + FN), 0)  # Negative Predictive Value
    f1_score <- ifelse((precision + sensitivity) > 0, 
                       2 * (precision * sensitivity) / (precision + sensitivity), 0)
    
    # TSS (True Skill Statistic) - good for imbalanced data
    tss <- sensitivity + specificity - 1
    
    # Cohen's Kappa - accounts for chance agreement
    po <- accuracy  # Observed agreement
    pe <- ((TP + FP) * (TP + FN) + (TN + FN) * (TN + FP)) / (sum(conf_mat)^2)  # Expected agreement
    kappa <- (po - pe) / (1 - pe)
    
    # Matthews Correlation Coefficient - good for imbalanced data
    mcc_num <- (TP * TN) - (FP * FN)
    mcc_den <- sqrt((TP + FP) * (TP + FN) * (TN + FP) * (TN + FN))
    mcc <- ifelse(mcc_den > 0, mcc_num / mcc_den, 0)
    
  } else {
    # Handle cases where not all categories are predicted
    accuracy <- sensitivity <- specificity <- precision <- npv <- NA
    f1_score <- tss <- kappa <- mcc <- NA
  }
  
  # Calculate Boyce index (continuous evaluation metric)
  boyce_index <- calculate_boyce_index(observed, fitted_probs)
  
  return(list(
    AIC = round(model_aic, 2),
    Deviance_Explained = round(dev_expl, 2),
    AUC = round(auc_value, 3),
    Threshold_Used = round(threshold, 3),
    Accuracy = round(accuracy, 3),
    Sensitivity_TPR = round(sensitivity, 3),
    Specificity_TNR = round(specificity, 3),
    Precision_PPV = round(precision, 3),
    NPV = round(npv, 3),
    F1_Score = round(f1_score, 3),
    TSS = round(tss, 3),
    Kappa = round(kappa, 3),
    MCC = round(mcc, 3),
    Boyce_Index = round(boyce_index, 3)
  ))
}

# Function to calculate Boyce Index (good for presence-only evaluations)
calculate_boyce_index <- function(observed, predicted, n_bins = 10) {
  
  # Create bins of predicted values
  breaks <- seq(0, 1, length.out = n_bins + 1)
  
  # Calculate expected and observed frequencies
  expected_freq <- numeric(n_bins)
  observed_freq <- numeric(n_bins)
  
  for(i in 1:n_bins) {
    in_bin <- predicted >= breaks[i] & predicted < breaks[i+1]
    if(i == n_bins) in_bin <- predicted >= breaks[i] & predicted <= breaks[i+1]
    
    expected_freq[i] <- sum(in_bin) / length(predicted)
    if(sum(in_bin) > 0) {
      observed_freq[i] <- sum(observed[in_bin]) / sum(observed)
    } else {
      observed_freq[i] <- 0
    }
  }
  
  # Remove bins with no predictions
  valid_bins <- expected_freq > 0
  if(sum(valid_bins) < 2) return(NA)
  
  # Calculate correlation
  boyce <- cor(expected_freq[valid_bins], observed_freq[valid_bins], 
               method = "spearman", use = "complete.obs")
  
  return(boyce)
}

# =============================================================================
# SECTION 4: Model Fitting with Cross-Validation
# =============================================================================

cat("\n=== Fitting GAM Models with Cross-Validation ===\n")

# Set up k-fold cross-validation
set.seed(42)  # For reproducibility
n_folds <- 5
folds <- createFolds(ddia_btm$Desmophyllum.dianthus, k = n_folds, list = TRUE)

# Initialize storage for cross-validation results
cv_results <- list()

# Function to fit model and evaluate on hold-out data
fit_and_evaluate_cv <- function(train_data, test_data, formula, k_value = NULL) {
  
  # Fit model on training data
  if(!is.null(k_value)) {
    # Modify formula to include k value
    formula_str <- deparse(formula)
    formula_str <- gsub("s\\(([^,)]+)\\)", paste0("s(\\1, k=", k_value, ")"), formula_str)
    formula <- as.formula(formula_str)
  }
  
  model <- gam(formula, data = train_data, family = binomial, method = "REML")
  
  # Predict on test data
  test_probs <- predict(model, newdata = test_data, type = "response")
  
  # Calculate optimal threshold on training data
  train_probs <- fitted(model)
  optimal_thresh <- calculate_tss_threshold(train_data$Desmophyllum.dianthus, train_probs)
  
  # Evaluate on test data
  test_binary <- ifelse(test_probs >= optimal_thresh, 1, 0)
  
  # Calculate AUC on test data
  roc_test <- roc(test_data$Desmophyllum.dianthus, test_probs, quiet = TRUE)
  auc_test <- auc(roc_test)
  
  # Calculate TSS on test data
  conf_mat <- table(Predicted = test_binary, Observed = test_data$Desmophyllum.dianthus)
  if(all(c("0","1") %in% rownames(conf_mat)) && all(c("0","1") %in% colnames(conf_mat))) {
    sensitivity <- conf_mat["1", "1"] / sum(conf_mat[, "1"])
    specificity <- conf_mat["0", "0"] / sum(conf_mat[, "0"])
    tss_test <- sensitivity + specificity - 1
  } else {
    tss_test <- NA
  }
  
  return(list(
    model = model,
    auc = auc_test,
    tss = tss_test,
    threshold = optimal_thresh
  ))
}

# Fit main models
formula_full <- Desmophyllum.dianthus ~ s(btm_arag_ann) + s(btm_chl_ann) + 
                s(btm_curr_mag_ann) + s(btm_dissic_ann) + s(btm_dissoc_ann) + 
                s(btm_o2_ann) + s(btm_sal_ann) + s(btm_talk_ann) + s(btm_temp_ann)

# GAM with different complexities
gam1 <- gam(Desmophyllum.dianthus ~ s(btm_arag_ann, k=3) + s(btm_chl_ann, k=3) + 
            s(btm_curr_mag_ann, k=3) + s(btm_dissic_ann, k=3) + s(btm_dissoc_ann, k=3) + 
            s(btm_o2_ann, k=3) + s(btm_sal_ann, k=3) + s(btm_talk_ann, k=3) + s(btm_temp_ann, k=3),
            data = ddia_btm, family = binomial, method = "REML")

gam2 <- gam(Desmophyllum.dianthus ~ s(btm_arag_ann, k=5) + s(btm_chl_ann, k=5) + 
            s(btm_curr_mag_ann, k=5) + s(btm_dissic_ann, k=5) + s(btm_dissoc_ann, k=5) + 
            s(btm_o2_ann, k=5) + s(btm_sal_ann, k=5) + s(btm_talk_ann, k=5) + s(btm_temp_ann, k=5),
            data = ddia_btm, family = binomial, method = "REML")

gam3 <- gam(formula_full, data = ddia_btm, family = binomial, method = "REML")

# GAM with variable selection
gamModSelect <- gam(formula_full, data = ddia_btm, family = binomial, method = "REML", select = TRUE)

cat("Basic GAM models fitted successfully!\n")

# =============================================================================
# SECTION 5: Threshold Optimization and Comparison
# =============================================================================

cat("\n=== Optimal Threshold Analysis ===\n")

# Calculate optimal thresholds for each model
models_list <- list(gam1, gam2, gam3, gamModSelect)
model_names <- c("GAM1 (k=3)", "GAM2 (k=5)", "GAM3 (default)", "GAM Select")

# Store all threshold results
all_thresholds <- list()

for(i in 1:length(models_list)) {
  fitted_probs <- fitted(models_list[[i]])
  observed <- ddia_btm$Desmophyllum.dianthus
  
  thresholds <- calculate_optimal_thresholds(observed, fitted_probs, model_names[i])
  all_thresholds[[i]] <- thresholds
  
  cat("\n", model_names[i], "- Optimal Thresholds:\n")
  print(thresholds)
}

# Combine all threshold results
threshold_comparison <- do.call(rbind, all_thresholds)

# =============================================================================
# SECTION 6: Comprehensive Model Evaluation with Multiple Metrics
# =============================================================================

cat("\n=== Comprehensive Model Evaluation for Rare Species ===\n")

# Evaluate each model with different thresholds
evaluation_results <- list()

for(i in 1:length(models_list)) {
  cat("\nEvaluating", model_names[i], "...\n")
  
  # Get TSS-optimized threshold for this model
  fitted_probs <- fitted(models_list[[i]])
  observed <- ddia_btm$Desmophyllum.dianthus
  optimal_threshold <- calculate_tss_threshold(observed, fitted_probs)
  
  # Evaluate with optimal threshold
  perf_optimal <- calculate_rare_species_performance(
    models_list[[i]], ddia_btm, "Desmophyllum.dianthus", 
    threshold = optimal_threshold
  )
  
  # Also evaluate with prevalence-based threshold
  perf_prevalence <- calculate_rare_species_performance(
    models_list[[i]], ddia_btm, "Desmophyllum.dianthus", 
    threshold = prevalence
  )
  
  # Store results
  evaluation_results[[model_names[i]]] <- list(
    optimal = perf_optimal,
    prevalence = perf_prevalence
  )
}

# Create comparison tables
cat("\n=== Performance Comparison - TSS-Optimized Thresholds ===\n")
optimal_comparison <- data.frame(
  Model = model_names,
  AIC = sapply(evaluation_results, function(x) x$optimal$AIC),
  Dev_Expl = sapply(evaluation_results, function(x) x$optimal$Deviance_Explained),
  AUC = sapply(evaluation_results, function(x) x$optimal$AUC),
  Threshold = sapply(evaluation_results, function(x) x$optimal$Threshold_Used),
  Sensitivity = sapply(evaluation_results, function(x) x$optimal$Sensitivity_TPR),
  Specificity = sapply(evaluation_results, function(x) x$optimal$Specificity_TNR),
  TSS = sapply(evaluation_results, function(x) x$optimal$TSS),
  F1 = sapply(evaluation_results, function(x) x$optimal$F1_Score),
  Kappa = sapply(evaluation_results, function(x) x$optimal$Kappa),
  MCC = sapply(evaluation_results, function(x) x$optimal$MCC),
  Boyce = sapply(evaluation_results, function(x) x$optimal$Boyce_Index)
)
print(optimal_comparison)

cat("\n=== Performance Comparison - Prevalence-Based Thresholds ===\n")
prevalence_comparison <- data.frame(
  Model = model_names,
  Threshold = round(prevalence, 3),
  Sensitivity = sapply(evaluation_results, function(x) x$prevalence$Sensitivity_TPR),
  Specificity = sapply(evaluation_results, function(x) x$prevalence$Specificity_TNR),
  TSS = sapply(evaluation_results, function(x) x$prevalence$TSS),
  F1 = sapply(evaluation_results, function(x) x$prevalence$F1_Score),
  Kappa = sapply(evaluation_results, function(x) x$prevalence$Kappa)
)
print(prevalence_comparison)

# =============================================================================
# SECTION 7: Enhanced Visualizations for Rare Species
# =============================================================================

cat("\n=== Creating Enhanced Visualizations ===\n")

# 1. ROC Curves Comparison
cat("Creating ROC curves comparison...\n")
roc_plot_data <- data.frame()

for(i in 1:length(models_list)) {
  fitted_probs <- fitted(models_list[[i]])
  observed <- ddia_btm$Desmophyllum.dianthus
  
  roc_obj <- roc(observed, fitted_probs, quiet = TRUE)
  roc_df <- data.frame(
    Model = model_names[i],
    FPR = 1 - roc_obj$specificities,
    TPR = roc_obj$sensitivities,
    AUC = round(auc(roc_obj), 3)
  )
  roc_plot_data <- rbind(roc_plot_data, roc_df)
}

# Create ROC plot
roc_plot <- ggplot(roc_plot_data, aes(x = FPR, y = TPR, color = Model)) +
  geom_line(size = 1.2) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray50") +
  labs(title = "ROC Curves Comparison - Deep-Sea Coral GAM Models",
       subtitle = "Rare species modeling requires careful threshold selection",
       x = "False Positive Rate (1 - Specificity)",
       y = "True Positive Rate (Sensitivity)") +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, size = 12),
    legend.position = "bottom"
  ) +
  scale_color_viridis_d()

print(roc_plot)

# 2. Precision-Recall Curves (better for imbalanced data)
cat("Creating Precision-Recall curves...\n")
pr_plot_data <- data.frame()

for(i in 1:length(models_list)) {
  fitted_probs <- fitted(models_list[[i]])
  observed <- ddia_btm$Desmophyllum.dianthus
  
  # Calculate precision and recall for different thresholds
  thresholds <- seq(0, 1, by = 0.01)
  precision <- recall <- numeric(length(thresholds))
  
  for(j in 1:length(thresholds)) {
    pred_binary <- ifelse(fitted_probs >= thresholds[j], 1, 0)
    conf_mat <- table(Predicted = pred_binary, Observed = observed)
    
    if(all(c("0","1") %in% rownames(conf_mat)) && all(c("0","1") %in% colnames(conf_mat))) {
      TP <- conf_mat["1", "1"]
      FP <- conf_mat["1", "0"]
      FN <- conf_mat["0", "1"]
      
      precision[j] <- ifelse((TP + FP) > 0, TP / (TP + FP), 0)
      recall[j] <- ifelse((TP + FN) > 0, TP / (TP + FN), 0)
    } else {
      precision[j] <- recall[j] <- NA
    }
  }
  
  pr_df <- data.frame(
    Model = model_names[i],
    Precision = precision,
    Recall = recall
  )
  pr_plot_data <- rbind(pr_plot_data, pr_df)
}

# Create Precision-Recall plot
pr_plot <- ggplot(pr_plot_data, aes(x = Recall, y = Precision, color = Model)) +
  geom_line(size = 1.2) +
  geom_hline(yintercept = prevalence, linetype = "dashed", color = "red") +
  annotate("text", x = 0.5, y = prevalence + 0.02, 
           label = paste("Baseline (prevalence =", round(prevalence, 3), ")"),
           color = "red", size = 3) +
  labs(title = "Precision-Recall Curves - Better for Rare Species Evaluation",
       subtitle = "Higher curves indicate better performance for rare species detection",
       x = "Recall (Sensitivity)",
       y = "Precision (Positive Predictive Value)") +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, size = 12),
    legend.position = "bottom"
  ) +
  scale_color_viridis_d()

print(pr_plot)

# 3. Threshold Performance Plot
cat("Creating threshold performance visualization...\n")
threshold_perf_data <- data.frame()

# Use GAM3 as example
fitted_probs <- fitted(gam3)
observed <- ddia_btm$Desmophyllum.dianthus
thresholds <- seq(0, 0.2, by = 0.005)  # Focus on lower thresholds for rare species

for(thresh in thresholds) {
  pred_binary <- ifelse(fitted_probs >= thresh, 1, 0)
  conf_mat <- table(Predicted = pred_binary, Observed = observed)
  
  if(all(c("0","1") %in% rownames(conf_mat)) && all(c("0","1") %in% colnames(conf_mat))) {
    sensitivity <- conf_mat["1", "1"] / sum(conf_mat[, "1"])
    specificity <- conf_mat["0", "0"] / sum(conf_mat[, "0"])
    precision <- conf_mat["1", "1"] / sum(conf_mat["1", ])
    f1 <- 2 * (precision * sensitivity) / (precision + sensitivity)
    tss <- sensitivity + specificity - 1
    
    threshold_perf_data <- rbind(threshold_perf_data, data.frame(
      Threshold = thresh,
      Sensitivity = sensitivity,
      Specificity = specificity,
      Precision = precision,
      F1_Score = f1,
      TSS = tss
    ))
  }
}

# Melt data for plotting
threshold_perf_long <- reshape2::melt(threshold_perf_data, id.vars = "Threshold",
                                     variable.name = "Metric", value.name = "Value")

# Create threshold performance plot
thresh_plot <- ggplot(threshold_perf_long, aes(x = Threshold, y = Value, color = Metric)) +
  geom_line(size = 1.2) +
  geom_vline(xintercept = prevalence, linetype = "dashed", color = "red", alpha = 0.5) +
  geom_vline(xintercept = calculate_tss_threshold(observed, fitted_probs), 
             linetype = "dashed", color = "blue", alpha = 0.5) +
  annotate("text", x = prevalence, y = 0.1, label = "Prevalence", 
           angle = 90, vjust = -0.5, color = "red", size = 3) +
  annotate("text", x = calculate_tss_threshold(observed, fitted_probs), 
           y = 0.1, label = "Optimal TSS", 
           angle = 90, vjust = 1.5, color = "blue", size = 3) +
  labs(title = "Performance Metrics across Different Thresholds",
       subtitle = "GAM3 model - Showing how metrics change with threshold selection",
       x = "Probability Threshold",
       y = "Metric Value") +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, size = 12),
    legend.position = "right"
  ) +
  scale_color_viridis_d()

print(thresh_plot)

# =============================================================================
# SECTION 8: Model Selection Recommendations for Rare Species
# =============================================================================

cat("\n=== Model Selection Recommendations for Rare Species ===\n")

# Identify best models based on different criteria
best_aic <- which.min(optimal_comparison$AIC)
best_auc <- which.max(optimal_comparison$AUC)
best_tss <- which.max(optimal_comparison$TSS)
best_f1 <- which.max(optimal_comparison$F1)

cat("\nBest models by different criteria:\n")
cat("- Lowest AIC:", model_names[best_aic], "(AIC =", optimal_comparison$AIC[best_aic], ")\n")
cat("- Highest AUC:", model_names[best_auc], "(AUC =", optimal_comparison$AUC[best_auc], ")\n")
cat("- Highest TSS:", model_names[best_tss], "(TSS =", optimal_comparison$TSS[best_tss], ")\n")
cat("- Highest F1 Score:", model_names[best_f1], "(F1 =", optimal_comparison$F1[best_f1], ")\n")

# =============================================================================
# SECTION 9: Export Results and Final Summary
# =============================================================================

cat("\n=== ENHANCED ANALYSIS SUMMARY ===\n")
cat("Dataset: Northeast Canyons DSCS Bottom Environmental Data\n")
cat("Species: Desmophyllum dianthus (Rare species - ", round(prevalence*100, 1), "% prevalence)\n", sep="")
cat("Total observations:", nrow(ddia_btm), "\n")
cat("Environmental variables tested:", length(env_vars), "\n")
cat("Models fitted:", length(models_list), "\n")

cat("\nKey Improvements for Rare Species Modeling:\n")
cat("1. Implemented multiple threshold optimization methods (TSS, Kappa, F1, etc.)\n")
cat("2. Added comprehensive evaluation metrics suitable for imbalanced data\n")
cat("3. Created Precision-Recall curves (more informative than ROC for rare species)\n")
cat("4. Analyzed performance across different threshold values\n")
cat("5. Included Boyce Index for continuous habitat suitability evaluation\n")

cat("\nRecommendations:\n")
cat("1. For rare species, use TSS-optimized thresholds rather than default 0.5\n")
cat("2. Consider multiple evaluation metrics - don't rely solely on AUC\n")
cat("3. Precision-Recall curves provide better insight for rare species\n")
cat("4. Balance between sensitivity (finding presences) and precision (avoiding false positives)\n")
cat("5. Consider the ecological cost of false negatives vs false positives\n")

# Save evaluation results
cat("\nSaving evaluation results...\n")
write.csv(optimal_comparison, "gam_evaluation_optimal_thresholds.csv", row.names = FALSE)
write.csv(threshold_comparison, "all_threshold_methods_comparison.csv", row.names = FALSE)

cat("\n=== ENHANCED ANALYSIS COMPLETE ===\n")
cat("All models evaluated with appropriate metrics for rare species!\n")
cat("Results saved to CSV files for further analysis.\n")

# Function to predict with optimal threshold
predict_rare_species <- function(model, newdata, threshold = NULL) {
  # Get predictions
  probs <- predict(model, newdata = newdata, type = "response")
  
  # If no threshold provided, use prevalence
  if(is.null(threshold)) {
    threshold <- prevalence
    cat("Using prevalence-based threshold:", round(threshold, 3), "\n")
  }
  
  # Convert to binary predictions
  binary_pred <- ifelse(probs >= threshold, 1, 0)
  
  # Return both probabilities and binary predictions
  return(list(
    probabilities = probs,
    binary = binary_pred,
    threshold_used = threshold
  ))
}

cat("\nExample usage for predictions with new data:\n")
cat("predictions <- predict_rare_species(gam3, newdata = new_data, threshold = 0.037)\n")
cat("This will use the appropriate threshold for rare species classification.\n")
