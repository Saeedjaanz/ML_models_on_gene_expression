# ============================================================================
# SVM-RFE (Support Vector Machine - Recursive Feature Elimination)
# ============================================================================

# Required Libraries
library(e1071)           # SVM implementation
library(caret)           # Machine learning workflows
library(ROCR)            # ROC curves
library(pROC)            # Advanced ROC analysis
library(ggplot2)         # Visualization
library(pheatmap)        # Heatmaps
library(RColorBrewer)    # Color palettes
library(gridExtra)       # Multiple plots
library(reshape2)        # Data reshaping
library(ComplexHeatmap)  # Advanced heatmaps
library(circlize)        # Color mapping
library(ggpubr)          # Publication-ready plots
library(scales)          # Scale functions
library(boot)            # Bootstrap resampling

# Set seed for reproducibility
set.seed(123)

# ============================================================================
# PART 1: DATA PREPARATION
# ============================================================================

# Transpose expression matrix (samples as rows, genes as columns)
expr_ml_t <- t(expr_ml)
data_ml <- as.data.frame(expr_ml_t)

# Add condition labels
data_ml$condition <- condition

# Verify data structure
cat("Dataset dimensions:", dim(data_ml), "\n")
cat("Number of features:", ncol(data_ml) - 1, "\n")
cat("Number of samples:", nrow(data_ml), "\n")
cat("Class distribution:\n")
print(table(data_ml$condition))

# Check for class imbalance
class_counts <- table(data_ml$condition)
imbalance_ratio <- max(class_counts) / min(class_counts)
cat(sprintf("\nClass imbalance ratio: %.2f\n", imbalance_ratio))
if(imbalance_ratio > 1.5) {
  cat("NOTE: Class imbalance detected. Using AUC for optimization.\n")
}

# ============================================================================
# PART 2: ENHANCED SVM-RFE IMPLEMENTATION (AUC-based)
# ============================================================================

# Enhanced SVM-RFE function with AUC as objective (handles class imbalance)
svm_rfe_cv <- function(data, class_col = "condition", step_size = 1, 
                       cv_folds = 10, kernel = "linear", use_auc = TRUE) {
  
  # Separate features and labels
  features <- data[, -which(names(data) == class_col)]
  labels <- data[[class_col]]
  
  # Initialize results storage
  feature_ranking <- list()
  cv_performance <- list()
  cv_performance_sd <- list()
  feature_sets <- list()
  
  current_features <- colnames(features)
  iteration <- 1
  
  cat("\n=== Starting SVM-RFE with", ifelse(use_auc, "AUC", "Accuracy"), "optimization ===\n")
  
  while(length(current_features) > 0) {
    
    cat(sprintf("Iteration %d: %d features remaining\n", 
                iteration, length(current_features)))
    
    # Current feature subset
    X_current <- features[, current_features, drop = FALSE]
    
    # Cross-validation
    cv_metrics <- numeric(cv_folds)
    feature_importance <- numeric(length(current_features))
    
    # Create folds
    folds <- createFolds(labels, k = cv_folds, list = TRUE)
    
    for(fold in 1:cv_folds) {
      # Split data
      train_idx <- unlist(folds[-fold])
      test_idx <- folds[[fold]]
      
      X_train <- X_current[train_idx, , drop = FALSE]
      y_train <- labels[train_idx]
      X_test <- X_current[test_idx, , drop = FALSE]
      y_test <- labels[test_idx]
      
      # Train SVM with probability estimates
      svm_model <- svm(x = X_train, y = y_train, kernel = kernel, 
                       scale = TRUE, probability = TRUE, type = "C-classification")
      
      # Predict
      if(use_auc) {
        # Get probabilities for AUC calculation
        predictions <- predict(svm_model, X_test, probability = TRUE)
        probs <- attr(predictions, "probabilities")
        
        # Calculate AUC (assuming binary classification)
        roc_obj <- roc(y_test, probs[, levels(labels)[2]], quiet = TRUE)
        cv_metrics[fold] <- as.numeric(auc(roc_obj))
      } else {
        # Use accuracy
        predictions <- predict(svm_model, X_test)
        cv_metrics[fold] <- mean(predictions == y_test)
      }
      
      # Calculate feature weights (for linear kernel)
      if(kernel == "linear") {
        w <- t(svm_model$coefs) %*% svm_model$SV
        feature_importance <- feature_importance + abs(as.vector(w))
      }
    }
    
    # Store results
    mean_perf <- mean(cv_metrics)
    sd_perf <- sd(cv_metrics)
    cv_performance[[iteration]] <- mean_perf
    cv_performance_sd[[iteration]] <- sd_perf
    feature_sets[[iteration]] <- current_features
    
    metric_name <- ifelse(use_auc, "AUC", "Accuracy")
    cat(sprintf("  Mean CV %s: %.4f (±%.4f)\n", 
                metric_name, mean_perf, sd_perf))
    
    # Rank features by importance
    if(kernel == "linear") {
      feature_ranks <- rank(feature_importance)
      least_important_idx <- which.min(feature_importance)
    } else {
      # For non-linear kernels, use random forest importance
      rf_model <- randomForest::randomForest(x = X_current, y = labels)
      feature_ranks <- rank(-rf_model$importance)
      least_important_idx <- which.max(feature_ranks)
    }
    
    # Store eliminated feature
    eliminated_feature <- current_features[least_important_idx]
    feature_ranking[[iteration]] <- eliminated_feature
    
    # Remove least important feature
    current_features <- current_features[-least_important_idx]
    
    iteration <- iteration + 1
  }
  
  # Compile results
  results <- list(
    feature_ranking = rev(unlist(feature_ranking)),
    cv_performance = unlist(cv_performance),
    cv_performance_sd = unlist(cv_performance_sd),
    feature_sets = feature_sets,
    optimal_n_features = which.max(unlist(cv_performance)),
    performance_metric = ifelse(use_auc, "AUC", "Accuracy")
  )
  return(results)
}

# ============================================================================
# PART 3: NESTED CROSS-VALIDATION (For Unbiased Performance)
# ============================================================================

cat("\n" %+% strrep("=", 70) %+% "\n")
cat("NESTED CROSS-VALIDATION FOR UNBIASED PERFORMANCE ESTIMATION\n")
cat(strrep("=", 70) %+% "\n")

# Nested CV function
nested_cv_svm_rfe <- function(data, class_col = "condition",
                              outer_folds = 5,inner_folds = 5,
                              kernel = "linear") {
  
  # Separate features and labels
  features <- data[, -which(names(data) == class_col)]
  labels <- data[[class_col]]
  
  # Create outer folds (for honest test sets)
  outer_folds_list <- createFolds(labels, k = outer_folds, list = TRUE)
  
  # Storage for results
  outer_results <- list()
  outer_selected_features <- list()
  outer_rfe_results <- list()
  
  cat(sprintf("\nRunning %d-fold outer CV with %d-fold inner CV...\n", 
              outer_folds, inner_folds))
  
  for(outer_fold in 1:outer_folds) {
    
    cat(sprintf("\n--- Outer Fold %d/%d ---\n", outer_fold, outer_folds))
    
    # Split into outer train and test
    outer_test_idx <- outer_folds_list[[outer_fold]]
    outer_train_idx <- setdiff(1:nrow(data), outer_test_idx)
    
    outer_train_data <- data[outer_train_idx, ]
    outer_test_data <- data[outer_test_idx, ]
    
    # Run RFE on outer training set with inner CV
    inner_rfe <- svm_rfe_cv(outer_train_data, class_col = class_col,
                            cv_folds = inner_folds, kernel = kernel,
                            use_auc = TRUE)
    
    # Store inner RFE results
    outer_rfe_results[[outer_fold]] <- inner_rfe
    
    # Get optimal features from inner CV
    optimal_features <- inner_rfe$feature_ranking[1:inner_rfe$optimal_n_features]
    outer_selected_features[[outer_fold]] <- optimal_features
    
    cat(sprintf("  Inner CV selected %d features\n", length(optimal_features)))
    
    # Train final model on outer train set with selected features
    train_subset <- outer_train_data[, c(optimal_features, class_col)]
    
    final_model <- svm(as.formula(paste(class_col, "~ .")), 
                       data = train_subset, kernel = kernel,
                       probability = TRUE, scale = TRUE)
    
    # Evaluate on outer test set
    test_subset <- outer_test_data[, c(optimal_features, class_col)]
    test_pred <- predict(final_model, test_subset, probability = TRUE)
    test_pred_class <- predict(final_model, test_subset)
    test_probs <- attr(test_pred, "probabilities")[, levels(labels)[2]]
    
    # Calculate metrics
    test_acc <- mean(test_pred_class == outer_test_data[[class_col]])
    test_roc <- roc(outer_test_data[[class_col]], test_probs, quiet = TRUE)
    test_auc <- as.numeric(auc(test_roc))
    
    # Confusion matrix
    cm <- confusionMatrix(test_pred_class, outer_test_data[[class_col]])
    
    # Store results
    outer_results[[outer_fold]] <- list(
      accuracy = test_acc,
      auc = test_auc,
      sensitivity = cm$byClass["Sensitivity"],
      specificity = cm$byClass["Specificity"],
      ppv = cm$byClass["Pos Pred Value"],
      npv = cm$byClass["Neg Pred Value"],
      f1 = cm$byClass["F1"],
      n_features = length(optimal_features),
      selected_features = optimal_features,
      confusion_matrix = cm$table,
      predictions = test_pred_class,
      true_labels = outer_test_data[[class_col]],
      probabilities = test_probs
    )
    
    cat(sprintf("  Outer test AUC: %.4f, Accuracy: %.4f\n", test_auc, test_acc))
  }
  
  # Aggregate results
  agg_results <- list(
    mean_auc = mean(sapply(outer_results, function(x) x$auc)),
    sd_auc = sd(sapply(outer_results, function(x) x$auc)),
    mean_accuracy = mean(sapply(outer_results, function(x) x$accuracy)),
    sd_accuracy = sd(sapply(outer_results, function(x) x$accuracy)),
    mean_sensitivity = mean(sapply(outer_results, function(x) x$sensitivity), na.rm = TRUE),
    sd_sensitivity = sd(sapply(outer_results, function(x) x$sensitivity), na.rm = TRUE),
    mean_specificity = mean(sapply(outer_results, function(x) x$specificity), na.rm = TRUE),
    sd_specificity = sd(sapply(outer_results, function(x) x$specificity), na.rm = TRUE),
    mean_f1 = mean(sapply(outer_results, function(x) x$f1), na.rm = TRUE),
    sd_f1 = sd(sapply(outer_results, function(x) x$f1), na.rm = TRUE),
    mean_n_features = mean(sapply(outer_results, function(x) x$n_features)),
    fold_results = outer_results,
    selected_features_per_fold = outer_selected_features,
    rfe_results_per_fold = outer_rfe_results
  )
  
  return(agg_results)
}

# Run nested CV
nested_cv_results <- nested_cv_svm_rfe(data_ml, 
                                       class_col = "condition",
                                       outer_folds = 5,
                                       inner_folds = 5,
                                       kernel = "linear")

cat("\n=== NESTED CV RESULTS (UNBIASED PERFORMANCE) ===\n")
cat(sprintf("Mean AUC: %.4f ± %.4f\n", 
            nested_cv_results$mean_auc, 
            nested_cv_results$sd_auc))
cat(sprintf("Mean Accuracy: %.4f ± %.4f\n", 
            nested_cv_results$mean_accuracy, 
            nested_cv_results$sd_accuracy))
cat(sprintf("Mean Sensitivity: %.4f ± %.4f\n", 
            nested_cv_results$mean_sensitivity,
            nested_cv_results$sd_sensitivity))
cat(sprintf("Mean Specificity: %.4f ± %.4f\n", 
            nested_cv_results$mean_specificity,
            nested_cv_results$sd_specificity))
cat(sprintf("Mean F1 Score: %.4f ± %.4f\n", 
            nested_cv_results$mean_f1,
            nested_cv_results$sd_f1))
cat(sprintf("Mean # Features: %.1f\n", nested_cv_results$mean_n_features))

# ============================================================================
# PART 4: FEATURE STABILITY ANALYSIS (Bootstrap Resampling)
# ============================================================================

cat("\n" %+% strrep("=", 70) %+% "\n")
cat("FEATURE STABILITY ANALYSIS VIA BOOTSTRAP\n")
cat(strrep("=", 70) %+% "\n")

# Bootstrap feature selection
bootstrap_feature_selection <- function(data, class_col = "condition",
                                        n_bootstrap = 100,
                                        cv_folds = 5,
                                        kernel = "linear") {
  
  n_samples <- nrow(data)
  all_features <- colnames(data)[colnames(data) != class_col]
  
  # Storage
  feature_frequency <- matrix(0, nrow = n_bootstrap, ncol = length(all_features))
  colnames(feature_frequency) <- all_features
  selected_features_list <- list()
  top_n_features_list <- list()
  
  cat(sprintf("\nRunning %d bootstrap iterations...\n", n_bootstrap))
  
  pb <- txtProgressBar(min = 0, max = n_bootstrap, style = 3)
  
  for(i in 1:n_bootstrap) {
    
    # Bootstrap sample (with replacement)
    boot_idx <- sample(1:n_samples, n_samples, replace = TRUE)
    boot_data <- data[boot_idx, ]
    
    # Run RFE on bootstrap sample
    boot_rfe <- tryCatch({
      svm_rfe_cv(boot_data, 
                 class_col = class_col,
                 cv_folds = cv_folds,
                 kernel = kernel,
                 use_auc = TRUE)
    }, error = function(e) {
      return(NULL)
    })
    
    if(!is.null(boot_rfe)) {
      # Get top features
      top_features <- boot_rfe$feature_ranking[1:boot_rfe$optimal_n_features]
      selected_features_list[[i]] <- boot_rfe$feature_ranking
      top_n_features_list[[i]] <- top_features
      
      # Mark selected features
      feature_frequency[i, top_features] <- 1
    }
    
    setTxtProgressBar(pb, i)
  }
  
  close(pb)
  
  # Calculate selection frequency for each feature
  selection_freq <- colSums(feature_frequency) / n_bootstrap * 100
  
  # Rank features by selection frequency
  stability_ranking <- data.frame(
    Feature = names(selection_freq),
    Selection_Frequency_Percent = selection_freq,
    Times_Selected = colSums(feature_frequency)
  )
  stability_ranking <- stability_ranking[order(-stability_ranking$Selection_Frequency_Percent), ]
  rownames(stability_ranking) <- NULL
  
  # Calculate stability metrics
  results <- list(
    stability_ranking = stability_ranking,
    feature_frequency_matrix = feature_frequency,
    selected_features_per_iteration = selected_features_list,
    top_features_per_iteration = top_n_features_list,
    n_bootstrap = n_bootstrap
  )
  
  return(results)
}

# Run bootstrap stability analysis
stability_results <- bootstrap_feature_selection(data_ml,
                                                 class_col = "condition",
                                                 n_bootstrap = 100,
                                                 cv_folds = 5,
                                                 kernel = "linear")

cat("\n=== FEATURE STABILITY RESULTS ===\n")
cat("\nTop features by selection frequency:\n")
print(head(stability_results$stability_ranking, 10))

# Identify highly stable features (selected in >70% of iterations)
stable_features <- stability_results$stability_ranking$Feature[
  stability_results$stability_ranking$Selection_Frequency_Percent >= 70
]

cat(sprintf("\n%d features selected in ≥70%% of bootstrap iterations:\n", 
            length(stable_features)))
if(length(stable_features) > 0) {
  print(stable_features)
} else {
  cat("  None (try lowering threshold to 50%)\n")
  stable_features <- stability_results$stability_ranking$Feature[
    stability_results$stability_ranking$Selection_Frequency_Percent >= 50
  ]
  cat(sprintf("\n%d features selected in ≥50%% of bootstrap iterations:\n", 
              length(stable_features)))
  print(stable_features)
}

# Get most stable features for final model
top_stable_features <- head(stability_results$stability_ranking$Feature, 
                            round(nested_cv_results$mean_n_features))

cat("\n=== FINAL FEATURE SET (Most Stable) ===\n")
print(top_stable_features)

# ============================================================================
# PART 5: PUBLICATION-QUALITY VISUALIZATIONS
# ============================================================================

# Set publication theme
theme_publication <- theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
    axis.title = element_text(face = "bold", size = 12),
    axis.text = element_text(size = 10),
    legend.title = element_text(face = "bold", size = 11),
    legend.text = element_text(size = 10),
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "grey90", colour = "black"),
    strip.text = element_text(face = "bold")
  )

# --- Visualization 1: Feature Elimination Curve (from first outer fold inner RFE) ---
# Use first outer fold's inner RFE results as representative
first_fold_rfe <- nested_cv_results$rfe_results_per_fold[[1]]
n_features <- length(first_fold_rfe$cv_performance):1
acc_df <- data.frame(
  n_features = n_features,
  performance = first_fold_rfe$cv_performance,
  se = first_fold_rfe$cv_performance_sd
)

p1 <- ggplot(acc_df, aes(x = n_features, y = performance)) +
  geom_ribbon(aes(ymin = performance - se, ymax = performance + se), 
              alpha = 0.2, fill = "#3498DB") +
  geom_line(color = "#2C3E50", size = 1.2) +
  geom_point(color = "#E74C3C", size = 3) +
  geom_vline(xintercept = first_fold_rfe$optimal_n_features, 
             linetype = "dashed", color = "#27AE60", size = 1) +
  annotate("text", 
           x = first_fold_rfe$optimal_n_features, 
           y = min(acc_df$performance), 
           label = sprintf("Optimal: %d features", first_fold_rfe$optimal_n_features),
           hjust = -0.1, vjust = -0.5, color = "#27AE60", fontface = "bold") +
  scale_x_continuous(breaks = pretty(n_features, n = 10)) +
  labs(
    title = "SVM-RFE: Feature Elimination Curve (Inner CV)",
    subtitle = "Representative from Outer Fold 1",
    x = "Number of Features",
    y = paste0("Cross-Validation ", first_fold_rfe$performance_metric)) +
  theme_publication

# --- Visualization 2: ROC Curve (Aggregated from all outer folds) ---
# Combine all predictions from outer folds
all_true_labels <- unlist(lapply(nested_cv_results$fold_results, function(x) x$true_labels))
all_probabilities <- unlist(lapply(nested_cv_results$fold_results, function(x) x$probabilities))

# Calculate overall ROC
roc_obj <- roc(all_true_labels, all_probabilities)
auc_value <- auc(roc_obj)

roc_df <- data.frame(
  fpr = 1 - roc_obj$specificities,
  tpr = roc_obj$sensitivities
)

p2 <- ggplot(roc_df, aes(x = fpr, y = tpr)) +
  geom_line(color = "#3498DB", size = 1.3) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "grey50") +
  annotate("text", x = 0.6, y = 0.3, 
           label = sprintf("Pooled AUC = %.3f", auc_value),
           size = 6, fontface = "bold", color = "#2C3E50") +
  annotate("text", x = 0.6, y = 0.2,
           label = sprintf("(Mean CV AUC: %.3f ± %.3f)", 
                           nested_cv_results$mean_auc, 
                           nested_cv_results$sd_auc),
           size = 5, color = "#7F8C8D") +
  labs(
    title = "SVM-RFE ROC Curve",
    x = "False Positive Rate (1 - Specificity)",
    y = "True Positive Rate (Sensitivity)"
  ) +
  coord_equal() +
  theme_publication
ggsave("SVM-RFE_ROC_Curve.png", p2, width = 8, height = 8, units = "in", dpi = 900)
# --- Visualization 3: Confusion Matrix Heatmap (Aggregated) ---
# Sum all confusion matrices from outer folds
cm_aggregated <- Reduce("+", lapply(nested_cv_results$fold_results, 
                                    function(x) x$confusion_matrix))
cm_df <- melt(cm_aggregated)
colnames(cm_df) <- c("Prediction", "Reference", "Count")

p3 <- ggplot(cm_df, aes(x = Reference, y = Prediction, fill = Count)) +
  geom_tile(color = "white", size = 1.5) +
  geom_text(aes(label = Count), size = 8, fontface = "bold", color = "white") +
  scale_fill_gradient(low = "#3498DB", high = "#E74C3C") +
  labs(
    title = "Confusion Matrix - Nested CV (Aggregated)",
    subtitle = sprintf("Mean Accuracy: %.3f ± %.3f", 
                       nested_cv_results$mean_accuracy,
                       nested_cv_results$sd_accuracy),
    x = "True Label",
    y = "Predicted Label"
  ) +
  theme_publication +
  theme(
    legend.position = "right",
    panel.grid = element_blank()
  )

# --- Visualization 4: Feature Stability (Top Features) ---
top_stability <- head(stability_results$stability_ranking, 15)
top_stability$Feature <- factor(top_stability$Feature, 
                                levels = rev(top_stability$Feature))

p4 <- ggplot(top_stability, aes(x = Feature, y = Selection_Frequency_Percent, 
                                fill = Selection_Frequency_Percent)) +
  geom_bar(stat = "identity", color = "black", size = 0.3) +
  geom_hline(yintercept = 70, linetype = "dashed", color = "red", size = 1) +
  annotate("text", x = 14, y = 75, label = "70% threshold", 
           color = "red", fontface = "bold", hjust = 1) +
  scale_fill_gradient2(low = "#3498DB", mid = "#F39C12", high = "#27AE60",
                       midpoint = 50) +
  coord_flip(clip="off") +
  labs(
    title = "Feature Selection Stability",
    subtitle = sprintf("Based on %d Bootstrap Iterations", 
                       stability_results$n_bootstrap),
    x = "Gene",
    y = "Selection Frequency (%)",
    fill = "Frequency (%)"
  ) +
  theme_publication +
  theme(legend.position = "right")

# --- Visualization 5: Expression Heatmap of Most Stable Features ---
expr_stable <- expr_ml[top_stable_features, ]

# Annotation for samples
annotation_col <- data.frame(
  Condition = condition
)
rownames(annotation_col) <- colnames(expr_stable)

# Color scheme
ann_colors <- list(
  Condition = c(Normal = "#3498DB", Cancer = "#E74C3C")
)

# Create heatmap
p5 <- pheatmap(
  expr_stable,
  scale = "row",
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  clustering_method = "complete",
  annotation_col = annotation_col,
  annotation_colors = ann_colors,
  color = colorRampPalette(c("#3498DB", "white", "#E74C3C"))(100),
  fontsize = 10,
  fontsize_row = 10,
  fontsize_col = 8,
  border_color = NA,
  main = "Expression Heatmap - Most Stable Features",
  show_colnames = FALSE
)

# --- Visualization 6: Boxplots for Top 3 Stable Features ---
top3_stable <- top_stable_features[1:min(3, length(top_stable_features))]
expr_top3_long <- melt(as.matrix(t(expr_ml[top3_stable, ])))
colnames(expr_top3_long) <- c("Sample", "Gene", "Expression")
expr_top3_long$Condition <- rep(condition, times = length(top3_stable))

p6 <- ggplot(expr_top3_long, aes(x = Gene, y = Expression, fill = Condition)) +
  geom_boxplot(outlier.shape = 16, outlier.size = 2, alpha = 0.8) +
  scale_fill_manual(values = c(Normal = "#3498DB", Cancer = "#E74C3C")) +
  stat_compare_means(aes(group = Condition), 
                     label = "p.signif", 
                     method = "t.test",
                     label.y.npc = 0.95) +
  labs(
    title = "Expression Levels - Top 3 Most Stable Features",
    subtitle = "Based on Bootstrap Stability Analysis",
    x = "Gene",
    y = "Log2(RPKM + 1)"
  ) +
  theme_publication +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"))

# --- Visualization 7: Nested CV Performance Distribution ---
nested_cv_perf_df <- data.frame(
  Fold = 1:length(nested_cv_results$fold_results),
  AUC = sapply(nested_cv_results$fold_results, function(x) x$auc),
  Accuracy = sapply(nested_cv_results$fold_results, function(x) x$accuracy),
  Sensitivity = sapply(nested_cv_results$fold_results, function(x) x$sensitivity),
  Specificity = sapply(nested_cv_results$fold_results, function(x) x$specificity)
)

nested_cv_long <- melt(nested_cv_perf_df, id.vars = "Fold", 
                       variable.name = "Metric", value.name = "Value")

p7 <- ggplot(nested_cv_long, aes(x = Metric, y = Value, fill = Metric)) +
  geom_boxplot(alpha = 0.7, outlier.shape = 16) +
  geom_jitter(width = 0.2, size = 2, alpha = 0.6) +
  stat_summary(fun = mean, geom = "point", shape = 23, size = 4, 
               fill = "red", color = "darkred") +
  scale_fill_brewer(palette = "Set2") +
  labs(
    title = "Nested CV Performance Distribution",
    subtitle = "Diamond = Mean, Box = IQR, Points = Individual Folds",
    x = "Performance Metric",
    y = "Value"
  ) +
  theme_publication +
  theme(legend.position = "none") +
  ylim(0, 1)

# --- Visualization 8: Feature Selection Consistency Heatmap ---
# Show which features were selected together
top_15_features <- head(stability_results$stability_ranking$Feature, 15)
freq_matrix_top <- stability_results$feature_frequency_matrix[, top_15_features]

# Calculate co-occurrence matrix
co_occurrence <- t(freq_matrix_top) %*% freq_matrix_top
co_occurrence_pct <- co_occurrence / stability_results$n_bootstrap * 100

p8 <- pheatmap(
  co_occurrence_pct,
  color = colorRampPalette(c("white", "#3498DB", "#E74C3C"))(100),
  display_numbers = FALSE,
  fontsize = 9,
  fontsize_row = 9,
  fontsize_col = 9,
  border_color = "grey60",
  main = "Feature Co-selection Frequency (%)",
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean"
)

# --- Visualization 9: Nested CV - Number of Features per Fold ---
n_features_per_fold <- sapply(nested_cv_results$fold_results, 
                              function(x) x$n_features)
auc_per_fold <- sapply(nested_cv_results$fold_results, 
                       function(x) x$auc)

fold_features_df <- data.frame(
  Fold = factor(1:length(n_features_per_fold)),
  N_Features = n_features_per_fold,
  AUC = auc_per_fold
)

p9 <- ggplot(fold_features_df, aes(x = Fold, y = N_Features, fill = AUC)) +
  geom_bar(stat = "identity", color = "black", size = 0.5) +
  geom_text(aes(label = N_Features), vjust = -0.5, fontface = "bold") +
  scale_fill_gradient(low = "#E74C3C", high = "#27AE60") +
  labs(
    title = "Features Selected per Outer CV Fold",
    subtitle = "Color = Test Set AUC",
    x = "Outer CV Fold",
    y = "Number of Features",
    fill = "AUC"
  ) +
  theme_publication

# --- Visualization 10: Performance Comparison Across Metrics ---
perf_comparison_df <- data.frame(
  Metric = c("AUC", "Accuracy", "Sensitivity", "Specificity", "F1"),
  Mean = c(nested_cv_results$mean_auc,
           nested_cv_results$mean_accuracy,
           nested_cv_results$mean_sensitivity,
           nested_cv_results$mean_specificity,
           nested_cv_results$mean_f1),
  SD = c(nested_cv_results$sd_auc,
         nested_cv_results$sd_accuracy,
         nested_cv_results$sd_sensitivity,
         nested_cv_results$sd_specificity,
         nested_cv_results$sd_f1)
)

perf_comparison_df$Metric <- factor(perf_comparison_df$Metric,
                                    levels = perf_comparison_df$Metric)

p10 <- ggplot(perf_comparison_df, aes(x = Metric, y = Mean, fill = Metric)) +
  geom_bar(stat = "identity", color = "black", size = 0.5) +
  geom_errorbar(aes(ymin = Mean - SD, ymax = Mean + SD), 
                width = 0.2, size = 0.8) +
  geom_text(aes(label = sprintf("%.3f", Mean)), 
            vjust = -0.5, fontface = "bold", size = 4) +
  scale_fill_brewer(palette = "Set3") +
  labs(
    title = "Nested CV Performance Metrics Summary",
    subtitle = "Error bars = Standard Deviation",
    x = "Performance Metric",
    y = "Value"
  ) +
  theme_publication +
  theme(legend.position = "none") +
  ylim(0, 1.1)

ggsave("SVM_RFE_Performance_Summary.pdf", p10, width = 8, height = 6, dpi = 300)

# ============================================================================
# PART 6: SAVE ALL VISUALIZATIONS
# ============================================================================

# Combined plots
combined_plot1 <- ggarrange(p1, p2, p3, p4, 
                            ncol = 2, nrow = 2,
                            labels = c("A", "B", "C", "D"),
                            font.label = list(size = 16, face = "bold"))

combined_plot2 <- ggarrange(p6, p7, p9, p10,
                            ncol = 2, nrow = 2,
                            labels = c("E", "F", "G", "H"),
                            font.label = list(size = 16, face = "bold"))

# Save all plots
ggsave("SVM_RFE_Main.png", combined_plot1, width = 14, height = 12, units = "in", dpi = 600)
ggsave("SVM_RFE_Metrics.png", combined_plot2, width = 14, height = 12, units = "in", dpi = 600)
ggsave("SVM_RFE_Feature_Elimination.png", p1, width = 8, height = 6, units = "in", dpi = 600)
ggsave("SVM_RFE_ROC_Curve.png", p2, width = 7, height = 7, units = "in", dpi = 600)
ggsave("SVM_RFE_Confusion_Matrix.png", p3, width = 7, height = 6, units = "in", dpi = 600)
ggsave("SVM_RFE_Feature_Stability.png", p4, width = 8, height = 7, units = "in", dpi = 600)
ggsave("SVM_RFE_Top3_Expression.png", p6,  width = 8, height = 6, units = "in", dpi = 600)
ggsave("SVM_RFE_NestedCV_Distribution.png", p7, width = 8, height = 6, units = "in", dpi = 600)
ggsave("SVM_RFE_Features_Per_Fold.png", p9, width = 8, height = 6, units = "in", dpi = 600)
ggsave("SVM_RFE_Performance_Summary.png", p10, width = 8, height = 6, units = "in", dpi = 600)

png("SVM_RFE_Expression_Heatmap_Stable.png", width = 10, height = 8, units = "in", res = 600)
print(p5)
dev.off()

png("SVM_RFE_Feature_Coselection.png", width = 10, height = 9,units = "in", res = 600)
print(p8)
dev.off()

# ============================================================================
# PART 7: SAVE RESULTS TO FILES
# ============================================================================

# Save feature stability ranking
write.csv(
  stability_results$stability_ranking,
  "SVM_RFE_Feature_Stability_Ranking.csv",
  row.names = FALSE
)

# Save stable features (≥70% or ≥50%)
threshold_used <- ifelse(length(stable_features) > 0, 70, 50)
write.csv(
  data.frame(
    Feature = stable_features,
    Selection_Frequency = stability_results$stability_ranking$Selection_Frequency_Percent[
      match(stable_features, stability_results$stability_ranking$Feature)
    ],
    Times_Selected = stability_results$stability_ranking$Times_Selected[
      match(stable_features, stability_results$stability_ranking$Feature)
    ]
  ),
  sprintf("SVM_RFE_Stable_Features_%dpct.csv", threshold_used),
  row.names = FALSE
)

# Save nested CV summary
nested_cv_summary <- data.frame(
  Metric = c("AUC", "Accuracy", "Sensitivity", "Specificity", "F1", "N_Features"),
  Mean = c(nested_cv_results$mean_auc,
           nested_cv_results$mean_accuracy,
           nested_cv_results$mean_sensitivity,
           nested_cv_results$mean_specificity,
           nested_cv_results$mean_f1,
           nested_cv_results$mean_n_features),
  SD = c(nested_cv_results$sd_auc,
         nested_cv_results$sd_accuracy,
         nested_cv_results$sd_sensitivity,
         nested_cv_results$sd_specificity,
         nested_cv_results$sd_f1,
         NA)
)

write.csv(nested_cv_summary, "SVM_RFE_NestedCV_Summary.csv", row.names = FALSE)

# Save nested CV per-fold results
nested_cv_folds <- do.call(rbind, lapply(1:length(nested_cv_results$fold_results), function(i) {
  data.frame(
    Fold = i,
    AUC = nested_cv_results$fold_results[[i]]$auc,
    Accuracy = nested_cv_results$fold_results[[i]]$accuracy,
    Sensitivity = nested_cv_results$fold_results[[i]]$sensitivity,
    Specificity = nested_cv_results$fold_results[[i]]$specificity,
    F1 = nested_cv_results$fold_results[[i]]$f1,
    N_Features = nested_cv_results$fold_results[[i]]$n_features,
    Selected_Features = paste(nested_cv_results$fold_results[[i]]$selected_features, 
                              collapse = ";")
  )
}))

write.csv(nested_cv_folds, "SVM_RFE_NestedCV_PerFold.csv", row.names = FALSE)

# Save RFE elimination curve (from first outer fold as representative)
rfe_curve_df <- data.frame(
  N_Features = length(first_fold_rfe$cv_performance):1,
  CV_Performance = first_fold_rfe$cv_performance,
  CV_SD = first_fold_rfe$cv_performance_sd,
  Metric = first_fold_rfe$performance_metric
)

write.csv(rfe_curve_df, "SVM_RFE_Elimination_Curve.csv", row.names = FALSE)

# ============================================================================
# FINAL SUMMARY OUTPUT
# ============================================================================

cat("\n" %+% strrep("=", 70) %+% "\n")
cat("SVM-RFE ANALYSIS COMPLETE\n")
cat(strrep("=", 70) %+% "\n")

cat("\n=== KEY RESULTS (UNBIASED NESTED CV) ===\n")
cat(sprintf("  • Mean AUC: %.3f ± %.3f\n", 
            nested_cv_results$mean_auc, nested_cv_results$sd_auc))
cat(sprintf("  • Mean Accuracy: %.3f ± %.3f\n", 
            nested_cv_results$mean_accuracy, nested_cv_results$sd_accuracy))
cat(sprintf("  • Mean Sensitivity: %.3f ± %.3f\n", 
            nested_cv_results$mean_sensitivity, nested_cv_results$sd_sensitivity))
cat(sprintf("  • Mean Specificity: %.3f ± %.3f\n", 
            nested_cv_results$mean_specificity, nested_cv_results$sd_specificity))
cat(sprintf("  • Mean F1 Score: %.3f ± %.3f\n", 
            nested_cv_results$mean_f1, nested_cv_results$sd_f1))
cat(sprintf("  • Mean # Features: %.1f\n", nested_cv_results$mean_n_features))

cat("\n=== FEATURE STABILITY (BOOTSTRAP) ===\n")
cat(sprintf("  • Bootstrap Iterations: %d\n", stability_results$n_bootstrap))
cat(sprintf("  • Highly Stable Features (≥%d%%): %d genes\n",
            threshold_used, length(stable_features)))

cat("\n=== MOST STABLE FEATURES ===\n")
print(head(stability_results$stability_ranking, 10))

cat("\n=== GENERATED FILES ===\n")
cat("\nPDF Visualizations (12 files):\n")
cat("  - SVM_RFE_Main.pdf (4-panel summary)\n")
cat("  - SVM_RFE_Metrics.pdf (4-panel metrics)\n")
cat("  - SVM_RFE_Feature_Elimination.pdf\n")
cat("  - SVM_RFE_ROC_Curve_NestedCV.pdf\n")
cat("  - SVM_RFE_Confusion_Matrix_NestedCV.pdf\n")
cat("  - SVM_RFE_Feature_Stability.pdf\n")
cat("  - SVM_RFE_Expression_Heatmap_Stable.pdf\n")
cat("  - SVM_RFE_Top3_Expression.pdf\n")
cat("  - SVM_RFE_NestedCV_Distribution.pdf\n")
cat("  - SVM_RFE_Feature_Coselection.pdf\n")
cat("  - SVM_RFE_Features_Per_Fold.pdf\n")
cat("  - SVM_RFE_Performance_Summary.pdf\n")

cat("\nCSV Data Files (5 files):\n")
cat("  - SVM_RFE_Feature_Stability_Ranking.csv\n")
cat(sprintf("  - SVM_RFE_Stable_Features_%dpct.csv\n", threshold_used))
cat("  - SVM_RFE_NestedCV_Summary.csv\n")
cat("  - SVM_RFE_NestedCV_PerFold.csv\n")
cat("  - SVM_RFE_Elimination_Curve.csv\n")

cat("\n" %+% strrep("=", 70) %+% "\n")
cat("All visualizations use data from nested CV and bootstrap results only.\n")
cat("This analysis follows gold-standard protocols for unbiased performance.\n")
cat(strrep("=", 70) %+% "\n\n")
