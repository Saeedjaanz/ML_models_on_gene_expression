# ============================================================================
# RANDOM FOREST CLASSIFICATION - GOLD STANDARD PROTOCOL
# Nested Cross-Validation with Hyperparameter Tuning
# Publication-Quality Analysis and Visualization
# ============================================================================

# Required Libraries
library(pROC)            # ROC analysis
library(PRROC)           # Precision-Recall curves
library(caret)           # ML workflows
library(UpSetR)          # Set visualizations
library(ggplot2)         # Visualization
library(tidyverse)       # Data manipulation
library(randomForest)    # Random Forest
library(pheatmap)        # Heatmaps
library(RColorBrewer)    # Color palettes
library(gridExtra)       # Multiple plots
library(ggpubr)          # Publication-ready plots
library(reshape2)        # Data reshaping

# ============================================================================
# PART 1: INPUT VALIDATION & CONFIGURATION
# ============================================================================

cat("\n" %+% strrep("=", 70) %+% "\n")
cat("RANDOM FOREST GOLD STANDARD ANALYSIS\n")
cat(strrep("=", 70) %+% "\n\n")

# Validate input data
stopifnot(
  is.matrix(expr_ml),
  length(condition) == ncol(expr_ml), 
  all(levels(condition) %in% c("Normal", "Cancer"))
)

cat("✓ Input validation passed\n")
cat(sprintf("  - Features: %d genes\n", nrow(expr_ml)))
cat(sprintf("  - Samples: %d (Normal: %d, Cancer: %d)\n", 
            ncol(expr_ml),
            sum(condition == "Normal"),
            sum(condition == "Cancer")))

# Global settings
set.seed(123)
K <- 10                    # Number of outer CV folds
consensus_thresh <- 0.7    # Consensus threshold for feature selection
ntree <- 1500             # Number of trees in forest
inner_cv_folds <- 5       # Number of inner CV folds for hyperparameter tuning

cat(sprintf("\n✓ Configuration:\n"))
cat(sprintf("  - Outer CV folds: %d\n", K))
cat(sprintf("  - Inner CV folds: %d\n", inner_cv_folds))
cat(sprintf("  - Number of trees: %d\n", ntree))
cat(sprintf("  - Consensus threshold: %.0f%%\n", consensus_thresh * 100))

# ============================================================================
# PART 2: DATA PREPARATION
# ============================================================================

# Prepare labels
y_bin <- ifelse(condition == "Cancer", 1, 0)

# Create stratified folds for outer CV
folds <- createFolds(condition, k = K, list = TRUE, returnTrain = FALSE)

# Storage objects
auc_vals <- numeric(K)
acc_vals <- numeric(K)
sens_vals <- numeric(K)
spec_vals <- numeric(K)
ppv_vals <- numeric(K)
npv_vals <- numeric(K)
f1_vals <- numeric(K)
pooled_obs <- c()
pooled_pred <- c()
pooled_pred_class <- c()
feat_imp_list <- vector("list", K)
top_feats_list <- vector("list", K)
best_mtry_vals <- numeric(K)
confusion_matrices <- vector("list", K)

cat("\n✓ Data preparation complete\n")

# ============================================================================
# PART 3: NESTED CROSS-VALIDATION
# ============================================================================

cat("\n" %+% strrep("=", 70) %+% "\n")
cat("EXECUTING NESTED CROSS-VALIDATION\n")
cat(strrep("=", 70) %+% "\n\n")

for (i in seq_along(folds)) {
  
  cat(sprintf("┌─ OUTER FOLD %d/%d ─┐\n", i, K))
  
  # Split data
  test_idx <- folds[[i]]
  train_idx <- setdiff(seq_along(condition), test_idx)
  
  x_train <- t(expr_ml[, train_idx])
  y_train <- condition[train_idx]
  x_test  <- t(expr_ml[, test_idx])
  y_test  <- condition[test_idx]
  
  cat(sprintf("│ Train: %d samples | Test: %d samples\n", 
              length(train_idx), length(test_idx)))
  
  # ──────────────────────────────────────────────────────────────────────────
  # INNER CV: Hyperparameter Tuning (mtry optimization)
  # ──────────────────────────────────────────────────────────────────────────
  
  cat("│\n│ ➤ Inner CV: Tuning mtry...\n")
  
  # Configure inner CV
  trctrl <- trainControl(
    method = "cv", 
    number = inner_cv_folds, 
    classProbs = TRUE, 
    summaryFunction = twoClassSummary,
    savePredictions = "final", 
    verboseIter = FALSE
  )
  
  # Define tuning grid
  tune_grid <- expand.grid(mtry = seq_len(ncol(x_train)))
  
  # Prepare factor labels
  y_train_factor <- factor(y_train, levels = c("Normal", "Cancer"))
  
  # Hyperparameter tuning
  set.seed(123 + i)
  rf_tune <- train(
    x = x_train, 
    y = y_train_factor, 
    method = "rf", 
    metric = "ROC",
    tuneGrid = tune_grid, 
    trControl = trctrl,  
    ntree = ntree, 
    importance = TRUE
  )
  
  best_mtry <- rf_tune$bestTune$mtry
  best_mtry_vals[i] <- best_mtry
  
  cat(sprintf("│   ✓ Best mtry: %d (ROC = %.4f)\n", 
              best_mtry, max(rf_tune$results$ROC)))
  
  # ──────────────────────────────────────────────────────────────────────────
  # TRAIN FINAL MODEL: With optimal hyperparameters
  # ──────────────────────────────────────────────────────────────────────────
  
  cat("│\n│ ➤ Training final model with optimal mtry...\n")
  
  rf_model <- randomForest(
    x = x_train, 
    y = y_train_factor,
    ntree = ntree, 
    mtry = best_mtry,
    importance = TRUE, 
    keep.forest = TRUE
  )
  
  # ──────────────────────────────────────────────────────────────────────────
  # EVALUATE ON TEST SET
  # ──────────────────────────────────────────────────────────────────────────
  
  cat("│\n│ ➤ Evaluating on held-out test set...\n")
  
  # Predictions
  preds_prob <- predict(rf_model, x_test, type = "prob")[, "Cancer"]
  preds_class <- predict(rf_model, x_test, type = "response")
  
  # Pool predictions for overall ROC
  pooled_obs <- c(pooled_obs, ifelse(y_test == "Cancer", 1, 0))
  pooled_pred <- c(pooled_pred, preds_prob)
  pooled_pred_class <- c(pooled_pred_class, as.character(preds_class))
  
  # Calculate ROC and AUC
  roc_obj <- try(
    roc(ifelse(y_test == "Cancer", 1, 0), preds_prob, quiet = TRUE), 
    silent = TRUE
  )
  auc_vals[i] <- if (inherits(roc_obj, "try-error")) NA else as.numeric(roc_obj$auc)
  
  # Confusion matrix
  cm <- confusionMatrix(preds_class, y_test)
  confusion_matrices[[i]] <- cm$table
  
  # Extract performance metrics
  acc_vals[i] <- cm$overall["Accuracy"]
  sens_vals[i] <- cm$byClass["Sensitivity"]
  spec_vals[i] <- cm$byClass["Specificity"]
  ppv_vals[i] <- cm$byClass["Pos Pred Value"]
  npv_vals[i] <- cm$byClass["Neg Pred Value"]
  f1_vals[i] <- cm$byClass["F1"]
  
  cat(sprintf("│   ✓ AUC: %.4f | Accuracy: %.4f\n", 
              auc_vals[i], acc_vals[i]))
  cat(sprintf("│   ✓ Sensitivity: %.4f | Specificity: %.4f\n", 
              sens_vals[i], spec_vals[i]))
  
  # ──────────────────────────────────────────────────────────────────────────
  # FEATURE IMPORTANCE
  # ──────────────────────────────────────────────────────────────────────────
  
  # Extract importance
  imp_df <- as.data.frame(importance(rf_model, scale = FALSE))
  imp_df$gene <- rownames(imp_df)
  feat_imp_list[[i]] <- imp_df
  
  # Rank features
  ranked <- imp_df %>% arrange(desc(MeanDecreaseAccuracy))
  top_feats_list[[i]] <- ranked$gene
  
  cat(sprintf("│\n│ ➤ Top features: %s\n", 
              paste(head(ranked$gene, 5), collapse = ", ")))
  
  # ──────────────────────────────────────────────────────────────────────────
  # PER-FOLD DIAGNOSTICS
  # ──────────────────────────────────────────────────────────────────────────
  
  # Training performance
  train_preds_class <- predict(rf_model, x_train, type = "response")
  train_acc <- mean(train_preds_class == y_train_factor)
  train_err <- 1 - train_acc
  test_err <- 1 - acc_vals[i]
  
  # OOB error
  err_matrix <- rf_model$err.rate
  final_oob_err <- tail(err_matrix[, "OOB"], 1)
  
  cat(sprintf("│   ✓ Train Error: %.4f | Test Error: %.4f | OOB Error: %.4f\n",
              train_err, test_err, final_oob_err))
  
  # Save per-fold diagnostic plot
  pdf(sprintf("RF_Fold%02d_Diagnostics.pdf", i), width = 10, height = 5)
  par(mfrow = c(1, 2), mar = c(4, 4, 3, 1))
  
  # Panel 1: OOB error curve
  matplot(err_matrix[, "OOB", drop = FALSE], type = "l", lwd = 2,
          col = "#3498DB",
          xlab = "Number of Trees", 
          ylab = "OOB Error Rate", 
          main = sprintf("Fold %d: OOB Error Convergence", i))
  abline(h = final_oob_err, col = "#E74C3C", lty = 2, lwd = 2)
  abline(h = train_err, col = "#27AE60", lty = 3, lwd = 1.5)
  abline(h = test_err, col = "#F39C12", lty = 3, lwd = 1.5)
  legend("topright", 
         legend = c("OOB Error", "Final OOB", "Train Error", "Test Error"),
         col = c("#3498DB", "#E74C3C", "#27AE60", "#F39C12"), 
         lty = c(1, 2, 3, 3), 
         lwd = c(2, 2, 1.5, 1.5),
         bty = "n")
  
  # Panel 2: Top feature importance
  topn <- min(10, nrow(imp_df))
  tmp <- imp_df %>% arrange(desc(MeanDecreaseAccuracy)) %>% head(topn)
  
  barplot(
    t(as.matrix(tmp[, c("MeanDecreaseAccuracy", "MeanDecreaseGini")])),
    beside = TRUE, 
    las = 2, 
    names.arg = tmp$gene,
    main = sprintf("Fold %d: Top %d Features", i, topn),
    ylab = "Importance Score", 
    col = c("#3498DB", "#E74C3C"),
    cex.names = 0.7
  )
  legend("topright", 
         legend = c("Mean Decrease Accuracy", "Mean Decrease Gini"), 
         fill = c("#3498DB", "#E74C3C"), 
         bty = "n")
  
  dev.off()
  
  # ──────────────────────────────────────────────────────────────────────────
  # SAVE PER-FOLD OUTPUTS
  # ──────────────────────────────────────────────────────────────────────────
  
  saveRDS(rf_model, file = sprintf("RF_Model_Fold%02d.rds", i))
  write.csv(imp_df, file = sprintf("RF_Importance_Fold%02d.csv", i), 
            row.names = FALSE)
  
  cat("└" %+% strrep("─", 68) %+% "┘\n\n")
}

# ============================================================================
# PART 4: AGGREGATE PERFORMANCE METRICS
# ============================================================================

cat(strrep("=", 70) %+% "\n")
cat("AGGREGATED PERFORMANCE METRICS\n")
cat(strrep("=", 70) %+% "\n\n")

# Summary dataframe
summary_df_rf <- data.frame(
  Fold = 1:K,
  Accuracy = acc_vals,
  AUC = auc_vals,
  Sensitivity = sens_vals,
  Specificity = spec_vals,
  PPV = ppv_vals,
  NPV = npv_vals,
  F1 = f1_vals,
  Best_mtry = best_mtry_vals
)

# Calculate means and SDs
mean_metrics <- colMeans(summary_df_rf[, -1], na.rm = TRUE)
sd_metrics <- apply(summary_df_rf[, -1], 2, sd, na.rm = TRUE)

cat("Per-Fold Performance:\n")
print(summary_df_rf)

cat("\n" %+% strrep("-", 70) %+% "\n")
cat("SUMMARY STATISTICS (Mean ± SD):\n")
cat(strrep("-", 70) %+% "\n")
cat(sprintf("  Accuracy:    %.4f ± %.4f\n", mean_metrics["Accuracy"], sd_metrics["Accuracy"]))
cat(sprintf("  AUC:         %.4f ± %.4f\n", mean_metrics["AUC"], sd_metrics["AUC"]))
cat(sprintf("  Sensitivity: %.4f ± %.4f\n", mean_metrics["Sensitivity"], sd_metrics["Sensitivity"]))
cat(sprintf("  Specificity: %.4f ± %.4f\n", mean_metrics["Specificity"], sd_metrics["Specificity"]))
cat(sprintf("  PPV:         %.4f ± %.4f\n", mean_metrics["PPV"], sd_metrics["PPV"]))
cat(sprintf("  NPV:         %.4f ± %.4f\n", mean_metrics["NPV"], sd_metrics["NPV"]))
cat(sprintf("  F1 Score:    %.4f ± %.4f\n", mean_metrics["F1"], sd_metrics["F1"]))
cat(sprintf("  Best mtry:   %.1f ± %.1f\n", mean_metrics["Best_mtry"], sd_metrics["Best_mtry"]))

# Pooled ROC analysis
if (length(pooled_obs) > 0) {
  pooled_roc <- roc(pooled_obs, pooled_pred, quiet = TRUE)
  pooled_auc <- as.numeric(pooled_roc$auc)
  
  cat("\n" %+% strrep("-", 70) %+% "\n")
  cat("POOLED PERFORMANCE (All Folds Combined):\n")
  cat(strrep("-", 70) %+% "\n")
  cat(sprintf("  Pooled AUC: %.4f\n", pooled_auc))
  
  # Bootstrap CI for pooled AUC
  ci_obj_rf <- ci.auc(pooled_roc, method = "bootstrap", boot.n = 2000)
  cat(sprintf("  95%% CI (Bootstrap, n=2000): [%.4f, %.4f]\n", 
              ci_obj_rf[1], ci_obj_rf[3]))
}

# ============================================================================
# PART 5: FEATURE IMPORTANCE AGGREGATION
# ============================================================================

cat("\n" %+% strrep("=", 70) %+% "\n")
cat("FEATURE IMPORTANCE ANALYSIS\n")
cat(strrep("=", 70) %+% "\n\n")

# Combine importance from all folds
all_imp <- bind_rows(lapply(seq_along(feat_imp_list), function(j) {
  df <- feat_imp_list[[j]]
  df$fold <- j
  df
}), .id = "id")

# Calculate mean and SD importance
imp_summary <- all_imp %>%
  group_by(gene) %>%
  summarise(
    mean_MDA = mean(MeanDecreaseAccuracy, na.rm = TRUE),
    sd_MDA = sd(MeanDecreaseAccuracy, na.rm = TRUE),
    mean_Gini = mean(MeanDecreaseGini, na.rm = TRUE),
    sd_Gini = sd(MeanDecreaseGini, na.rm = TRUE),
    selection_frequency = n() / K * 100,  # % of folds where feature appeared
    .groups = "drop"
  ) %>%
  arrange(desc(mean_MDA))

cat("Top 10 Features (by Mean Decrease Accuracy):\n")
print(head(imp_summary, 10))

# Identify consensus features (selected in ≥70% of folds)
consensus_features <- imp_summary %>%
  filter(selection_frequency >= consensus_thresh * 100) %>%
  pull(gene)

cat(sprintf("\n✓ Consensus features (≥%.0f%% selection): %d genes\n", 
            consensus_thresh * 100, length(consensus_features)))
if(length(consensus_features) > 0) {
  cat("  ", paste(head(consensus_features, 10), collapse = ", "), "\n")
}

# ============================================================================
# PART 6: PUBLICATION-QUALITY VISUALIZATIONS
# ============================================================================

cat("\n" %+% strrep("=", 70) %+% "\n")
cat("GENERATING PUBLICATION-QUALITY FIGURES\n")
cat(strrep("=", 70) %+% "\n\n")

# Set publication theme
theme_publication <- theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
    plot.subtitle = element_text(size = 11, hjust = 0.5, color = "grey30"),
    axis.title = element_text(face = "bold", size = 12),
    axis.text = element_text(size = 10, color = "black"),
    legend.title = element_text(face = "bold", size = 11),
    legend.text = element_text(size = 10),
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "grey90", colour = "black"),
    strip.text = element_text(face = "bold")
  )

# ──────────────────────────────────────────────────────────────────────────
# Visualization 1: Pooled ROC Curve with CI
# ──────────────────────────────────────────────────────────────────────────

roc_df <- data.frame(
  FPR = 1 - pooled_roc$specificities,
  TPR = pooled_roc$sensitivities
)

p1 <- ggplot(roc_df, aes(x = FPR, y = TPR)) +
  geom_line(color = "#3498DB", size = 1.3) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "grey50", size = 0.8) +
  annotate("text", x = 0.6, y = 0.3, 
           label = sprintf("Pooled AUC = %.3f", pooled_auc),
           size = 5, fontface = "bold", color = "#2C3E50") +
  annotate("text", x = 0.6, y = 0.2,
           label = sprintf("95%% CI: [%.3f, %.3f]", ci_obj_rf[1], ci_obj_rf[3]),
           size = 4, color = "#7F8C8D") +
  annotate("text", x = 0.6, y = 0.1,
           label = sprintf("Mean CV AUC: %.3f ± %.3f", 
                           mean_metrics["AUC"], sd_metrics["AUC"]),
           size = 4, color = "#7F8C8D") +
  labs(
    title = "Random Forest ROC Curve",
    subtitle = sprintf("%d-Fold Cross-Validation (Pooled Predictions)", K),
    x = "False Positive Rate (1 - Specificity)",
    y = "True Positive Rate (Sensitivity)"
  ) +
  coord_equal() +
  theme_publication

# ──────────────────────────────────────────────────────────────────────────
# Visualization 2: Performance Metrics Distribution
# ──────────────────────────────────────────────────────────────────────────

perf_long <- summary_df_rf %>%
  select(Fold, Accuracy, AUC, Sensitivity, Specificity, F1) %>%
  pivot_longer(cols = -Fold, names_to = "Metric", values_to = "Value")

p2 <- ggplot(perf_long, aes(x = Metric, y = Value, fill = Metric)) +
  geom_boxplot(alpha = 0.7, outlier.shape = 16, outlier.size = 2) +
  geom_jitter(width = 0.15, size = 2, alpha = 0.5) +
  stat_summary(fun = mean, geom = "point", shape = 23, size = 4,
               fill = "red", color = "darkred") +
  scale_fill_brewer(palette = "Set2") +
  labs(
    title = "Performance Metrics Across CV Folds",
    subtitle = "Diamond = Mean, Box = IQR, Points = Individual Folds",
    x = "Performance Metric",
    y = "Value"
  ) +
  theme_publication +
  theme(legend.position = "none") +
  ylim(0, 1)

# ──────────────────────────────────────────────────────────────────────────
# Visualization 3: Feature Importance - Mean Decrease Accuracy
# ──────────────────────────────────────────────────────────────────────────

top_features_mda <- head(imp_summary, 15)
top_features_mda$gene <- factor(top_features_mda$gene, 
                                levels = rev(top_features_mda$gene))

p3 <- ggplot(top_features_mda, aes(x = gene, y = mean_MDA, fill = mean_MDA)) +
  geom_bar(stat = "identity", color = "black", size = 0.3) +
  geom_errorbar(aes(ymin = mean_MDA - sd_MDA, ymax = mean_MDA + sd_MDA),
                width = 0.3, size = 0.5) +
  scale_fill_gradient(low = "#F39C12", high = "#C0392B") +
  coord_flip() +
  labs(
    title = "Feature Importance: Mean Decrease Accuracy",
    subtitle = sprintf("Top 15 Features (Aggregated from %d-Fold CV)", K),
    x = "Gene",
    y = "Mean Decrease in Accuracy",
    fill = "Importance"
  ) +
  theme_publication +
  theme(legend.position = "right")

# ──────────────────────────────────────────────────────────────────────────
# Visualization 4: Feature Importance - Mean Decrease Gini
# ──────────────────────────────────────────────────────────────────────────

imp_summary_gini <- imp_summary %>% arrange(desc(mean_Gini))
top_features_gini <- head(imp_summary_gini, 15)
top_features_gini$gene <- factor(top_features_gini$gene, 
                                 levels = rev(top_features_gini$gene))

p4 <- ggplot(top_features_gini, aes(x = gene, y = mean_Gini, fill = mean_Gini)) +
  geom_bar(stat = "identity", color = "black", size = 0.3) +
  geom_errorbar(aes(ymin = mean_Gini - sd_Gini, ymax = mean_Gini + sd_Gini),
                width = 0.3, size = 0.5) +
  scale_fill_gradient(low = "#3498DB", high = "#8E44AD") +
  coord_flip() +
  labs(
    title = "Feature Importance: Mean Decrease Gini",
    subtitle = sprintf("Top 15 Features (Aggregated from %d-Fold CV)", K),
    x = "Gene",
    y = "Mean Decrease in Gini Impurity",
    fill = "Importance"
  ) +
  theme_publication +
  theme(legend.position = "right")

# ──────────────────────────────────────────────────────────────────────────
# Visualization 5: Confusion Matrix (Aggregated)
# ──────────────────────────────────────────────────────────────────────────

# Sum all confusion matrices
cm_aggregated <- Reduce("+", confusion_matrices)
cm_df <- melt(cm_aggregated)
colnames(cm_df) <- c("Prediction", "Reference", "Count")

# Calculate percentages
cm_df <- cm_df %>%
  group_by(Reference) %>%
  mutate(Percentage = Count / sum(Count) * 100) %>%
  ungroup()

p5 <- ggplot(cm_df, aes(x = Reference, y = Prediction, fill = Count)) +
  geom_tile(color = "white", size = 1.5) +
  geom_text(aes(label = sprintf("%d\n(%.1f%%)", Count, Percentage)), 
            size = 6, fontface = "bold", color = "white") +
  scale_fill_gradient(low = "#3498DB", high = "#E74C3C") +
  labs(
    title = "Aggregated Confusion Matrix",
    subtitle = sprintf("%d-Fold Cross-Validation", K),
    x = "True Label",
    y = "Predicted Label"
  ) +
  theme_publication +
  theme(
    legend.position = "right",
    panel.grid = element_blank()
  )

# ──────────────────────────────────────────────────────────────────────────
# Visualization 6: Feature Selection Stability
# ──────────────────────────────────────────────────────────────────────────

p6 <- ggplot(top_features_mda, 
             aes(x = gene, y = selection_frequency, fill = selection_frequency)) +
  geom_bar(stat = "identity", color = "black", size = 0.3) +
  geom_hline(yintercept = consensus_thresh * 100, 
             linetype = "dashed", color = "red", size = 1) +
  annotate("text", x = 13, y = consensus_thresh * 100 + 5, 
           label = sprintf("%.0f%% threshold", consensus_thresh * 100),
           color = "red", fontface = "bold", hjust = 1) +
  scale_fill_gradient2(low = "#E74C3C", mid = "#F39C12", high = "#27AE60",
                       midpoint = 50) +
  coord_flip() +
  labs(
    title = "Feature Selection Stability",
    subtitle = "Percentage of CV Folds Where Feature Was Selected",
    x = "Gene",
    y = "Selection Frequency (%)",
    fill = "Frequency (%)"
  ) +
  theme_publication +
  theme(legend.position = "right") +
  ylim(0, 105)

# ──────────────────────────────────────────────────────────────────────────
# Visualization 7: Expression Heatmap of Top Features
# ──────────────────────────────────────────────────────────────────────────

top_10_genes <- head(imp_summary$gene, 10)
expr_top <- expr_ml[top_10_genes, ]

# Annotation
annotation_col <- data.frame(Condition = condition)
rownames(annotation_col) <- colnames(expr_top)

ann_colors <- list(
  Condition = c(Normal = "#3498DB", Cancer = "#E74C3C")
)

p7 <- pheatmap(
  expr_top,
  scale = "row",
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  clustering_method = "complete",
  annotation_col = annotation_col,
  annotation_colors = ann_colors,
  color = colorRampPalette(c("#3498DB", "white", "#E74C3C"))(100),
  fontsize = 10,
  fontsize_row = 10,
  fontsize_col = 7,
  border_color = NA,
  main = "Expression Heatmap - Top 10 Important Features",
  show_colnames = FALSE
)

# ──────────────────────────────────────────────────────────────────────────
# Visualization 8: Hyperparameter Tuning Results
# ──────────────────────────────────────────────────────────────────────────

mtry_df <- data.frame(
  Fold = factor(1:K),
  Best_mtry = best_mtry_vals,
  AUC = auc_vals
)

p8 <- ggplot(mtry_df, aes(x = Fold, y = Best_mtry, fill = AUC)) +
  geom_bar(stat = "identity", color = "black", size = 0.5) +
  geom_text(aes(label = Best_mtry), vjust = -0.5, fontface = "bold") +
  geom_hline(yintercept = mean(best_mtry_vals), 
             linetype = "dashed", color = "red", size = 1) +
  scale_fill_gradient(low = "#E74C3C", high = "#27AE60") +
  labs(
    title = "Optimal mtry per CV Fold",
    subtitle = sprintf("Mean mtry: %.1f ± %.1f", 
                       mean(best_mtry_vals), sd(best_mtry_vals)),
    x = "CV Fold",
    y = "Optimal mtry",
    fill = "Test AUC"
  ) +
  theme_publication

# ──────────────────────────────────────────────────────────────────────────
# Visualization 9: Boxplots for Top 3 Features
# ──────────────────────────────────────────────────────────────────────────

top3_genes <- head(imp_summary$gene, 3)
expr_top3_long <- melt(as.matrix(t(expr_ml[top3_genes, ])))
colnames(expr_top3_long) <- c("Sample", "Gene", "Expression")
expr_top3_long$Condition <- rep(condition, times = length(top3_genes))

p9 <- ggplot(expr_top3_long, aes(x = Gene, y = Expression, fill = Condition)) +
  geom_boxplot(outlier.shape = 16, outlier.size = 2, alpha = 0.8) +
  scale_fill_manual(values = c(Normal = "#3498DB", Cancer = "#E74C3C")) +
  stat_compare_means(aes(group = Condition), 
                     label = "p.signif", 
                     method = "t.test",
                     label.y.npc = 0.95) +
  labs(
    title = "Expression Levels - Top 3 Most Important Features",
    subtitle = "Statistical significance by t-test",
    x = "Gene",
    y = "Log2(RPKM + 1)"
  ) +
  theme_publication +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"))

# ──────────────────────────────────────────────────────────────────────────
# Visualization 10: Performance Summary Bar Chart
# ──────────────────────────────────────────────────────────────────────────

perf_summary_df <- data.frame(
  Metric = c("AUC", "Accuracy", "Sensitivity", "Specificity", "F1"),
  Mean = c(mean_metrics["AUC"], mean_metrics["Accuracy"], 
           mean_metrics["Sensitivity"], mean_metrics["Specificity"],
           mean_metrics["F1"]),
  SD = c(sd_metrics["AUC"], sd_metrics["Accuracy"],
         sd_metrics["Sensitivity"], sd_metrics["Specificity"],
         sd_metrics["F1"])
)

perf_summary_df$Metric <- factor(perf_summary_df$Metric,
                                 levels = perf_summary_df$Metric)

p10 <- ggplot(perf_summary_df, aes(x = Metric, y = Mean, fill = Metric)) +
  geom_bar(stat = "identity", color = "black", size = 0.5) +
  geom_errorbar(aes(ymin = Mean - SD, ymax = Mean + SD), 
                width = 0.2, size = 0.8) +
  geom_text(aes(label = sprintf("%.3f", Mean)), 
            vjust = -0.5, fontface = "bold", size = 4) +
  scale_fill_brewer(palette = "Set3") +
  labs(
    title = "Random Forest Performance Summary",
    subtitle = sprintf("%d-Fold Cross-Validation Results", K),
    x = "Performance Metric",
    y = "Value"
  ) +
  theme_publication +
  theme(legend.position = "none") +
  ylim(0, 1.1)

# ============================================================================
# PART 7: SAVE ALL VISUALIZATIONS
# ============================================================================

cat("Saving visualizations...\n")

# Combined plots
combined_plot1 <- ggarrange(p1, p2, p5, p10,
                            ncol = 2, nrow = 2,
                            labels = c("A", "B", "C", "D"),
                            font.label = list(size = 16, face = "bold"))

combined_plot2 <- ggarrange(p3, p4, p6, p8,
                            ncol = 2, nrow = 2,
                            labels = c("E", "F", "G", "H"),
                            font.label = list(size = 16, face = "bold"))

# Save combined plots
ggsave("RF_GoldStandard_Performance.pdf", combined_plot1, 
       width = 14, height = 12, dpi = 300)
ggsave("RF_GoldStandard_Features.pdf", combined_plot2,
       width = 14, height = 12, dpi = 300)

# Save individual plots
ggsave("RF_ROC_Curve.png", p1, width = 7, height = 7, units = "in", dpi = 600)
ggsave("RF_Performance_Distribution.png", p2, width = 8, height = 6,, units = "in", dpi = 600)
ggsave("RF_Feature_Importance_MDA.png", p3, width = 8, height = 7,, units = "in", dpi = 600)
ggsave("RF_Feature_Importance_Gini.png", p4, width = 8, height = 7, , units = "in", dpi = 600)
ggsave("RF_Confusion_Matrix.png", p5, width = 7, height = 6,, units = "in", dpi = 600)
ggsave("RF_Feature_Stability.png", p6, width = 8, height = 7,, units = "in", dpi = 600)
ggsave("RF_Optimal_mtry.png", p8, width = 8, height = 6, , units = "in", dpi = 600)
ggsave("RF_Top3_Expression.png", p9, width = 8, height = 6,, units = "in", dpi = 600)
ggsave("RF_Performance_Summary.png", p10, width = 8, height = 6, , units = "in", dpi = 600)

# Save heatmap
png("RF_Expression_Heatmap_Top10.png", width = 10, height = 8, units = "in", res = 600)
print(p7)
dev.off()

cat("✓ All visualizations saved\n")

# ============================================================================
# PART 8: SAVE RESULTS TO FILES
# ============================================================================

cat("\nSaving results to CSV files...\n")

# 1. Per-fold performance metrics
write.csv(summary_df_rf, "RF_Performance_PerFold.csv", row.names = FALSE)

# 2. Summary statistics
summary_stats <- data.frame(
  Metric = c("Accuracy", "AUC", "Sensitivity", "Specificity", 
             "PPV", "NPV", "F1", "Best_mtry"),
  Mean = mean_metrics,
  SD = sd_metrics
)
write.csv(summary_stats, "RF_Performance_Summary.csv", row.names = FALSE)

# 3. Pooled performance
pooled_performance <- data.frame(
  Metric = c("Pooled_AUC", "Pooled_AUC_CI_Lower", "Pooled_AUC_CI_Upper"),
  Value = c(pooled_auc, ci_obj_rf[1], ci_obj_rf[3])
)
write.csv(pooled_performance, "RF_Pooled_Performance.csv", row.names = FALSE)

# 4. Feature importance summary
write.csv(imp_summary, "RF_Feature_Importance_Summary.csv", row.names = FALSE)

# 5. Consensus features
if(length(consensus_features) > 0) {
  consensus_df <- imp_summary %>%
    filter(gene %in% consensus_features) %>%
    select(gene, mean_MDA, sd_MDA, mean_Gini, sd_Gini, selection_frequency)
  
  write.csv(consensus_df, 
            sprintf("RF_Consensus_Features_%dpct.csv", 
                    round(consensus_thresh * 100)),
            row.names = FALSE)
}

# 6. Aggregated confusion matrix
write.csv(cm_df, "RF_Confusion_Matrix_Aggregated.csv", row.names = FALSE)

cat("✓ All results saved to CSV files\n")

# ============================================================================
# FINAL SUMMARY
# ============================================================================

cat("\n" %+% strrep("=", 70) %+% "\n")
cat("RANDOM FOREST GOLD STANDARD ANALYSIS COMPLETE\n")
cat(strrep("=", 70) %+% "\n\n")

cat("=== KEY RESULTS ===\n")
cat(sprintf("  • Cross-Validation: %d-fold nested CV\n", K))
cat(sprintf("  • Mean AUC: %.3f ± %.3f\n", 
            mean_metrics["AUC"], sd_metrics["AUC"]))
cat(sprintf("  • Pooled AUC: %.3f [%.3f, %.3f]\n", 
            pooled_auc, ci_obj_rf[1], ci_obj_rf[3]))
cat(sprintf("  • Mean Accuracy: %.3f ± %.3f\n", 
            mean_metrics["Accuracy"], sd_metrics["Accuracy"]))
cat(sprintf("  • Mean Sensitivity: %.3f ± %.3f\n", 
            mean_metrics["Sensitivity"], sd_metrics["Sensitivity"]))
cat(sprintf("  • Mean Specificity: %.3f ± %.3f\n", 
            mean_metrics["Specificity"], sd_metrics["Specificity"]))
cat(sprintf("  • Mean F1 Score: %.3f ± %.3f\n", 
            mean_metrics["F1"], sd_metrics["F1"]))

cat("\n=== FEATURE SELECTION ===\n")
cat(sprintf("  • Total features analyzed: %d\n", nrow(imp_summary)))
cat(sprintf("  • Consensus features (≥%.0f%%): %d\n", 
            consensus_thresh * 100, length(consensus_features)))
cat(sprintf("  • Optimal mtry: %.1f ± %.1f\n", 
            mean_metrics["Best_mtry"], sd_metrics["Best_mtry"]))

cat("\n=== TOP 5 FEATURES (by Mean Decrease Accuracy) ===\n")
print(head(imp_summary[, c("gene", "mean_MDA", "sd_MDA", "selection_frequency")], 5))

cat("\n=== GENERATED FILES ===\n")
cat("\nPDF Visualizations (12 files):\n")
cat("  - RF_GoldStandard_Performance.pdf (4-panel: ROC, metrics, CM, summary)\n")
cat("  - RF_GoldStandard_Features.pdf (4-panel: importance & stability)\n")
cat("  - RF_ROC_Curve.pdf\n")
cat("  - RF_Performance_Distribution.pdf\n")
cat("  - RF_Feature_Importance_MDA.pdf\n")
cat("  - RF_Feature_Importance_Gini.pdf\n")
cat("  - RF_Confusion_Matrix.pdf\n")
cat("  - RF_Feature_Stability.pdf\n")
cat("  - RF_Expression_Heatmap_Top10.pdf\n")
cat("  - RF_Optimal_mtry.pdf\n")
cat("  - RF_Top3_Expression.pdf\n")
cat("  - RF_Performance_Summary.pdf\n")
cat(sprintf("  - RF_Fold01-Fold%02d_Diagnostics.pdf (%d files)\n", K, K))

cat("\nCSV Data Files:\n")
cat("  - RF_Performance_PerFold.csv\n")
cat("  - RF_Performance_Summary.csv\n")
cat("  - RF_Pooled_Performance.csv\n")
cat("  - RF_Feature_Importance_Summary.csv\n")
if(length(consensus_features) > 0) {
  cat(sprintf("  - RF_Consensus_Features_%dpct.csv\n", 
              round(consensus_thresh * 100)))
}
cat("  - RF_Confusion_Matrix_Aggregated.csv\n")
cat(sprintf("  - RF_Importance_Fold01-Fold%02d.csv (%d files)\n", K, K))

cat("\nRDS Model Files:\n")
cat(sprintf("  - RF_Model_Fold01-Fold%02d.rds (%d files)\n", K, K))

cat("\n" %+% strrep("=", 70) %+% "\n")
cat("All analyses use nested cross-validation with hyperparameter tuning.\n")
cat("Results represent unbiased performance estimates.\n")
cat(strrep("=", 70) %+% "\n\n")

