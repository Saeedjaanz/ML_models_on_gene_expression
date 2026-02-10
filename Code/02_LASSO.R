# ============================================================================
# LASSO LOGISTIC REGRESSION
# ============================================================================

# Required Libraries
library(tidyverse)       # Data manipulation
library(glmnet)          # LASSO regression
library(caret)           # ML workflows
library(pROC)            # ROC analysis
library(UpSetR)          # Set visualizations
library(ggplot2)         # Visualization
library(pheatmap)        # Heatmaps
library(RColorBrewer)    # Color palettes
library(gridExtra)       # Multiple plots
library(ggpubr)          # Publication-ready plots
library(reshape2)        # Data reshaping

# ============================================================================
# PART 1: INPUT VALIDATION & CONFIGURATION
# ============================================================================

cat("\n" %+% strrep("=", 70) %+% "\n")
cat("LASSO LOGISTIC REGRESSION - GOLD STANDARD ANALYSIS\n")
cat(strrep("=", 70) %+% "\n\n")

# Validate input data
stopifnot(is.matrix(expr_ml),
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
K <- 10                       # Outer CV folds
consensus_thresh <- 0.7       # Gene must appear in ≥70% folds
use_lambda <- "lambda.min"    # Options: "lambda.min" or "lambda.1se"
perform_permutation <- FALSE  # Set TRUE for permutation test (computationally expensive)
n_perm <- 200                 # Number of permutations (1000 recommended)

cat(sprintf("\n✓ Configuration:\n"))
cat(sprintf("  - Outer CV folds: %d\n", K))
cat(sprintf("  - Lambda selection: %s\n", use_lambda))
cat(sprintf("  - Consensus threshold: %.0f%%\n", consensus_thresh * 100))
cat(sprintf("  - Permutation test: %s\n", ifelse(perform_permutation, "Yes", "No")))

# ============================================================================
# PART 2: DATA PREPARATION
# ============================================================================

# Binary labels (0 = Normal, 1 = Cancer)
y_bin <- ifelse(condition == "Cancer", 1, 0)

# Create stratified outer folds
folds <- createFolds(condition, k = K, list = TRUE, returnTrain = FALSE)

# Storage objects for comprehensive metrics
mse_vals <- numeric(K)
acc_vals <- numeric(K)
auc_vals <- numeric(K)
sens_vals <- numeric(K)
spec_vals <- numeric(K)
ppv_vals <- numeric(K)
npv_vals <- numeric(K)
f1_vals <- numeric(K)
lambda_vals <- numeric(K)
n_features_vals <- numeric(K)

conf_mats <- vector("list", K)
selected_genes_list <- vector("list", K)
coef_matrix_list <- vector("list", K)

# For pooled ROC
pooled_obs <- c()
pooled_pred <- c()
pooled_pred_class <- c()

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
  
  # Prepare matrices: glmnet expects [samples × features]
  x_train <- t(expr_ml[, train_idx, drop = FALSE])
  x_test  <- t(expr_ml[, test_idx, drop = FALSE])
  
  y_train <- y_bin[train_idx]
  y_test  <- y_bin[test_idx]
  
  cat(sprintf("│ Train: %d samples | Test: %d samples\n", 
              length(train_idx), length(test_idx)))
  
  # ──────────────────────────────────────────────────────────────────────────
  # INNER CV: Lambda Selection (Hyperparameter Tuning)
  # ──────────────────────────────────────────────────────────────────────────
  
  cat("│\n│ ➤ Inner CV: Tuning lambda via AUC optimization...\n")
  
  # Determine inner folds
  inner_nfolds <- min(10, max(3, length(train_idx)))
  
  # Fit cv.glmnet (LASSO logistic regression)
  cvfit <- cv.glmnet(
    x = x_train, y = y_train,
    family = "binomial",
    alpha = 1,              # LASSO (L1 penalty)
    type.measure = "auc",   # Optimize AUC
    nfolds = inner_nfolds, standardize = TRUE)
  
  # Select lambda
  chosen_lambda <- if (use_lambda == "lambda.1se") {
    cvfit$lambda.1se
  } else {
    cvfit$lambda.min
  }
  lambda_vals[i] <- chosen_lambda
  
  cat(sprintf("│   ✓ Chosen lambda (%s): %.6f\n", use_lambda, chosen_lambda))
  cat(sprintf("│   ✓ Inner CV AUC: %.4f\n", max(cvfit$cvm)))
  
  # ──────────────────────────────────────────────────────────────────────────
  # PREDICT ON HELD-OUT TEST SET
  # ──────────────────────────────────────────────────────────────────────────
  
  cat("│\n│ ➤ Evaluating on held-out test set...\n")
  
  # Predictions
  preds_prob <- as.numeric(predict(cvfit, newx = x_test, s = chosen_lambda, 
                                   type = "response"))
  preds_class <- factor(ifelse(preds_prob >= 0.5, "Cancer", "Normal"), 
                        levels = levels(condition))
  
  # Pool predictions
  pooled_obs <- c(pooled_obs, y_test)
  pooled_pred <- c(pooled_pred, preds_prob)
  pooled_pred_class <- c(pooled_pred_class, as.character(preds_class))
  
  # Calculate metrics
  mse_vals[i] <- mean((preds_prob - y_test)^2)
  acc_vals[i] <- mean(preds_class == condition[test_idx])
  
  # ROC and AUC
  roc_obj_i <- try(roc(y_test, preds_prob, quiet = TRUE), silent = TRUE)
  auc_vals[i] <- if (inherits(roc_obj_i, "try-error")) NA else as.numeric(roc_obj_i$auc)
  
  # Confusion matrix
  cm <- confusionMatrix(preds_class, condition[test_idx])
  conf_mats[[i]] <- cm$table
  
  # Extract additional metrics
  sens_vals[i] <- cm$byClass["Sensitivity"]
  spec_vals[i] <- cm$byClass["Specificity"]
  ppv_vals[i] <- cm$byClass["Pos Pred Value"]
  npv_vals[i] <- cm$byClass["Neg Pred Value"]
  f1_vals[i] <- cm$byClass["F1"]
  
  cat(sprintf("│   ✓ AUC: %.4f | Accuracy: %.4f\n", auc_vals[i], acc_vals[i]))
  cat(sprintf("│   ✓ Sensitivity: %.4f | Specificity: %.4f\n", 
              sens_vals[i], spec_vals[i]))
  
  # ──────────────────────────────────────────────────────────────────────────
  # FEATURE SELECTION: Extract non-zero coefficients
  # ──────────────────────────────────────────────────────────────────────────
  
  coef_sparse <- as.matrix(coef(cvfit, s = chosen_lambda))
  coef_df <- data.frame(
    gene = rownames(coef_sparse),
    coef = as.numeric(coef_sparse),
    stringsAsFactors = FALSE
  )
  coef_df <- coef_df[coef_df$gene != "(Intercept)", ]
  nz <- coef_df[coef_df$coef != 0, ]
  
  selected_genes_list[[i]] <- nz$gene
  coef_matrix_list[[i]] <- nz
  n_features_vals[i] <- nrow(nz)
  
  cat(sprintf("│\n│ ➤ Selected features: %d genes\n", nrow(nz)))
  if(nrow(nz) > 0) {
    cat(sprintf("│   %s\n", paste(head(nz$gene, 5), collapse = ", ")))
    if(nrow(nz) > 5) cat("│   ...\n")
  }
  
  cat("└" %+% strrep("─", 68) %+% "┘\n\n")
}

# ============================================================================
# PART 4: AGGREGATE PERFORMANCE METRICS
# ============================================================================

cat(strrep("=", 70) %+% "\n")
cat("AGGREGATED PERFORMANCE METRICS\n")
cat(strrep("=", 70) %+% "\n\n")

# Summary dataframe
summary_df <- data.frame(
  Fold = seq_len(K), MSE = mse_vals,
  Accuracy = acc_vals, AUC = auc_vals,
  Sensitivity = sens_vals,
  Specificity = spec_vals,
  PPV = ppv_vals, NPV = npv_vals,
  F1 = f1_vals, Lambda = lambda_vals,
  N_Features = n_features_vals)

cat("Per-Fold Performance:\n")
print(summary_df)

# Calculate means and SDs
mean_metrics <- colMeans(summary_df[, -1], na.rm = TRUE)
sd_metrics <- apply(summary_df[, -1], 2, sd, na.rm = TRUE)

cat("\n" %+% strrep("-", 70) %+% "\n")
cat("SUMMARY STATISTICS (Mean ± SD):\n")
cat(strrep("-", 70) %+% "\n")
cat(sprintf("  MSE:         %.4f ± %.4f\n", mean_metrics["MSE"], sd_metrics["MSE"]))
cat(sprintf("  Accuracy:    %.4f ± %.4f\n", mean_metrics["Accuracy"], sd_metrics["Accuracy"]))
cat(sprintf("  AUC:         %.4f ± %.4f\n", mean_metrics["AUC"], sd_metrics["AUC"]))
cat(sprintf("  Sensitivity: %.4f ± %.4f\n", mean_metrics["Sensitivity"], sd_metrics["Sensitivity"]))
cat(sprintf("  Specificity: %.4f ± %.4f\n", mean_metrics["Specificity"], sd_metrics["Specificity"]))
cat(sprintf("  PPV:         %.4f ± %.4f\n", mean_metrics["PPV"], sd_metrics["PPV"]))
cat(sprintf("  NPV:         %.4f ± %.4f\n", mean_metrics["NPV"], sd_metrics["NPV"]))
cat(sprintf("  F1 Score:    %.4f ± %.4f\n", mean_metrics["F1"], sd_metrics["F1"]))
cat(sprintf("  Lambda:      %.6f ± %.6f\n", mean_metrics["Lambda"], sd_metrics["Lambda"]))
cat(sprintf("  N Features:  %.1f ± %.1f\n", mean_metrics["N_Features"], sd_metrics["N_Features"]))

# Pooled ROC analysis
if (length(pooled_obs) > 0) {
  pooled_roc <- roc(pooled_obs, pooled_pred, quiet = TRUE)
  pooled_auc <- as.numeric(pooled_roc$auc)
  
  cat("\n" %+% strrep("-", 70) %+% "\n")
  cat("POOLED PERFORMANCE (All Folds Combined):\n")
  cat(strrep("-", 70) %+% "\n")
  cat(sprintf("  Pooled AUC: %.4f\n", pooled_auc))
  
  # Bootstrap CI
  ci_obj_lasso <- ci.auc(pooled_roc, method = "bootstrap", boot.n = 2000)
  cat(sprintf("  95%% CI (Bootstrap, n=2000): [%.4f, %.4f]\n", 
              ci_obj_lasso[1], ci_obj_lasso[3]))
} else {
  cat("\nNo pooled predictions available.\n")
}

# ============================================================================
# PART 5: FEATURE STABILITY ANALYSIS
# ============================================================================

cat("\n" %+% strrep("=", 70) %+% "\n")
cat("FEATURE SELECTION STABILITY\n")
cat(strrep("=", 70) %+% "\n\n")

all_sel <- unlist(selected_genes_list)

if (length(all_sel) == 0) {
  cat("⚠ No genes selected in any outer fold.\n")
  cons <- data.frame(gene = character(0), Freq = integer(0), 
                     stringsAsFactors = FALSE)
  consensus_genes <- character(0)
} else {
  # Feature frequency table
  cons <- as.data.frame(table(all_sel), stringsAsFactors = FALSE)
  colnames(cons) <- c("gene", "Freq")
  cons <- cons %>% arrange(desc(Freq))
  cons$FreqPercent <- cons$Freq / K * 100
  
  # Calculate mean coefficient across folds
  cons$mean_coef <- sapply(cons$gene, function(g) {
    coefs <- sapply(coef_matrix_list, function(cm) {
      if(g %in% cm$gene) cm$coef[cm$gene == g] else 0
    })
    mean(coefs[coefs != 0])
  })
  
  cons$sd_coef <- sapply(cons$gene, function(g) {
    coefs <- sapply(coef_matrix_list, function(cm) {
      if(g %in% cm$gene) cm$coef[cm$gene == g] else 0
    })
    sd(coefs[coefs != 0])
  })
  
  cat("Feature Selection Frequency (Top 20):\n")
  print(head(cons, 20))
  
  # Consensus features
  consensus_count <- ceiling(K * consensus_thresh)
  consensus_genes <- cons$gene[cons$Freq >= consensus_count]
  
  cat(sprintf("\n✓ Consensus features (≥%.0f%% selection): %d genes\n", 
              consensus_thresh * 100, length(consensus_genes)))
  if(length(consensus_genes) > 0) {
    cat("  ", paste(consensus_genes, collapse = ", "), "\n")
  }
}

# ============================================================================
# PART 6: FINAL MODEL ON FULL DATA
# ============================================================================

cat("\n" %+% strrep("=", 70) %+% "\n")
cat("FINAL MODEL FIT ON FULL DATASET\n")
cat(strrep("=", 70) %+% "\n\n")

# Prepare full data
full_x <- t(expr_ml)
full_y <- y_bin

# Inner CV folds for full data
nfolds_full <- min(10, max(3, ncol(expr_ml)))

# Fit cv.glmnet on full dataset
cv_full <- cv.glmnet(
  x = full_x, y = full_y,
  family = "binomial",
  alpha = 1, type.measure = "auc",
  nfolds = nfolds_full, standardize = TRUE)

cat(sprintf("Full-data lambda.min: %.6f\n", cv_full$lambda.min))
cat(sprintf("Full-data lambda.1se: %.6f\n", cv_full$lambda.1se))

# Extract coefficients at both lambdas
final_lambda_min <- cv_full$lambda.min
final_lambda_1se <- cv_full$lambda.1se

coefs_min <- as.matrix(coef(cv_full, s = final_lambda_min))
coefs_1se <- as.matrix(coef(cv_full, s = final_lambda_1se))

# Process lambda.min coefficients
coefs_min_df <- data.frame(
  gene = rownames(coefs_min),
  coef = as.numeric(coefs_min),
  stringsAsFactors = FALSE
)
coefs_min_df <- coefs_min_df[coefs_min_df$gene != "(Intercept)", ]
nonzero_min <- coefs_min_df[coefs_min_df$coef != 0, ]

# Process lambda.1se coefficients
coefs_1se_df <- data.frame(
  gene = rownames(coefs_1se),
  coef = as.numeric(coefs_1se),
  stringsAsFactors = FALSE
)
coefs_1se_df <- coefs_1se_df[coefs_1se_df$gene != "(Intercept)", ]
nonzero_1se <- coefs_1se_df[coefs_1se_df$coef != 0, ]

cat(sprintf("\nNon-zero coefficients at lambda.min: %d genes\n", nrow(nonzero_min)))
cat(sprintf("Non-zero coefficients at lambda.1se: %d genes\n", nrow(nonzero_1se)))

# Also fit full path for visualization
fit_glmnet_full <- glmnet(x = full_x, y = full_y, family = "binomial", 
                          alpha = 1, standardize = TRUE)

# ============================================================================
# PART 7: VISUALIZATIONS
# ============================================================================

cat("\n" %+% strrep("=", 70) %+% "\n")
cat("GENERATING FIGURES\n")
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
# Visualization 1: Pooled ROC Curve
# ──────────────────────────────────────────────────────────────────────────

roc_df <- data.frame(
  FPR = 1 - pooled_roc$specificities,
  TPR = pooled_roc$sensitivities
)

p1 <- ggplot(roc_df, aes(x = FPR, y = TPR)) +
  geom_line(color = "#3498DB", size = 1.3) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "grey50", size = 0.8) +
  annotate("text", x = 0.6, y = 0.3, label = sprintf("Pooled AUC = %.3f", pooled_auc),
           size = 5, fontface = "bold", color = "#2C3E50") +
  annotate("text", x = 0.6, y = 0.2, label = sprintf("95%% CI: [%.3f, %.3f]", ci_obj_lasso[1], 
           ci_obj_lasso[3]), size = 4, color = "#7F8C8D") +
  annotate("text", x = 0.6, y = 0.1, label = sprintf("Mean CV AUC: %.3f ± %.3f", 
           mean_metrics["AUC"], sd_metrics["AUC"]), size = 4, color = "#7F8C8D") +
  labs(
    title = "LASSO ROC Curve",
    subtitle = sprintf("%d-Fold Cross-Validation (Pooled Predictions)", K),
    x = "False Positive Rate (1 - Specificity)",
    y = "True Positive Rate (Sensitivity)") +
  coord_equal() +
  theme_publication

# ──────────────────────────────────────────────────────────────────────────
# Visualization 2: CV Plot (Lambda Selection)
# ──────────────────────────────────────────────────────────────────────────

# Extract CV curve data
cv_data <- data.frame(
  lambda = cv_full$lambda,
  cvm = cv_full$cvm, cvsd = cv_full$cvsd,
  cvup = cv_full$cvup, cvlo = cv_full$cvlo,
  nzero = cv_full$nzero)

p2 <- ggplot(cv_data, aes(x = log(lambda), y = cvm)) +
  geom_ribbon(aes(ymin = cvlo, ymax = cvup), alpha = 0.2, fill = "#3498DB") +
  geom_line(color = "#2C3E50", size = 1.2) +
  geom_point(color = "#E74C3C", size = 2) +
  geom_vline(xintercept = log(cv_full$lambda.min), linetype = "dashed", 
             color = "#27AE60", size = 1) +
  geom_vline(xintercept = log(cv_full$lambda.1se), linetype = "dashed", 
             color = "#F39C12", size = 1) +
  annotate("text", x = log(cv_full$lambda.min), y = min(cv_data$cvm), 
           label = "lambda.min", hjust = -0.1, color = "#27AE60", fontface = "bold") +
  annotate("text", x = log(cv_full$lambda.1se), y = min(cv_data$cvm), 
           label = "lambda.1se", hjust = 1.1,vjust = 3.5, color = "#F39C12", fontface = "bold") +
  labs(title = "LASSO Cross-Validation Curve",
    subtitle = sprintf("%d-Fold CV on Full Dataset", nfolds_full),
    x = "log(Lambda)", y = "Cross-Validated AUC") +
  theme_publication

# ──────────────────────────────────────────────────────────────────────────
# Visualization 3: Coefficient Paths
# ──────────────────────────────────────────────────────────────────────────

# Extract coefficient paths
coef_path <- as.matrix(fit_glmnet_full$beta)
lambda_seq <- fit_glmnet_full$lambda

# Select top genes by maximum absolute coefficient
max_coef <- apply(abs(coef_path), 1, max)
top_genes_path <- names(sort(max_coef, decreasing = TRUE)[1:min(20, nrow(coef_path))])

# Create long format
path_df <- data.frame(t(coef_path[top_genes_path, ]))
path_df$lambda <- lambda_seq
path_df$l1_norm <- fit_glmnet_full$dev.ratio

path_long <- path_df %>%
  pivot_longer(cols = -c(lambda, l1_norm), names_to = "gene", values_to = "coefficient")

p3 <- ggplot(path_long, aes(x = log(lambda), y = coefficient, color = gene)) +
  geom_line(size = 0.8, alpha = 0.7) +
  geom_vline(xintercept = log(cv_full$lambda.min), linetype = "dashed", 
             color = "#27AE60", size = 1) +
  geom_vline(xintercept = log(cv_full$lambda.1se), linetype = "dashed", 
             color = "#F39C12", size = 1) +
  labs(title = "LASSO Coefficient Regularization Paths",
    subtitle = "Top 20 Genes by Maximum |Coefficient|",
    x = "log(Lambda)", y = "Coefficient Value", color = "Gene" ) +
  theme_publication +
  theme(legend.position = "right", legend.key.size = unit(0.4, "cm"))

# ──────────────────────────────────────────────────────────────────────────
# Visualization 4: Performance Metrics Distribution
# ──────────────────────────────────────────────────────────────────────────

perf_long <- summary_df %>%
  select(Fold, Accuracy, AUC, Sensitivity, Specificity, F1) %>%
  pivot_longer(cols = -Fold, names_to = "Metric", values_to = "Value")

p4 <- ggplot(perf_long, aes(x = Metric, y = Value, fill = Metric)) +
  geom_boxplot(alpha = 0.7, outlier.shape = 16, outlier.size = 2) +
  geom_jitter(width = 0.15, size = 2, alpha = 0.5) +
  stat_summary(fun = mean, geom = "point", shape = 23, size = 4,
               fill = "red", color = "darkred") +
  scale_fill_brewer(palette = "Set2") +
  labs(title = "Performance Metrics Across CV Folds",
    subtitle = "Diamond = Mean, Box = IQR, Points = Individual Folds",
    x = "Performance Metric", y = "Value") +
  theme_publication +
  theme(legend.position = "none") +
  ylim(0, 1)

# ──────────────────────────────────────────────────────────────────────────
# Visualization 5: Feature Selection Stability
# ──────────────────────────────────────────────────────────────────────────

if(nrow(cons) > 0) {
  top_features_stability <- head(cons, 15)
  top_features_stability$gene <- factor(top_features_stability$gene,
                                        levels = rev(top_features_stability$gene))
  
  p5 <- ggplot(top_features_stability, 
               aes(x = gene, y = FreqPercent, fill = FreqPercent)) +
    geom_bar(stat = "identity", color = "black", size = 0.3) +
    geom_hline(yintercept = consensus_thresh * 100, 
               linetype = "dashed", color = "red", size = 1) +
    annotate("text", x = 13, y = consensus_thresh * 100 + 5, 
             label = sprintf("%.0f%% threshold", consensus_thresh * 100),
             color = "red", fontface = "bold", hjust = 1) +
    scale_fill_gradient2(low = "#E74C3C", mid = "#F39C12", high = "#27AE60",
                         midpoint = 50) +
    coord_flip(clip="off") +
    labs(title = "Feature Selection Stability",
      subtitle = sprintf("Top 15 Features (from %d-Fold CV)", K),
      x = "Gene", y = "Selection Frequency (%)", fill = "Frequency (%)") +
    theme_publication +
    theme(legend.position = "right") +
    ylim(0, 105)
} else {
  p5 <- ggplot() + theme_void() + 
    annotate("text", x = 0.5, y = 0.5, label = "No features selected", 
             size = 6, fontface = "bold")
}

# ──────────────────────────────────────────────────────────────────────────
# Visualization 6: Confusion Matrix (Aggregated)
# ──────────────────────────────────────────────────────────────────────────

cm_aggregated <- Reduce("+", conf_mats)
cm_df <- melt(cm_aggregated)
colnames(cm_df) <- c("Prediction", "Reference", "Count")

cm_df <- cm_df %>%
  group_by(Reference) %>%
  mutate(Percentage = Count / sum(Count) * 100) %>%
  ungroup()

p6 <- ggplot(cm_df, aes(x = Reference, y = Prediction, fill = Count)) +
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
# Visualization 7: Coefficient Heatmap (Consensus Features)
# ──────────────────────────────────────────────────────────────────────────

if(length(consensus_genes) > 0) {
  # Extract coefficients for consensus features across folds
  coef_matrix <- matrix(0, nrow = length(consensus_genes), ncol = K)
  rownames(coef_matrix) <- consensus_genes
  colnames(coef_matrix) <- paste0("Fold", 1:K)
  
  for(i in 1:K) {
    fold_coefs <- coef_matrix_list[[i]]
    for(g in consensus_genes) {
      if(g %in% fold_coefs$gene) {
        coef_matrix[g, i] <- fold_coefs$coef[fold_coefs$gene == g]
      }
    }
  }
  
  p7 <- pheatmap(
    coef_matrix,
    scale = "none",
    cluster_rows = TRUE,
    cluster_cols = FALSE,
    color = colorRampPalette(c("#3498DB", "white", "#E74C3C"))(100),
    breaks = seq(-max(abs(coef_matrix)), max(abs(coef_matrix)), length.out = 101),
    fontsize = 10,
    fontsize_row = 9,
    fontsize_col = 9,
    border_color = "grey60",
    main = "LASSO Coefficients - Consensus Features Across Folds",
    display_numbers = FALSE
  )
} else {
  p7 <- NULL
}

# ──────────────────────────────────────────────────────────────────────────
# Visualization 8: Number of Features Selected per Fold
# ──────────────────────────────────────────────────────────────────────────

features_df <- data.frame(
  Fold = factor(1:K),
  N_Features = n_features_vals,
  AUC = auc_vals
)

p8 <- ggplot(features_df, aes(x = Fold, y = N_Features, fill = AUC)) +
  geom_bar(stat = "identity", color = "black", size = 0.5) +
  geom_text(aes(label = N_Features), vjust = -0.5, fontface = "bold") +
  geom_hline(yintercept = mean(n_features_vals), 
             linetype = "dashed", color = "red", size = 1) +
  scale_fill_gradient(low = "#E74C3C", high = "#27AE60") +
  labs(
    title = "Features Selected per CV Fold",
    subtitle = sprintf("Mean: %.1f ± %.1f features", 
                       mean_metrics["N_Features"], sd_metrics["N_Features"]),
    x = "CV Fold",
    y = "Number of Selected Features",
    fill = "Test AUC"
  ) +
  theme_publication

# ──────────────────────────────────────────────────────────────────────────
# Visualization 9: Expression Heatmap of Consensus Features
# ──────────────────────────────────────────────────────────────────────────

if(length(consensus_genes) > 0) {
  expr_consensus <- expr_ml[consensus_genes, , drop = FALSE]
  
  # Annotation
  annotation_col <- data.frame(Condition = condition)
  rownames(annotation_col) <- colnames(expr_consensus)
  
  ann_colors <- list(
    Condition = c(Normal = "#3498DB", Cancer = "#E74C3C")
  )
  
  p9 <- pheatmap(
    expr_consensus, scale = "row",
    clustering_distance_rows = "euclidean",
    clustering_distance_cols = "euclidean",
    clustering_method = "complete",
    annotation_col = annotation_col,
    annotation_colors = ann_colors,
    color = colorRampPalette(c("#3498DB", "white", "#E74C3C"))(100),
    fontsize = 10, fontsize_row = 9, fontsize_col = 7,
    border_color = NA,
    main = "Expression Heatmap - Consensus Features",
    show_colnames = FALSE)
} else {
  p9 <- NULL
}

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
  geom_errorbar(aes(ymin = Mean - SD, ymax = Mean + SD), width = 0.2, size = 0.8) +
  geom_text(aes(label = sprintf("%.3f", Mean)),  vjust = -0.5, fontface = "bold", size = 4) +
  scale_fill_brewer(palette = "Set3") +
  labs(title = "LASSO Performance Summary",
    subtitle = sprintf("%d-Fold Cross-Validation Results", K),
    x = "Performance Metric", y = "Value") +
  theme_publication +
  theme(legend.position = "none") + ylim(0, 1.1)

# ──────────────────────────────────────────────────────────────────────────
# Visualization 11: Boxplots for Top 3 Consensus Features
# ──────────────────────────────────────────────────────────────────────────

if(length(consensus_genes) >= 3) {
  top3_genes <- consensus_genes[1:3]
  expr_top3_long <- melt(as.matrix(t(expr_ml[top3_genes, ])))
  colnames(expr_top3_long) <- c("Sample", "Gene", "Expression")
  expr_top3_long$Condition <- rep(condition, times = length(top3_genes))
  
  p11 <- ggplot(expr_top3_long, aes(x = Gene, y = Expression, fill = Condition)) +
    geom_boxplot(outlier.shape = 16, outlier.size = 2, alpha = 0.8) +
    scale_fill_manual(values = c(Normal = "#3498DB", Cancer = "#E74C3C")) +
    stat_compare_means(aes(group = Condition), 
                       label = "p.signif", 
                       method = "t.test",
                       label.y.npc = 0.95) +
    labs(
      title = "Expression Levels - Top 3 Consensus Features",
      subtitle = "Statistical significance by t-test",
      x = "Gene", y = "Log2(RPKM + 1)") +
    theme_publication +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"))
} else {
  p11 <- NULL
}

# ──────────────────────────────────────────────────────────────────────────
# Visualization 12: Lambda Selection Across Folds
# ──────────────────────────────────────────────────────────────────────────

lambda_df <- data.frame(Fold = factor(1:K), Lambda = lambda_vals, AUC = auc_vals)

p12 <- ggplot(lambda_df, aes(x = Fold, y = log(Lambda), fill = AUC)) +
  geom_bar(stat = "identity", color = "black", size = 0.5) +
  geom_hline(yintercept = log(mean(lambda_vals)), 
             linetype = "dashed", color = "red", size = 1) +
  scale_fill_gradient(low = "#E74C3C", high = "#27AE60") +
  labs(
    title = "Selected Lambda per CV Fold",
    subtitle = sprintf("Mean log(lambda): %.3f ± %.3f", 
                       log(mean(lambda_vals)), log(sd(lambda_vals))),
    x = "CV Fold", y = "log(Lambda)", fill = "Test AUC") +
  theme_publication

# ============================================================================
# PART 8: SAVE ALL VISUALIZATIONS
# ============================================================================

cat("Saving visualizations...\n")

# Combined plots
combined_plot1 <- ggarrange(p1, p2, p4, p6,
                            ncol = 2, nrow = 2,
                            labels = c("A", "B", "C", "D"),
                            font.label = list(size = 16, face = "bold"))

combined_plot2 <- ggarrange(p3, p5, p8, p10,
                            ncol = 2, nrow = 2,
                            labels = c("E", "F", "G", "H"),
                            font.label = list(size = 16, face = "bold"))

# Save combined plots
ggsave("LASSO_Performance.pdf", combined_plot1, width = 14, height = 12, dpi = 300)
ggsave("LASSO_Features.pdf", combined_plot2, width = 14, height = 12, dpi = 300)

# Save individual plots
ggsave("LASSO_ROC_Curve.png", p1, width = 7, height = 7, units = "in", dpi = 600)
ggsave("LASSO_CV_Lambda_Selection.png", p2, width = 8, height = 6, units = "in", dpi = 600)
ggsave("LASSO_Coefficient_Paths.png", p3, width = 10, height = 6, units = "in", dpi = 600)
ggsave("LASSO_Performance_Distribution.png", p4, width = 8, height = 6, units = "in", dpi = 600)
ggsave("LASSO_Feature_Stability.png", p5, width = 8, height = 7, units = "in", dpi = 600)
ggsave("LASSO_Confusion_Matrix.png", p6, width = 7, height = 6, units = "in", dpi = 600)
ggsave("LASSO_Features_Per_Fold.png", p8, width = 8, height = 6, units = "in", dpi = 600)
ggsave("LASSO_Performance_Summary.png", p10, width = 8, height = 6, units = "in", dpi = 600)
if(!is.null(p11)) ggsave("LASSO_Top3_Expression.png", p11, width = 8, height = 6, units = "in", dpi = 600)
ggsave("LASSO_Lambda_Per_Fold.png", p12, width = 8, height = 6, units = "in", dpi = 600)

# Save heatmaps
if(!is.null(p7)) {
  png("LASSO_Coefficient_Heatmap_Consensus.png", width = 10, height = 8, units = "in", res = 600)
  print(p7)
  dev.off()
}

if(!is.null(p9)) {
  png("LASSO_Expression_Heatmap_Consensus.png", width = 10, height = 8, units = "in", res = 600)
  print(p9)
  dev.off()
}
cat("✓ All visualizations saved\n")

# ============================================================================
# PART 9: SAVE RESULTS TO FILES
# ============================================================================

cat("\nSaving results to CSV files...\n")

# 1. Per-fold performance metrics
write.csv(summary_df, "LASSO_Performance_PerFold.csv", row.names = FALSE)

# 2. Summary statistics
summary_stats <- data.frame(Metric = names(mean_metrics), Mean = mean_metrics, 
                            SD = sd_metrics)
write.csv(summary_stats, "LASSO_Performance_Summary.csv", row.names = FALSE)

# 3. Pooled performance
pooled_performance <- data.frame(
  Metric = c("Pooled_AUC", "Pooled_AUC_CI_Lower", "Pooled_AUC_CI_Upper"),
  Value = c(pooled_auc, ci_obj_lasso[1], ci_obj_lasso[3]))

write.csv(pooled_performance, "LASSO_Pooled_Performance.csv", row.names = FALSE)

# 4. Feature selection frequency/stability
if(nrow(cons) > 0) {
  write.csv(cons, "LASSO_Feature_Selection_Stability.csv", row.names = FALSE)
}

# 5. Consensus features
if(length(consensus_genes) > 0) {
  consensus_df <- cons %>%
    filter(gene %in% consensus_genes) %>%
    select(gene, Freq, FreqPercent, mean_coef, sd_coef)
  
  write.csv(consensus_df, sprintf("LASSO_Consensus_Features_%dpct.csv", 
                    round(consensus_thresh * 100)), row.names = FALSE)
}

# 6. Final model coefficients (lambda.min)
write.csv(nonzero_min, "LASSO_Final_Coefficients_LambdaMin.csv", row.names = FALSE)

# 7. Final model coefficients (lambda.1se)
write.csv(nonzero_1se, "LASSO_Final_Coefficients_Lambda1se.csv", row.names = FALSE)

# 8. Aggregated confusion matrix
write.csv(cm_df, "LASSO_Confusion_Matrix_Aggregated.csv", row.names = FALSE)

# 9. Per-fold selected features
for(i in 1:K) {
  if(nrow(coef_matrix_list[[i]]) > 0) {
    write.csv(coef_matrix_list[[i]], 
              sprintf("LASSO_Selected_Features_Fold%02d.csv", i),
              row.names = FALSE)
  }
}

cat("✓ All results saved to CSV files\n")

# ============================================================================
# FINAL SUMMARY
# ============================================================================

cat("\n" %+% strrep("=", 70) %+% "\n")
cat("LASSO LOGISTIC REGRESSION - GOLD STANDARD ANALYSIS COMPLETE\n")
cat(strrep("=", 70) %+% "\n\n")

cat("=== KEY RESULTS ===\n")
cat(sprintf("  • Cross-Validation: %d-fold nested CV\n", K))
cat(sprintf("  • Lambda selection: %s\n", use_lambda))
cat(sprintf("  • Mean AUC: %.3f ± %.3f\n", mean_metrics["AUC"], sd_metrics["AUC"]))
cat(sprintf("  • Pooled AUC: %.3f [%.3f, %.3f]\n", pooled_auc, ci_obj_lasso[1], ci_obj_lasso[3]))
cat(sprintf("  • Mean Accuracy: %.3f ± %.3f\n", mean_metrics["Accuracy"], sd_metrics["Accuracy"]))
cat(sprintf("  • Mean Sensitivity: %.3f ± %.3f\n", mean_metrics["Sensitivity"], sd_metrics["Sensitivity"]))
cat(sprintf("  • Mean Specificity: %.3f ± %.3f\n", mean_metrics["Specificity"], sd_metrics["Specificity"]))
cat(sprintf("  • Mean F1 Score: %.3f ± %.3f\n", mean_metrics["F1"], sd_metrics["F1"]))

cat("\n=== FEATURE SELECTION ===\n")
cat(sprintf("  • Mean features selected: %.1f ± %.1f per fold\n", 
            mean_metrics["N_Features"], sd_metrics["N_Features"]))
cat(sprintf("  • Consensus features (≥%.0f%%): %d genes\n", 
            consensus_thresh * 100, length(consensus_genes)))
cat(sprintf("  • Lambda.min features (full data): %d\n", nrow(nonzero_min)))
cat(sprintf("  • Lambda.1se features (full data): %d\n", nrow(nonzero_1se)))

cat("\n=== REGULARIZATION PARAMETERS ===\n")
cat(sprintf("  • Mean selected lambda: %.6f ± %.6f\n", 
            mean_metrics["Lambda"], sd_metrics["Lambda"]))
cat(sprintf("  • Full-data lambda.min: %.6f\n", cv_full$lambda.min))
cat(sprintf("  • Full-data lambda.1se: %.6f\n", cv_full$lambda.1se))

if(length(consensus_genes) > 0) {
  cat("\n=== CONSENSUS FEATURES ===\n")
  consensus_table <- cons %>%
    filter(gene %in% consensus_genes) %>%
    select(gene, Freq, FreqPercent, mean_coef)
  print(consensus_table)
}

cat("\n=== GENERATED FILES ===\n")
cat("\nPDF Visualizations:\n")
cat("  - LASSO_GoldStandard_Performance.pdf (4-panel: ROC, CV, metrics, CM)\n")
cat("  - LASSO_GoldStandard_Features.pdf (4-panel: paths, stability, features)\n")
cat("  - LASSO_ROC_Curve.pdf\n")
cat("  - LASSO_CV_Lambda_Selection.pdf\n")
cat("  - LASSO_Coefficient_Paths.pdf\n")
cat("  - LASSO_Performance_Distribution.pdf\n")
cat("  - LASSO_Feature_Stability.pdf\n")
cat("  - LASSO_Confusion_Matrix.pdf\n")
if(!is.null(p7)) cat("  - LASSO_Coefficient_Heatmap_Consensus.pdf\n")
cat("  - LASSO_Features_Per_Fold.pdf\n")
if(!is.null(p9)) cat("  - LASSO_Expression_Heatmap_Consensus.pdf\n")
cat("  - LASSO_Performance_Summary.pdf\n")
if(!is.null(p11)) cat("  - LASSO_Top3_Expression.pdf\n")
cat("  - LASSO_Lambda_Per_Fold.pdf\n")

cat("\nCSV Data Files:\n")
cat("  - LASSO_Performance_PerFold.csv\n")
cat("  - LASSO_Performance_Summary.csv\n")
cat("  - LASSO_Pooled_Performance.csv\n")
if(nrow(cons) > 0) cat("  - LASSO_Feature_Selection_Stability.csv\n")
if(length(consensus_genes) > 0) {
  cat(sprintf("  - LASSO_Consensus_Features_%dpct.csv\n", 
              round(consensus_thresh * 100)))
}
cat("  - LASSO_Final_Coefficients_LambdaMin.csv\n")
cat("  - LASSO_Final_Coefficients_Lambda1se.csv\n")
cat("  - LASSO_Confusion_Matrix_Aggregated.csv\n")
cat(sprintf("  - LASSO_Selected_Features_Fold01-Fold%02d.csv\n", K))

cat("\n" %+% strrep("=", 70) %+% "\n")
cat("All analyses use nested cross-validation with hyperparameter tuning.\n")
cat("Feature selection via L1 regularization (LASSO penalty).\n")
cat("Results represent unbiased performance estimates.\n")
cat(strrep("=", 70) %+% "\n\n")
