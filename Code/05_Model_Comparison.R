# ============================================================================
# VENN DIAGRAM & UPSET PLOT - MODEL COMPARISON
# Compare Selected Features Across LASSO, Random Forest, and SVM-RFE
# Publication-Quality Visualization
# ============================================================================

# Required Libraries
library(ggplot2)
library(VennDiagram)
library(UpSetR)
library(ggpubr)
library(RColorBrewer)
library(gridExtra)
library(grid)

# ============================================================================
# PART 1: EXTRACT SELECTED GENES FROM EACH MODEL
# ============================================================================

cat("\n" %+% strrep("=", 70) %+% "\n")
cat("MULTI-MODEL FEATURE COMPARISON\n")
cat(strrep("=", 70) %+% "\n\n")

# ──────────────────────────────────────────────────────────────────────────
# LASSO: Extract consensus features (≥70% selection)
# ──────────────────────────────────────────────────────────────────────────

# From your LASSO code, consensus_genes are already defined
# If not, extract from stability analysis:
if(!exists("consensus_genes") || length(consensus_genes) == 0) {
  # Alternative: use features selected in ≥50% of folds
  lasso_threshold <- 0.5  # Adjust if needed
  lasso_consensus_count <- ceiling(K * lasso_threshold)
  lasso_genes <- cons$gene[cons$Freq >= lasso_consensus_count]
} else {
  lasso_genes <- consensus_genes
}

cat(sprintf("LASSO selected features: %d genes\n", length(lasso_genes)))
if(length(lasso_genes) > 0) {
  cat("  ", paste(head(lasso_genes, 10), collapse = ", "), 
      ifelse(length(lasso_genes) > 10, "...\n", "\n"))
}

# ──────────────────────────────────────────────────────────────────────────
# Random Forest: Extract consensus features (≥70% selection)
# ──────────────────────────────────────────────────────────────────────────

# From RF code, use consensus features or top features by importance
if(exists("consensus_features") && length(consensus_features) > 0) {
  rf_genes <- consensus_features
} else {
  # Alternative: use top N features by mean importance
  # Adjust N as needed (e.g., top 10-20 features)
  rf_top_n <- min(15, nrow(imp_summary))
  rf_genes <- head(imp_summary$gene, rf_top_n)
}

cat(sprintf("Random Forest selected features: %d genes\n", length(rf_genes)))
if(length(rf_genes) > 0) {
  cat("  ", paste(head(rf_genes, 10), collapse = ", "), 
      ifelse(length(rf_genes) > 10, "...\n", "\n"))
}

# ──────────────────────────────────────────────────────────────────────────
# SVM-RFE: Extract stable features (≥70% or ≥50% from bootstrap)
# ──────────────────────────────────────────────────────────────────────────

# From SVM-RFE code, stable_features are already defined
if(exists("stable_features") && length(stable_features) > 0) {
  svm_genes <- stable_features
} else {
  # Alternative: use top N features by selection frequency
  svm_top_n <- min(15, nrow(stability_results$stability_ranking))
  svm_genes <- head(stability_results$stability_ranking$Feature, svm_top_n)
}

cat(sprintf("SVM-RFE selected features: %d genes\n", length(svm_genes)))
if(length(svm_genes) > 0) {
  cat("  ", paste(head(svm_genes, 10), collapse = ", "), 
      ifelse(length(svm_genes) > 10, "...\n", "\n"))
}

# ──────────────────────────────────────────────────────────────────────────
# Create gene list for visualization
# ──────────────────────────────────────────────────────────────────────────

gene_list <- list(
  LASSO = lasso_genes,
  RandomForest = rf_genes,
  SVM_RFE = svm_genes
)

# Remove empty sets
gene_list <- gene_list[sapply(gene_list, length) > 0]

cat(sprintf("\n✓ Extracted features from %d models\n", length(gene_list)))

# ──────────────────────────────────────────────────────────────────────────
# Calculate overlaps
# ──────────────────────────────────────────────────────────────────────────

all_genes <- unique(unlist(gene_list))
cat(sprintf("✓ Total unique genes across models: %d\n", length(all_genes)))

# Pairwise overlaps
if(length(gene_list) >= 2) {
  model_names <- names(gene_list)
  cat("\nPairwise Overlaps:\n")
  for(i in 1:(length(gene_list)-1)) {
    for(j in (i+1):length(gene_list)) {
      overlap <- intersect(gene_list[[i]], gene_list[[j]])
      cat(sprintf("  %s ∩ %s: %d genes\n", 
                  model_names[i], model_names[j], length(overlap)))
    }
  }
}

# Core consensus (genes in all models)
if(length(gene_list) == 3) {
  core_genes <- Reduce(intersect, gene_list)
  cat(sprintf("\nCore Consensus (all 3 models): %d genes\n", length(core_genes)))
  if(length(core_genes) > 0) {
    cat("  ", paste(core_genes, collapse = ", "), "\n")
  }
} else if(length(gene_list) == 2) {
  core_genes <- intersect(gene_list[[1]], gene_list[[2]])
  cat(sprintf("\nCore Consensus (both models): %d genes\n", length(core_genes)))
  if(length(core_genes) > 0) {
    cat("  ", paste(core_genes, collapse = ", "), "\n")
  }
}

# ============================================================================
# PART 2: VENN DIAGRAM (EQUAL CIRCLE SIZE)
# ============================================================================

cat("\n" %+% strrep("=", 70) %+% "\n")
cat("GENERATING VENN DIAGRAM\n")
cat(strrep("=", 70) %+% "\n\n")

# Set colors
venn_colors <- c("#3498DB", "#E74C3C", "#27AE60")
names(venn_colors) <- names(gene_list)

# Create Venn diagram based on number of models
if(length(gene_list) == 2) {
  
  # Two-way Venn
  venn_plot <- venn.diagram(
    x = gene_list,
    category.names = names(gene_list),
    filename = NULL,
    output = TRUE,
    
    # Circles - EQUAL SIZE
    lwd = 2,
    lty = 'solid',
    col = venn_colors[1:2],
    fill = venn_colors[1:2],
    alpha = 0.3,
    
    # Circle sizes - force equal
    scaled = FALSE,
    euler.d = FALSE,
    
    # Numbers
    cex = 1.5,
    fontface = "bold",
    fontfamily = "sans",
    
    # Category names
    cat.cex = 1.3,
    cat.fontface = "bold",
    cat.default.pos = "outer",
    cat.pos = c(-27, 27),
    cat.dist = c(0.055, 0.055),
    cat.fontfamily = "sans",
    cat.col = venn_colors[1:2]
  )
  
} else if(length(gene_list) == 3) {
  
  # Three-way Venn
  venn_plot <- venn.diagram(
    x = gene_list,
    category.names = names(gene_list),
    filename = NULL,
    output = TRUE,
    
    # Circles - EQUAL SIZE
    lwd = 2,
    lty = 'solid',
    col = venn_colors,
    fill = venn_colors,
    alpha = 0.3,
    
    # Circle sizes - force equal
    scaled = FALSE,
    euler.d = FALSE,
    
    # Numbers
    cex = 1.3,
    fontface = "bold",
    fontfamily = "sans",
    
    # Category names
    cat.cex = 1.2,
    cat.fontface = "bold",
    cat.default.pos = "outer",
    cat.pos = c(-27, 27, 135),
    cat.dist = c(0.055, 0.055, 0.085),
    cat.fontfamily = "sans",
    cat.col = venn_colors
  )
  
} else {
  stop("Venn diagram requires 2 or 3 models. Currently have: ", length(gene_list))
}

# Save Venn diagram
png("Model_Comparison_Venn.png", width = 8, height = 8, units="in", res=600)
grid.newpage()
grid.draw(venn_plot)
dev.off()

png("Model_Comparison_Venn.png", width = 8, height = 8, units = "in", res = 600)
grid.newpage()
grid.draw(venn_plot)
dev.off()

cat("✓ Venn diagram saved\n")

# ============================================================================
# PART 3: UPSET PLOT
# ============================================================================

cat("\n" %+% strrep("=", 70) %+% "\n")
cat("GENERATING UPSET PLOT\n")
cat(strrep("=", 70) %+% "\n\n")

# Prepare data for UpSet plot
# Create binary matrix: genes × models
upset_df <- data.frame(
  Gene = all_genes,
  stringsAsFactors = FALSE
)

for(model in names(gene_list)) {
  upset_df[[model]] <- as.integer(upset_df$Gene %in% gene_list[[model]])
}

# Save UpSet plot with publication quality
pdf("Model_Comparison_UpSet.pdf", width = 10, height = 6)
upset(
  upset_df,
  sets = names(gene_list),
  order.by = "freq",
  decreasing = TRUE,
  
  # Main bar plot
  main.bar.color = "#2C3E50",
  
  # Matrix
  matrix.color = "#3498DB",
  matrix.dot.alpha = 0.5,
  
  # Set size bar
  sets.bar.color = venn_colors[names(gene_list)],
  
  # Text
  text.scale = c(1.5, 1.3, 1.2, 1.2, 1.5, 1.2),
  
  # Keep empty intersections
  empty.intersections = "on",
  
  # Number of intersections to show
  nintersects = NA,  # Show all
  
  # Point size
  point.size = 3.5,
  line.size = 1.2,
  
  # Set order
  #sets = rev(names(gene_list))  # Reverse for better visualization
)
dev.off()

png("Model_Comparison_UpSet.png", width = 10, height = 6, units = "in", res = 600)
upset(
  upset_df,
  sets = names(gene_list),
  order.by = "freq",
  decreasing = TRUE,
  main.bar.color = "#2C3E50",
  matrix.color = "#3498DB",
  matrix.dot.alpha = 0.5,
  sets.bar.color = venn_colors[names(gene_list)],
  text.scale = c(1.5, 1.3, 1.2, 1.2, 1.5, 1.2),
  empty.intersections = "on",
  nintersects = NA,
  point.size = 3.5,
  line.size = 1.2,
  #sets = rev(names(gene_list))
)
dev.off()

cat("✓ UpSet plot saved\n")

# ============================================================================
# PART 4: ADDITIONAL VISUALIZATIONS
# ============================================================================

cat("\n" %+% strrep("=", 70) %+% "\n")
cat("GENERATING ADDITIONAL COMPARISON PLOTS\n")
cat(strrep("=", 70) %+% "\n\n")

# Set publication theme
theme_publication <- theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
    plot.subtitle = element_text(size = 13, hjust = 0.5, color = "grey30"),
    axis.title = element_text(face = "bold", size = 13),
    axis.text = element_text(size = 11, color = "black"),
    legend.title = element_text(face = "bold", size = 11),
    legend.text = element_text(size = 11),
    panel.grid.minor = element_blank()
  )

# ──────────────────────────────────────────────────────────────────────────
# Plot 1: Bar Chart of Feature Counts
# ──────────────────────────────────────────────────────────────────────────

feature_counts <- data.frame(
  Model = names(gene_list),
  Count = sapply(gene_list, length),
  stringsAsFactors = FALSE
)

feature_counts$Model <- factor(feature_counts$Model, levels = feature_counts$Model)

p_counts <- ggplot(feature_counts, aes(x = Model, y = Count, fill = Model)) +
  geom_bar(stat = "identity", color = "black", size = 0.5, width = 0.6) +
  geom_text(aes(label = Count), vjust = -0.5, fontface = "bold", size = 5) +
  scale_fill_manual(values = venn_colors[feature_counts$Model]) +
  labs(title = "Number of Selected Features per Model",
    subtitle = sprintf("Total unique genes: %d", length(all_genes)),
    x = "Model", y = "Number of Features") +
  theme_publication +
  theme(legend.position = "none") +
  ylim(0, max(feature_counts$Count) * 1.15)

ggsave("Model_Comparison_Feature_Counts.png", p_counts, width = 7, height = 6, units="in", dpi = 1200)

# ──────────────────────────────────────────────────────────────────────────
# Plot 2: Jaccard Similarity Heatmap (if 3 models)
# ──────────────────────────────────────────────────────────────────────────

if(length(gene_list) >= 2) {
  
  # Calculate Jaccard similarity matrix
  n_models <- length(gene_list)
  jaccard_matrix <- matrix(0, nrow = n_models, ncol = n_models)
  rownames(jaccard_matrix) <- names(gene_list)
  colnames(jaccard_matrix) <- names(gene_list)
  
  for(i in 1:n_models) {
    for(j in 1:n_models) {
      set1 <- gene_list[[i]]
      set2 <- gene_list[[j]]
      intersection <- length(intersect(set1, set2))
      union <- length(union(set1, set2))
      jaccard_matrix[i, j] <- if(union > 0) intersection / union else 0
    }
  }
  
  # Convert to long format for ggplot
  jaccard_df <- as.data.frame(as.table(jaccard_matrix))
  colnames(jaccard_df) <- c("Model1", "Model2", "Jaccard")
  
  p_jaccard <- ggplot(jaccard_df, aes(x = Model1, y = Model2, fill = Jaccard)) +
    geom_tile(color = "white", size = 1.5) +
    geom_text(aes(label = sprintf("%.2f", Jaccard)), 
              size = 6, fontface = "bold", color = "white") +
    scale_fill_gradient2(low = "#3498DB", mid = "#F39C12", high = "#E74C3C",
                         midpoint = 0.5, limits = c(0, 1)) +
    labs(
      title = "Jaccard Similarity Between Models",
      subtitle = "Intersection / Union of Selected Features",
      x = "Model",
      y = "Model",
      fill = "Jaccard\nIndex"
    ) +
    theme_publication +
    theme(
      panel.grid = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1)
    ) +
    coord_equal()
  
  ggsave("Model_Comparison_Jaccard_Similarity.png", p_jaccard,
         width = 7, height = 6, unit="in",dpi = 300)
}

# ──────────────────────────────────────────────────────────────────────────
# Plot 3: Expression Heatmap of Core Consensus Genes
# ──────────────────────────────────────────────────────────────────────────

if(exists("core_genes") && length(core_genes) > 0) {
  
  library(pheatmap)
  
  expr_core <- expr_ml[core_genes, , drop = FALSE]
  
  # Annotation
  annotation_col <- data.frame(Condition = condition)
  rownames(annotation_col) <- colnames(expr_core)
  
  ann_colors <- list(
    Condition = c(Normal = "#3498DB", Cancer = "#E74C3C")
  )
  
  p_heatmap <- pheatmap(
    expr_core,
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
    main = sprintf("Expression Heatmap - Core Consensus Genes (n=%d)", 
                   length(core_genes)),
    show_colnames = FALSE
  )
  
  png("Model_Comparison_Core_Genes_Heatmap.png", width = 10, height = 8, unit="in",res=600)
  print(p_heatmap)
  dev.off()
  
  cat("✓ Core consensus genes heatmap saved\n")
}

cat("✓ Additional plots saved\n")

# ============================================================================
# PART 5: SAVE RESULTS TO FILES
# ============================================================================

cat("\n" %+% strrep("=", 70) %+% "\n")
cat("SAVING COMPARISON RESULTS\n")
cat(strrep("=", 70) %+% "\n\n")

# 1. Save selected genes per model
for(model_name in names(gene_list)) {
  gene_df <- data.frame(
    Gene = gene_list[[model_name]],
    Model = model_name,
    stringsAsFactors = FALSE
  )
  write.csv(gene_df, 
            sprintf("Selected_Genes_%s.csv", model_name),
            row.names = FALSE)
}

# 2. Save all unique genes
all_genes_df <- data.frame(
  Gene = all_genes,
  stringsAsFactors = FALSE
)

# Add model membership columns
for(model_name in names(gene_list)) {
  all_genes_df[[model_name]] <- as.integer(all_genes_df$Gene %in% gene_list[[model_name]])
}

# Add total count
all_genes_df$N_Models <- rowSums(all_genes_df[, names(gene_list), drop = FALSE])

# Sort by number of models
all_genes_df <- all_genes_df[order(-all_genes_df$N_Models, all_genes_df$Gene), ]

write.csv(all_genes_df, "All_Selected_Genes_Comparison.csv", row.names = FALSE)

# 3. Save core consensus genes
if(exists("core_genes") && length(core_genes) > 0) {
  core_genes_df <- data.frame(
    Gene = core_genes,
    Selected_By = "All Models",
    stringsAsFactors = FALSE
  )
  write.csv(core_genes_df, "Core_Consensus_Genes.csv", row.names = FALSE)
}

# 4. Save overlap summary
overlap_summary <- data.frame(
  Category = character(),
  Count = integer(),
  Genes = character(),
  stringsAsFactors = FALSE
)

# Individual models
for(model_name in names(gene_list)) {
  overlap_summary <- rbind(
    overlap_summary,
    data.frame(
      Category = model_name,
      Count = length(gene_list[[model_name]]),
      Genes = paste(gene_list[[model_name]], collapse = "; "),
      stringsAsFactors = FALSE
    )
  )
}

# Core consensus
if(exists("core_genes")) {
  overlap_summary <- rbind(
    overlap_summary,
    data.frame(
      Category = "Core_Consensus",
      Count = length(core_genes),
      Genes = paste(core_genes, collapse = "; "),
      stringsAsFactors = FALSE
    )
  )
}

# Pairwise overlaps
if(length(gene_list) >= 2) {
  model_names <- names(gene_list)
  for(i in 1:(length(gene_list)-1)) {
    for(j in (i+1):length(gene_list)) {
      overlap <- intersect(gene_list[[i]], gene_list[[j]])
      overlap_summary <- rbind(
        overlap_summary,
        data.frame(
          Category = sprintf("%s_AND_%s", model_names[i], model_names[j]),
          Count = length(overlap),
          Genes = paste(overlap, collapse = "; "),
          stringsAsFactors = FALSE
        )
      )
    }
  }
}

write.csv(overlap_summary, "Model_Overlap_Summary.csv", row.names = FALSE)

# 5. Save Jaccard similarity matrix
if(exists("jaccard_matrix")) {
  write.csv(jaccard_matrix, "Model_Jaccard_Similarity.csv", row.names = TRUE)
}

cat("✓ All comparison results saved\n")

# ============================================================================
# FINAL SUMMARY
# ============================================================================

cat("\n" %+% strrep("=", 70) %+% "\n")
cat("MULTI-MODEL COMPARISON COMPLETE\n")
cat(strrep("=", 70) %+% "\n\n")

cat("=== FEATURE SELECTION SUMMARY ===\n")
for(model_name in names(gene_list)) {
  cat(sprintf("  %s: %d genes\n", model_name, length(gene_list[[model_name]])))
}
cat(sprintf("  Total unique: %d genes\n", length(all_genes)))
if(exists("core_genes")) {
  cat(sprintf("  Core consensus: %d genes\n", length(core_genes)))
}

if(length(gene_list) >= 2) {
  cat("\n=== PAIRWISE OVERLAPS ===\n")
  model_names <- names(gene_list)
  for(i in 1:(length(gene_list)-1)) {
    for(j in (i+1):length(gene_list)) {
      overlap <- intersect(gene_list[[i]], gene_list[[j]])
      cat(sprintf("  %s ∩ %s: %d genes\n", 
                  model_names[i], model_names[j], length(overlap)))
    }
  }
}

cat("\n=== GENERATED FILES ===\n")
cat("\nVisualization Files:\n")
cat("  - Model_Comparison_Venn.pdf\n")
cat("  - Model_Comparison_Venn.png\n")
cat("  - Model_Comparison_UpSet.pdf\n")
cat("  - Model_Comparison_UpSet.png\n")
cat("  - Model_Comparison_Feature_Counts.pdf\n")
if(exists("jaccard_matrix")) {
  cat("  - Model_Comparison_Jaccard_Similarity.pdf\n")
}
if(exists("core_genes") && length(core_genes) > 0) {
  cat("  - Model_Comparison_Core_Genes_Heatmap.pdf\n")
}

cat("\nData Files:\n")
for(model_name in names(gene_list)) {
  cat(sprintf("  - Selected_Genes_%s.csv\n", model_name))
}
cat("  - All_Selected_Genes_Comparison.csv\n")
if(exists("core_genes") && length(core_genes) > 0) {
  cat("  - Core_Consensus_Genes.csv\n")
}
cat("  - Model_Overlap_Summary.csv\n")
if(exists("jaccard_matrix")) {
  cat("  - Model_Jaccard_Similarity.csv\n")
}

cat("\n" %+% strrep("=", 70) %+% "\n")
cat("Venn diagram uses EQUAL circle sizes for visual consistency.\n")
cat("UpSet plot shows all intersection patterns.\n")
cat(strrep("=", 70) %+% "\n\n")