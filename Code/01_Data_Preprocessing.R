# ============================================================================
# GSE144269 DATA PREPARATION FOR MACHINE LEARNING
# Hepatocellular Carcinoma (HCC) Paired Samples
# Prepares expr_ml matrix and condition factor for ML models
# ============================================================================

# Required Libraries
library(GEOquery)
library(edgeR)
library(limma)

# ============================================================================
# PART 1: DATA DOWNLOAD AND LOADING
# ============================================================================

cat("\n" %+% strrep("=", 70) %+% "\n")
cat("DOWNLOADING AND LOADING GSE144269 DATA\n")
cat(strrep("=", 70) %+% "\n\n")

# Download supplementary files
getGEOSuppFiles("GSE144269")

# Read count data
count_file <- "GSE144269/GSE144269_RSEM_GeneCounts.txt.gz"
counts <- read.delim(gzfile(count_file), row.names = 1, check.names = FALSE)

cat(sprintf("✓ Loaded count matrix: %d genes × %d samples\n", 
            nrow(counts), ncol(counts)))

# ============================================================================
# PART 2: METADATA PREPARATION
# ============================================================================

cat("\n" %+% strrep("=", 70) %+% "\n")
cat("PREPARING METADATA\n")
cat(strrep("=", 70) %+% "\n\n")

# Download metadata from GEO
gset <- getGEO("GSE144269", GSEMatrix = TRUE, AnnotGPL = FALSE)
meta <- pData(phenoData(gset[[1]]))

# Extract sample information
meta$sample_id <- meta$geo_accession
meta$condition <- sub("tumor/non-tumor: ", "", meta$`tumor/non-tumor:ch1`)
meta$patient_id <- sub("pat ", "", sapply(strsplit(meta$title, " "), `[`, 2))

# Extract sample codes from title [brackets]
meta$sample_code <- ifelse(
  grepl("\\[.*\\]", meta$title),
  sub(".*\\[(.*)\\].*", "\\1", meta$title),
  as.character(meta$geo_accession)
)

cat(sprintf("✓ Metadata prepared: %d samples\n", nrow(meta)))
cat(sprintf("  - Condition distribution:\n"))
print(table(meta$condition))

# ============================================================================
# PART 3: GENE NAME CLEANING AND AGGREGATION
# ============================================================================

cat("\n" %+% strrep("=", 70) %+% "\n")
cat("CLEANING GENE NAMES\n")
cat(strrep("=", 70) %+% "\n\n")

# Extract clean gene names (remove ENSG IDs)
orig_rownames <- rownames(counts)
clean_genes <- sub("^[^|]*\\|", "", orig_rownames)

# Convert to numeric matrix
counts_mat <- as.matrix(counts)
if(!is.numeric(counts_mat)) {
  counts_mat <- apply(counts_mat, 2, as.numeric)
}

# Aggregate duplicate gene names by summing counts
if(any(duplicated(clean_genes))) {
  cat("⚠ Found duplicated gene names. Aggregating by summing counts...\n")
  counts_agg_mat <- rowsum(counts_mat, group = clean_genes)
  counts <- as.data.frame(counts_agg_mat)
  rownames(counts) <- rownames(counts_agg_mat)
} else {
  counts <- as.data.frame(counts_mat)
  rownames(counts) <- clean_genes
}

cat(sprintf("✓ Gene names cleaned: %d unique genes\n", nrow(counts)))

# ============================================================================
# PART 4: ALIGN COUNTS AND METADATA
# ============================================================================

cat("\n" %+% strrep("=", 70) %+% "\n")
cat("ALIGNING COUNTS AND METADATA\n")
cat(strrep("=", 70) %+% "\n\n")

# Find common samples
counts_cols <- colnames(counts)
meta_codes <- meta$sample_code

common_samples <- intersect(counts_cols, meta_codes)

if(length(common_samples) == 0) {
  stop("ERROR: No overlapping sample names found between counts and metadata!")
}

cat(sprintf("  - Samples in counts: %d\n", length(counts_cols)))
cat(sprintf("  - Samples in metadata: %d\n", length(meta_codes)))
cat(sprintf("  - Common samples: %d\n", length(common_samples)))

# Align metadata to counts
meta_aligned <- meta[match(common_samples, meta$sample_code), , drop = FALSE]

# Align counts to metadata
counts_aligned <- counts[, meta_aligned$sample_code, drop = FALSE]

# Verify alignment
if(!all(colnames(counts_aligned) == meta_aligned$sample_code)) {
  stop("ERROR: Sample alignment failed!")
}

# Update objects
counts <- counts_aligned
meta <- meta_aligned

cat(sprintf("✓ Alignment complete: %d matched samples\n", ncol(counts)))

# ============================================================================
# PART 5: DIFFERENTIAL EXPRESSION ANALYSIS
# ============================================================================

cat("\n" %+% strrep("=", 70) %+% "\n")
cat("DIFFERENTIAL EXPRESSION ANALYSIS\n")
cat(strrep("=", 70) %+% "\n\n")

# Create DGEList object
dge <- DGEList(counts = counts)
dge <- calcNormFactors(dge)

# Filter low-expression genes
meta$condition_factor <- factor(meta$condition, levels = c("non-tumor", "tumor"))
keep <- filterByExpr(dge, group = meta$condition_factor)
dge <- dge[keep, , keep.lib.sizes = FALSE]

cat(sprintf("✓ Filtered genes: %d genes retained\n", nrow(dge)))

# Design matrix (paired samples: patient + condition)
meta$patient_factor <- factor(meta$patient_id)
design <- model.matrix(~ patient_factor + condition_factor, data = meta)
rownames(design) <- meta$sample_code

# Estimate dispersion and voom transformation
dge <- estimateDisp(dge, design)
v <- voom(dge, design, plot = FALSE)

# Fit linear model
fit <- lmFit(v, design)
fit <- eBayes(fit)

# Extract differential expression results
de_results <- topTable(fit, coef = "condition_factortumor", 
                       number = Inf, sort.by = "P")

de_results$gene <- rownames(de_results)

cat(sprintf("✓ Differential expression complete: %d genes tested\n", 
            nrow(de_results)))

# Significant genes
de_sig <- de_results[de_results$adj.P.Val < 0.05 & abs(de_results$logFC) > 1, ]
cat(sprintf("✓ Significant DE genes (FDR<0.05, |logFC|>1): %d genes\n", 
            nrow(de_sig)))

# ============================================================================
# PART 6: PREPARE ML INPUT DATA
# ============================================================================

cat("\n" %+% strrep("=", 70) %+% "\n")
cat("PREPARING DATA FOR MACHINE LEARNING MODELS\n")
cat(strrep("=", 70) %+% "\n\n")

# Define hub genes (same as GSE40419 analysis)
hub_genes <- c("AKT1", "CASP3", "CTNNB1", "MMP9", 
               "ESR1", "SRC", "EGFR", "BCL2", "HSP90AA1")

# Get normalized expression matrix (log2-CPM from voom)
expr_normalized <- v$E

# Check which hub genes are present
genes_in_expr <- rownames(expr_normalized)
hub_genes_present <- intersect(hub_genes, genes_in_expr)
missing_hubs <- setdiff(hub_genes, genes_in_expr)

if(length(missing_hubs) > 0) {
  cat("⚠ Missing hub genes:\n")
  cat("  ", paste(missing_hubs, collapse = ", "), "\n")
}

cat(sprintf("✓ Hub genes present: %d/%d\n", 
            length(hub_genes_present), length(hub_genes)))
cat("  ", paste(hub_genes_present, collapse = ", "), "\n")

# Subset expression matrix to hub genes
expr_ml <- expr_normalized[hub_genes_present, , drop = FALSE]

# Create condition factor (MUST match ML model requirements)
# ML models expect levels: c("Normal", "Cancer")
condition <- factor(
  ifelse(meta$condition_factor == "tumor", "Cancer", "Normal"),
  levels = c("Normal", "Cancer")
)

# Verify alignment
if(!all(colnames(expr_ml) == meta$sample_code)) {
  stop("ERROR: Expression matrix and metadata are not aligned!")
}

if(length(condition) != ncol(expr_ml)) {
  stop("ERROR: Condition vector length does not match expression matrix!")
}

# ============================================================================
# FINAL VERIFICATION AND SUMMARY
# ============================================================================

cat("\n" %+% strrep("=", 70) %+% "\n")
cat("FINAL DATA SUMMARY\n")
cat(strrep("=", 70) %+% "\n\n")

cat("Expression Matrix (expr_ml):\n")
cat(sprintf("  - Dimensions: %d genes × %d samples\n", 
            nrow(expr_ml), ncol(expr_ml)))
cat(sprintf("  - Genes: %s\n", paste(rownames(expr_ml), collapse = ", ")))
cat(sprintf("  - Data type: %s\n", class(expr_ml)))
cat(sprintf("  - Value range: [%.2f, %.2f]\n", 
            min(expr_ml), max(expr_ml)))

cat("\nCondition Factor (condition):\n")
cat(sprintf("  - Length: %d\n", length(condition)))
cat(sprintf("  - Levels: %s\n", paste(levels(condition), collapse = ", ")))
cat(sprintf("  - Distribution:\n"))
print(table(condition))

cat("\n" %+% strrep("=", 70) %+% "\n")
cat("✓ DATA PREPARATION COMPLETE\n")
cat("✓ Ready for LASSO, Random Forest, and SVM-RFE analysis\n")
cat(strrep("=", 70) %+% "\n\n")

# ============================================================================
# VERIFY ML MODEL REQUIREMENTS
# ============================================================================

# Check all requirements are met
stopifnot(
  is.matrix(expr_ml),
  length(condition) == ncol(expr_ml),
  all(levels(condition) %in% c("Normal", "Cancer"))
)

cat("✓ All ML model requirements verified:\n")
cat("  [✓] expr_ml is a matrix\n")
cat("  [✓] condition length matches expr_ml columns\n")
cat("  [✓] condition levels are 'Normal' and 'Cancer'\n")
cat("\nYou can now run LASSO, Random Forest, and SVM-RFE models!\n\n")
