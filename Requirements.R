# requirements.R
# Required R packages for ML-based hub gene identification from gene expression data
# Author: Saeed
# Date: February 2026

cat("Installing and loading required packages...\n\n")

# Complete list of required packages
required_packages <- c(
  # Data acquisition and preprocessing
  "GEOquery",        # Download GEO datasets
  "edgeR",           # RNA-seq differential expression
  "limma",           # Microarray analysis
  
  # Machine Learning
  "caret",           # ML framework and workflows
  "e1071",           # SVM implementation
  "randomForest",    # Random Forest
  "glmnet",          # LASSO regression
  
  # Performance evaluation
  "pROC",            # ROC analysis
  "PRROC",           # Precision-Recall curves
  "ROCR",            # ROC curves and performance
  "boot",            # Bootstrap resampling
  
  # Visualization - Core
  "ggplot2",         # Main plotting framework
  "ggpubr",          # Publication-ready plots
  "pheatmap",        # Simple heatmaps
  "ComplexHeatmap",  # Advanced heatmaps
  
  # Visualization - Specialized
  "VennDiagram",     # Venn diagrams
  "UpSetR",          # UpSet plots
  "RColorBrewer",    # Color palettes
  "circlize",        # Color mapping for heatmaps
  "scales",          # Scale functions
  
  # Plot arrangement
  "gridExtra",       # Arrange multiple plots
  "grid",            # Base grid graphics
  
  # Data manipulation
  "tidyverse",       # Complete tidyverse (includes dplyr, tidyr, etc.)
  "reshape2"         # Data reshaping
)

# Function to install missing packages
install_if_missing <- function(packages) {
  # Check BiocManager packages separately
  bioc_packages <- c("GEOquery", "edgeR", "limma", "ComplexHeatmap")
  cran_packages <- setdiff(packages, bioc_packages)
  
  # Install BiocManager if needed
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
  }
  
  # Install missing CRAN packages
  new_cran <- cran_packages[!(cran_packages %in% installed.packages()[,"Package"])]
  if(length(new_cran) > 0) {
    cat("Installing CRAN packages:", paste(new_cran, collapse=", "), "\n\n")
    install.packages(new_cran, dependencies = TRUE)
  }
  
  # Install missing Bioconductor packages
  new_bioc <- bioc_packages[!(bioc_packages %in% installed.packages()[,"Package"])]
  if(length(new_bioc) > 0) {
    cat("Installing Bioconductor packages:", paste(new_bioc, collapse=", "), "\n\n")
    BiocManager::install(new_bioc, update = FALSE)
  }
  
  if(length(new_cran) == 0 && length(new_bioc) == 0) {
    cat("✓ All required packages are already installed.\n\n")
  }
}

# Install missing packages
install_if_missing(required_packages)

# Load all packages and check for errors
cat("Loading packages...\n")
load_status <- sapply(required_packages, function(pkg) {
  suppressPackageStartupMessages(
    tryCatch({
      library(pkg, character.only = TRUE)
      return(TRUE)
    }, error = function(e) {
      cat("✗ Failed to load:", pkg, "\n")
      return(FALSE)
    })
  )
})

# Summary
if(all(load_status)) {
  cat("\n✓ All packages loaded successfully!\n\n")
} else {
  cat("\n⚠ Some packages failed to load. Please check the errors above.\n\n")
}

# Display session info
cat("R version:", R.version.string, "\n")
cat("Platform:", R.version$platform, "\n")
cat("\nPackage versions:\n")
for(pkg in required_packages) {
  if(pkg %in% installed.packages()[,"Package"]) {
    cat(sprintf("  %-20s %s\n", pkg, packageVersion(pkg)))
  }
}