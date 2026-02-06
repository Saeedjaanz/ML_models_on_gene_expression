# Machine Learning Models for Hub Gene Identification from Gene Expression Data

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
![R Version](https://img.shields.io/badge/R-%E2%89%A5%204.0.0-blue)

## 📋 Overview

This repository implements three **gold-standard machine learning methods** for identifying hub/marker genes from gene expression data:

- **🔵 LASSO** - L1-regularized regression with coefficient path analysis
- **🟢 Random Forest** - Ensemble learning with feature importance ranking  
- **🔴 SVM-RFE** - Support Vector Machine with Recursive Feature Elimination

The analysis follows rigorous protocols including nested cross-validation, stability selection, and consensus feature identification across all methods.

---

## 🗂️ Repository Structure
```
ML_models_on_gene_expression/
├── README.md                         # This file
├── .gitignore                        # Git ignore rules
├── requirements.R                    # R package dependencies
├── LICENSE                           # MIT License
│
├── Code/                             # Analysis scripts
│   ├── 01_Data_Preprocessing.R       # Data loading and preprocessing
│   ├── 02_LASSO.R                    # LASSO with CV lambda selection
│   ├── 03_Random_Forest.R            # Random Forest with mtry tuning
│   ├── 04_SVM_RFE.R                  # SVM-RFE with nested CV
│   └── 05_Model_Comparison.R         # Venn/UpSet plots and consensus
│
├── Data/                             # Gene expression data
│   └── README.md                     # Data format and source info
│
└── Results/                          # All outputs
    ├── Figures/                      # Visualizations
    │   ├── 02_LASSO/                 # LASSO plots
    │   ├── 03_Random_Forest/         # RF plots
    │   ├── 04_SVM_RFE/               # SVM-RFE plots
    │   └── 05_Model_Comparison/      # Comparison plots
    │
    └── Results and Metrices/         # CSV outputs and metrics
        ├── 02_LASSO/                 # Selected features, coefficients
        ├── 03_Random_Forest/         # Feature importance, models
        ├── 04_SVM_RFE/               # Stable features, rankings
        └── 05_Model_Comparison/      # Consensus genes, overlap
```

---

## 🔬 Methods & Protocols

### 1️⃣ LASSO (Least Absolute Shrinkage and Selection Operator)

**Key Features:**
- 10-fold cross-validation for lambda selection
- Coefficient path analysis across regularization strengths
- Feature stability assessment across folds
- Selection of features with non-zero coefficients at lambda.min and lambda.1se

**Outputs:**
- Selected gene lists
- Coefficient heatmaps
- Cross-validation curves
- Feature stability metrics
- ROC curves and performance metrics

### 2️⃣ Random Forest

**Key Features:**
- Feature importance via Mean Decrease Accuracy (MDA) and Gini index
- 10-fold cross-validation with mtry optimization
- Consensus features (≥70% fold appearance threshold)
- Per-fold diagnostic plots

**Outputs:**
- Feature importance rankings
- Stability analysis across folds
- Confusion matrices
- ROC curves
- Expression heatmaps for top genes

### 3️⃣ SVM-RFE (Support Vector Machine - Recursive Feature Elimination)

**Key Features:**
- Nested cross-validation for unbiased evaluation
- Recursive feature elimination to identify minimal gene sets
- Stability selection (70% threshold across folds)
- Feature co-selection analysis

**Outputs:**
- Feature elimination curves
- Stable feature rankings
- Nested CV performance distributions
- Feature co-selection matrices
- ROC curves and confusion matrices

### 4️⃣ Model Comparison & Consensus Analysis

**Key Features:**
- Venn diagrams (3-way overlap visualization)
- UpSet plots (multi-set intersection analysis)
- Jaccard similarity coefficients
- Core consensus genes (selected by all 3 methods)

**Outputs:**
- Comprehensive comparison tables
- Overlap visualizations
- Core consensus gene list
- Expression heatmaps for consensus genes

---

## 🚀 Getting Started

### Prerequisites

- **R** version ≥ 4.0.0
- **RStudio** (recommended)
- **Git** for version control

### Installation

1. **Clone this repository:**
```bash
git clone https://github.com/Saeedjaanz/ML_models_on_gene_expression.git
cd ML_models_on_gene_expression
```

2. **Install required R packages:**

Open R/RStudio and run:
```R
source("requirements.R")
```

This will automatically install all required packages from CRAN and Bioconductor.

3. **Prepare your data:**

Place your gene expression data in the `Data/` folder following the format specified in `Data/README.md`.

---

## 📊 Usage

### Running the Complete Pipeline

Execute scripts in sequence:
```R
# 1. Preprocess data
source("Code/01_Data_Preprocessing.R")

# 2. Run LASSO
source("Code/02_LASSO.R")

# 3. Run Random Forest
source("Code/03_Random_Forest.R")

# 4. Run SVM-RFE
source("Code/04_SVM_RFE.R")

# 5. Compare models and find consensus
source("Code/05_Model_Comparison.R")
```

### Running Individual Methods

Each script can be run independently after data preprocessing:
```R
# Run only LASSO
source("Code/01_Data_Preprocessing.R")
source("Code/02_LASSO.R")
```

**Note:** The model comparison script (05) requires outputs from all three methods (02, 03, 04).

---

## 📈 Key Results

### Model Performance Summary

| Method         | Accuracy | Sensitivity | Specificity | AUC   | Selected Features |
|----------------|----------|-------------|-------------|-------|-------------------|
| **LASSO**      | 0.873    | 0.911       | 0.842       | 0.951 | 8 genes           |
| **RandomForest** | 0.897  | 0.898       | 0.899       | 0.964 | 8 genes           |
| **SVM-RFE**    | 0.842    | 0.898       | 0.793       | 0.907 | 8 genes           |

### Output Files Guide

#### Selected Features

**LASSO:**
- `Selected_Genes_LASSO.csv` - Final gene list
- `LASSO_Consensus_Features_70pct.csv` - Stable features across folds
- `LASSO_Final_Coefficients_*.csv` - Coefficient values

**Random Forest:**
- `Selected_Genes_RandomForest.csv` - Final gene list
- `RF_Consensus_Features_70pct.csv` - Consensus features
- `RF_Feature_Importance_Summary.csv` - Importance scores

**SVM-RFE:**
- `Selected_Genes_SVM_RFE.csv` - Final gene list
- `SVM_RFE_Stable_Features_70pct.csv` - Stable features
- `SVM_RFE_Feature_Stability_Ranking.csv` - Stability scores

#### Consensus Analysis
- `Core_Consensus_Genes.csv` - **Genes selected by all 3 methods** ⭐
- `All_Selected_Genes_Comparison.csv` - Complete comparison table
- `Model_Overlap_Summary.csv` - Method overlap statistics

---

## 📁 Important Figures

### LASSO Visualizations
- `LASSO_CV_Lambda_Selection.png` - Cross-validation curve
- `LASSO_Coefficient_Paths.png` - Coefficient paths
- `LASSO_Feature_Stability.png` - Stability across folds
- `LASSO_ROC_Curve.png` - ROC curve
- `LASSO_Coefficient_Heatmap_Consensus.png` - Consensus gene heatmap

### Random Forest Visualizations
- `RF_Feature_Importance_MDA.png` - Mean Decrease Accuracy
- `RF_Feature_Importance_Gini.png` - Gini importance
- `RF_Feature_Stability.png` - Feature stability
- `RF_ROC_Curve.png` - ROC curve
- `RF_Expression_Heatmap_Top10.png` - Top 10 genes

### SVM-RFE Visualizations
- `SVM_RFE_Feature_Elimination.png` - Elimination curve
- `SVM_RFE_Feature_Stability.png` - Stability heatmap
- `SVM_RFE_NestedCV_Distribution.png` - Performance distribution
- `SVM_RFE_ROC_Curve.png` - ROC curve
- `SVM_RFE_Expression_Heatmap_Stable.png` - Stable genes

### Model Comparison
- `Model_Comparison_Venn.png` - 3-way Venn diagram
- `Model_Comparison_UpSet.png` - UpSet plot
- `Model_Comparison_Jaccard_Similarity.png` - Similarity heatmap
- `Model_Comparison_Core_Genes_Heatmap.png` - Consensus genes

---

## 🔄 Reproducibility

All analyses ensure reproducibility through:

- **Fixed random seeds**: `set.seed(123)` in all scripts
- **Cross-validation**: 10-fold CV for all methods
- **Nested CV**: For SVM-RFE hyperparameter tuning
- **Stability thresholds**: 70% consensus for feature selection
- **Documented parameters**: All hyperparameters logged

---

## 📦 Dependencies

### Core R Packages

**Data Processing:**
- GEOquery, edgeR, limma, tidyverse

**Machine Learning:**
- caret, e1071, randomForest, glmnet

**Evaluation:**
- pROC, PRROC, ROCR, boot

**Visualization:**
- ggplot2, ggpubr, pheatmap, ComplexHeatmap, VennDiagram, UpSetR

See `requirements.R` for complete list with versions.

---

## 📝 Citation

If you use this code in your research, please cite:
```bibtex
@software{Saeed2026MLGenes,
  author = {Saeed},
  title = {Machine Learning Models for Hub/Marker Gene Identification from Gene Expression Data},
  year = {2026},
  url = {https://github.com/Saeedjaanz/ML_models_on_gene_expression}
}
```

---

## 📄 License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

---

## 👤 Author

**Abdullah**

- GitHub: [@Saeedjaanz](https://github.com/Saeedjaanz)
- Email: Saeedjaan.smd.fsd@gmail.com

---

## 🙏 Acknowledgments

- Gold-standard protocols for feature selection in high-dimensional biological data
- Nested cross-validation methodologies for unbiased ML evaluation
- Stability-based consensus approaches for robust feature selection

---

## 📞 Support

For questions or issues:
1. Open an [Issue](https://github.com/Saeedjaanz/ML_models_on_gene_expression/issues)
2. Email: Saeedjaan.smd.fsd@gmail.com

---

**Last Updated:** February 2026
