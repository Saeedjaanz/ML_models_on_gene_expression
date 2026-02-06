# Data Directory

## Overview
This directory contains gene expression data used for hub gene identification.

## Data Source
- **Dataset**: [Specify GEO accession number, e.g., GSE12345 or your data source]
- **Platform**: [e.g., GPL570 - Affymetrix Human Genome U133 Plus 2.0 Array]
- **Samples**: [Number of disease vs control samples]

## Required Files

### Expression Matrix
- **File**: `expression_matrix.csv`
- **Format**: Genes as rows, samples as columns
- **Content**: Normalized gene expression values
```
Gene,Sample1,Sample2,Sample3,...
GENE1,5.234,6.112,4.891,...
GENE2,3.456,2.987,3.678,...
```

### Sample Labels
- **File**: `sample_labels.csv`
- **Format**: Sample IDs with class labels
```
Sample,Class
Sample1,Disease
Sample2,Control
Sample3,Disease
```

## Data Preprocessing

Data preprocessing includes:
- Quality control and normalization
- Batch effect correction (if applicable)
- Filtering of low-expressed genes
- Log2 transformation

See `Code/01_Data_Preprocessing.R` for details.

## Note
Gene expression data files are not included in this repository due to:
- File size constraints (GitHub limit: 100MB per file)
- Data privacy/sharing agreements

**To obtain data:**
1. Contact: [Your email]
2. Or download from GEO: [Link if public]
3. Or use your own gene expression dataset following the format above