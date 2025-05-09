# AUCell Scoring

This script calculates cell type scores using the AUCell algorithm based on marker gene sets.

## Prerequisites

Required R packages:
- Seurat
- AUCell
- GSEABase

Install required packages:
```r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("AUCell")
install.packages("Seurat")
BiocManager::install("GSEABase")
```

## Input files

1. **Seurat RDS file**: Contains the Seurat object for analysis
2. **Geneset CSV file**: Contains marker genes for each cell type with the following format:

```
cell_type,marker_gene
PC S,Gene1
PC S,Gene2
Neuronal,Gene3
...
```

## Usage

```bash
Rscript aucell_score.r <seurat_rds_file> <geneset_csv_file>
```

Example:
```bash
Rscript aucell_score.r data/seurat_object.rds data/marker_genes.csv
```

## Output

For each cell type in the CSV file:
- Adds a score to the Seurat object as `score_[cell_type]`
- Generates a feature plot saved as `score_[cell_type].png`
