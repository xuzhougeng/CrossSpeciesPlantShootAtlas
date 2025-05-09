> The pipeline code was modified from https://github.com/snap-stanford/SATURN

# SATURN Cross-Species Integration

This module provides a pipeline for cross-species integration of single-cell transcriptomics data using SATURN.

## Overview

SATURN enables integration of single-cell data across species by leveraging protein sequence embeddings to create a unified cross-species cell atlas.

## Pipeline Components

* **data_loader.ipynb**: Prepares protein embeddings and organizes input data for integration
* **extract.py**: Extracts protein embeddings using ESM2 models
* **post-analysis.ipynb**: Processes integration results with visualization and clustering

## Usage

1. **Prepare Protein Embeddings**:
   - Process protein sequences for each species using ESM2 embeddings
   - Run `extract.py` to generate embeddings for each species

2. **Data Integration**:
   - Use the prepared embeddings to integrate cross-species data
   - Run the integration pipeline to generate a unified atlas

3. **Post-Analysis**:
   - Perform dimensionality reduction (PCA, UMAP)
   - Analyze cell clusters with Leiden algorithm
   - Visualize results by species and cell types

## Note

1. The pipeline requires a GPU with approximately 60GB of VRAM
2. Create a conda environment named SATURN and configure it according to the official documentation