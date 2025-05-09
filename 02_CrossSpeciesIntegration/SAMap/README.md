> The pipeline code was modified from https://github.com/atarashansky/SAMap

# SAMap Cross-Species Integration

This module provides a pipeline for cross-species integration of single-cell transcriptomics data using SAMap.

## Overview

SAMap enables integration of single-cell data across species by identifying homologous genes and aligning datasets.

## Pipeline Components

* **SAMap_prepare.ipynb**: Data preprocessing
* **SAMap_integration.ipynb**: Cross-species integration
* **pairwise_blastp.sh**: Script for generating gene homology maps between species

## Usage

1. **Data Preprocessing**:
   - Prepare input data in appropriate format
   - Process gene homology information

2. **Homology Mapping**:
   - Place protein sequence files (*.pep) in the `processed` directory
   - Run the pairwise BLAST script to generate homology maps:
     ```
     bash pairwise_blastp.sh
     ```
   - This script:
     - Creates BLAST databases for each protein file
     - Performs pairwise BLAST searches between all species
     - Outputs homology maps in `maps/` directory with E-value threshold 1e-6

3. **Integration**:
   - Align datasets across species using the generated homology maps
   - Generate integrated embedding

For more analysis please read https://github.com/atarashansky/SAMap/blob/main/SAMap_vignette.ipynb

