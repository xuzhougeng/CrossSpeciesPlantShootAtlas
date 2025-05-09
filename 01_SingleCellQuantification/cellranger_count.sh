#!/bin/bash

# CellRanger pipeline for single-cell RNA-seq data processing

if [ $# -ne 4 ]; then
    echo "Usage: $0 <fastq_dir> <transcriptome> <output_dir> <sample_name>"
    exit 1
fi

# check cellranger command
if ! command -v cellranger &> /dev/null; then
    echo "cellranger could not be found"
    exit 1
fi

# Set variables
FASTQ_DIR=$1
TRANSCRIPTOME=$2
OUTPUT_DIR=$3
SAMPLE_NAME=$4

# set cores and memory
CORES=64
MEMORY=128

# Run cellranger count for single sample
cellranger count \
    --id=$SAMPLE_NAME \
    --transcriptome=$TRANSCRIPTOME \
    --fastqs=$FASTQ_DIR \
    --localcores=$CORES \
    --localmem=$MEMORY \
    --expect-cells=10000
