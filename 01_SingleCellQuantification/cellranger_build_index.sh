#!/bin/bash
set -e
set -u

if [ $# -ne 3 ]; then
    echo "Usage: $0 <fasta> <gtf> <output_dir>"
    exit 1
fi

if ! command -v cellranger &> /dev/null; then
    echo "cellranger could not be found"
    exit 1
fi

# set variables
FASTA=$1
GTF=$2
OUTPUT_DIR=$3

MEMGB=100

# build index
cellranger mkref --genome=$OUTPUT_DIR --fasta=$FASTA --genes=$GTF --memgb=$MEMGB
