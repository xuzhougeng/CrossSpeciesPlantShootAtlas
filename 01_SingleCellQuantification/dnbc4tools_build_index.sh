#!/bin/bash
set -e
set -u

if [ $# -ne 3 ]; then
    echo "Usage: $0 <fasta> <gtf> <species>"
    exit 1
fi

if ! command -v dnbc4tools &> /dev/null; then
    echo "dnbc4tools could not be found"
    exit 1
fi

# set variables
FASTA=$1
GTF=$2
SPECIES=$3

# build index
dnbc4tools  rna mkref --ingtf $GTF  --fasta $FASTA  --threads 30 --species $SPECIES
