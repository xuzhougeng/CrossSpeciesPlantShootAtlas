#!/bin/bash
set -e
set -u

if [ $# -ne 7 ]; then
    echo "Usage: $0 <cDNAfastq1> <cDNAfastq2> <oligofastq1> <oligofastq2> <genomeDir> <sampleName> <threads>"
    exit 1
fi

if ! command -v dnbc4tools &> /dev/null; then
    echo "dnbc4tools could not be found"
    exit 1
fi

# set variables
cDNAfastq1=$1
cDNAfastq2=$2
oligofastq1=$3
oligofastq2=$4

genomeDir=$5
sampleName=$6
threads=$7


# count
dnbc4tools rna run \
  --cDNAfastq1 $cDNAfastq1 \
  --cDNAfastq2 $cDNAfastq2 \
  --oligofastq1 $oligofastq1 \
  --oligofastq2 $oligofastq2 \
  --genomeDir $genomeDir  \
  --name $sampleName --threads $threads