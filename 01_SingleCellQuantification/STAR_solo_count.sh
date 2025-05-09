#!/bin/bash
set -e
set -u

if [ $# -ne 3 ]; then
    echo "Usage: $0 <ref> <r2> <r1>"
    exit 1
fi
# set variables
REF=$1
R2=$2
R1=$3

# check STAR command
if ! command -v STAR &> /dev/null; then
    echo "STAR could not be found"
    exit 1
fi

BAM="--outSAMtype BAM SortedByCoordinate --outBAMsortingBinsN 500 --limitBAMsortRAM 60000000000 --outMultimapperOrder Random --runRNGseed 1 --outSAMattributes NH HI AS nM CB UB CR CY UR UY GX GN"
#BAM="--outSAMtype None"
GZIP="--readFilesCommand zcat"


BC1=BD/CLS1.txt
BC2=BD/CLS2.txt
BC3=BD/CLS3.txt
CPUS=30


# 按需调整
#ulimit -n 4096
ulimit -n 20480

STAR --runThreadN $CPUS \
     --genomeDir $REF \
     --outFileNamePrefix ./star_out \
     --readFilesIn $R2 $R1 \
     $GZIP $BAM \
     --soloType CB_UMI_Complex \
     --soloAdapterSequence GTGANNNNNNNNNGACA \
     --soloCBposition 2_-9_2_-1 2_4_2_12 2_17_2_25 \
     --soloUMIposition 3_10_3_17 \
     --soloCBwhitelist $BC1 $BC2 $BC3 \
     --soloCBmatchWLtype 1MM \
     --soloCellFilter EmptyDrops_CR