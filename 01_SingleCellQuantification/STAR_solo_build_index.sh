
#!/bin/bash
set -e
set -u

if [ $# -ne 3 ]; then
    echo "Usage: $0 <fasta> <gtf> <output_dir>"
    exit 1
fi

if ! command -v STAR &> /dev/null; then
    echo "STAR could not be found"
    exit 1
fi

# set variables
FASTA=$1
GTF=$2
OUTPUT_DIR=$3

# build index
STAR --runMode genomeGenerate \
     --genomeDir $OUTPUT_DIR \
     --genomeFastaFiles $FASTA \
     --sjdbGTFfile $GTF \
     --runThreadN 30 \
     --sjdbOverhang 100
