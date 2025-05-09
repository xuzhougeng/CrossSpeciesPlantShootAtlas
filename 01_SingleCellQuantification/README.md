# Quantification of Transcriptome Data on Different Platforms

## BGI Platform (DNBSEQ)

**`dnb4tools_build_index.sh`** – Builds the reference genome index
**`dnbc4tools_count.sh`** – Performs single-cell quantification using DNBSEQ platform data

**Usage example:**

```bash
# Build index
./dnb4tools_build_index.sh /path/to/genome.fa /path/to/gtf species_name

# Quantification
./dnbc4tools_count.sh /path/to/fastq /path/to/index /path/to/output sample_name
```

## BD Platform (BD Rhapsody)

**`STAR_solo_build_index.sh`** – Builds STAR reference genome index
**`STAR_solo_count.sh`** – Performs single-cell quantification using STAR-solo, supporting BD Rhapsody's unique barcode design

**Parameter details:**

* Input: Reference genome directory, R2 file, R1 file
* Output: Sorted BAM file and quantification results
* Features: Supports multi-barcode design, uses EmptyDrops_CR for cell filtering

**Usage example:**

```bash
# Build index
./STAR_solo_build_index.sh /path/to/genome.fa /path/to/gtf /path/to/output_dir

# Quantification
./STAR_solo_count.sh /path/to/star_index sample_R2.fastq.gz sample_R1.fastq.gz
```

## 10X Platform (10X Genomics)

**`cellranger_build_index.sh`** – Builds CellRanger reference genome index
**`cellranger_count.sh`** – Performs single-cell quantification using CellRanger

**Parameter details:**

* Input: fastq directory, reference genome, output directory, sample name
* Output: Quantification results and feature matrix
* Configuration: Defaults to 64 cores and 128 GB memory

**Usage example:**

```bash
# Build index
./cellranger_build_index.sh /path/to/genome.fa /path/to/gtf

# Quantification
./cellranger_count.sh /path/to/fastq /path/to/transcriptome /path/to/output sample_name
```

# WRT Pipeline (Transcriptome-based)

Compared to the conventional genome-based scRNA-seq pipeline (such as WRG), the transcriptome-based WRT pipeline adopts a different data processing strategy. The FASTA file obtained from transcriptome assembly consists entirely of transcript sequences, which are made up of exonic regions only. Therefore, to build an index for the WRT pipeline, it is sufficient to create a simplified GTF file in which each transcript is annotated as a single exon, using the transcriptome FASTA as input.

We provide a convenient script, `create_WRT_gtf.py`, which automatically generates the corresponding GTF file from a transcriptome FASTA file.

Prerequisites

* Python 3.x
* BioPython package (`pip install biopython`)

Generate a Simplified GTF from Transcriptome FASTA

```bash
python create_WRT_gtf.py transcriptome.fasta transcriptome.gtf
```

The resulting GTF file will follow this structure:

```
transcript_id    WRT    exon    1    length    .    +    .    gene_id "gene_transcript_id"; transcript_id "transcript_transcript_id";
```

